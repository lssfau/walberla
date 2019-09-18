#include "field/GhostLayerField.h"
#include "core/mpi/MPIManager.h"
#include "Fft.h"

namespace walberla {
namespace fft {

template <typename Field_T>
FourierTransform<Field_T>::FourierTransform( shared_ptr< StructuredBlockForest > & blocks, BlockDataID fieldId,
                                             const std::function<real_t(uint_t,uint_t,uint_t)>& greens )
: blocks_(blocks), fieldId_(fieldId), greens_()
{
#ifdef WALBERLA_USE_PFFT
   pfft_init();
#else
#ifdef WALBERLA_BUILD_WITH_MPI
   WALBERLA_CHECK_EQUAL(MPIManager::instance()->numProcesses(), 1, "To use multiple MPI processes, you need PFFT");
#endif
#endif

   WALBERLA_CHECK_EQUAL(blocks_->size(), 1, "FFT needs one block per process");
   
   ptrdiff_t n[3];
   n[0] = int_c( blocks->getNumberOfXCells() );
   n[1] = int_c( blocks->getNumberOfYCells() );
   n[2] = int_c( blocks->getNumberOfZCells() );
   
   WALBERLA_ASSERT(blocks_->isXPeriodic());
   WALBERLA_ASSERT(blocks_->isYPeriodic());
   WALBERLA_ASSERT(blocks_->isZPeriodic());
   
#if !defined(WALBERLA_USE_PFFT) || !defined(NDEBUG)
   auto block = blocks_->begin();
   Field_T *data = block->getData< Field_T >( fieldId_ );
#endif
   
#ifdef WALBERLA_USE_PFFT
   ptrdiff_t local_ni[3];
   ptrdiff_t local_i_start[3];
   ptrdiff_t local_o_start[3];
   MPI_Comm comm = MPIManager::instance()->comm();
   ptrdiff_t alloc_local = pfft_local_size_dft_r2c_3d(n, comm, PFFT_TRANSPOSED_NONE,
                                                      local_ni, local_i_start, local_no, local_o_start);
   
   // make sure that PFFT and Walberla use the same domain decomposition
   WALBERLA_ASSERT_EQUAL(local_ni[0], data->xSize());
   WALBERLA_ASSERT_EQUAL(local_ni[1], data->ySize());
   WALBERLA_ASSERT_EQUAL(local_ni[2], data->zSize());
   WALBERLA_ASSERT_EQUAL(local_i_start[0], int_c(block->getAABB().xMin()));
   WALBERLA_ASSERT_EQUAL(local_i_start[1], int_c(block->getAABB().yMin()));
   WALBERLA_ASSERT_EQUAL(local_i_start[2], int_c(block->getAABB().zMin()));
   
   in_ = FFTReal(pfft_alloc_real(2 * uint_c(alloc_local)), pfft_free);
   out_ = FFTComplex(pfft_alloc_complex(uint_c(alloc_local)), pfft_free);
   if (greens)
      greens_ = FFTReal(pfft_alloc_real(uint_c(alloc_local)+1), pfft_free);
   
   plan_forward_  = FFTPlan(pfft_plan_dft_r2c_3d(n, in_.get(), out_.get(), comm,
                                         PFFT_FORWARD, PFFT_TRANSPOSED_NONE | PFFT_DESTROY_INPUT), pfft_destroy_plan);
   plan_backward_ = FFTPlan(pfft_plan_dft_c2r_3d(n, out_.get(), in_.get(), comm,
                                         PFFT_BACKWARD, PFFT_TRANSPOSED_NONE | PFFT_DESTROY_INPUT), pfft_destroy_plan);
#else
   ptrdiff_t local_o_start[3] = {0,0,0};
   local_no[0] = int_c(data->xSize());
   local_no[1] = int_c(data->ySize());
   local_no[2] = int_c(data->zSize())/2+1;
   ptrdiff_t alloc_local = local_no[0]*local_no[1]*local_no[2];
   
   in_ = FFTReal(fftw_alloc_real(2 * uint_c(alloc_local)), fftw_free);
   out_ = FFTComplex(fftw_alloc_complex(uint_c(alloc_local)), fftw_free);
   if (greens)
      greens_ = FFTReal(fftw_alloc_real(uint_c(alloc_local)+1), fftw_free);
   
   plan_forward_ = FFTPlan(fftw_plan_dft_r2c_3d(int_c(n[0]),int_c(n[1]),int_c(n[2]), in_.get(), out_.get(), FFTW_DESTROY_INPUT), fftw_destroy_plan);
   plan_backward_ = FFTPlan(fftw_plan_dft_c2r_3d(int_c(n[0]),int_c(n[1]),int_c(n[2]), out_.get(), in_.get(), FFTW_DESTROY_INPUT), fftw_destroy_plan);
#endif
   
   if (greens)
      for (ptrdiff_t x = 0; x < local_no[0]; ++x)
      for (ptrdiff_t y = 0; y < local_no[1]; ++y)
      for (ptrdiff_t z = 0; z < local_no[2]; ++z)
         greens_[size_t(z + local_no[2] * (y + local_no[1] * x))] = greens(uint_c(x+local_o_start[0]),uint_c(y+local_o_start[1]),uint_c(z+local_o_start[2]));
   
   WALBERLA_ASSERT_NOT_NULLPTR(plan_forward_.get());
   WALBERLA_ASSERT_NOT_NULLPTR(plan_backward_.get());
}

template <typename Field_T>
void FourierTransform<Field_T>::operator() ()
{
   auto block = blocks_->begin();
   Field_T *data = block->getData< Field_T >( fieldId_ );
   
   WALBERLA_FOR_ALL_CELLS_XYZ(data, {
      in_[uint_c(z) + data->zSize() * (uint_c(y) + data->ySize() * uint_c(x))] = data->get(x,y,z);
   });
   
#ifdef WALBERLA_USE_PFFT
   pfft_execute(plan_forward_.get());
#else
   fftw_execute(plan_forward_.get());
#endif
   
   if (greens_)
   {
      for (size_t i = 0; i < size_t(local_no[0]*local_no[1]*local_no[2]); ++i)
      {
         out_[i][0] *= greens_[i];
         out_[i][1] *= greens_[i];
      }
   }

#ifdef WALBERLA_USE_PFFT
   pfft_execute(plan_backward_.get());
#else
   fftw_execute(plan_backward_.get());
#endif
   
   WALBERLA_FOR_ALL_CELLS_XYZ(data, {
      data->get(x,y,z) = (typename Field_T::value_type)( in_[uint_c(z) + data->zSize() * (uint_c(y) + data->ySize() * uint_c(x))] );
   });
}

// explicit instantiation
template class FourierTransform<GhostLayerField< double, 1 > >;
template class FourierTransform<GhostLayerField< float, 1 > >;
template class FourierTransform<Field< double, 1 > >;
template class FourierTransform<Field< float, 1 > >;

} // namespace fft
} // namespace walberla
