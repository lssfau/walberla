#pragma once

#include "blockforest/StructuredBlockForest.h"

#include "waLBerlaDefinitions.h"

#ifdef WALBERLA_BUILD_WITH_MPI
#define WALBERLA_USE_PFFT
#endif

#ifdef WALBERLA_USE_PFFT
#include "pfft.h"
#else
#include "fftw3.h"
#endif

namespace walberla {
namespace fft {

template <typename Field_T>
class FourierTransform
{
   public:
      FourierTransform( shared_ptr< StructuredBlockForest > & blocks, BlockDataID fieldId,
                        const std::function<real_t(uint_t,uint_t,uint_t)>& greens = std::function<real_t(uint_t,uint_t,uint_t)>() );
      void operator() ();
      
   private:
      shared_ptr< StructuredBlockForest > & blocks_;
      BlockDataID fieldId_;
   
      ptrdiff_t local_no[3];
#ifdef WALBERLA_USE_PFFT
      typedef std::unique_ptr<double[], std::function<void(double*)> > FFTReal;
      typedef std::unique_ptr<pfft_complex[], std::function<void(pfft_complex*)> > FFTComplex;
      typedef std::unique_ptr<std::remove_pointer<pfft_plan>::type, std::function<void(pfft_plan)> > FFTPlan;
#else
      using FFTReal = std::unique_ptr<double [], std::function<void (double *)>>;
      typedef std::unique_ptr<fftw_complex[], std::function<void(fftw_complex*)> > FFTComplex;
      typedef std::unique_ptr<std::remove_pointer<fftw_plan>::type, std::function<void(fftw_plan)> > FFTPlan;
#endif
      FFTReal in_, greens_;
      FFTComplex out_;
      FFTPlan plan_forward_, plan_backward_;
};

} // namespace fft
} // namespace walberla
