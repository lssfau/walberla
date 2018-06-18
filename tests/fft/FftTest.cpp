#include "blockforest/Initialization.h"
#include "core/Abort.h"
#include "core/Environment.h"
#include "core/debug/TestSubsystem.h"
#include "core/logging/Logging.h"
#include "field/AddToStorage.h"

#include "fft/Fft.h"
#include <random>

using namespace walberla;
typedef GhostLayerField< real_t, 1 > Field_T;

const real_t factor = real_t(2.745652);

real_t greens(uint_t /* x */, uint_t /* y */, uint_t /* z */)
{
   return factor;
}

int main (int argc, char** argv)
{
   debug::enterTestMode();
   walberla::Environment env(argc, argv);
   WALBERLA_LOG_INFO_ON_ROOT("Test starting on " << (MPIManager::instance()->numProcesses()) << " processes."  );
   
   uint_t L = 16;
   uint_t processes = uint_c(cbrt(MPIManager::instance()->numProcesses()));
   WALBERLA_ASSERT_EQUAL(processes*processes*processes, MPIManager::instance()->numProcesses(), "Number of processes must be a cubic number.");
   
   uint_t num_blocks = processes;
   uint_t cells_per_block = L/processes;
   WALBERLA_ASSERT_EQUAL(cells_per_block*processes, L, "Number of processes per direction must evenly divide " << L);
   
   auto blocks = blockforest::createUniformBlockGrid(num_blocks,num_blocks,num_blocks, cells_per_block,cells_per_block,cells_per_block, 1.0, processes,processes,processes, true,true,true);
   BlockDataID originalFieldId = field::addToStorage<Field_T>( blocks, "original" );
   BlockDataID fftFieldId = field::addToStorage<Field_T>( blocks, "result" );
   
   // set input data
   for( auto block = blocks->begin(); block != blocks->end(); ++block )
   {
      Field_T *data_in = block->getData< Field_T >( originalFieldId );
      Field_T *data_out = block->getData< Field_T >( fftFieldId );
      WALBERLA_FOR_ALL_CELLS_XYZ_OMP(data_in, omp critical, {
         Vector3<real_t> point( real_c(x), real_c(y), real_c(z) );
         blocks->transformBlockLocalToGlobal(point, *block);
         data_in->get(x,y,z) = real_c(std::ranlux48_base(uint_c(point[0])+(uint_c(point[1])*L+uint_c(point[2]))*L)())*real_c(std::pow(2,-48));
      });
      
      data_out->set(data_in);
   }
   
   fft::FourierTransform<Field_T> ft(blocks, fftFieldId, greens);
   ft();
   
   for( auto block = blocks->begin(); block != blocks->end(); ++block )
   {
      Field_T *data_in = block->getData< Field_T >( originalFieldId );
      Field_T *data_out = block->getData< Field_T >( fftFieldId );
      WALBERLA_FOR_ALL_CELLS_XYZ(data_in, {
         Vector3<real_t> point( real_c(x), real_c(y), real_c(z) );
         WALBERLA_CHECK_FLOAT_EQUAL(data_in->get(x,y,z)*real_c(L*L*L)*factor, data_out->get(x,y,z), "Incorrect value at " << point);
      });
   }
   
   return 0;
}
