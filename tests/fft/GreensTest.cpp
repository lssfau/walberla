#include "blockforest/Initialization.h"
#include "core/Abort.h"
#include "core/Environment.h"
#include "core/debug/TestSubsystem.h"
#include "core/logging/Logging.h"
#include "core/math/Constants.h"
#include "field/AddToStorage.h"
#include "blockforest/communication/UniformBufferedScheme.h"
#include "stencil/D3Q7.h"
#include "field/communication/PackInfo.h"

#include "fft/Fft.h"
#include <random>

using namespace walberla;
typedef GhostLayerField< real_t, 1 > Field_T;

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
   BlockDataID originalFieldId = field::addToStorage<Field_T>( blocks, "original", real_t(0), field::zyxf, 1 );
   BlockDataID fftFieldId = field::addToStorage<Field_T>( blocks, "result", real_t(0), field::zyxf, 1 );
   
   auto comm = blockforest::communication::UniformBufferedScheme< stencil::D3Q7 >( blocks );
   comm.addPackInfo(make_shared< field::communication::PackInfo< Field_T > >( fftFieldId ));
   
   // set input data
   Vector3<real_t> position(7,5,3);
   for( auto block = blocks->begin(); block != blocks->end(); ++block )
   {
      Field_T *data_in = block->getData< Field_T >( originalFieldId );
      Field_T *data_out = block->getData< Field_T >( fftFieldId );
      WALBERLA_FOR_ALL_CELLS_XYZ(data_in, {
         Vector3<real_t> point( real_c(x), real_c(y), real_c(z) );
         blocks->transformBlockLocalToGlobal(point, *block);
         if (point == position)
            data_in->get(x,y,z) = 1;
      });
      
      data_out->set(data_in);
   }
	
   Vector3<uint_t> dim(blocks->getNumberOfXCells(), blocks->getNumberOfYCells(), blocks->getNumberOfZCells());
   auto greens = [&dim] (uint_t x, uint_t y, uint_t z) -> real_t {
      if (x == 0 && y == 0 && z == 0)
         return 0;
      return real_c(0.5) / ( std::cos( real_c(2) * real_c(math::pi) * real_c(x) / real_c(dim[0])) +
                             std::cos( real_c(2) * real_c(math::pi) * real_c(y) / real_c(dim[1])) +
                             std::cos( real_c(2) * real_c(math::pi) * real_c(z) / real_c(dim[2])) -
                             real_c(3) ) / real_c(dim[0]*dim[1]*dim[2]);
   };
   
   fft::FourierTransform<Field_T> ft(blocks, fftFieldId, greens);
   ft();
   
   comm();
   for( auto block = blocks->begin(); block != blocks->end(); ++block )
   {
      Field_T *data_in = block->getData< Field_T >( originalFieldId );
      Field_T *data_out = block->getData< Field_T >( fftFieldId );
      WALBERLA_FOR_ALL_CELLS_XYZ(data_in, {
         real_t charge = data_out->get(x+1,y  ,z  ) +
                         data_out->get(x-1,y  ,z  ) +
                         data_out->get(x  ,y+1,z  ) +
                         data_out->get(x  ,y-1,z  ) +
                         data_out->get(x  ,y  ,z+1) +
                         data_out->get(x  ,y  ,z-1) -
                         data_out->get(x  ,y  ,z  ) * 6;
         Vector3<real_t> point( real_c(x), real_c(y), real_c(z) );
         WALBERLA_CHECK_LESS(std::fabs(charge-data_in->get(x,y,z)), 1e-3, "Incorrect charge at " << point);
      });
   }
   
   return 0;
}
