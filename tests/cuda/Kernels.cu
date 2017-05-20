
#include <iostream>

#include "cuda/FieldAccessor.h"
#include "cuda/FieldIndexing.h"

namespace walberla {


namespace cuda {
   template<typename T>
   class GPUField;
}

__global__ void kernel_double( cuda::FieldAccessor<double> f )
{
   f.set( blockIdx, threadIdx );
   f.get() *= 2.0;
}

void kernel_double_field( const cuda::GPUField<double> & field )
{
   using namespace std;
   cuda::FieldIndexing<double> iter = cuda::FieldIndexing<double>::sliceBeforeGhostLayerXYZ( field, 1, stencil::E, true );
   std::cout << "Kernel call dims "
             << iter.blockDim().x << ","
             << iter.gridDim().x << ","
             << iter.gridDim().y << ","
             << iter.gridDim().z << endl;
   kernel_double<<< iter.gridDim(), iter.blockDim() >>> ( iter.gpuAccess() );
}

} // namespace walberla
