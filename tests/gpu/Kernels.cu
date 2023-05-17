#include "gpu/FieldAccessor.h"

namespace walberla {


namespace gpu
{
   template<typename T>
   class GPUField;
}

__global__ void kernel_double(gpu::FieldAccessor<double> f )
{
   f.set( blockIdx, threadIdx );
   f.get() *= 2.0;
}


} // namespace walberla
