#include "cuda/FieldAccessor.h"

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


} // namespace walberla
