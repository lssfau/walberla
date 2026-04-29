#include "TestAllocators.h"

namespace test_allocators {

template<typename T>
__global__ void checkVectorKernel(T* ptr, size_t size, bool* flag){
    if (blockIdx.x != 0 || threadIdx.x != 0)
        return;

    *flag = true;
    for(size_t i = 0; i < size; ++i){
        if(ptr[i] != static_cast<T>(i)) 
            *flag = false;
    }

}

template<typename T>
void checkVector(T* ptr, size_t size, bool* flag){
    checkVectorKernel<<<1,1>>>(ptr, size, flag);
    WALBERLA_GPU_CHECK(gpuDeviceSynchronize());
}

template void checkVector<int>(int*, size_t, bool*);
template void checkVector<float>(float*, size_t, bool*);
template void checkVector<double>(double*, size_t, bool*);

} // namespace test_allocators
