#pragma once

#include "core/all.h"
#include "gpu/GPURAII.h"


__global__ void kernel(double * a, double * b) {
    uint32_t idx = threadIdx.x;

    a[idx] = 2.0 * b[idx];
}

using namespace walberla::gpu;

constexpr uint32_t N { 64 };

void run(double * a, double * b, StreamRAII stream) {
    kernel<<< 1, N, 0, stream >>>(a, b);
    stream.synchronize();
}

int main(void) {
    double * a;
    WALBERLA_GPU_CHECK(gpuMallocManaged(&a, N * sizeof(double)));

    for(int i = 0; i < N; ++i) {
        a[i] = 2.0 * double(i) + 0.75;
    }

    double * b;
    WALBERLA_GPU_CHECK(gpuMallocManaged(&b, N * sizeof(double)));;

    //  Default stream
    auto defaultStream = StreamRAII::defaultStream();

    kernel<<< 1, N, 0, defaultStream >>>(a, b);
    defaultStream.synchronize();

    for(int i = 0; i < N; ++i) {
        WALBERLA_CHECK_EQUAL(a[i], 2.0 * b[i]);
    }

    [[maybe_unused]] gpuStream_t underlyingStream;

    {
        std::swap(a, b);
        
        auto s1 = StreamRAII::newStream();
        kernel<<< 1, N, 0, s1 >>>(a, b);
        s1.synchronize();

        for(int i = 0; i < N; ++i) {
            WALBERLA_CHECK_EQUAL(a[i], 2.0 * b[i]);
        }

        underlyingStream = s1;
    }
    
#if defined(WALBERLA_BUILD_WITH_HIP)
    //  Underlying stream must be destroyed at this point, so any operations on it should fail
    //  (HIP raises a proper error here, CUDA just segfaults.)
    WALBERLA_CHECK(gpuStreamSynchronize(underlyingStream) != gpuSuccess);
#endif

    {
        std::swap(a, b);
        
        auto s2 = StreamRAII::newStream();
        underlyingStream = s2;

        //  Move constructor test
        run(a, b, std::move(s2));

        for(int i = 0; i < N; ++i) {
            WALBERLA_CHECK_EQUAL(a[i], 2.0 * b[i]);
        }

#if defined(WALBERLA_BUILD_WITH_HIP)
        WALBERLA_CHECK(gpuStreamSynchronize(underlyingStream) != gpuSuccess);
#endif
    }

    {
        //  Move assignment test
        auto raiiA = StreamRAII::newStream();
        gpuStream_t streamA { raiiA };
        
        auto raiiB = StreamRAII::newStream();
        gpuStream_t streamB { raiiB };

        raiiA = std::move(raiiB);
        raiiA.synchronize();
        WALBERLA_GPU_CHECK(gpuStreamSynchronize(streamB));

#if defined(WALBERLA_BUILD_WITH_HIP)
        WALBERLA_CHECK(gpuStreamSynchronize(streamA) != gpuSuccess);
#endif
    }

    return EXIT_SUCCESS;
}
