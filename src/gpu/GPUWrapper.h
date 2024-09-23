//======================================================================================================================
//
//  This file is part of waLBerla. waLBerla is free software: you can
//  redistribute it and/or modify it under the terms of the GNU General Public
//  License as published by the Free Software Foundation, either version 3 of
//  the License, or (at your option) any later version.
//
//  waLBerla is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
//  for more details.
//
//  You should have received a copy of the GNU General Public License along
//  with waLBerla (see COPYING.txt). If not, see <http://www.gnu.org/licenses/>.
//
//! \file GPUWrapper.h
//! \ingroup gpu
//! \author Markus Holzer <markus.holzer@fau.de>
//
//======================================================================================================================

#pragma once

// https://rocmdocs.amd.com/en/latest/Programming_Guides/CUDAAPIHIPTEXTURE.html
#if defined(WALBERLA_BUILD_WITH_CUDA)
    #include <cuda_runtime.h>


    using gpuError_t = cudaError_t;
    #define gpuSuccess cudaSuccess
    #define gpuGetErrorName cudaGetErrorName
    #define gpuGetErrorString cudaGetErrorString
    #define gpuPeekAtLastError cudaPeekAtLastError
    #define gpuGetLastError cudaGetLastError

    #define gpuMalloc cudaMalloc
    #define gpuMallocHost cudaMallocHost
    #define gpuHostAllocDefault cudaHostAllocDefault
    #define gpuHostAlloc cudaHostAlloc
    #define gpuMemcpyHostToDevice cudaMemcpyHostToDevice
    #define gpuMemcpyDeviceToHost cudaMemcpyDeviceToHost
    #define gpuMemcpyDeviceToDevice cudaMemcpyDeviceToDevice
    #define gpuMemcpy cudaMemcpy
    #define gpuMemcpyAsync cudaMemcpyAsync
    #define gpuMemcpy3D cudaMemcpy3D
    #define gpuMemcpy3DParms cudaMemcpy3DParms
    #define gpuMemcpy3DAsync cudaMemcpy3DAsync

    #define gpuMemset cudaMemset
    #define gpuMemsetAsync cudaMemsetAsync
    #define gpuMemset2D cudaMemset2D
    #define gpuMemset2DAsync cudaMemset2DAsync
    #define gpuMemset3D cudaMemset3D
    #define gpuMemset3DAsync cudaMemset3DAsync

    #define make_gpuPos make_cudaPos
    #define make_gpuPitchedPtr make_cudaPitchedPtr
    #define gpuPitchedPtr cudaPitchedPtr
    #define make_gpuExtent make_cudaExtent
    using gpuExtent = cudaExtent;

    #define gpuFree cudaFree
    #define gpuFreeHost cudaFreeHost

    using gpuStream_t = cudaStream_t;
    #define gpuStreamDestroy cudaStreamDestroy
    #define gpuStreamCreateWithPriority cudaStreamCreateWithPriority
    #define gpuDeviceGetStreamPriorityRange cudaDeviceGetStreamPriorityRange
    #define gpuStreamCreate cudaStreamCreate
    #define gpuStreamSynchronize cudaStreamSynchronize
    #define gpuDeviceSynchronize cudaDeviceSynchronize

    using gpuEvent_t = cudaEvent_t;
    #define gpuEventCreate cudaEventCreate
    #define gpuEventRecord cudaEventRecord
    #define gpuEventDestroy cudaEventDestroy
    #define gpuStreamWaitEvent cudaStreamWaitEvent
    #define gpuStreamDefault cudaStreamDefault

    #define gpuGetDeviceCount cudaGetDeviceCount
    #define gpuSetDevice cudaSetDevice
    #define gpuDeviceProp cudaDeviceProp
    #define gpuGetDeviceProperties cudaGetDeviceProperties

    #define gpuLaunchKernel cudaLaunchKernel
#endif


#ifdef WALBERLA_BUILD_WITH_HIP
    #include <hip/hip_runtime.h>


    using gpuError_t = hipError_t;
    #define gpuSuccess hipSuccess
    #define gpuGetErrorName hipGetErrorName
    #define gpuGetErrorString hipGetErrorString
    #define gpuPeekAtLastError hipPeekAtLastError
    #define gpuGetLastError hipGetLastError

    #define gpuMalloc hipMalloc
    #define gpuMallocHost hipHostMalloc
    #define gpuHostAllocDefault hipHostMallocDefault
    // warning: 'hipHostAlloc' is deprecated: use hipHostMalloc insteadwarning: 'hipHostAlloc' is deprecated: use hipHostMalloc instead
    #define gpuHostAlloc hipHostMalloc
    #define gpuMemcpyHostToDevice hipMemcpyHostToDevice
    #define gpuMemcpyDeviceToHost hipMemcpyDeviceToHost
    #define gpuMemcpyDeviceToDevice hipMemcpyDeviceToDevice
    #define gpuMemcpy hipMemcpy
    #define gpuMemcpyAsync hipMemcpyAsync
    #define gpuMemcpy3D hipMemcpy3D
    #define gpuMemcpy3DParms hipMemcpy3DParms
    #define gpuMemcpy3DAsync hipMemcpy3DAsync

    #define gpuMemset hipMemset
    #define gpuMemsetAsync hipMemsetAsync
    #define gpuMemset2D hipMemset2D
    #define gpuMemset2DAsync hipMemset2DAsync
    #define gpuMemset3D hipMemset3D
    #define gpuMemset3DAsync hipMemset3DAsync

    #define make_gpuPitchedPtr make_hipPitchedPtr
    #define make_gpuPos make_hipPos
    using gpuPitchedPtr = hipPitchedPtr;
    #define make_gpuExtent make_hipExtent
    using gpuExtent = hipExtent;

    #define gpuFree hipFree
    #define gpuFreeHost hipHostFree

    using gpuStream_t = hipStream_t;
    #define gpuStreamDestroy hipStreamDestroy
    #define gpuStreamCreateWithPriority hipStreamCreateWithPriority
    #define gpuDeviceGetStreamPriorityRange hipDeviceGetStreamPriorityRange
    #define gpuStreamCreate hipStreamCreate
    #define gpuStreamSynchronize hipStreamSynchronize
    #define gpuDeviceSynchronize hipDeviceSynchronize

    using gpuEvent_t = hipEvent_t;
    #define gpuEventCreate hipEventCreate
    #define gpuEventRecord hipEventRecord
    #define gpuEventDestroy hipEventDestroy
    #define gpuStreamWaitEvent hipStreamWaitEvent
    #define gpuStreamDefault hipStreamDefault

    #define gpuGetDeviceCount hipGetDeviceCount
    #define gpuSetDevice hipSetDevice
    #define gpuDeviceProp hipDeviceProp
    #define gpuGetDeviceProperties hipGetDeviceProperties

    #define gpuLaunchKernel hipLaunchKernel
#endif
