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
//! \file DeviceWrapper.h
//! \ingroup gpu
//! \author Richard Angersbach <richard.angersbach@fau.de>
//
//======================================================================================================================

#pragma once

/// \cond internal

#include <sstream>
#include "core/Abort.h"

// CMake generated header
#include "waLBerlaDefinitions.h"

// DEVICE SECTION //

#if defined(WALBERLA_BUILD_WITH_GPU_SUPPORT)

#define WALBERLA_DEVICE_SECTION() if (true)
#define WALBERLA_NON_DEVICE_SECTION() if (false)

#else

#define WALBERLA_DEVICE_SECTION() if (false)
#define WALBERLA_NON_DEVICE_SECTION() if (true)

#endif

namespace walberla {
namespace gpustubs {
   // empty namespace which can be used
} // namespace gpustubs
} // namespace walberla

#if defined(WALBERLA_BUILD_WITH_GPU_SUPPORT)

// include runtime header
#include "gpu/GPUWrapper.h"

#else // WALBERLA_BUILD_WITH_GPU_SUPPORT

namespace walberla {
namespace gpustubs {

// dummy definitions for CUDA/HIP data types and functions in order to guarantee successful compilation without CUDA/HIP enabled

#define WALBERLA_DEVICE_FUNCTION_ERROR \
   WALBERLA_ABORT("Invalid device function call! In case of compiling without CUDA/HIP, functions are not " \
                  "available and shouldn't be called!");

#ifndef __CUDACC__
   #define __device__
   #define __global__
   #define __host__
   #define __forceinline__
#endif

using gpuError_t = int;
const gpuError_t gpuSuccess = 0;

#define gpuHostAllocDefault 0x00
#define gpuHostAllocMapped 0x02
#define gpuHostAllocPortable 0x01
#define gpuHostAllocWriteCombined 0x04

using gpuMemcpyKind                          = int;
const gpuMemcpyKind gpuMemcpyHostToHost     = 0;
const gpuMemcpyKind gpuMemcpyHostToDevice   = 1;
const gpuMemcpyKind gpuMemcpyDeviceToHost   = 2;
const gpuMemcpyKind gpuMemcpyDeviceToDevice = 3;
const gpuMemcpyKind gpuMemcpyDefault        = 4;

inline const char* gpuGetErrorName(gpuError_t /*code*/) { WALBERLA_DEVICE_FUNCTION_ERROR }
inline const char* gpuGetErrorString(gpuError_t /*code*/) { WALBERLA_DEVICE_FUNCTION_ERROR }
inline gpuError_t gpuGetLastError(void) { WALBERLA_DEVICE_FUNCTION_ERROR }
inline gpuError_t gpuPeekAtLastError(void) { WALBERLA_DEVICE_FUNCTION_ERROR }

inline gpuError_t gpuMalloc(void** /*devPtr*/, size_t /*size*/) { WALBERLA_DEVICE_FUNCTION_ERROR }
inline gpuError_t gpuMallocHost(void** /*ptr*/, size_t /*size*/) { WALBERLA_DEVICE_FUNCTION_ERROR }
inline gpuError_t gpuHostAlloc(void** /*pHost*/, size_t /*size*/, unsigned int /*flags*/) { WALBERLA_DEVICE_FUNCTION_ERROR }

struct gpuPos
{
   size_t x, y, z;
};

struct gpuPitchedPtr
{
   size_t pitch;
   void* ptr;
   size_t xsize;
   size_t ysize;
};

struct gpuExtent
{
   size_t depth;
   size_t height;
   size_t width;
};

struct gpuArray;
typedef struct gpuArray* gpuArray_t;
typedef struct gpuArray* gpuArray_const_t;

struct CUstream_st;
typedef struct CUstream_st* gpuStream_t;
inline gpuError_t gpuStreamDestroy(gpuStream_t /*stream*/) { WALBERLA_DEVICE_FUNCTION_ERROR }
inline gpuError_t gpuStreamCreateWithPriority(gpuStream_t* /*pStream*/, unsigned int /*flags*/, int /*priority*/) { WALBERLA_DEVICE_FUNCTION_ERROR }
inline gpuError_t gpuStreamCreateWithFlags(gpuStream_t* /*pStream*/, unsigned int /*flags*/) { WALBERLA_DEVICE_FUNCTION_ERROR }
inline gpuError_t gpuDeviceGetStreamPriorityRange(int* /*leastPriority*/, int* /*greatestPriority*/) { WALBERLA_DEVICE_FUNCTION_ERROR }
inline gpuError_t gpuStreamCreate(gpuStream_t* /*pStream*/) { WALBERLA_DEVICE_FUNCTION_ERROR }
inline gpuError_t gpuStreamSynchronize(gpuStream_t /*stream*/) { WALBERLA_DEVICE_FUNCTION_ERROR }

struct gpuMemcpy3DParms
{
   gpuArray_t dstArray;
   gpuPos dstPos;
   gpuPitchedPtr dstPtr;
   gpuExtent extent;
   gpuMemcpyKind kind;
   gpuArray_t srcArray;
   gpuPos srcPos;
   gpuPitchedPtr srcPtr;
};

inline gpuError_t gpuMemcpy(void* /*dst*/, const void* /*src*/, size_t /*count*/, gpuMemcpyKind /*kind*/) { WALBERLA_DEVICE_FUNCTION_ERROR }
inline gpuError_t gpuMemcpyAsync(void* /*dst*/, const void* /*src*/, size_t /*count*/, gpuMemcpyKind /*kind*/, gpuStream_t /*stream*/) { WALBERLA_DEVICE_FUNCTION_ERROR }
inline gpuError_t gpuMemcpy3D(const gpuMemcpy3DParms* /*p*/) { WALBERLA_DEVICE_FUNCTION_ERROR }
inline gpuError_t gpuMemcpy3DAsync(const gpuMemcpy3DParms* /*p*/, gpuStream_t /*stream*/) { WALBERLA_DEVICE_FUNCTION_ERROR }

inline gpuPos make_gpuPos(size_t /*x*/, size_t /*y*/, size_t /*z*/) { WALBERLA_DEVICE_FUNCTION_ERROR }
inline gpuPitchedPtr make_gpuPitchedPtr (void* /*d*/, size_t /*p*/, size_t /*xsz*/, size_t /*ysz*/) { WALBERLA_DEVICE_FUNCTION_ERROR }
inline gpuExtent make_gpuExtent(size_t /*w*/, size_t /*h*/, size_t /*d*/) { WALBERLA_DEVICE_FUNCTION_ERROR }

inline gpuError_t gpuFree(void* /*devPtr*/) { WALBERLA_DEVICE_FUNCTION_ERROR }
inline gpuError_t gpuFreeHost(void* /*ptr*/) { WALBERLA_DEVICE_FUNCTION_ERROR }

inline gpuError_t gpuDeviceSynchronize(void) { WALBERLA_DEVICE_FUNCTION_ERROR }

struct CUevent_st;
typedef struct CUevent_st* gpuEvent_t;
inline gpuError_t gpuEventCreate(gpuEvent_t* /*event*/) { WALBERLA_DEVICE_FUNCTION_ERROR }
inline gpuError_t gpuEventCreateWithFlags(gpuEvent_t* /*event*/, unsigned int /*flags*/) { WALBERLA_DEVICE_FUNCTION_ERROR }
inline gpuError_t gpuEventRecord(gpuEvent_t /*event*/, gpuStream_t /*stream*/) { WALBERLA_DEVICE_FUNCTION_ERROR }
inline gpuError_t gpuEventDestroy(gpuEvent_t /*event*/) { WALBERLA_DEVICE_FUNCTION_ERROR }
inline gpuError_t gpuStreamWaitEvent (gpuStream_t /*stream*/, gpuEvent_t /*event*/, unsigned int /*flags*/) { WALBERLA_DEVICE_FUNCTION_ERROR }

#define gpuStreamDefault 0x00

inline gpuError_t gpuGetDeviceCount(int* /*count*/) { WALBERLA_DEVICE_FUNCTION_ERROR }
inline gpuError_t gpuSetDevice(int /*device*/) { WALBERLA_DEVICE_FUNCTION_ERROR }

struct CUuuid_st
{
   char bytes;
};
typedef struct CUuuid_st gpuUUID_t;

struct gpuDeviceProp
{
   char name[256];
   gpuUUID_t uuid;
   size_t totalGlobalMem;
   size_t sharedMemPerBlock;
   int regsPerBlock;
   int warpSize;
   size_t memPitch;
   int maxThreadsPerBlock;
   int maxThreadsDim[3];
   int maxGridSize[3];
   int clockRate;
   size_t totalConstMem;
   int major;
   int minor;
   size_t textureAlignment;
   size_t texturePitchAlignment;
   int deviceOverlap;
   int multiProcessorCount;
   int kernelExecTimeoutEnabled;
   int integrated;
   int canMapHostMemory;
   int computeMode;
   int maxTexture1D;
   int maxTexture1DMipmap;
   int maxTexture1DLinear;
   int maxTexture2D[2];
   int maxTexture2DMipmap[2];
   int maxTexture2DLinear[3];
   int maxTexture2DGather[2];
   int maxTexture3D[3];
   int maxTexture3DAlt[3];
   int maxTextureCubemap;
   int maxTexture1DLayered[2];
   int maxTexture2DLayered[3];
   int maxTextureCubemapLayered[2];
   int maxSurface1D;
   int maxSurface2D[2];
   int maxSurface3D[3];
   int maxSurface1DLayered[2];
   int maxSurface2DLayered[3];
   int maxSurfaceCubemap;
   int maxSurfaceCubemapLayered[2];
   size_t surfaceAlignment;
   int concurrentKernels;
   int ECCEnabled;
   int pciBusID;
   int pciDeviceID;
   int pciDomainID;
   int tccDriver;
   int asyncEngineCount;
   int unifiedAddressing;
   int memoryClockRate;
   int memoryBusWidth;
   int l2CacheSize;
   int persistingL2CacheMaxSize;
   int maxThreadsPerMultiProcessor;
   int streamPrioritiesSupported;
   int globalL1CacheSupported;
   int localL1CacheSupported;
   size_t sharedMemPerMultiprocessor;
   int regsPerMultiprocessor;
   int managedMemory;
   int isMultiGpuBoard;
   int multiGpuBoardGroupID;
   int singleToDoublePrecisionPerfRatio;
   int pageableMemoryAccess;
   int concurrentManagedAccess;
   int computePreemptionSupported;
   int canUseHostPointerForRegisteredMem;
   int cooperativeLaunch;
   int cooperativeMultiDeviceLaunch;
   int pageableMemoryAccessUsesHostPageTables;
   int directManagedMemAccessFromHost;
   int accessPolicyMaxWindowSize;
};
inline gpuError_t gpuGetDeviceProperties(gpuDeviceProp* /*prop*/, int /*device*/) { WALBERLA_DEVICE_FUNCTION_ERROR }

struct uint3
{
   unsigned int x, y, z;
};
typedef struct uint3 uint3;

struct dim3
{
   unsigned int x, y, z;
   dim3(unsigned int vx = 1, unsigned int vy = 1, unsigned int vz = 1) : x(vx), y(vy), z(vz) {}
   dim3(uint3 v) : x(v.x), y(v.y), z(v.z) {}
   operator uint3(void) { uint3 t; t.x = x; t.y = y; t.z = z; return t; }
};
typedef struct dim3 dim3;

inline gpuError_t gpuLaunchKernel(const void* /*func*/, dim3 /*gridDim*/, dim3 /*blockDim*/, void** /*args*/, size_t /*sharedMem*/, gpuStream_t /*stream*/) { WALBERLA_DEVICE_FUNCTION_ERROR }

#ifdef _WIN32
#define GPURT_CB __stdcall
#else
#define GPURT_CB
#endif

typedef void(GPURT_CB* gpuHostFn_t)(void* /*userData*/);
inline gpuError_t gpuLaunchHostFunc(gpuStream_t /*stream*/, gpuHostFn_t /*fn*/, void* /*userData*/) { WALBERLA_DEVICE_FUNCTION_ERROR }

#undef WALBERLA_DEVICE_FUNCTION_ERROR

} // namespace gpustubs
using namespace gpustubs;

} // namespace walberla


#endif // WALBERLA_BUILD_WITH_GPU_SUPPORT

/// \endcond
