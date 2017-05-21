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
//! \file GPUCopy.cpp
//! \ingroup cuda
//! \author Paulo Carvalho <prcjunior@inf.ufpr.br>
//! \brief Copy routines of 4D intervals involving GPU buffers.
//
//======================================================================================================================

#include "core/debug/Debug.h"

#include "GPUCopy.h"
#include "ErrorChecking.h"


namespace walberla {
namespace cuda {


void copyDevToDevFZYXRestricted( const cudaPitchedPtr& dst, const cudaPitchedPtr& src,
                                 uint_t typeSz, uint_t dstAllocZ, uint_t srcAllocZ,
                                 uint_t dstX, uint_t dstY, uint_t dstZ, uint_t dstF,
                                 uint_t srcX, uint_t srcY, uint_t srcZ, uint_t srcF,
                                 uint_t Nx, uint_t Ny, uint_t Nz, uint_t Nf )
{
   WALBERLA_ASSERT( Nf == 1 || ( Nz == dstAllocZ && Nz == srcAllocZ ) );

   cudaMemcpy3DParms p;
   memset( &p, 0, sizeof(p) );

   p.srcPos.x = srcX * typeSz;
   p.srcPos.y = srcY;
   p.srcPos.z = srcF * srcAllocZ + srcZ;
   p.srcPtr.ptr = src.ptr;
   p.srcPtr.pitch = src.pitch;
   p.srcPtr.xsize = src.xsize;
   p.srcPtr.ysize = src.ysize;

   p.dstPos.x = dstX * typeSz;
   p.dstPos.y = dstY;
   p.dstPos.z = dstF * dstAllocZ + dstZ;
   p.dstPtr.ptr = dst.ptr;
   p.dstPtr.pitch = dst.pitch;
   p.dstPtr.xsize = dst.xsize;
   p.dstPtr.ysize = dst.ysize;

   p.extent.width = Nx * typeSz;
   p.extent.height = Ny;
   p.extent.depth = Nz * Nf;

   p.kind = cudaMemcpyDeviceToDevice;

   WALBERLA_CUDA_CHECK( cudaMemcpy3D(&p) );
}


void copyHostToDevFZYXRestricted( const cudaPitchedPtr& dst, unsigned char* src,
                                  uint_t typeSz, uint_t dstAllocZ, uint_t srcAllocZ,
                                  uint_t dstX, uint_t dstY, uint_t dstZ, uint_t dstF,
                                  uint_t srcX, uint_t srcY, uint_t srcZ, uint_t srcF,
                                  uint_t Nx, uint_t Ny, uint_t Nz, uint_t Nf )
{
   WALBERLA_ASSERT( Nf == 1 || ( Nz == dstAllocZ && Nz == srcAllocZ ) );

   cudaMemcpy3DParms p;
   memset( &p, 0, sizeof(p) );

   p.srcPos.x = srcX * typeSz;
   p.srcPos.y = srcY;
   p.srcPos.z = srcF * srcAllocZ + srcZ;
   p.srcPtr.ptr = src;
   p.srcPtr.pitch = Nx * typeSz;
   p.srcPtr.xsize = Nx * typeSz;
   p.srcPtr.ysize = Ny;

   p.dstPos.x = dstX * typeSz;
   p.dstPos.y = dstY;
   p.dstPos.z = dstF * dstAllocZ + dstZ;
   p.dstPtr.ptr = dst.ptr;
   p.dstPtr.pitch = dst.pitch;
   p.dstPtr.xsize = dst.xsize;
   p.dstPtr.ysize = dst.ysize;

   p.extent.width = Nx * typeSz;
   p.extent.height = Ny;
   p.extent.depth = Nz * Nf;

   p.kind = cudaMemcpyHostToDevice;

   WALBERLA_CUDA_CHECK( cudaMemcpy3D(&p) );
}


void copyDevToHostFZYXRestricted( unsigned char* dst, const cudaPitchedPtr& src,
                                  uint_t typeSz, uint_t dstAllocZ, uint_t srcAllocZ,
                                  uint_t dstX, uint_t dstY, uint_t dstZ, uint_t dstF,
                                  uint_t srcX, uint_t srcY, uint_t srcZ, uint_t srcF,
                                  uint_t Nx, uint_t Ny, uint_t Nz, uint_t Nf )
{
   cudaMemcpy3DParms p;
   memset( &p, 0, sizeof(p) );

   p.srcPos.x = srcX * typeSz;
   p.srcPos.y = srcY;
   p.srcPos.z = srcF * srcAllocZ + srcZ;
   p.srcPtr.ptr = src.ptr;
   p.srcPtr.pitch = src.pitch;
   p.srcPtr.xsize = src.xsize;
   p.srcPtr.ysize = src.ysize;

   p.dstPos.x = dstX * typeSz;
   p.dstPos.y = dstY;
   p.dstPos.z = dstF * dstAllocZ + dstZ;
   p.dstPtr.ptr = dst;
   p.dstPtr.pitch = Nx * typeSz;
   p.dstPtr.xsize = Nx * typeSz;
   p.dstPtr.ysize = Ny;

   p.extent.width = Nx * typeSz;
   p.extent.height = Ny;
   p.extent.depth = Nz * Nf;

   p.kind = cudaMemcpyDeviceToHost;

   WALBERLA_CUDA_CHECK( cudaMemcpy3D(&p) );
}



} // namespace cuda
} // namespace walberla
