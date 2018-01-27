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
//! \author João Victor Tozatti Risso <jvtrisso@inf.ufpr.br>
//! \brief Copy routines of 4D intervals involving GPU buffers.
//
//======================================================================================================================

#include "core/debug/Debug.h"

#include "GPUCopy.h"
#include "ErrorChecking.h"

#include <cstring>


namespace walberla {
namespace cuda {

void copyDevToDevFZYX( const cudaPitchedPtr& dst, const cudaPitchedPtr& src,
                       std::tuple< uint_t, uint_t, uint_t, uint_t > & dstOffset,
                       std::tuple< uint_t, uint_t, uint_t, uint_t > & srcOffset,
                       uint_t dstAllocSizeZ, uint_t srcAllocSizeZ, uint_t typeSize,
                       std::tuple< uint_t, uint_t, uint_t, uint_t > & intervalSize,
                       cudaStream_t copyStream )
{
   const uint_t & Nx = std::get<0>(intervalSize),
                & Ny = std::get<1>(intervalSize),
                & Nz = std::get<2>(intervalSize),
                & Nf = std::get<3>(intervalSize);

   const uint_t & srcX = std::get<0>(srcOffset),
                & srcY = std::get<1>(srcOffset),
                & srcZ = std::get<2>(srcOffset),
                & srcF = std::get<3>(srcOffset);

   const uint_t & dstX = std::get<0>(dstOffset),
                & dstY = std::get<1>(dstOffset),
                & dstZ = std::get<2>(dstOffset),
                & dstF = std::get<3>(dstOffset);

   auto copyFunctor = [&](uint_t dstCoordF, uint_t srcCoordF, uint_t fIntervalSize) {
      WALBERLA_ASSERT( fIntervalSize == 1 || ( Nz == dstAllocSizeZ && Nz == srcAllocSizeZ ) );

      cudaMemcpy3DParms p;
      std::memset( &p, 0, sizeof(p) );

      p.srcPos = make_cudaPos( srcX * typeSize, srcY, srcCoordF * srcAllocSizeZ + srcZ );
      p.srcPtr = make_cudaPitchedPtr( src.ptr, src.pitch, src.xsize, src.ysize );

      p.dstPos = make_cudaPos( dstX * typeSize, dstY, dstCoordF * dstAllocSizeZ + dstZ );
      p.dstPtr = make_cudaPitchedPtr( dst.ptr, dst.pitch, dst.xsize, dst.ysize );

      p.extent = make_cudaExtent( Nx * typeSize, Ny, Nz * fIntervalSize );
      p.kind = cudaMemcpyDeviceToDevice;

      if ( copyStream == 0 )
      {
         WALBERLA_CUDA_CHECK( cudaMemcpy3D(&p) );
      }
      else
      {
         WALBERLA_CUDA_CHECK( cudaMemcpy3DAsync(&p, copyStream) );
      }
   };

   if( Nf == 1 || ( Nz == dstAllocSizeZ && Nz == srcAllocSizeZ ) )
   {
      copyFunctor( dstF, srcF, Nf );
   }
   else
   {
      for( uint_t f = 0; f < Nf; ++f )
      {
         copyFunctor( dstF + f, srcF + f, uint_c(1) );
      }
   }
}


void copyDevToDevZYXF( const cudaPitchedPtr& dst, const cudaPitchedPtr& src,
                       std::tuple< uint_t, uint_t, uint_t, uint_t > & dstOffset,
                       std::tuple< uint_t, uint_t, uint_t, uint_t > & srcOffset,
                       uint_t dstAllocSizeY, uint_t srcAllocSizeY, uint_t typeSize,
                       std::tuple< uint_t, uint_t, uint_t, uint_t > & intervalSize,
                       cudaStream_t copyStream )
{
   const uint_t & Nx = std::get<0>(intervalSize),
                & Ny = std::get<1>(intervalSize),
                & Nz = std::get<2>(intervalSize),
                & Nf = std::get<3>(intervalSize);

   const uint_t & srcX = std::get<0>(srcOffset),
                & srcY = std::get<1>(srcOffset),
                & srcZ = std::get<2>(srcOffset),
                & srcF = std::get<3>(srcOffset);

   const uint_t & dstX = std::get<0>(dstOffset),
                & dstY = std::get<1>(dstOffset),
                & dstZ = std::get<2>(dstOffset),
                & dstF = std::get<3>(dstOffset);

   auto copyFunctor = [&](uint_t dstCoordZ, uint_t srcCoordZ, uint_t zIntervalSize) {
      cudaMemcpy3DParms p;
      std::memset( &p, 0, sizeof(p) );

      p.srcPos = make_cudaPos( srcF * typeSize, srcX, srcCoordZ * srcAllocSizeY + srcY );
      p.srcPtr = make_cudaPitchedPtr( src.ptr, src.pitch, src.xsize, src.ysize );

      p.dstPos = make_cudaPos( dstF * typeSize, dstX, dstCoordZ * dstAllocSizeY + dstY );
      p.dstPtr = make_cudaPitchedPtr( dst.ptr, dst.pitch, dst.xsize, dst.ysize );

      p.extent = make_cudaExtent( Nf * typeSize, Nx, Ny * zIntervalSize );
      p.kind = cudaMemcpyDeviceToDevice;

      if ( copyStream == 0 )
      {
         WALBERLA_CUDA_CHECK( cudaMemcpy3D(&p) );
      }
      else
      {
         WALBERLA_CUDA_CHECK( cudaMemcpy3DAsync(&p, copyStream) );
      }
   };

   if ( Nz == 1 || ( Ny == dstAllocSizeY && Ny == srcAllocSizeY ) )
   {
      copyFunctor( dstZ, srcZ, Nz );
   }
   else
   {
      for( uint_t z = 0; z < Nz; ++z )
      {
         copyFunctor( dstZ + z, srcZ + z, 1 );
      }
   }
}


void copyHostToDevFZYX( const cudaPitchedPtr& dst, unsigned char* src,
                        std::tuple< uint_t, uint_t, uint_t, uint_t > & dstOffset,
                        std::tuple< uint_t, uint_t, uint_t, uint_t > & srcOffset,
                        uint_t dstAllocSizeZ, uint_t srcAllocSizeZ, uint_t typeSize,
                        std::tuple< uint_t, uint_t, uint_t, uint_t > & intervalSize,
                        cudaStream_t copyStream )
{
   const uint_t & Nx = std::get<0>(intervalSize),
                & Ny = std::get<1>(intervalSize),
                & Nz = std::get<2>(intervalSize),
                & Nf = std::get<3>(intervalSize);

   const uint_t & srcX = std::get<0>(srcOffset),
                & srcY = std::get<1>(srcOffset),
                & srcZ = std::get<2>(srcOffset),
                & srcF = std::get<3>(srcOffset);

   const uint_t & dstX = std::get<0>(dstOffset),
                & dstY = std::get<1>(dstOffset),
                & dstZ = std::get<2>(dstOffset),
                & dstF = std::get<3>(dstOffset);

   auto copyFunctor = [&](uint_t dstCoordF, uint_t srcCoordF, uint_t fIntervalSize) {
      cudaMemcpy3DParms p;
      std::memset( &p, 0, sizeof(p) );

      p.srcPos = make_cudaPos( srcX * typeSize, srcY, srcCoordF * srcAllocSizeZ + srcZ );
      p.srcPtr = make_cudaPitchedPtr( src, Nx * typeSize, Nx * typeSize, Ny );

      p.dstPos = make_cudaPos( dstX * typeSize, dstY, dstCoordF * dstAllocSizeZ + dstZ );
      p.dstPtr = make_cudaPitchedPtr( dst.ptr, dst.pitch, dst.xsize, dst.ysize );

      p.extent = make_cudaExtent( Nx * typeSize, Ny, Nz * fIntervalSize );
      p.kind = cudaMemcpyHostToDevice;

      if (copyStream == 0)
      {
         WALBERLA_CUDA_CHECK( cudaMemcpy3D(&p) );
      }
      else
      {
         // Using cudaMemcpy3DAsync requires page-locked memory on the host!
         WALBERLA_CUDA_CHECK( cudaMemcpy3DAsync(&p, copyStream) );
      }
   };

   if ( Nf == 1 || ( Nz == dstAllocSizeZ ) )
   {
      copyFunctor( dstF, srcF, Nf );
   }
   else
   {
      for( uint_t f = 0; f < Nf; ++f )
      {
         copyFunctor( dstF + f, srcF + f, uint_c(1) );
      }
   }
}

void copyHostToDevZYXF( const cudaPitchedPtr& dst, unsigned char* src,
                        std::tuple< uint_t, uint_t, uint_t, uint_t > & dstOffset,
                        std::tuple< uint_t, uint_t, uint_t, uint_t > & srcOffset,
                        uint_t dstAllocSizeY, uint_t srcAllocSizeY, uint_t typeSize,
                        std::tuple< uint_t, uint_t, uint_t, uint_t > & intervalSize,
                        cudaStream_t copyStream )
{
   const uint_t & Nx = std::get<0>(intervalSize),
                & Ny = std::get<1>(intervalSize),
                & Nz = std::get<2>(intervalSize),
                & Nf = std::get<3>(intervalSize);

   const uint_t & srcX = std::get<0>(srcOffset),
                & srcY = std::get<1>(srcOffset),
                & srcZ = std::get<2>(srcOffset),
                & srcF = std::get<3>(srcOffset);

   const uint_t & dstX = std::get<0>(dstOffset),
                & dstY = std::get<1>(dstOffset),
                & dstZ = std::get<2>(dstOffset),
                & dstF = std::get<3>(dstOffset);

   auto copyFunctor = [&](uint_t dstCoordZ, uint_t srcCoordZ, uint_t zIntervalSize) {
         cudaMemcpy3DParms p;
         std::memset( &p, 0, sizeof(p) );

         p.srcPos = make_cudaPos( srcF * typeSize, srcX, srcCoordZ * srcAllocSizeY + srcY );
         p.srcPtr = make_cudaPitchedPtr( src, Nf * typeSize, Nf * typeSize, Nx );

         p.dstPos = make_cudaPos( dstF * typeSize, dstX, dstCoordZ * dstAllocSizeY + dstY );
         p.dstPtr = make_cudaPitchedPtr( dst.ptr, dst.pitch, dst.xsize, dst.ysize );

         p.extent = make_cudaExtent( Nf * typeSize, Nx, Ny * zIntervalSize );
         p.kind = cudaMemcpyHostToDevice;

         if ( copyStream == 0 )
         {
            WALBERLA_CUDA_CHECK( cudaMemcpy3D(&p) );
         }
         else
         {
            // Using cudaMemcpy3DAsync requires page-locked memory on the host!
            WALBERLA_CUDA_CHECK( cudaMemcpy3DAsync(&p, copyStream) );
         }
   };

   if ( Nz == 1 || ( Ny == dstAllocSizeY && Ny == srcAllocSizeY ) )
   {
      copyFunctor( dstZ, srcZ, Nz );
   }
   else
   {
      for( uint_t z = 0; z < Nz; ++z )
      {
         copyFunctor( dstZ + z, srcZ + z, 1 );
      }
   }
}


void copyDevToHostFZYX( unsigned char* dst, const cudaPitchedPtr& src,
                        std::tuple< uint_t, uint_t, uint_t, uint_t > & dstOffset,
                        std::tuple< uint_t, uint_t, uint_t, uint_t > & srcOffset,
                        uint_t dstAllocSizeZ, uint_t srcAllocSizeZ, uint_t typeSize,
                        std::tuple< uint_t, uint_t, uint_t, uint_t > & intervalSize,
                        cudaStream_t copyStream )
{
   const uint_t & Nx = std::get<0>(intervalSize),
                & Ny = std::get<1>(intervalSize),
                & Nz = std::get<2>(intervalSize),
                & Nf = std::get<3>(intervalSize);

   const uint_t & srcX = std::get<0>(srcOffset),
                & srcY = std::get<1>(srcOffset),
                & srcZ = std::get<2>(srcOffset),
                & srcF = std::get<3>(srcOffset);

   const uint_t & dstX = std::get<0>(dstOffset),
                & dstY = std::get<1>(dstOffset),
                & dstZ = std::get<2>(dstOffset),
                & dstF = std::get<3>(dstOffset);

   auto copyFunctor = [&](uint_t dstCoordF, uint_t srcCoordF, uint_t fIntervalSize) {
      cudaMemcpy3DParms p;
      std::memset( &p, 0, sizeof(p) );

      p.srcPos = make_cudaPos( srcX * typeSize, srcY, srcCoordF * srcAllocSizeZ + srcZ );
      p.srcPtr = make_cudaPitchedPtr( src.ptr, src.pitch, src.xsize, src.ysize );

      p.dstPos = make_cudaPos( dstX * typeSize, dstY, dstCoordF * dstAllocSizeZ + dstZ );
      p.dstPtr = make_cudaPitchedPtr( dst, Nx * typeSize, Nx * typeSize, Ny );

      p.extent = make_cudaExtent( Nx * typeSize, Ny, Nz * fIntervalSize );
      p.kind = cudaMemcpyDeviceToHost;

      if ( copyStream == 0 )
      {
         WALBERLA_CUDA_CHECK( cudaMemcpy3D(&p) );
      }
      else
      {
         // Using cudaMemcpy3DAsync requires page-locked memory on the host!
         WALBERLA_CUDA_CHECK( cudaMemcpy3DAsync(&p, copyStream) );
      }
   };

   if( Nf == 1 || ( Nz == dstAllocSizeZ && Nz == srcAllocSizeZ ) )
   {
      copyFunctor( dstF, srcF, Nf );
   }
   else
   {
      for( uint_t f = 0; f < Nf; ++f )
      {
         copyFunctor( dstF + f, srcF + f, 1 );
      }
   }
}


void copyDevToHostZYXF( unsigned char* dst, const cudaPitchedPtr& src,
                        std::tuple< uint_t, uint_t, uint_t, uint_t > & dstOffset,
                        std::tuple< uint_t, uint_t, uint_t, uint_t > & srcOffset,
                        uint_t dstAllocSizeY, uint_t srcAllocSizeY, uint_t typeSize,
                        std::tuple< uint_t, uint_t, uint_t, uint_t > & intervalSize,
                        cudaStream_t copyStream )
{
   const uint_t & Nx = std::get<0>(intervalSize),
                & Ny = std::get<1>(intervalSize),
                & Nz = std::get<2>(intervalSize),
                & Nf = std::get<3>(intervalSize);

   const uint_t & srcX = std::get<0>(srcOffset),
                & srcY = std::get<1>(srcOffset),
                & srcZ = std::get<2>(srcOffset),
                & srcF = std::get<3>(srcOffset);

   const uint_t & dstX = std::get<0>(dstOffset),
                & dstY = std::get<1>(dstOffset),
                & dstZ = std::get<2>(dstOffset),
                & dstF = std::get<3>(dstOffset);

   auto copyFunctor = [&](uint_t dstCoordZ, uint_t srcCoordZ, uint_t zIntervalSize) {
      cudaMemcpy3DParms p;
      std::memset( &p, 0, sizeof(p) );

      p.srcPos = make_cudaPos( srcF * typeSize, srcX, srcCoordZ * srcAllocSizeY + srcY );
      p.srcPtr = make_cudaPitchedPtr( src.ptr, src.pitch, src.xsize, src.ysize );

      p.dstPos = make_cudaPos( dstF * typeSize, dstX, dstCoordZ * dstAllocSizeY + dstY );
      p.dstPtr = make_cudaPitchedPtr( dst, Nf * typeSize, Nf * typeSize, Nx );

      p.extent = make_cudaExtent( Nf * typeSize, Nx, Ny * zIntervalSize );

      p.kind = cudaMemcpyDeviceToHost;

      if ( copyStream == 0 )
      {
         WALBERLA_CUDA_CHECK( cudaMemcpy3D(&p) );
      }
      else
      {
         // Using cudaMemcpy3DAsync requires page-locked memory on the host!
         WALBERLA_CUDA_CHECK( cudaMemcpy3DAsync(&p, copyStream) );
      }
   };


   if ( Nz == 1 || ( Ny == dstAllocSizeY && Ny == srcAllocSizeY ) )
   {
      copyFunctor( dstZ, srcZ, Nz );
   }
   else
   {
      for( uint_t z = 0; z < Nz; ++z )
      {
         copyFunctor( dstZ + z, srcZ + z, 1 );
      }
   }
}

} // namespace cuda
} // namespace walberla
