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
//! \ingroup gpu
//! \author Paulo Carvalho <prcjunior@inf.ufpr.br>
//! \author João Victor Tozatti Risso <jvtrisso@inf.ufpr.br>
//! \brief Copy routines of 4D intervals involving GPU buffers.
//
//======================================================================================================================

#include "core/debug/Debug.h"

#include "GPUCopy.h"

#include <cstring>


namespace walberla {
namespace gpu
{

void copyDevToDevFZYX( const gpuPitchedPtr& dst, const gpuPitchedPtr& src,
                       std::tuple< uint_t, uint_t, uint_t, uint_t > & dstOffset,
                       std::tuple< uint_t, uint_t, uint_t, uint_t > & srcOffset,
                       uint_t dstAllocSizeZ, uint_t srcAllocSizeZ, uint_t typeSize,
                       std::tuple< uint_t, uint_t, uint_t, uint_t > & intervalSize,
                       gpuStream_t copyStream )
{
   const uint_t & Nx = std::get<0>(intervalSize);
   const uint_t & Ny = std::get<1>(intervalSize);
   const uint_t & Nz = std::get<2>(intervalSize);
   const uint_t & Nf = std::get<3>(intervalSize);

   const uint_t & srcX = std::get<0>(srcOffset);
   const uint_t & srcY = std::get<1>(srcOffset);
   const uint_t & srcZ = std::get<2>(srcOffset);
   const uint_t & srcF = std::get<3>(srcOffset);

   const uint_t & dstX = std::get<0>(dstOffset);
   const uint_t & dstY = std::get<1>(dstOffset);
   const uint_t & dstZ = std::get<2>(dstOffset);
   const uint_t & dstF = std::get<3>(dstOffset);

   auto copyFunctor = [&](uint_t dstCoordF, uint_t srcCoordF, uint_t fIntervalSize) {
      WALBERLA_ASSERT( fIntervalSize == 1 || ( Nz == dstAllocSizeZ && Nz == srcAllocSizeZ ) );

      gpuMemcpy3DParms p;
      std::memset( &p, 0, sizeof(p) );

      p.srcPos = make_gpuPos( srcX * typeSize, srcY, srcCoordF * srcAllocSizeZ + srcZ );
      p.srcPtr = make_gpuPitchedPtr( src.ptr, src.pitch, src.xsize, src.ysize );

      p.dstPos = make_gpuPos( dstX * typeSize, dstY, dstCoordF * dstAllocSizeZ + dstZ );
      p.dstPtr = make_gpuPitchedPtr( dst.ptr, dst.pitch, dst.xsize, dst.ysize );

      p.extent = make_gpuExtent( Nx * typeSize, Ny, Nz * fIntervalSize );
      p.kind = gpuMemcpyDeviceToDevice;

      if ( copyStream == nullptr )
      {
         WALBERLA_GPU_CHECK( gpuMemcpy3D(&p) )
      }
      else
      {
         // Using hipMemcpy3DAsync requires page-locked memory on the host!
         WALBERLA_GPU_CHECK( gpuMemcpy3DAsync(&p, copyStream) )
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


void copyDevToDevZYXF( const gpuPitchedPtr& dst, const gpuPitchedPtr& src,
                       std::tuple< uint_t, uint_t, uint_t, uint_t > & dstOffset,
                       std::tuple< uint_t, uint_t, uint_t, uint_t > & srcOffset,
                       uint_t dstAllocSizeY, uint_t srcAllocSizeY, uint_t typeSize,
                       std::tuple< uint_t, uint_t, uint_t, uint_t > & intervalSize,
                       gpuStream_t copyStream )
{
   const uint_t & Nx = std::get<0>(intervalSize);
   const uint_t & Ny = std::get<1>(intervalSize);
   const uint_t & Nz = std::get<2>(intervalSize);
   const uint_t & Nf = std::get<3>(intervalSize);

   const uint_t & srcX = std::get<0>(srcOffset);
   const uint_t & srcY = std::get<1>(srcOffset);
   const uint_t & srcZ = std::get<2>(srcOffset);
   const uint_t & srcF = std::get<3>(srcOffset);

   const uint_t & dstX = std::get<0>(dstOffset);
   const uint_t & dstY = std::get<1>(dstOffset);
   const uint_t & dstZ = std::get<2>(dstOffset);
   const uint_t & dstF = std::get<3>(dstOffset);

   auto copyFunctor = [&](uint_t dstCoordZ, uint_t srcCoordZ, uint_t zIntervalSize) {
      gpuMemcpy3DParms p;
      std::memset( &p, 0, sizeof(p) );

      p.srcPos = make_gpuPos( srcF * typeSize, srcX, srcCoordZ * srcAllocSizeY + srcY );
      p.srcPtr = make_gpuPitchedPtr( src.ptr, src.pitch, src.xsize, src.ysize );

      p.dstPos = make_gpuPos( dstF * typeSize, dstX, dstCoordZ * dstAllocSizeY + dstY );
      p.dstPtr = make_gpuPitchedPtr( dst.ptr, dst.pitch, dst.xsize, dst.ysize );

      p.extent = make_gpuExtent( Nf * typeSize, Nx, Ny * zIntervalSize );
      p.kind = gpuMemcpyDeviceToDevice;

      if ( copyStream == nullptr )
      {
         WALBERLA_GPU_CHECK( gpuMemcpy3D(&p) )
      }
      else
      {
         WALBERLA_GPU_CHECK( gpuMemcpy3DAsync(&p, copyStream) )
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


void copyHostToDevFZYX( const gpuPitchedPtr& dst, unsigned char* src,
                        std::tuple< uint_t, uint_t, uint_t, uint_t > & dstOffset,
                        std::tuple< uint_t, uint_t, uint_t, uint_t > & srcOffset,
                        uint_t dstAllocSizeZ, uint_t srcAllocSizeZ, uint_t typeSize,
                        std::tuple< uint_t, uint_t, uint_t, uint_t > & intervalSize,
                        gpuStream_t copyStream )
{
   const uint_t & Nx = std::get<0>(intervalSize);
   const uint_t & Ny = std::get<1>(intervalSize);
   const uint_t & Nz = std::get<2>(intervalSize);
   const uint_t & Nf = std::get<3>(intervalSize);

   const uint_t & srcX = std::get<0>(srcOffset);
   const uint_t & srcY = std::get<1>(srcOffset);
   const uint_t & srcZ = std::get<2>(srcOffset);
   const uint_t & srcF = std::get<3>(srcOffset);

   const uint_t & dstX = std::get<0>(dstOffset);
   const uint_t & dstY = std::get<1>(dstOffset);
   const uint_t & dstZ = std::get<2>(dstOffset);
   const uint_t & dstF = std::get<3>(dstOffset);

   auto copyFunctor = [&](uint_t dstCoordF, uint_t srcCoordF, uint_t fIntervalSize) {
      gpuMemcpy3DParms p;
      std::memset( &p, 0, sizeof(p) );

      p.srcPos = make_gpuPos( srcX * typeSize, srcY, srcCoordF * srcAllocSizeZ + srcZ );
      p.srcPtr = make_gpuPitchedPtr( src, Nx * typeSize, Nx * typeSize, Ny );

      p.dstPos = make_gpuPos( dstX * typeSize, dstY, dstCoordF * dstAllocSizeZ + dstZ );
      p.dstPtr = make_gpuPitchedPtr( dst.ptr, dst.pitch, dst.xsize, dst.ysize );

      p.extent = make_gpuExtent( Nx * typeSize, Ny, Nz * fIntervalSize );
      p.kind = gpuMemcpyHostToDevice;

      if (copyStream == nullptr)
      {
         WALBERLA_GPU_CHECK( gpuMemcpy3D(&p) )
      }
      else
      {
         // Using gpuMemcpy3DAsync requires page-locked memory on the host!
         WALBERLA_GPU_CHECK( gpuMemcpy3DAsync(&p, copyStream) )
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

void copyHostToDevZYXF( const gpuPitchedPtr& dst, unsigned char* src,
                        std::tuple< uint_t, uint_t, uint_t, uint_t > & dstOffset,
                        std::tuple< uint_t, uint_t, uint_t, uint_t > & srcOffset,
                        uint_t dstAllocSizeY, uint_t srcAllocSizeY, uint_t typeSize,
                        std::tuple< uint_t, uint_t, uint_t, uint_t > & intervalSize,
                        gpuStream_t copyStream )
{
   const uint_t & Nx = std::get<0>(intervalSize);
   const uint_t & Ny = std::get<1>(intervalSize);
   const uint_t & Nz = std::get<2>(intervalSize);
   const uint_t & Nf = std::get<3>(intervalSize);

   const uint_t & srcX = std::get<0>(srcOffset);
   const uint_t & srcY = std::get<1>(srcOffset);
   const uint_t & srcZ = std::get<2>(srcOffset);
   const uint_t & srcF = std::get<3>(srcOffset);

   const uint_t & dstX = std::get<0>(dstOffset);
   const uint_t & dstY = std::get<1>(dstOffset);
   const uint_t & dstZ = std::get<2>(dstOffset);
   const uint_t & dstF = std::get<3>(dstOffset);

   auto copyFunctor = [&](uint_t dstCoordZ, uint_t srcCoordZ, uint_t zIntervalSize) {
         gpuMemcpy3DParms p;
         std::memset( &p, 0, sizeof(p) );

         p.srcPos = make_gpuPos( srcF * typeSize, srcX, srcCoordZ * srcAllocSizeY + srcY );
         p.srcPtr = make_gpuPitchedPtr( src, Nf * typeSize, Nf * typeSize, Nx );

         p.dstPos = make_gpuPos( dstF * typeSize, dstX, dstCoordZ * dstAllocSizeY + dstY );
         p.dstPtr = make_gpuPitchedPtr( dst.ptr, dst.pitch, dst.xsize, dst.ysize );

         p.extent = make_gpuExtent( Nf * typeSize, Nx, Ny * zIntervalSize );
         p.kind = gpuMemcpyHostToDevice;

         if ( copyStream == nullptr )
         {
            WALBERLA_GPU_CHECK( gpuMemcpy3D(&p) )
         }
         else
         {
            // Using gpuMemcpy3DAsync requires page-locked memory on the host!
            WALBERLA_GPU_CHECK( gpuMemcpy3DAsync(&p, copyStream) )
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


void copyDevToHostFZYX( unsigned char* dst, const gpuPitchedPtr& src,
                        std::tuple< uint_t, uint_t, uint_t, uint_t > & dstOffset,
                        std::tuple< uint_t, uint_t, uint_t, uint_t > & srcOffset,
                        uint_t dstAllocSizeZ, uint_t srcAllocSizeZ, uint_t typeSize,
                        std::tuple< uint_t, uint_t, uint_t, uint_t > & intervalSize,
                        gpuStream_t copyStream )
{
   const uint_t & Nx = std::get<0>(intervalSize);
   const uint_t & Ny = std::get<1>(intervalSize);
   const uint_t & Nz = std::get<2>(intervalSize);
   const uint_t & Nf = std::get<3>(intervalSize);

   const uint_t & srcX = std::get<0>(srcOffset);
   const uint_t & srcY = std::get<1>(srcOffset);
   const uint_t & srcZ = std::get<2>(srcOffset);
   const uint_t & srcF = std::get<3>(srcOffset);

   const uint_t & dstX = std::get<0>(dstOffset);
   const uint_t & dstY = std::get<1>(dstOffset);
   const uint_t & dstZ = std::get<2>(dstOffset);
   const uint_t & dstF = std::get<3>(dstOffset);

   auto copyFunctor = [&](uint_t dstCoordF, uint_t srcCoordF, uint_t fIntervalSize) {
      gpuMemcpy3DParms p;
      std::memset( &p, 0, sizeof(p) );

      p.srcPos = make_gpuPos( srcX * typeSize, srcY, srcCoordF * srcAllocSizeZ + srcZ );
      p.srcPtr = make_gpuPitchedPtr( src.ptr, src.pitch, src.xsize, src.ysize );

      p.dstPos = make_gpuPos( dstX * typeSize, dstY, dstCoordF * dstAllocSizeZ + dstZ );
      p.dstPtr = make_gpuPitchedPtr( dst, Nx * typeSize, Nx * typeSize, Ny );

      p.extent = make_gpuExtent( Nx * typeSize, Ny, Nz * fIntervalSize );
      p.kind = gpuMemcpyDeviceToHost;

      if ( copyStream == nullptr )
      {
         WALBERLA_GPU_CHECK( gpuMemcpy3D(&p) );
      }
      else
      {
         // Using gpuMemcpy3DAsync requires page-locked memory on the host!
         WALBERLA_GPU_CHECK( gpuMemcpy3DAsync(&p, copyStream) )
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


void copyDevToHostZYXF( unsigned char* dst, const gpuPitchedPtr& src,
                        std::tuple< uint_t, uint_t, uint_t, uint_t > & dstOffset,
                        std::tuple< uint_t, uint_t, uint_t, uint_t > & srcOffset,
                        uint_t dstAllocSizeY, uint_t srcAllocSizeY, uint_t typeSize,
                        std::tuple< uint_t, uint_t, uint_t, uint_t > & intervalSize,
                        gpuStream_t copyStream )
{
   const uint_t & Nx = std::get<0>(intervalSize);
   const uint_t & Ny = std::get<1>(intervalSize);
   const uint_t & Nz = std::get<2>(intervalSize);
   const uint_t & Nf = std::get<3>(intervalSize);

   const uint_t & srcX = std::get<0>(srcOffset);
   const uint_t & srcY = std::get<1>(srcOffset);
   const uint_t & srcZ = std::get<2>(srcOffset);
   const uint_t & srcF = std::get<3>(srcOffset);

   const uint_t & dstX = std::get<0>(dstOffset);
   const uint_t & dstY = std::get<1>(dstOffset);
   const uint_t & dstZ = std::get<2>(dstOffset);
   const uint_t & dstF = std::get<3>(dstOffset);

   auto copyFunctor = [&](uint_t dstCoordZ, uint_t srcCoordZ, uint_t zIntervalSize) {
      gpuMemcpy3DParms p;
      std::memset( &p, 0, sizeof(p) );

      p.srcPos = make_gpuPos( srcF * typeSize, srcX, srcCoordZ * srcAllocSizeY + srcY );
      p.srcPtr = make_gpuPitchedPtr( src.ptr, src.pitch, src.xsize, src.ysize );

      p.dstPos = make_gpuPos( dstF * typeSize, dstX, dstCoordZ * dstAllocSizeY + dstY );
      p.dstPtr = make_gpuPitchedPtr( dst, Nf * typeSize, Nf * typeSize, Nx );

      p.extent = make_gpuExtent( Nf * typeSize, Nx, Ny * zIntervalSize );

      p.kind = gpuMemcpyDeviceToHost;

      if ( copyStream == nullptr )
      {
         WALBERLA_GPU_CHECK( gpuMemcpy3D(&p) )
      }
      else
      {
         // Using gpuMemcpy3DAsync requires page-locked memory on the host!
         WALBERLA_GPU_CHECK( gpuMemcpy3DAsync(&p, copyStream) )
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

} // namespace gpu
} // namespace walberla
