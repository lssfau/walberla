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
//! \file GPUCopy.h
//! \ingroup cuda
//! \author Paulo Carvalho <prcjunior@inf.ufpr.br>
//! \brief Copy routines of 4D intervals involving GPU buffers.
//
//======================================================================================================================

#pragma once

#include "core/DataTypes.h"

#include <cuda_runtime.h>


namespace walberla {
namespace cuda {


//*******************************************************************************************************************
/*! Restricted version of copyDevToDevFZYX() that requires Nf == 1 || ( Nz == dstAllocZ && Nz == srcAllocZ ).
 * See copyDevToDevFZYX() for more details.
 *******************************************************************************************************************/
void copyDevToDevFZYXRestricted( const cudaPitchedPtr& dst, const cudaPitchedPtr& src,
                                 uint_t typeSz, uint_t dstAllocZ, uint_t srcAllocZ,
                                 uint_t dstX, uint_t dstY, uint_t dstZ, uint_t dstF,
                                 uint_t srcX, uint_t srcY, uint_t srcZ, uint_t srcF,
                                 uint_t Nx, uint_t Ny, uint_t Nz, uint_t Nf );


//*******************************************************************************************************************
/*! Copy a 4D interval of a device buffer to another device buffer with fzyx memory layout.
 *
 * \param dst        destination buffer
 * \param src        source buffer
 * \param typeSz     size of an f element
 * \param dstAllocZ  allocation size in z direction of the destination buffer
 * \param srcAllocZ  allocation size in z direction of the source buffer
 * \param dstX       x coordinate of the interval start point in the destination buffer
 * \param dstY       y coordinate of the interval start point in the destination buffer
 * \param dstZ       z coordinate of the interval start point in the destination buffer
 * \param dstF       f coordinate of the interval start point in the destination buffer
 * \param srcX       x coordinate of the interval start point in the source buffer
 * \param srcY       y coordinate of the interval start point in the source buffer
 * \param srcZ       z coordinate of the interval start point in the source buffer
 * \param srcF       f coordinate of the interval start point in the source buffer
 * \param Nx         interval size in x direction
 * \param Ny         interval size in y direction
 * \param Nz         interval size in z direction
 * \param Nf         interval size in f direction
 *******************************************************************************************************************/
inline void copyDevToDevFZYX( const cudaPitchedPtr& dst, const cudaPitchedPtr& src,
                              uint_t typeSz, uint_t dstAllocZ, uint_t srcAllocZ,
                              uint_t dstX, uint_t dstY, uint_t dstZ, uint_t dstF,
                              uint_t srcX, uint_t srcY, uint_t srcZ, uint_t srcF,
                              uint_t Nx, uint_t Ny, uint_t Nz, uint_t Nf )
{
   if( Nf == 1 || ( Nz == dstAllocZ && Nz == srcAllocZ ) )
   {
      copyDevToDevFZYXRestricted( dst, src,
                                  typeSz, dstAllocZ, srcAllocZ,
                                  dstX, dstY, dstZ, dstF,
                                  srcX, srcY, srcZ, srcF,
                                  Nx, Ny, Nz, Nf );
   }
   else
   {
      for( uint_t f = 0; f < Nf; ++f )
      {
         copyDevToDevFZYXRestricted( dst, src,
                                     typeSz, dstAllocZ, srcAllocZ,
                                     dstX, dstY, dstZ, dstF + f,
                                     srcX, srcY, srcZ, srcF + f,
                                     Nx, Ny, Nz, 1 );
      }
   }
}


//*******************************************************************************************************************
/*! Restricted version of copyHostToDevFZYX() that requires Nf == 1 || ( Nz == dstAllocZ && Nz == srcAllocZ ).
 * See copyHostToDevFZYX() for more details.
 *******************************************************************************************************************/
void copyHostToDevFZYXRestricted( const cudaPitchedPtr& dst, unsigned char* src,
                                  uint_t typeSz, uint_t dstAllocZ, uint_t srcAllocZ,
                                  uint_t dstX, uint_t dstY, uint_t dstZ, uint_t dstF,
                                  uint_t srcX, uint_t srcY, uint_t srcZ, uint_t srcF,
                                  uint_t Nx, uint_t Ny, uint_t Nz, uint_t Nf );


//*******************************************************************************************************************
/*! Copy a 4D interval of a host buffer to a device buffer with fzyx memory layout. See copyDevToDevFZYX() for
 * parameter information.
 *******************************************************************************************************************/
inline void copyHostToDevFZYX( const cudaPitchedPtr& dst, unsigned char* src,
                               uint_t typeSz, uint_t dstAllocZ, uint_t srcAllocZ,
                               uint_t dstX, uint_t dstY, uint_t dstZ, uint_t dstF,
                               uint_t srcX, uint_t srcY, uint_t srcZ, uint_t srcF,
                               uint_t Nx, uint_t Ny, uint_t Nz, uint_t Nf )
{
   if( Nf == 1 || ( Nz == dstAllocZ && Nz == srcAllocZ ) )
   {
      copyHostToDevFZYXRestricted( dst, src,
                                   typeSz, dstAllocZ, srcAllocZ,
                                   dstX, dstY, dstZ, dstF,
                                   srcX, srcY, srcZ, srcF,
                                   Nx, Ny, Nz, Nf );
   }
   else
   {
      for( uint_t f = 0; f < Nf; ++f )
      {
         copyHostToDevFZYXRestricted( dst, src,
                                      typeSz, dstAllocZ, srcAllocZ,
                                      dstX, dstY, dstZ, dstF + f,
                                      srcX, srcY, srcZ, srcF + f,
                                      Nx, Ny, Nz, 1 );
      }
   }
}


//*******************************************************************************************************************
/*! Restricted version of copyDevToHostFZYX() that requires Nf == 1 || ( Nz == dstAllocZ && Nz == srcAllocZ ).
 * See copyDevToHostFZYX() for more details.
 *******************************************************************************************************************/
void copyDevToHostFZYXRestricted( unsigned char* dst, const cudaPitchedPtr& src,
                                  uint_t typeSz, uint_t dstAllocZ, uint_t srcAllocZ,
                                  uint_t dstX, uint_t dstY, uint_t dstZ, uint_t dstF,
                                  uint_t srcX, uint_t srcY, uint_t srcZ, uint_t srcF,
                                  uint_t Nx, uint_t Ny, uint_t Nz, uint_t Nf );


//*******************************************************************************************************************
/*! Copy a 4D interval of a device buffer to a host buffer with fzyx memory layout. See copyDevToDevFZYX() for
 * parameter information.
 *******************************************************************************************************************/
inline void copyDevToHostFZYX( unsigned char* dst, const cudaPitchedPtr& src,
                               uint_t typeSz, uint_t dstAllocZ, uint_t srcAllocZ,
                               uint_t dstX, uint_t dstY, uint_t dstZ, uint_t dstF,
                               uint_t srcX, uint_t srcY, uint_t srcZ, uint_t srcF,
                               uint_t Nx, uint_t Ny, uint_t Nz, uint_t Nf )
{
   if( Nf == 1 || ( Nz == dstAllocZ && Nz == srcAllocZ ) )
   {
      copyDevToHostFZYXRestricted( dst, src,
                                   typeSz, dstAllocZ, srcAllocZ,
                                   dstX, dstY, dstZ, dstF,
                                   srcX, srcY, srcZ, srcF,
                                   Nx, Ny, Nz, Nf );
   }
   else
   {
      for( uint_t f = 0; f < Nf; ++f )
      {
         copyDevToHostFZYXRestricted( dst, src,
                                      typeSz, dstAllocZ, srcAllocZ,
                                      dstX, dstY, dstZ, dstF + f,
                                      srcX, srcY, srcZ, srcF + f,
                                      Nx, Ny, Nz, 1 );
      }
   }
}



} // namespace cuda
} // namespace walberla
