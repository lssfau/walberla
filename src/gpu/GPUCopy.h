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
//! \ingroup gpu
//! \author Paulo Carvalho <prcjunior@inf.ufpr.br>
//! \author João Victor Tozatti Risso <jvtrisso@inf.ufpr.br>
//! \brief Copy routines of 4D intervals involving GPU buffers.
//
//======================================================================================================================

#pragma once

#include "core/DataTypes.h"

#include "gpu/ErrorChecking.h"
#include "gpu/DeviceWrapper.h"

#include <tuple>


namespace walberla {
namespace gpu
{


//****************************************************************************************************************************
/*! Copy a 4D interval of a device buffer to another device buffer with fzyx memory layout.
 *
 * \param dst           destination buffer
 * \param src           source buffer
 * \param dstOffset     (x, y, z, f)-tuple containing the coordinate of the interval start point in the destination buffer
 * \param srcOffset     (x, y, z, f)-tuple containing the coordinate of the interval start point in the source buffer
 * \param dstAllocSizeZ allocation size in z direction of the destination buffer
 * \param srcAllocSizeZ allocation size in z direction of the source buffer
 * \param typeSize      size of an f element
 * \param intervalSize  interval size
 * \param copyStream    CUDA/HIP stream, if not NULL copy operations will be performed asynchronously
 *****************************************************************************************************************************/
void copyDevToDevFZYX( const gpuPitchedPtr& dst, const gpuPitchedPtr& src,
                       std::tuple< uint_t, uint_t, uint_t, uint_t > & dstOffset,
                       std::tuple< uint_t, uint_t, uint_t, uint_t > & srcOffset,
                       uint_t dstAllocSizeZ, uint_t srcAllocSizeZ, uint_t typeSize,
                       std::tuple< uint_t, uint_t, uint_t, uint_t > & intervalSize,
                       gpuStream_t copyStream );

//****************************************************************************************************************************
/*! Copy a 4D interval of a device buffer to another device buffer with zyxf memory layout.
 *
 * \param dst           destination buffer
 * \param src           source buffer
 * \param dstOffset     (x, y, z, f)-tuple containing the coordinate of the interval start point in the destination buffer
 * \param srcOffset     (x, y, z, f)-tuple containing the coordinate of the interval start point in the source buffer
 * \param dstAllocSizeY allocation size in y direction of the destination buffer
 * \param srcAllocSizeY allocation size in y direction of the source buffer
 * \param typeSize      size of an f element
 * \param intervalSize  interval size
 * \param copyStream    CUDA/HIP stream, if not NULL copy operations will be performed asynchronously
 *****************************************************************************************************************************/
void copyDevToDevZYXF( const gpuPitchedPtr& dst, const gpuPitchedPtr& src,
                       std::tuple< uint_t, uint_t, uint_t, uint_t > & dstOffset,
                       std::tuple< uint_t, uint_t, uint_t, uint_t > & srcOffset,
                       uint_t dstAllocSizeY, uint_t srcAllocSizeY, uint_t typeSize,
                       std::tuple< uint_t, uint_t, uint_t, uint_t > & intervalSize,
                       gpuStream_t copyStream );

//*******************************************************************************************************************
/*! Copy a 4D interval of a host buffer to a device buffer with fzyx memory layout. See \ref copyDevToDevFZYX() for
 * parameter information.
 *******************************************************************************************************************/
void copyHostToDevFZYX( const gpuPitchedPtr& dst, unsigned char* src,
                        std::tuple< uint_t, uint_t, uint_t, uint_t > & dstOffset,
                        std::tuple< uint_t, uint_t, uint_t, uint_t > & srcOffset,
                        uint_t dstAllocSizeZ, uint_t srcAllocSizeZ, uint_t typeSize,
                        std::tuple< uint_t, uint_t, uint_t, uint_t > & intervalSize,
                        gpuStream_t copyStream );

//*******************************************************************************************************************
/*! Copy a 4D interval of a host buffer to a device buffer with zyxf memory layout. See \ref copyDevToDevZYXF() for
 * parameter information.
 *******************************************************************************************************************/
void copyHostToDevZYXF( const gpuPitchedPtr& dst, unsigned char* src,
                        std::tuple< uint_t, uint_t, uint_t, uint_t > & dstOffset,
                        std::tuple< uint_t, uint_t, uint_t, uint_t > & srcOffset,
                        uint_t dstAllocSizeY, uint_t srcAllocSizeY, uint_t typeSize,
                        std::tuple< uint_t, uint_t, uint_t, uint_t > & intervalSize,
                        gpuStream_t copyStream );

//*******************************************************************************************************************
/*! Copy a 4D interval of a device buffer to a host buffer with fzyx memory layout. See \ref copyDevToDevFZYX() for
 * parameter information.
 *******************************************************************************************************************/
void copyDevToHostFZYX( unsigned char* dst, const gpuPitchedPtr& src,
                        std::tuple< uint_t, uint_t, uint_t, uint_t > & dstOffset,
                        std::tuple< uint_t, uint_t, uint_t, uint_t > & srcOffset,
                        uint_t dstAllocSizeZ, uint_t srcAllocSizeZ, uint_t typeSize,
                        std::tuple< uint_t, uint_t, uint_t, uint_t > & intervalSize,
                        gpuStream_t copyStream );

//*******************************************************************************************************************
/*! Copy a 4D interval of a device buffer to a host buffer with zyxf memory layout. See \ref copyDevToDevZYXF() for
 * parameter information.
 *******************************************************************************************************************/
void copyDevToHostZYXF( unsigned char* dst, const gpuPitchedPtr& src,
                        std::tuple< uint_t, uint_t, uint_t, uint_t > & dstOffset,
                        std::tuple< uint_t, uint_t, uint_t, uint_t > & srcOffset,
                        uint_t dstAllocSizeY, uint_t srcAllocSizeY, uint_t typeSize,
                        std::tuple< uint_t, uint_t, uint_t, uint_t > & intervalSize,
                        gpuStream_t copyStream );

} // namespace gpu
} // namespace walberla
