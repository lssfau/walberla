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
//! \file GeneratedGPUPackInfo.h
//! \ingroup gpu
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#pragma once
#include "domain_decomposition/IBlock.h"

#include "gpu/GPUWrapper.h"

#include "stencil/Directions.h"

namespace walberla::gpu {

/**
 * \brief Data packing/unpacking for ghost layer based communication of a \ref GPUField.
 *
 * Encapsulate information on how to extract data from blocks that should be
 * communicated to neighboring blocks (see \ref pack())
 * and how to inject this data in a receiving block (see \ref unpack()).
 * This involves a memory buffer and two memory copy operations.
 *
 * A special method exists for communication between two blocks which are
 * allocated on the same process (see \ref communicateLocal()).
 * In this case the data does not have be communicated via a buffer,
 * but can be copied directly.
 *
 * Data that is packed in direction "dir" at one block is unpacked in
 * direction "stencil::inverseDir[dir]" at the neighboring block. This
 * behavior must be implemented in \ref communicateLocal()!
 *
 * \ingroup gpu
 */
class GeneratedGPUPackInfo
{
public:
  GeneratedGPUPackInfo() = default;
  virtual ~GeneratedGPUPackInfo() = default;

   /**
    * \brief Pack data from a block into a send buffer.
    *
    * \param dir        pack data for neighbor in this direction
    * \param buffer     buffer for writing the data into
    * \param block      the block whose data should be packed into a buffer
    * \param stream     GPU stream
    */
   virtual void pack  ( stencil::Direction dir, unsigned char *buffer, IBlock *block, gpuStream_t stream ) = 0;
   /**
    * \brief Copy data from one local block to another local block.
    *
    * Both blocks are allocated on the same MPI rank.
    *
    * \param dir       the direction of the communication (from sender to receiver)
    * \param sender    id of block where the data should be copied from
    * \param receiver  id of block where the data should be copied to
    * \param stream     GPU stream
    */
   virtual void communicateLocal  ( stencil::Direction dir, const IBlock *sender, IBlock *receiver, gpuStream_t stream ) = 0;
   /**
    * \brief Unpack data from a receive buffer into a block.
    *
    * \param dir        receive data from neighbor in this direction
    * \param buffer     buffer for reading the data from
    * \param block      the block where the unpacked data should be stored into
    * \param stream     GPU stream
    */
   virtual void unpack( stencil::Direction dir, unsigned char *buffer, IBlock *block, gpuStream_t stream ) = 0;
   virtual uint_t size( stencil::Direction dir, IBlock *block ) = 0;
};

} //namespace walberla::gpu