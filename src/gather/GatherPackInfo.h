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
//! \file GatherPackInfo.h
//! \ingroup gather
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/mpi/RecvBuffer.h"
#include "core/mpi/SendBuffer.h"
#include "domain_decomposition/IBlock.h"


namespace walberla {
namespace gather {


   /****************************************************************************************************************//**
   * Interface for packing & unpacking of data to be gathered
   *
   * Implement this interface, and register your object at MPIGatherScheme
   ********************************************************************************************************************/
   class GatherPackInfo
   {
      public:
         virtual ~GatherPackInfo() = default;

         /**********************************************************************************************************//**
         * Packs all data to be gathered into the given buffer
         *
         * Restrictions:
         * - The packed amount of data for a given block has to stay constant
         * - The packed amount of data for a given block has to be constant for all timesteps
         *    (i.e. for all calls of this function)
         *
         **************************************************************************************************************/
         virtual void packData  ( const IBlock * sender,
                                  mpi::SendBuffer & outBuffer ) = 0;


         /**********************************************************************************************************//**
         * Unpacks the data, packed by the packData() function
         *
         * - function is called on gathering process multiple times per timestep, once for each block
         *   that packed a non-zero message. It has to read the same amount as packData() has written.
         *   Therefore is usually necessary that packData inserts a header describing how much data the block packed
         **************************************************************************************************************/
         virtual void unpackData( mpi::RecvBuffer & buffer ) = 0;


         /**********************************************************************************************************//**
         * Called after a timestep has finished, and unpackData was called for each sending block
         *
         **************************************************************************************************************/
         virtual void gatherFinished() {}
   };


} // namespace gather
} // namespace walberla


