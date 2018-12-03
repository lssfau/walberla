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
//! \file BlockInfo.h
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//
//======================================================================================================================

#pragma once

#include <core/mpi/RecvBuffer.h>
#include <core/mpi/SendBuffer.h>

#include <ostream>

namespace walberla {
namespace pe_coupling {

struct BlockInfo {
   // lbm quantities
   uint_t numberOfCells;
   uint_t numberOfFluidCells;
   uint_t numberOfNearBoundaryCells;
   // pe quantities
   uint_t numberOfLocalBodies;
   uint_t numberOfShadowBodies;
   uint_t numberOfContacts;
   // coupling quantities
   uint_t numberOfPeSubCycles;

   BlockInfo()
         : numberOfCells(0), numberOfFluidCells(0), numberOfNearBoundaryCells(0),
           numberOfLocalBodies(0), numberOfShadowBodies(0), numberOfContacts(0),
           numberOfPeSubCycles(0) {}

   BlockInfo(const uint_t numCells, const uint_t numFluidCells, const uint_t numNearBoundaryCells,
             const uint_t numLocalBodies, const uint_t numShadowBodies, const uint_t numContacts,
             const uint_t numPeSubCycles)
         : numberOfCells(numCells), numberOfFluidCells(numFluidCells), numberOfNearBoundaryCells(numNearBoundaryCells),
           numberOfLocalBodies(numLocalBodies), numberOfShadowBodies(numShadowBodies), numberOfContacts(numContacts),
           numberOfPeSubCycles(numPeSubCycles) {}
};


inline
std::ostream& operator<<( std::ostream& os, const BlockInfo& bi )
{
   os << bi.numberOfCells << " / " << bi.numberOfFluidCells << " / " << bi.numberOfNearBoundaryCells << " / "
      << bi.numberOfLocalBodies << " / " << bi.numberOfShadowBodies << " / " << bi.numberOfContacts << " / "
      << bi.numberOfPeSubCycles;
   return os;
}

template< typename T,    // Element type of SendBuffer
          typename G>    // Growth policy of SendBuffer
mpi::GenericSendBuffer<T,G>& operator<<( mpi::GenericSendBuffer<T,G> & buf, const BlockInfo& info )
{
   buf.addDebugMarker( "pca" );
   buf << info.numberOfCells << info.numberOfFluidCells << info.numberOfNearBoundaryCells
       << info.numberOfLocalBodies << info.numberOfShadowBodies << info.numberOfContacts
       << info.numberOfPeSubCycles;
   return buf;
}

template< typename T>    // Element type of SendBuffer
mpi::GenericRecvBuffer<T>& operator>>( mpi::GenericRecvBuffer<T> & buf, BlockInfo& info )
{
   buf.readDebugMarker( "pca" );
   buf >> info.numberOfCells >> info.numberOfFluidCells >> info.numberOfNearBoundaryCells
       >> info.numberOfLocalBodies >> info.numberOfShadowBodies >> info.numberOfContacts
       >> info.numberOfPeSubCycles;
   return buf;
}

} // namespace pe_coupling
} // namespace walberla
