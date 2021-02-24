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
//! \author Lukas Werner <lks.werner@fau.de>
//! //
//======================================================================================================================

#pragma once

#include <core/mpi/RecvBuffer.h>
#include <core/mpi/SendBuffer.h>

#include <ostream>

namespace walberla {
namespace lbm_mesapd_coupling {
namespace amr {

struct BlockInfo {
   // lbm quantities
   uint_t numberOfCells;
   uint_t numberOfFluidCells;
   uint_t numberOfNearBoundaryCells;
   // rpd quantities
   uint_t numberOfLocalParticles;
   uint_t numberOfGhostParticles;
   // coupling quantities
   uint_t numberOfRPDSubCycles;

   BlockInfo()
         : numberOfCells(0), numberOfFluidCells(0), numberOfNearBoundaryCells(0),
           numberOfLocalParticles(0), numberOfGhostParticles(0), numberOfRPDSubCycles(0) {}

   BlockInfo(const uint_t numCells, const uint_t numFluidCells, const uint_t numNearBoundaryCells,
             const uint_t numLocalParticles, const uint_t numGhostParticles, const uint_t numRPDSubCycles)
         : numberOfCells(numCells), numberOfFluidCells(numFluidCells), numberOfNearBoundaryCells(numNearBoundaryCells),
           numberOfLocalParticles(numLocalParticles), numberOfGhostParticles(numGhostParticles), numberOfRPDSubCycles(numRPDSubCycles) {}
};


inline
std::ostream& operator<<( std::ostream& os, const BlockInfo& bi )
{
   os << bi.numberOfCells << " / " << bi.numberOfFluidCells << " / " << bi.numberOfNearBoundaryCells << " / "
      << bi.numberOfLocalParticles << " / "<< bi.numberOfGhostParticles << " / " << bi.numberOfRPDSubCycles;
   return os;
}

template< typename T,    // Element type of SendBuffer
      typename G>    // Growth policy of SendBuffer
mpi::GenericSendBuffer<T,G>& operator<<( mpi::GenericSendBuffer<T,G> & buf, const BlockInfo& info )
{
   buf.addDebugMarker( "pca" );
   buf << info.numberOfCells << info.numberOfFluidCells << info.numberOfNearBoundaryCells
       << info.numberOfLocalParticles << info.numberOfGhostParticles << info.numberOfRPDSubCycles;
   return buf;
}

template< typename T>    // Element type of SendBuffer
mpi::GenericRecvBuffer<T>& operator>>( mpi::GenericRecvBuffer<T> & buf, BlockInfo& info )
{
   buf.readDebugMarker( "pca" );
   buf >> info.numberOfCells >> info.numberOfFluidCells >> info.numberOfNearBoundaryCells
       >> info.numberOfLocalParticles >> info.numberOfGhostParticles >> info.numberOfRPDSubCycles;
   return buf;
}

} // namespace amr
} // namespace lbm_mesapd_coupling
} // namespace walberla
