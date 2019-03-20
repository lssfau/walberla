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
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#pragma once

#include <core/mpi/RecvBuffer.h>
#include <core/mpi/SendBuffer.h>

#include <ostream>

namespace walberla {
namespace pe {

struct BlockInfo
{
   uint_t computationalWeight;
   uint_t communicationWeight;

   BlockInfo()
      : computationalWeight(0)
      , communicationWeight(0)
   {}

   BlockInfo(const uint_t compWeight, const uint_t commWeight)
      : computationalWeight(compWeight)
      , communicationWeight(commWeight)
   {}

   BlockInfo&      operator+=( const BlockInfo& rhs )
   {
      computationalWeight += rhs.computationalWeight;
      communicationWeight += rhs.communicationWeight;
      return *this;
   }
};

inline
BlockInfo operator+ ( const BlockInfo& lhs, const BlockInfo& rhs )
{
   BlockInfo result(lhs);
   result += rhs;
   return result;
}

inline
std::ostream& operator<<( std::ostream& os, const BlockInfo& bi )
{
   os << bi.computationalWeight << " / " << bi.communicationWeight;
   return os;
}

template< typename T,    // Element type of SendBuffer
          typename G>    // Growth policy of SendBuffer
mpi::GenericSendBuffer<T,G>& operator<<( mpi::GenericSendBuffer<T,G> & buf, const BlockInfo& info )
{
   buf.addDebugMarker( "pa" );
   buf << info.computationalWeight << info.communicationWeight;
   return buf;
}

template< typename T>    // Element type of SendBuffer
mpi::GenericRecvBuffer<T>& operator>>( mpi::GenericRecvBuffer<T> & buf, BlockInfo& info )
{
   buf.readDebugMarker( "pa" );
   buf >> info.computationalWeight >> info.communicationWeight;
   return buf;
}

}
}
