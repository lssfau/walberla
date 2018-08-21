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
//! \file Owner.h
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#pragma once

#include "domain_decomposition/IBlockID.h"
#include "core/mpi/RecvBuffer.h"
#include "core/mpi/SendBuffer.h"

#include <iosfwd>
#include <vector>

namespace walberla{
namespace pe{

//=================================================================================================
//
//  STRUCT DEFINITION
//
//=================================================================================================
struct Owner
{
   int                  rank_;     //< rank of the owner of the shadow copy
   IBlockID::IDType     blockID_;  //< block id of the block the shadow copy is located in

   Owner(): rank_(-1), blockID_(0) {}
   Owner(const int rank, const IBlockID::IDType blockID) : rank_(rank), blockID_(blockID) {}
};
//*************************************************************************************************

inline
std::ostream& operator<<(std::ostream& lhs, const Owner& rhs)
{
   lhs << "( " << rhs.rank_ << ", " << rhs.blockID_  << " )";
   return lhs;
}

inline
bool operator==(const Owner& lhs, const Owner& rhs)
{
   return (lhs.rank_ == rhs.rank_) && (lhs.blockID_ == rhs.blockID_);
}

}  // namespace pe
}  // namespace walberla

//======================================================================================================================
//
//  Send/Recv Buffer Serialization Specialization
//
//======================================================================================================================

namespace walberla {
namespace mpi {

   template< typename T,    // Element type of SendBuffer
             typename G >   // Growth policy of SendBuffer
   inline mpi::GenericSendBuffer<T,G>& operator<<( mpi::GenericSendBuffer<T,G> & buf, const pe::Owner& owner )
   {
      buf.addDebugMarker( "ow" );
      buf << owner.rank_ << owner.blockID_;
      return buf;
   }

   template< typename T >    // Element type  of RecvBuffer
   inline mpi::GenericRecvBuffer<T>& operator>>( mpi::GenericRecvBuffer<T> & buf, pe::Owner& owner )
   {
      buf.readDebugMarker( "ow" );
      buf >> owner.rank_ >> owner.blockID_;
      return buf;
   }

   template < >
   struct BufferSizeTrait< pe::Owner > {
      static const bool constantSize = false;
      static const uint_t size = BufferSizeTrait<int>::size + BufferSizeTrait<IBlockID::IDType>::size + mpi::BUFFER_DEBUG_OVERHEAD;
   };
}  // namespace mpi
}  // namespace walberla
