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
//! \file BodyData.h
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#pragma once


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <pe/Types.h>

#include <core/math/Quaternion.h>
#include <core/math/Vector3.h>

namespace walberla {
namespace pe {
namespace debug {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

class BodyData
{
public:
   BodyData();
   explicit BodyData( ConstBodyID rb );

   id_t uid;
   id_t sid;
   Vec3 pos;
   Quat rot;
   Vec3 v;
   Vec3 w;

   static
   bool sortUID(const BodyData& a, const BodyData& b) {return a.uid < b.uid;};

   static
   bool sortSID(const BodyData& a, const BodyData& b) {return a.sid < b.sid;};
};

bool checkEqual(const BodyData& bd1, const BodyData& bd2);

} // debug
} // namespace pe
} // namespace walberla

//======================================================================================================================
//
//  Send/Recv Buffer Serialization Specialization
//
//======================================================================================================================

namespace walberla {
namespace mpi {

   template< typename T,    // Element type of SendBuffer
             typename G>    // Growth policy of SendBuffer
   mpi::GenericSendBuffer<T,G>& operator<<( mpi::GenericSendBuffer<T,G> & buf, const pe::debug::BodyData & bd )
   {
      buf.addDebugMarker( "bd" );
      buf << bd.uid << bd.sid << bd.pos << bd.rot << bd.v << bd.w;
      return buf;
   }

   template< typename T >    // Element type  of RecvBuffer
   mpi::GenericRecvBuffer<T>& operator>>( mpi::GenericRecvBuffer<T> & buf, pe::debug::BodyData & bd )
   {
      buf.readDebugMarker( "bd" );
      buf >> bd.uid >> bd.sid >> bd.pos >> bd.rot >> bd.v >> bd.w;
      return buf;
   }

   template <>
   struct BufferSizeTrait< walberla::pe::debug::BodyData >
   {
      static const bool constantSize = true;
      static const uint_t size = 2 * BufferSizeTrait<id_t>::size + 3 * BufferSizeTrait<pe::Vec3>::size + 3 * BufferSizeTrait<pe::Quat>::size + mpi::BUFFER_DEBUG_OVERHEAD;
   };
}
}
