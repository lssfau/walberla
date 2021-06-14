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
//! \file
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

//======================================================================================================================
//
//  THIS FILE IS GENERATED - PLEASE CHANGE THE TEMPLATE !!!
//
//======================================================================================================================

#pragma once

#include <mesa_pd/data/ContactHistory.h>
#include <mesa_pd/data/DataTypes.h>
#include <mesa_pd/mpi/notifications/reset.h>

#include <core/mpi/BufferDataTypeExtensions.h>
#include <core/mpi/Datatype.h>
#include <core/mpi/RecvBuffer.h>
#include <core/mpi/SendBuffer.h>

#include <map>

namespace walberla {
namespace mesa_pd {

/**
 * Transmits the contact history
 */
class ContactHistoryNotification
{
public:
   struct Parameters
   {
      id_t uid_;
      std::map<walberla::id_t, walberla::mesa_pd::data::ContactHistory> contactHistory_;
   };

   inline explicit ContactHistoryNotification( const data::Particle& p ) : p_(p) {}

   const data::Particle& p_;
};

template <>
void reset<ContactHistoryNotification>(data::Particle& p)
{
   p.setNewContactHistory(std::map<walberla::id_t, walberla::mesa_pd::data::ContactHistory>());
}

void reduce(data::Particle&& p, const ContactHistoryNotification::Parameters& objparam)
{
   auto& ch = p.getNewContactHistoryRef();
   for (auto& entry : objparam.contactHistory_)
   {
      auto ret = ch.insert(entry);
      WALBERLA_CHECK(ret.second, "entry already present: " << entry.first << " / " << entry.second);
   }
}

}  // namespace mesa_pd
}  // namespace walberla

//======================================================================================================================
//
//  Send/Recv Buffer Serialization Specialization
//
//======================================================================================================================

namespace walberla {
namespace mpi {

template< typename T,    // Element type of SendBuffer
          typename G>    // Growth policy of SendBuffer
mpi::GenericSendBuffer<T,G>& operator<<( mpi::GenericSendBuffer<T,G> & buf, const mesa_pd::ContactHistoryNotification& obj )
{
   buf.addDebugMarker( "ft" );
   buf << obj.p_.getUid();
   buf << obj.p_.getNewContactHistory();
   return buf;
}

template< typename T>    // Element type  of RecvBuffer
mpi::GenericRecvBuffer<T>& operator>>( mpi::GenericRecvBuffer<T> & buf, mesa_pd::ContactHistoryNotification::Parameters& objparam )
{
   buf.readDebugMarker( "ft" );
   buf >> objparam.uid_;
   buf >> objparam.contactHistory_;
   return buf;
}

} // mpi
} // walberla
