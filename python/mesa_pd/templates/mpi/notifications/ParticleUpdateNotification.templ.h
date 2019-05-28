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
//! \file ParticleUpdateNotification.h
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

//======================================================================================================================
//
//  THIS FILE IS GENERATED - PLEASE CHANGE THE TEMPLATE !!!
//
//======================================================================================================================

#pragma once

#include <mesa_pd/data/DataTypes.h>
#include <mesa_pd/data/ParticleStorage.h>
#include <mesa_pd/mpi/notifications/NotificationType.h>

#include <core/mpi/Datatype.h>
#include <core/mpi/RecvBuffer.h>
#include <core/mpi/SendBuffer.h>

namespace walberla {
namespace mesa_pd {

/**
 * Updates a ghost particle.
 *
 * Sends all properties marked as ALWAYS.
 */
class ParticleUpdateNotification {
public:
   struct Parameters {
   {%- for prop in properties %}
   {%- if prop.syncMode in ["ALWAYS"] %}
   {{prop.type}} {{prop.name}} {{'{'}}{{prop.defValue}}{{'}'}};
   {%- endif %}
   {%- endfor %}
   };

   inline explicit ParticleUpdateNotification( const data::Particle& particle ) : particle_(particle) {}
   const data::Particle& particle_;
};

template<>
struct NotificationTrait<ParticleUpdateNotification>
{
   static const NotificationType id = PARTICLE_UPDATE_NOTIFICATION;
};

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
mpi::GenericSendBuffer<T,G>& operator<<( mpi::GenericSendBuffer<T,G> & buf, const mesa_pd::ParticleUpdateNotification& obj )
{
   buf.addDebugMarker( "un" );
   {%- for prop in properties %}
   {%- if prop.syncMode in ["ALWAYS"] %}
   buf << obj.particle_.get{{prop.name | capFirst}}();
   {%- endif %}
   {%- endfor %}
   return buf;
}

template< typename T>    // Element type  of RecvBuffer
mpi::GenericRecvBuffer<T>& operator>>( mpi::GenericRecvBuffer<T> & buf, mesa_pd::ParticleUpdateNotification::Parameters& objparam )
{
   buf.readDebugMarker( "un" );
   {%- for prop in properties %}
   {%- if prop.syncMode in ["ALWAYS"] %}
   buf >> objparam.{{prop.name}};
   {%- endif %}
   {%- endfor %}
   return buf;
}

} // mpi
} // walberla
