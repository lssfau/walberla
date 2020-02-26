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
//! \file ParticleCopyNotification.h
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
#include <mesa_pd/mpi/ShapePackUnpack.h>
#include <mesa_pd/mpi/notifications/NotificationType.h>

#include <core/mpi/Datatype.h>
#include <core/mpi/RecvBuffer.h>
#include <core/mpi/SendBuffer.h>

namespace walberla {
namespace mesa_pd {

/**
 * A complete particle copy.
 *
 * Copies all properties that are not marked NEVER.
 */
class ParticleCopyNotification
{
public:
   struct Parameters
   {
      {%- for prop in particle.properties %}
      {%- if not prop.syncMode == "NEVER" %}
      {{prop.type}} {{prop.name}} {{'{'}}{{prop.defValue}}{{'}'}};
      {%- endif %}
      {%- endfor %}
   };

   inline explicit ParticleCopyNotification( const data::Particle& particle ) : particle_(particle) {}
   const data::Particle& particle_;
};

inline data::ParticleStorage::iterator createNewParticle(data::ParticleStorage& ps, const ParticleCopyNotification::Parameters& data)
{
   WALBERLA_ASSERT_EQUAL(ps.find(data.uid), ps.end(), "Particle with same uid already existent!");

   auto pIt = ps.create(data.uid);
   {%- for prop in particle.properties %}
   {%- if not prop.syncMode == "NEVER" %}
   pIt->set{{prop.name | capFirst}}(data.{{prop.name}});
   {%- endif %}
   {%- endfor %}
   return pIt;
}

template<>
struct NotificationTrait<ParticleCopyNotification>
{
   static const NotificationType id = PARTICLE_COPY_NOTIFICATION;
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
mpi::GenericSendBuffer<T,G>& operator<<( mpi::GenericSendBuffer<T,G> & buf, const mesa_pd::ParticleCopyNotification& obj )
{
   buf.addDebugMarker( "cn" );
   {%- for prop in particle.properties %}
   {%- if not prop.syncMode == "NEVER" %}
   buf << obj.particle_.get{{prop.name | capFirst}}();
   {%- endif %}
   {%- endfor %}
   return buf;
}

template< typename T>    // Element type  of RecvBuffer
mpi::GenericRecvBuffer<T>& operator>>( mpi::GenericRecvBuffer<T> & buf, mesa_pd::ParticleCopyNotification::Parameters& objparam )
{
   buf.readDebugMarker( "cn" );
   {%- for prop in particle.properties %}
   {%- if not prop.syncMode == "NEVER" %}
   buf >> objparam.{{prop.name}};
   {%- endif %}
   {%- endfor %}
   return buf;
}

} // mpi
} // walberla
