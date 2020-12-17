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

#include <mesa_pd/data/DataTypes.h>
#include <mesa_pd/data/ParticleStorage.h>
#include <mesa_pd/mpi/notifications/NotificationType.h>
#include <mesa_pd/mpi/notifications/reset.h>

#include <core/mpi/Datatype.h>
#include <core/mpi/RecvBuffer.h>
#include <core/mpi/SendBuffer.h>

namespace walberla {
namespace mesa_pd {

/**
 * Trasmits force and torque information.
 */
class {{name}}
{
public:
   struct Parameters
   {
      id_t uid_;
      {%- for prop in properties %}
      {{prop.type}} {{prop.name}}_;
      {%- endfor %}
   };

   inline explicit {{name}}( const data::Particle& p ) : p_(p) {}

   const data::Particle& p_;
};

template <>
void reset<{{name}}>(data::Particle& p)
{
   {%- for prop in properties %}
   p.set{{prop.name | capFirst}}( {{prop.resetValue}} );
   {%- endfor %}
}

void reduce(data::Particle&& p, const {{name}}::Parameters& objparam)
{
   {%- for prop in properties %}
   p.get{{prop.name | capFirst}}Ref() += objparam.{{prop.name}}_;
   {%- endfor %}
}

void update(data::Particle&& p, const {{name}}::Parameters& objparam)
{
   {%- for prop in properties %}
   p.set{{prop.name | capFirst}}( objparam.{{prop.name}}_ );
   {%- endfor %}
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
mpi::GenericSendBuffer<T,G>& operator<<( mpi::GenericSendBuffer<T,G> & buf, const mesa_pd::{{name}}& obj )
{
   buf.addDebugMarker( "pn" );
   buf << obj.p_.getUid();
   {%- for prop in properties %}
   buf << obj.p_.get{{prop.name | capFirst}}();
   {%- endfor %}
   return buf;
}

template< typename T>    // Element type  of RecvBuffer
mpi::GenericRecvBuffer<T>& operator>>( mpi::GenericRecvBuffer<T> & buf, mesa_pd::{{name}}::Parameters& objparam )
{
   buf.readDebugMarker( "pn" );
   buf >> objparam.uid_;
   {%- for prop in properties %}
   buf >> objparam.{{prop.name}}_;
   {%- endfor %}
   return buf;
}

template< >
struct BufferSizeTrait< mesa_pd::{{name}} > {
   static const bool constantSize = true;
   static const uint_t size = BufferSizeTrait<id_t>::size +
   {%- for prop in properties %}
                              BufferSizeTrait<{{prop.type}}>::size +
   {%- endfor %}
                              mpi::BUFFER_DEBUG_OVERHEAD;
};

} // mpi
} // walberla
