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
 * Transmits force and torque information.
 */
class ForceTorqueNotification
{
public:
   struct Parameters
   {
      id_t uid_;
      mesa_pd::Vec3 force_;
      mesa_pd::Vec3 torque_;
   };

   inline explicit ForceTorqueNotification( const data::Particle& p ) : p_(p) {}

   const data::Particle& p_;
};

template <>
inline void reset<ForceTorqueNotification>(data::Particle& p)
{
   p.setForce( Vec3(real_t(0)) );
   p.setTorque( Vec3(real_t(0)) );
}

inline void reduce(data::Particle&& p, const ForceTorqueNotification::Parameters& objparam)
{
   p.getForceRef() += objparam.force_;
   p.getTorqueRef() += objparam.torque_;
}

inline void update(data::Particle&& p, const ForceTorqueNotification::Parameters& objparam)
{
   p.setForce( objparam.force_ );
   p.setTorque( objparam.torque_ );
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
mpi::GenericSendBuffer<T,G>& operator<<( mpi::GenericSendBuffer<T,G> & buf, const mesa_pd::ForceTorqueNotification& obj )
{
   buf.addDebugMarker( "pn" );
   buf << obj.p_.getUid();
   buf << obj.p_.getForce();
   buf << obj.p_.getTorque();
   return buf;
}

template< typename T>    // Element type  of RecvBuffer
mpi::GenericRecvBuffer<T>& operator>>( mpi::GenericRecvBuffer<T> & buf, mesa_pd::ForceTorqueNotification::Parameters& objparam )
{
   buf.readDebugMarker( "pn" );
   buf >> objparam.uid_;
   buf >> objparam.force_;
   buf >> objparam.torque_;
   return buf;
}

template< >
struct BufferSizeTrait< mesa_pd::ForceTorqueNotification > {
   static const bool constantSize = true;
   static const uint_t size = BufferSizeTrait<id_t>::size +
                              BufferSizeTrait<mesa_pd::Vec3>::size +
                              BufferSizeTrait<mesa_pd::Vec3>::size +
                              mpi::BUFFER_DEBUG_OVERHEAD;
};

} // mpi
} // walberla