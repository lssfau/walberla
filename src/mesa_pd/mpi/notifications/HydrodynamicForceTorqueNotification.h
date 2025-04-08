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
class HydrodynamicForceTorqueNotification
{
public:
   struct Parameters
   {
      id_t uid_;
      mesa_pd::Vec3 hydrodynamicForce_;
      mesa_pd::Vec3 hydrodynamicTorque_;
   };

   inline explicit HydrodynamicForceTorqueNotification( const data::Particle& p ) : p_(p) {}

   const data::Particle& p_;
};

template <>
inline void reset<HydrodynamicForceTorqueNotification>(data::Particle& p)
{
   p.setHydrodynamicForce( Vec3(real_t(0)) );
   p.setHydrodynamicTorque( Vec3(real_t(0)) );
}

inline void reduce(data::Particle&& p, const HydrodynamicForceTorqueNotification::Parameters& objparam)
{
   p.getHydrodynamicForceRef() += objparam.hydrodynamicForce_;
   p.getHydrodynamicTorqueRef() += objparam.hydrodynamicTorque_;
}

inline void update(data::Particle&& p, const HydrodynamicForceTorqueNotification::Parameters& objparam)
{
   p.setHydrodynamicForce( objparam.hydrodynamicForce_ );
   p.setHydrodynamicTorque( objparam.hydrodynamicTorque_ );
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
mpi::GenericSendBuffer<T,G>& operator<<( mpi::GenericSendBuffer<T,G> & buf, const mesa_pd::HydrodynamicForceTorqueNotification& obj )
{
   buf.addDebugMarker( "pn" );
   buf << obj.p_.getUid();
   buf << obj.p_.getHydrodynamicForce();
   buf << obj.p_.getHydrodynamicTorque();
   return buf;
}

template< typename T>    // Element type  of RecvBuffer
mpi::GenericRecvBuffer<T>& operator>>( mpi::GenericRecvBuffer<T> & buf, mesa_pd::HydrodynamicForceTorqueNotification::Parameters& objparam )
{
   buf.readDebugMarker( "pn" );
   buf >> objparam.uid_;
   buf >> objparam.hydrodynamicForce_;
   buf >> objparam.hydrodynamicTorque_;
   return buf;
}

template< >
struct BufferSizeTrait< mesa_pd::HydrodynamicForceTorqueNotification > {
   static const bool constantSize = true;
   static const uint_t size = BufferSizeTrait<id_t>::size +
                              BufferSizeTrait<mesa_pd::Vec3>::size +
                              BufferSizeTrait<mesa_pd::Vec3>::size +
                              mpi::BUFFER_DEBUG_OVERHEAD;
};

} // mpi
} // walberla