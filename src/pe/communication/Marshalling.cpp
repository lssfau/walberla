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
//! \file Marshalling.cpp
//! \author Tobias Preclik
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//! \brief Marshalling of objects for data transmission or storage.
//
//======================================================================================================================

#include "Marshalling.h"

namespace walberla {
namespace pe {
namespace communication {

//*************************************************************************************************
void marshal( mpi::SendBuffer& buffer, const MPIRigidBodyTrait& obj ) {
   buffer << obj.getOwner();
}
//*************************************************************************************************

//*************************************************************************************************
void unmarshal( mpi::RecvBuffer& buffer, MPIRigidBodyTraitParameter& objparam ) {
   buffer >> objparam.owner_;
}

//*************************************************************************************************
void marshal( mpi::SendBuffer& buffer, const RigidBody& obj ) {
   marshal( buffer, obj.MPITrait );

   buffer << obj.getSystemID();
   buffer << obj.getID();
   buffer << obj.isCommunicating();
   buffer << obj.hasInfiniteMass();
   buffer << obj.getPosition();
   buffer << obj.hasSuperBody();
   buffer << obj.getQuaternion();
   if( !obj.hasSuperBody() )
   {
      buffer << obj.getLinearVel();
      buffer << obj.getAngularVel();
   }
}
//*************************************************************************************************

//*************************************************************************************************
void unmarshal( mpi::RecvBuffer& buffer, RigidBodyParameters& objparam ) {
   unmarshal( buffer, objparam.mpiTrait_ );

   buffer >> objparam.sid_;
   buffer >> objparam.uid_;
   buffer >> objparam.communicating_;
   buffer >> objparam.infiniteMass_;
   buffer >> objparam.gpos_;
   buffer >> objparam.hasSuperBody_;
   buffer >> objparam.q_;

   if( !objparam.hasSuperBody_ )
   {
      buffer >> objparam.v_;
      buffer >> objparam.w_;
   }
}

//*************************************************************************************************
//*************************************************************************************************
//*************************************************************************************************
//*************************************************************************************************
//*************************************************************************************************
//*************************************************************************************************

//*************************************************************************************************
void marshal( mpi::SendBuffer& buffer, const GeomPrimitive& obj ) {
   marshal( buffer, static_cast<const RigidBody&>( obj ) );
   buffer << obj.getMaterial();
}
//*************************************************************************************************


//*************************************************************************************************
void unmarshal( mpi::RecvBuffer& buffer, GeomPrimitiveParameters& objparam ) {
   unmarshal( buffer, static_cast<RigidBodyParameters&>( objparam ) );
   buffer >> objparam.material_;
}

}  // namespace communication
}  // namespace pe
}  // namespace walberla
