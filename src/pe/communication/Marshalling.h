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
//! \file Marshalling.h
//! \author Tobias Preclik
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//! \brief Marshalling of objects for data transmission or storage.
//
//======================================================================================================================

#pragma once


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <stdexcept>
#include "pe/rigidbody/GeomPrimitive.h"

#include "core/math/Vector3.h"
#include "core/debug/Debug.h"
#include "core/mpi/RecvBuffer.h"
#include "core/mpi/SendBuffer.h"


namespace walberla {
namespace pe {
namespace communication {

struct MPIRigidBodyTraitParameter {
   Owner owner_;
};

//*************************************************************************************************
/*!\brief Marshalling rigid body parameters.
 *
 * \param buffer The buffer to be filled.
 * \param obj The object to be marshalled.
 * \return void
 */
void marshal( mpi::SendBuffer& buffer, const MPIRigidBodyTrait& obj );
//*************************************************************************************************

//*************************************************************************************************
/*!\brief Unmarshalling rigid body parameters.
 *
 * \param buffer The buffer from where to read.
 * \param objparam The object to be reconstructed.
 * \param hasSuperBody False if body is not part of a union. Subordinate bodies in unions do not encode velocities but encode relative positions.
 * \return void
 */
void unmarshal( mpi::RecvBuffer& buffer, MPIRigidBodyTraitParameter& objparam );

//*************************************************************************************************
//*************************************************************************************************
//*************************************************************************************************
//*************************************************************************************************
//*************************************************************************************************
//*************************************************************************************************

struct RigidBodyParameters {
   MPIRigidBodyTraitParameter mpiTrait_;
   bool communicating_, infiniteMass_;
   id_t sid_, uid_;
   Vec3 gpos_, v_, w_;
   bool hasSuperBody_;
   Quat q_;
};

//*************************************************************************************************
/*!\brief Marshalling rigid body parameters.
 *
 * \param buffer The buffer to be filled.
 * \param obj The object to be marshalled.
 * \return void
 */
void marshal( mpi::SendBuffer& buffer, const RigidBody& obj );
//*************************************************************************************************

//*************************************************************************************************
/*!\brief Unmarshalling rigid body parameters.
 *
 * \param buffer The buffer from where to read.
 * \param objparam The object to be reconstructed.
 * \param hasSuperBody False if body is not part of a union. Subordinate bodies in unions do not encode velocities but encode relative positions.
 * \return void
 */
void unmarshal( mpi::RecvBuffer& buffer, RigidBodyParameters& objparam );

//*************************************************************************************************
//*************************************************************************************************
//*************************************************************************************************
//*************************************************************************************************
//*************************************************************************************************
//*************************************************************************************************

struct GeomPrimitiveParameters : public RigidBodyParameters {
   MaterialID material_;
};

//*************************************************************************************************
/*!\brief Marshalling parameters of a geometric primitive.
 *
 * \param buffer The buffer to be filled.
 * \param obj The object to be marshalled.
 * \return void
 */
void marshal( mpi::SendBuffer& buffer, const GeomPrimitive& obj );
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Unmarshalling parameters of a geometric primitive.
 *
 * \param buffer The buffer to be filled.
 * \param obj The object to be marshalled.
 * \param hasSuperBody False if body is not part of a union. Passed on to rigid body unmarshalling.
 * \return void
 */
void unmarshal( mpi::RecvBuffer& buffer, GeomPrimitiveParameters& objparam );

}  // namespace communication
}  // namespace pe
}  // namespace walberla
