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
//! \file Polyhedron.h
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//! \brief Marshalling of objects for data transmission or storage.
//
//======================================================================================================================

#pragma once


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include "mesh/pe/rigid_body/ConvexPolyhedron.h"
#include "pe/communication/Instantiate.h"
#include "pe/communication/Marshalling.h"

namespace walberla {
namespace mesh {
namespace pe {

using namespace walberla::pe::communication;

struct ConvexPolyhedronParameters : public GeomPrimitiveParameters {
   TriangleMesh mesh_;
};

//*************************************************************************************************
/*!\brief Marshalling a polyhedron primitive.
 *
 * \param buffer The buffer to be filled.
 * \param obj The object to be marshalled.
 * \return void
 */
void marshal( mpi::SendBuffer& buffer, const ConvexPolyhedron& obj );
//*************************************************************************************************

//*************************************************************************************************
/*!\brief Unmarshalling a polyhedron primitive.
 *
 * \param buffer The buffer from where to read.
 * \param objparam The object to be reconstructed.
 * \param hasSuperBody False if body is not part of a union. Passed on to rigid body unmarshalling.
 * \return void
 */
void unmarshal( mpi::RecvBuffer& buffer, ConvexPolyhedronParameters& objparam );
//*************************************************************************************************

}  // namespace pe
}  // namespace mesh


namespace pe {
namespace communication {

template<>
inline mesh::pe::ConvexPolyhedronPtr instantiate( mpi::RecvBuffer& buffer, const math::AABB& domain, const math::AABB& block, mesh::pe::ConvexPolyhedronID& newBody )
{
   mesh::pe::ConvexPolyhedronParameters subobjparam;
   unmarshal( buffer, subobjparam );
   correctBodyPosition(domain, block.center(), subobjparam.gpos_);
   auto cp = std::make_unique<mesh::pe::ConvexPolyhedron>( subobjparam.sid_, subobjparam.uid_, subobjparam.gpos_, subobjparam.q_, subobjparam.mesh_, subobjparam.material_, false, subobjparam.communicating_, subobjparam.infiniteMass_ );
   cp->setLinearVel( subobjparam.v_ );
   cp->setAngularVel( subobjparam.w_ );
   cp->MPITrait.setOwner( subobjparam.mpiTrait_.owner_ );
   newBody = static_cast<mesh::pe::ConvexPolyhedronID>(cp.get());
   return cp;
}

}  // namespace communication
}  // namespace pe
}  // namespace walberla
