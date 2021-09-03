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
//! \file Ellipsoid.h
//! \author Tobias Preclik
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//! \brief Marshalling of objects for data transmission or storage.
//
//======================================================================================================================

#pragma once


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include "pe/communication/Instantiate.h"
#include "pe/communication/Marshalling.h"
#include "pe/rigidbody/Ellipsoid.h"

namespace walberla {
namespace pe {
namespace communication {

struct EllipsoidParameters : public GeomPrimitiveParameters {
   Vec3 semiAxes_;
};

//*************************************************************************************************
/*!\brief Marshalling a Ellipsoid primitive.
 *
 * \param buffer The buffer to be filled.
 * \param obj The object to be marshalled.
 * \return void
 */
void marshal( mpi::SendBuffer& buffer, const Ellipsoid& obj );
//*************************************************************************************************

//*************************************************************************************************
/*!\brief Unmarshalling a Ellipsoid primitive.
 *
 * \param buffer The buffer from where to read.
 * \param objparam The object to be reconstructed.
 * \return void
 */
void unmarshal( mpi::RecvBuffer& buffer, EllipsoidParameters& objparam );
//*************************************************************************************************


inline EllipsoidPtr instantiate( mpi::RecvBuffer& buffer, const math::AABB& domain, const math::AABB& block, EllipsoidID& newBody )
{
   EllipsoidParameters subobjparam;
   unmarshal( buffer, subobjparam );
   correctBodyPosition(domain, block.center(), subobjparam.gpos_);
   auto el = std::make_unique<Ellipsoid>( subobjparam.sid_, subobjparam.uid_, subobjparam.gpos_,  subobjparam.q_, subobjparam.semiAxes_, subobjparam.material_, false, subobjparam.communicating_, subobjparam.infiniteMass_ );
   el->setLinearVel( subobjparam.v_ );
   el->setAngularVel( subobjparam.w_ );
   el->MPITrait.setOwner( subobjparam.mpiTrait_.owner_ );
   newBody = el.get();
   return el;
}

}  // namespace communication
}  // namespace pe
}  // namespace walberla
