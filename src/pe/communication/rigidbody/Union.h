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
//! \file Union.h
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
#include "pe/rigidbody/Union.h"

namespace walberla {
namespace pe {
namespace communication {

//forward declaration
template < typename BodyTypeTuple >
class MarshalDynamically;

//forward declaration
template < typename BodyTypeTuple >
class UnmarshalDynamically;

struct UnionParameters : public GeomPrimitiveParameters {
   real_t         m_;
   Mat3           I_;
   math::AABB     aabb_;
   size_t         size_;
};

//*************************************************************************************************
/*!\brief Marshalling a box primitive.
 *
 * \param buffer The buffer to be filled.
 * \param obj The object to be marshalled.
 * \return void
 */
template < typename BodyTypeTuple >
void marshal( mpi::SendBuffer& buffer, const Union<BodyTypeTuple>& obj )
{
   // Material of union is not relevant for reconstruction thus marshal RigidBody instead of GeomPrimitive.
   marshal( buffer, static_cast<const RigidBody&>( obj ) );

   buffer << obj.getMass();               // Encoding the total mass
   buffer << obj.getBodyInertia();        // Encoding the moment of inertia

   // Preparing the axis-aligned bounding box
   buffer << obj.getAABB();               // Encoding the axis-aligned bounding box

   buffer << static_cast<size_t> (obj.size());                  // Encoding the number of contained bodies

   // Encoding the contained primitives
   const typename Union<BodyTypeTuple>::const_iterator begin( obj.begin() );
   const typename Union<BodyTypeTuple>::const_iterator end  ( obj.end()   );
   for(  typename Union<BodyTypeTuple>::const_iterator body=begin; body!=end; ++body )
   {
      buffer << body->getTypeID();
      MarshalDynamically<BodyTypeTuple>::execute( buffer, *body );
   }
}

//*************************************************************************************************

//*************************************************************************************************
/*!\brief Unmarshalling a box primitive.
 *
 * \param buffer The buffer from where to read.
 * \param objparam The object to be reconstructed.
 * \param hasSuperBody False if body is not part of a union. Passed on to rigid body unmarshalling.
 * \return void
 */
inline
void unmarshal( mpi::RecvBuffer& buffer, UnionParameters& objparam )
{
   // Material of union is not relevant for reconstruction thus marshal RigidBody instead of GeomPrimitive.
   unmarshal( buffer, static_cast<RigidBodyParameters&>( objparam ) );

   buffer >> objparam.m_;
   buffer >> objparam.I_;
   buffer >> objparam.aabb_;
   //std::cout << "mass of union: " << objparam.m_ << std::endl;

   // Decode the number of contained bodies
   buffer >> objparam.size_;
   //std::cout << "size of union: " << size << std::endl;
}
//*************************************************************************************************


template <typename BodyTypeTuple>
inline std::unique_ptr<Union<BodyTypeTuple>> instantiate( mpi::RecvBuffer& buffer, const math::AABB& domain, const math::AABB& block, Union<BodyTypeTuple>*& newBody )
{
   UnionParameters subobjparam;
   unmarshal( buffer, subobjparam );
   correctBodyPosition(domain, block.center(), subobjparam.gpos_);
   auto un = std::make_unique<Union<BodyTypeTuple>>( subobjparam.sid_,
                                                     subobjparam.uid_,
                                                     subobjparam.gpos_,
                                                     subobjparam.rpos_,
                                                     subobjparam.q_,
                                                     false,
                                                     subobjparam.communicating_,
                                                     subobjparam.infiniteMass_ );

   un->MPITrait.setOwner( subobjparam.mpiTrait_.owner_ );

   // Decoding the contained primitives
   for( size_t i = 0; i < subobjparam.size_; ++i )
   {
      decltype ( static_cast<BodyID>(nullptr)->getTypeID() ) type;
      buffer >> type;
      BodyPtr obj( UnmarshalDynamically<BodyTypeTuple>::execute(buffer, type, domain, block) );
      obj->setRemote( true );
      un->add(std::move(obj));
   }
   un->setLinearVel( subobjparam.v_ );
   un->setAngularVel( subobjparam.w_ );
   newBody = un.get();
   return un;
}

}  // namespace communication
}  // namespace pe
}  // namespace walberla
