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
#include "pe/communication/Marshalling.h"
#include "pe/communication/rigidbody/Box.h"
#include "pe/communication/rigidbody/Capsule.h"
#include "pe/communication/rigidbody/Sphere.h"
#include "pe/communication/rigidbody/Union.h"
#include "pe/communication/rigidbody/Ellipsoid.h"
#include "pe/utility/BodyCast.h"

#include "core/Abort.h"


namespace walberla {
namespace pe {
namespace communication {

template < typename BodyTypeTuple >
class MarshalDynamically{
private:
   struct MarshalFunctor
   {
      mpi::SendBuffer& buffer_;

      MarshalFunctor(mpi::SendBuffer& buffer) : buffer_(buffer) {}

      template< typename BodyType >
      void operator()( const BodyType* bd) { marshal( buffer_, *bd); }
   };

public:
   //*************************************************************************************************
   /*!\brief Marshalling rigid body parameters dynamically.
    *
    * \param buffer The buffer to be filled.
    * \param obj The object to be marshalled dynamically.
    * \return void
    *
    * The rigid body is casted dynamically to its original type and then marshalled. For recognition
    * an identifying tag is prepended.
    */
   static void execute(mpi::SendBuffer& buffer, const RigidBody& b)
   {
      MarshalFunctor func(buffer);
      return SingleCast<BodyTypeTuple, MarshalFunctor, void>::execute (&b, func);
   }
};

template < typename BodyTypeTuple >
class UnmarshalDynamically{
private:
   struct UnmarshalFunctor
   {
      mpi::RecvBuffer& buffer_;
      const math::AABB& domain_;
      const math::AABB& block_;

      UnmarshalFunctor(mpi::RecvBuffer& buffer, const math::AABB& domain, const math::AABB& block)
         : buffer_(buffer)
         , domain_(domain)
         , block_(block) {}

      template< typename BodyType >
      BodyPtr operator()( BodyType* bd) { return instantiate( buffer_, domain_, block_, bd ); }
   };

public:
   //*************************************************************************************************
   /*!\brief Marshalling rigid body parameters dynamically.
    *
    * \param buffer The buffer to be filled.
    * \param obj The object to be marshalled dynamically.
    * \return void
    *
    * The rigid body is casted dynamically to its original type and then marshalled. For recognition
    * an identifying tag is prepended.
    */
   static BodyPtr execute(mpi::RecvBuffer& buffer, const id_t typeID, const math::AABB& domain, const math::AABB& block)
   {
      UnmarshalFunctor func(buffer, domain, block);
      return SingleCast<BodyTypeTuple, UnmarshalFunctor, BodyPtr>::execute (typeID, func);
   }
};

}  // namespace communication
}  // namespace pe
}  // namespace walberla
