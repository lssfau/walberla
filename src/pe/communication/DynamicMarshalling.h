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
#include "pe/communication/rigidbody/Squirmer.h"
#include "pe/communication/rigidbody/Union.h"

#include "core/Abort.h"

#include <boost/tuple/tuple.hpp>


namespace walberla {
namespace pe {
namespace communication {

template < typename BodyTypeTuple >
struct MarshalDynamically{
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
      static_assert(boost::is_base_of<RigidBody, typename BodyTypeTuple::head_type>::value, "only derived types from RigidBody are allowed!");

      if (BodyTypeTuple::head_type::getStaticTypeID() == b.getTypeID())
      {
         typedef typename BodyTypeTuple::head_type BodyT;
         marshal( buffer, static_cast<const BodyT&>(b));
      } else
      {
        MarshalDynamically<typename BodyTypeTuple::tail_type>::execute(buffer, b);
      }
   }
};

template < >
struct MarshalDynamically< boost::tuples::null_type>{
    static void execute(mpi::SendBuffer& /*buffer*/, const RigidBody& b)
    {
       WALBERLA_ABORT("MarshalDynamically: BodyTypeID " << b.getTypeID() << " could not be resolved!");
    }
};

template < typename BodyTypeTuple >
struct UnmarshalDynamically{
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
   static BodyID execute(mpi::RecvBuffer& buffer, const id_t typeID, const math::AABB& domain, const math::AABB& block)
   {
      if (BodyTypeTuple::head_type::getStaticTypeID() == typeID)
      {
         typedef typename BodyTypeTuple::head_type BodyT;
         BodyT* newBody;
         return instantiate( buffer, domain, block, newBody );
      } else
      {
         return UnmarshalDynamically<typename BodyTypeTuple::tail_type>::execute(buffer, typeID, domain, block);
      }
   }
};

template < >
struct UnmarshalDynamically< boost::tuples::null_type>{
    static BodyID execute(mpi::RecvBuffer& /*buffer*/, const id_t typeID, const math::AABB& /*domain*/, const math::AABB& /*block*/)
    {
       WALBERLA_ABORT("UnmarshalDynamically: BodyTypeID " << typeID << " could not be resolved!");
    }
};

}  // namespace communication
}  // namespace pe
}  // namespace walberla
