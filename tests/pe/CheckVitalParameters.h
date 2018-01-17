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
//! \file CheckVitalParameters.h
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#pragma once

#include "pe/rigidbody/RigidBody.h"
#include "pe/rigidbody/Sphere.h"

#include "core/debug/TestSubsystem.h"

namespace walberla {
namespace pe {

template <typename RigidBodyID>
void checkVitalParameters(RigidBodyID /*bd1*/, RigidBodyID /*bd2*/) { WALBERLA_CHECK_EQUAL(1,2);}

template <>
inline void checkVitalParameters(SphereID bd1, SphereID bd2)
{
   WALBERLA_CHECK_FLOAT_EQUAL(bd1->getRadius(), bd2->getRadius());

   WALBERLA_CHECK_FLOAT_EQUAL(bd1->getPosition(), bd2->getPosition());
   WALBERLA_CHECK_FLOAT_EQUAL(bd1->getLinearVel(), bd2->getLinearVel());
   WALBERLA_CHECK_FLOAT_EQUAL(bd1->getRotation(), bd2->getRotation());
   WALBERLA_CHECK_FLOAT_EQUAL(bd1->getAngularVel(), bd2->getAngularVel());

   if (std::isinf(bd1->getMass()))
   {
      WALBERLA_CHECK            ( std::isinf(bd2->getMass()) );
      WALBERLA_CHECK_FLOAT_EQUAL(bd1->getInvMass(), bd2->getInvMass());
      WALBERLA_CHECK            ( math::isinf(bd1->getBodyInertia()) );
      WALBERLA_CHECK            ( math::isinf(bd2->getBodyInertia()) );
      WALBERLA_CHECK_FLOAT_EQUAL(bd1->getInvInertia(), bd2->getInvInertia());
   } else {
      WALBERLA_CHECK_FLOAT_EQUAL(bd1->getMass(), bd2->getMass());
      WALBERLA_CHECK_FLOAT_EQUAL(bd1->getInvMass(), bd2->getInvMass());
      WALBERLA_CHECK_FLOAT_EQUAL(bd1->getInertia(), bd2->getInertia());
      WALBERLA_CHECK_FLOAT_EQUAL(bd1->getInvInertia(), bd2->getInvInertia());
   }

   WALBERLA_CHECK_EQUAL(bd1->isGlobal(), bd2->isGlobal());
   WALBERLA_CHECK_EQUAL(bd1->hasInfiniteMass(), bd2->hasInfiniteMass());
   WALBERLA_CHECK_EQUAL(bd1->isCommunicating(), bd2->isCommunicating());
   WALBERLA_CHECK_EQUAL(bd1->isFinite(), bd2->isFinite());
}

}
}
