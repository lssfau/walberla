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
//! \file BodyData.cpp
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include "BodyData.h"

#include <pe/rigidbody/RigidBody.h>

namespace walberla {
namespace pe {
namespace debug {

//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

BodyData::BodyData() = default;

BodyData::BodyData(ConstBodyID rb)
   : uid(rb->getID())
   , sid(rb->getSystemID())
   , pos(rb->getPosition())
   , rot(rb->getQuaternion())
   , v(rb->getLinearVel())
   , w(rb->getAngularVel())
{

}

bool checkEqual(const BodyData& bd1, const BodyData& bd2)
{
   //std::cout << std::setprecision(20) << bd1.w << "\n" << bd2.w;
   WALBERLA_CHECK_EQUAL(bd1.uid, bd2.uid);
   WALBERLA_CHECK_FLOAT_EQUAL(bd1.pos, bd2.pos);
   WALBERLA_CHECK_FLOAT_EQUAL(bd1.rot, bd2.rot);
   WALBERLA_CHECK_FLOAT_EQUAL(bd1.v, bd2.v);
   WALBERLA_CHECK_FLOAT_EQUAL(bd1.w, bd2.w, std::setprecision(20) << bd1.w << "\n" << bd2.w);

   return true;
}

} // namespace debug
} // namespace pe
} // namespace walberla
