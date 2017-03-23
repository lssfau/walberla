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
//! \file ParseMessage.cpp
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include "pe/communication/ParseMessage.h"

#include "core/math/AABB.h"

#include "core/debug/TestSubsystem.h"

using namespace walberla;
using namespace walberla::pe;

void correctBodyPositionTest()
{
   //********************************************************************************************************
   //ParseMessage.h correctBodyPosition Test
   //********************************************************************************************************
   math::AABB domain( Vec3(0,0,0), Vec3(3,3,3) );
   Vec3 center(0.5, 0.5, 0.5);
   Vec3 pos;

   pos = Vec3(1.5, 0.0, 0.0);
   pe::communication::correctBodyPosition(domain, center, pos);
   WALBERLA_CHECK_FLOAT_EQUAL( Vec3(1.5, 0.0, 0.0) , pos);

   pos = Vec3(1.5, 1.5, 0.0);
   pe::communication::correctBodyPosition(domain, center, pos);
   WALBERLA_CHECK_FLOAT_EQUAL( Vec3(1.5, 1.5, 0.0) , pos);

   pos = Vec3(0.0, 1.5, 0.0);
   pe::communication::correctBodyPosition(domain, center, pos);
   WALBERLA_CHECK_FLOAT_EQUAL( Vec3(0.0, 1.5, 0.0) , pos);

   pos = Vec3(0.0, 0.0, 1.5);
   pe::communication::correctBodyPosition(domain, center, pos);
   WALBERLA_CHECK_FLOAT_EQUAL( Vec3(0.0, 0.0, 1.5) , pos);

   pos = Vec3(1.5, 0.0, 1.5);
   pe::communication::correctBodyPosition(domain, center, pos);
   WALBERLA_CHECK_FLOAT_EQUAL( Vec3(1.5, 0.0, 1.5) , pos);

   pos = Vec3(1.5, 1.5, 1.5);
   pe::communication::correctBodyPosition(domain, center, pos);
   WALBERLA_CHECK_FLOAT_EQUAL( Vec3(1.5, 1.5, 1.5) , pos);

   pos = Vec3(0.0, 1.5, 1.5);
   pe::communication::correctBodyPosition(domain, center, pos);
   WALBERLA_CHECK_FLOAT_EQUAL( Vec3(0.0, 1.5, 1.5) , pos);

   // periodic tests

   pos = Vec3(2.5, 0.5, 0.5);
   pe::communication::correctBodyPosition(domain, center, pos);
   WALBERLA_CHECK_FLOAT_EQUAL( Vec3(-0.5, 0.5, 0.5) , pos);

   pos = Vec3(2.5, 2.5, 0.5);
   pe::communication::correctBodyPosition(domain, center, pos);
   WALBERLA_CHECK_FLOAT_EQUAL( Vec3(-0.5, -0.5, 0.5) , pos);

   pos = Vec3(0.5, 2.5, 0.5);
   pe::communication::correctBodyPosition(domain, center, pos);
   WALBERLA_CHECK_FLOAT_EQUAL( Vec3(0.5, -0.5, 0.5) , pos);

   pos = Vec3(0.5, 0.5, 2.5);
   pe::communication::correctBodyPosition(domain, center, pos);
   WALBERLA_CHECK_FLOAT_EQUAL( Vec3(0.5, 0.5, -0.5) , pos);

   pos = Vec3(2.5, 0.5, 2.5);
   pe::communication::correctBodyPosition(domain, center, pos);
   WALBERLA_CHECK_FLOAT_EQUAL( Vec3(-0.5, 0.5, -0.5) , pos);

   pos = Vec3(2.5, 2.5, 2.5);
   pe::communication::correctBodyPosition(domain, center, pos);
   WALBERLA_CHECK_FLOAT_EQUAL( Vec3(-0.5, -0.5, -0.5) , pos);

   pos = Vec3(0.5, 2.5, 2.5);
   pe::communication::correctBodyPosition(domain, center, pos);
   WALBERLA_CHECK_FLOAT_EQUAL( Vec3(0.5, -0.5, -0.5) , pos);

   // upper corner

   center = Vec3(2.5, 2.5, 2.5);

   pos = Vec3(1.5, 2.5, 2.5);
   pe::communication::correctBodyPosition(domain, center, pos);
   WALBERLA_CHECK_FLOAT_EQUAL( Vec3(1.5, 2.5, 2.5) , pos);

   pos = Vec3(1.5, 1.5, 2.5);
   pe::communication::correctBodyPosition(domain, center, pos);
   WALBERLA_CHECK_FLOAT_EQUAL( Vec3(1.5, 1.5, 2.5) , pos);

   pos = Vec3(2.5, 1.5, 2.5);
   pe::communication::correctBodyPosition(domain, center, pos);
   WALBERLA_CHECK_FLOAT_EQUAL( Vec3(2.5, 1.5, 2.5) , pos);

   pos = Vec3(2.5, 2.5, 1.5);
   pe::communication::correctBodyPosition(domain, center, pos);
   WALBERLA_CHECK_FLOAT_EQUAL( Vec3(2.5, 2.5, 1.5) , pos);

   pos = Vec3(1.5, 2.5, 1.5);
   pe::communication::correctBodyPosition(domain, center, pos);
   WALBERLA_CHECK_FLOAT_EQUAL( Vec3(1.5, 2.5, 1.5) , pos);

   pos = Vec3(1.5, 1.5, 1.5);
   pe::communication::correctBodyPosition(domain, center, pos);
   WALBERLA_CHECK_FLOAT_EQUAL( Vec3(1.5, 1.5, 1.5) , pos);

   pos = Vec3(2.5, 1.5, 1.5);
   pe::communication::correctBodyPosition(domain, center, pos);
   WALBERLA_CHECK_FLOAT_EQUAL( Vec3(2.5, 1.5, 1.5) , pos);

   // periodic tests

   pos = Vec3(0.5, 2.5, 2.5);
   pe::communication::correctBodyPosition(domain, center, pos);
   WALBERLA_CHECK_FLOAT_EQUAL( Vec3(3.5, 2.5, 2.5) , pos);

   pos = Vec3(0.5, 0.5, 2.5);
   pe::communication::correctBodyPosition(domain, center, pos);
   WALBERLA_CHECK_FLOAT_EQUAL( Vec3(3.5, 3.5, 2.5) , pos);

   pos = Vec3(2.5, 0.5, 2.5);
   pe::communication::correctBodyPosition(domain, center, pos);
   WALBERLA_CHECK_FLOAT_EQUAL( Vec3(2.5, 3.5, 2.5) , pos);

   pos = Vec3(2.5, 2.5, 0.5);
   pe::communication::correctBodyPosition(domain, center, pos);
   WALBERLA_CHECK_FLOAT_EQUAL( Vec3(2.5, 2.5, 3.5) , pos);

   pos = Vec3(0.5, 2.5, 0.5);
   pe::communication::correctBodyPosition(domain, center, pos);
   WALBERLA_CHECK_FLOAT_EQUAL( Vec3(3.5, 2.5, 3.5) , pos);

   pos = Vec3(0.5, 0.5, 0.5);
   pe::communication::correctBodyPosition(domain, center, pos);
   WALBERLA_CHECK_FLOAT_EQUAL( Vec3(3.5, 3.5, 3.5) , pos);

   pos = Vec3(2.5, 0.5, 0.5);
   pe::communication::correctBodyPosition(domain, center, pos);
   WALBERLA_CHECK_FLOAT_EQUAL( Vec3(2.5, 3.5, 3.5) , pos);

}

int main( int /*argc*/, char ** /*argv*/ )
{
   walberla::debug::enterTestMode();

   correctBodyPositionTest();

   return EXIT_SUCCESS;
}
