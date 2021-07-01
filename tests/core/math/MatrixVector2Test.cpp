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
//! \file MatrixVector2Test.cpp
//! \ingroup core
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#include "core/debug/TestSubsystem.h"
#include "core/math/Matrix2.h"
#include "core/math/Vector2.h"

#include <iostream>

using namespace walberla;
using walberla::uint8_t;

void matrix2Test()
{
   Matrix2< real_t > m1(2, 1, -1, 5);
   Matrix2< real_t > m2(4, 2, -2, 1);

   WALBERLA_CHECK_EQUAL(m1 * m2 * m2.getInverse(), m1);
   WALBERLA_CHECK_EQUAL(m1 * m1.getInverse() * m2, m2);

   Vector2< real_t > v1(5, 7);

   WALBERLA_CHECK_EQUAL(m2 * v1, Vector2< real_t >(34, -3));

   const real_t s = real_c(5);

   WALBERLA_CHECK_EQUAL(s * m1, Matrix2< real_t >(10, 5, -5, 25));
   WALBERLA_CHECK_EQUAL(m1 * s, Matrix2< real_t >(10, 5, -5, 25));
}

void vector2Test()
{
   Vector2< real_t > v1(1, 2);
   Vector2< real_t > v2(3, 4);
   Vector2< uint_t > v3(4, 3);

   auto sum = v1 + v2;
   WALBERLA_CHECK_EQUAL(sum, Vector2< real_t >(4, 6));

   auto product = v1 * v2;
   WALBERLA_CHECK_FLOAT_EQUAL(product, 1.0 * 3 + 2.0 * 4);

   WALBERLA_CHECK_FLOAT_EQUAL(v2.length(), 5.0);
   WALBERLA_CHECK_FLOAT_EQUAL(v2.sqrLength(), 25.0);

   WALBERLA_CHECK_EQUAL(v3.indexOfMin(), 1);
   WALBERLA_CHECK_EQUAL(v3.indexOfMax(), 0);
}

int main()
{
   vector2Test();
   matrix2Test();

   return 0;
}
