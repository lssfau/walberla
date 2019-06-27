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
//! \file   EllipsoidSupport.cpp
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include <mesa_pd/collision_detection/Support.h>
#include <mesa_pd/data/DataTypes.h>
#include <mesa_pd/data/shape/Ellipsoid.h>

#include <core/Environment.h>
#include <core/logging/Logging.h>

#include <cmath>
#include <iostream>

namespace walberla {
namespace mesa_pd {

inline
real_t contour(const Vec3& point, const Vec3& semiAxes)
{
   return (point[0] * point[0]) / (semiAxes[0] * semiAxes[0]) +
          (point[1] * point[1]) / (semiAxes[1] * semiAxes[1]) +
          (point[2] * point[2]) / (semiAxes[2] * semiAxes[2]);
}

void check( )
{
   using namespace walberla::mesa_pd::collision_detection;

   Vec3 pos      = Vec3(real_t(1),real_t(2),real_t(3));
   Vec3 semiAxes = Vec3(real_t(2),real_t(3),real_t(1));

   auto el = data::Ellipsoid(semiAxes);
   Support e0(pos, Rot3(), el);
   WALBERLA_CHECK_FLOAT_EQUAL(e0.support(Vec3(real_t(1),real_t(0),real_t(0))), Vec3(real_t(3),real_t(2),real_t(3)));
   WALBERLA_CHECK_FLOAT_EQUAL(e0.support(Vec3(real_t(0),real_t(1),real_t(0))), Vec3(real_t(1),real_t(5),real_t(3)));
   WALBERLA_CHECK_FLOAT_EQUAL(e0.support(Vec3(real_t(0),real_t(0),real_t(1))), Vec3(real_t(1),real_t(2),real_t(4)));
   WALBERLA_CHECK_FLOAT_EQUAL(e0.support(Vec3(real_t(-1),real_t(0),real_t(0))), Vec3(real_t(-1),real_t(2),real_t(3)));
   WALBERLA_CHECK_FLOAT_EQUAL(e0.support(Vec3(real_t(0),real_t(-1),real_t(0))), Vec3(real_t(1),real_t(-1),real_t(3)));
   WALBERLA_CHECK_FLOAT_EQUAL(e0.support(Vec3(real_t(0),real_t(0),real_t(-1))), Vec3(real_t(1),real_t(2),real_t(2)));

   WALBERLA_CHECK_FLOAT_EQUAL( contour(e0.support(Vec3(real_t(-1),real_t(-2),real_t(-3)).getNormalized()) - pos, semiAxes), real_t(1) );
   WALBERLA_CHECK_FLOAT_EQUAL( contour(e0.support(Vec3(real_t(-2),real_t(-3),real_t(-1)).getNormalized()) - pos, semiAxes), real_t(1) );
   WALBERLA_CHECK_FLOAT_EQUAL( contour(e0.support(Vec3(real_t(-3),real_t(-1),real_t(-2)).getNormalized()) - pos, semiAxes), real_t(1) );

   Rot3 rot = Rot3(Vec3(1,3,2).getNormalized(), real_t(1.56));
   Support e1(pos, rot, el);
   WALBERLA_CHECK_FLOAT_EQUAL( contour( rot.getMatrix().getTranspose() * (e1.support(Vec3(real_t(-1),real_t(-2),real_t(-3)).getNormalized()) - pos), semiAxes), real_t(1) );
   WALBERLA_CHECK_FLOAT_EQUAL( contour( rot.getMatrix().getTranspose() * (e1.support(Vec3(real_t(-2),real_t(-3),real_t(-1)).getNormalized()) - pos), semiAxes), real_t(1) );
   WALBERLA_CHECK_FLOAT_EQUAL( contour( rot.getMatrix().getTranspose() * (e1.support(Vec3(real_t(-3),real_t(-1),real_t(-2)).getNormalized()) - pos), semiAxes), real_t(1) );
}

} //namespace mesa_pd
} //namespace walberla

int main( int argc, char ** argv )
{
   walberla::Environment env(argc, argv);
   WALBERLA_UNUSED(env);
   walberla::mpi::MPIManager::instance()->useWorldComm();

   walberla::mesa_pd::check();

   return EXIT_SUCCESS;
}
