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
//! \file   AnalyticCollisionFunctions.cpp
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include <mesa_pd/collision_detection/AnalyticCollisionFunctions.h>
#include <mesa_pd/data/DataTypes.h>

#include <core/Environment.h>
#include <core/logging/Logging.h>

#include <cmath>
#include <iostream>

namespace walberla {
namespace mesa_pd {

void checkSphereSphereCollision( )
{
   using namespace walberla::mesa_pd::collision_detection::analytic;

   const real_t radius  = real_t(0.5);
   auto dir = Vec3(real_t(1), real_t(2), real_t(3)).getNormalized();
   real_t shift = real_c(0.75);

   Vec3   contactPoint;
   Vec3   contactNormal;
   real_t penetrationDepth;

   //check two spheres in contact
   WALBERLA_CHECK( detectSphereSphereCollision( Vec3(real_t(0), real_t(0), real_t(0)),
                                                radius,
                                                dir * shift,
                                                radius,
                                                contactPoint,
                                                contactNormal,
                                                penetrationDepth,
                                                real_t(0) ) );
   WALBERLA_CHECK_FLOAT_EQUAL( contactPoint, dir * shift * real_t(0.5));
   WALBERLA_CHECK_FLOAT_EQUAL( contactNormal, -dir );
   WALBERLA_CHECK_FLOAT_EQUAL( penetrationDepth, shift - real_t(1) );

   //no collision
   WALBERLA_CHECK( !detectSphereSphereCollision( Vec3(real_t(0), real_t(0), real_t(0)),
                                                 radius,
                                                 dir * real_t(1.1),
                                                 radius,
                                                 contactPoint,
                                                 contactNormal,
                                                 penetrationDepth,
                                                 real_t(0) ) );
}

void checkSphereCylindricalBoundaryCollision( )
{
   using namespace walberla::mesa_pd::collision_detection::analytic;

   const real_t radius  = real_t(0.5);
   auto pos = Vec3(real_t(6), real_t(-2), real_t(2));
   auto dir = Vec3(real_t(1), real_t(1), real_t(1)).getNormalized();

   Vec3   contactPoint;
   Vec3   contactNormal;
   real_t penetrationDepth;

   //check two spheres in contact
   WALBERLA_CHECK( detectSphereCylindricalBoundaryCollision( pos,
                                                             radius,
                                                             Vec3(real_t(0), real_t(0), real_t(0)),
                                                             real_t(5),
                                                             dir,
                                                             contactPoint,
                                                             contactNormal,
                                                             penetrationDepth,
                                                             real_t(0) ) );
   auto nearestPoint = dot(pos, dir) * dir;
   WALBERLA_CHECK_FLOAT_EQUAL( contactPoint, (pos - nearestPoint).getNormalized() * real_t(5) + nearestPoint);
   WALBERLA_CHECK_FLOAT_EQUAL( contactNormal, (nearestPoint - pos).getNormalized() );
   WALBERLA_CHECK_FLOAT_EQUAL( penetrationDepth, -((pos - nearestPoint).length() + radius - real_t(5)) );

   //no collision
   WALBERLA_CHECK( !detectSphereCylindricalBoundaryCollision( Vec3(real_t(0), real_t(0), real_t(0)),
                                                              radius,
                                                              Vec3(real_t(0), real_t(0), real_t(0)),
                                                              real_t(5),
                                                              Vec3(real_t(1), real_t(1), real_t(1)).getNormalized(),
                                                              contactPoint,
                                                              contactNormal,
                                                              penetrationDepth,
                                                              real_t(0) ) );
}

void checkSphereBoxCollision( )
{
   using namespace walberla::mesa_pd::collision_detection::analytic;

   Vec3   contactPoint;
   Vec3   contactNormal;
   real_t penetrationDepth;

   //check two spheres in contact
   WALBERLA_CHECK( detectSphereBoxCollision( Vec3(real_t(2), real_t(2), real_t(2)),
                                             real_c(std::sqrt(3)) + real_t(0.1),
                                             Vec3(real_t(0), real_t(0), real_t(0)),
                                             Vec3(real_t(2), real_t(2), real_t(2)),
                                             Rot3(),
                                             contactPoint,
                                             contactNormal,
                                             penetrationDepth,
                                             real_t(0) ) );
   WALBERLA_CHECK_FLOAT_EQUAL( contactPoint, Vec3(real_t(1), real_t(1), real_t(1)));
   WALBERLA_CHECK_FLOAT_EQUAL( contactNormal, Vec3(real_t(2), real_t(2), real_t(2)).getNormalized() );
   WALBERLA_CHECK_FLOAT_EQUAL( penetrationDepth, -real_t(0.1) );
}

void checkSphereHalfSpaceCollision( )
{
   using namespace walberla::mesa_pd::collision_detection::analytic;

   const real_t radius  = real_t(1.0);
   auto dir = Vec3(real_t(1), real_t(2), real_t(3)).getNormalized();
   real_t shift = real_c(0.75);

   Vec3   contactPoint;
   Vec3   contactNormal;
   real_t penetrationDepth;

   //check sphere - half space contact
   WALBERLA_CHECK( detectSphereHalfSpaceCollision( Vec3(real_t(0), real_t(0), real_t(0)),
                                                   radius,
                                                   dir * shift,
                                                   -dir,
                                                   contactPoint,
                                                   contactNormal,
                                                   penetrationDepth,
                                                   real_t(0) ) );
   WALBERLA_CHECK_FLOAT_EQUAL( contactPoint, dir * shift * real_t(1.0));
   WALBERLA_CHECK_FLOAT_EQUAL( contactNormal, -dir );
   WALBERLA_CHECK_FLOAT_EQUAL( penetrationDepth, shift - real_t(1) );


   auto pos = Vec3(shift, real_t(0), real_t(0));
   WALBERLA_CHECK( detectSphereHalfSpaceCollision( Vec3(real_t(0), real_t(0), real_t(0)),
                                                   radius,
                                                   pos,
                                                   -Vec3(real_t(1), real_t(0), real_t(0)),
                                                   contactPoint,
                                                   contactNormal,
                                                   penetrationDepth,
                                                   real_t(0) ) );
   WALBERLA_CHECK_FLOAT_EQUAL( contactPoint, pos);
   WALBERLA_CHECK_FLOAT_EQUAL( contactNormal, -Vec3(real_t(1), real_t(0), real_t(0)) );
   WALBERLA_CHECK_FLOAT_EQUAL( penetrationDepth, pos[0] - real_t(1) );
}

} //namespace mesa_pd
} //namespace walberla

int main( int argc, char ** argv )
{
   walberla::Environment env(argc, argv);
   WALBERLA_UNUSED(env);
   walberla::mpi::MPIManager::instance()->useWorldComm();

   walberla::mesa_pd::checkSphereSphereCollision();
   walberla::mesa_pd::checkSphereCylindricalBoundaryCollision();
   walberla::mesa_pd::checkSphereBoxCollision();
   walberla::mesa_pd::checkSphereHalfSpaceCollision();

   return EXIT_SUCCESS;
}
