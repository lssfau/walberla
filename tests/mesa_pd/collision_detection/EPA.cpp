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
//! \file
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include <mesa_pd/collision_detection/AnalyticCollisionFunctions.h>
#include <mesa_pd/collision_detection/EPA.h>
#include <mesa_pd/collision_detection/GJK.h>
#include <mesa_pd/collision_detection/Support.h>
#include <mesa_pd/data/shape/Sphere.h>

#include <core/Environment.h>
#include <core/logging/Logging.h>

#include <cmath>
#include <iostream>
#include <memory>

namespace walberla {
namespace mesa_pd {

void check( )
{
   using namespace walberla::mesa_pd::collision_detection;

   auto pos0    = Vec3(real_t(1),real_t(2),real_t(3));
   auto radius0 = real_t(0.6);
   auto pos1    = Vec3(real_t(2),real_t(2),real_t(3));
   auto radius1 = real_t(0.6);

   auto sp0 = data::Sphere(real_t(0.6));
   auto geom0 = Support(pos0, Rot3(), sp0);
   auto geom1 = Support(pos1, Rot3(), sp0);

   Vec3   epa_normal;
   Vec3   epa_contactPoint;
   real_t epa_penetrationDepth;

   real_t margin = real_t(1e-4);
   GJK gjk;
   WALBERLA_CHECK(gjk.doGJKmargin(geom0, geom1, margin));
   EPA epa;
   epa.useSphereOptimization(true);
   WALBERLA_CHECK(epa.doEPAmargin(geom0, geom1, gjk, epa_normal, epa_contactPoint, epa_penetrationDepth, margin));

   Vec3   analytic_normal;
   Vec3   analytic_contactPoint;
   real_t analytic_penetrationDepth;
   WALBERLA_CHECK(analytic::detectSphereSphereCollision(pos0,
                                                        radius0,
                                                        pos1,
                                                        radius1,
                                                        analytic_contactPoint,
                                                        analytic_normal,
                                                        analytic_penetrationDepth,
                                                        margin));

   WALBERLA_CHECK_FLOAT_EQUAL(epa_contactPoint, analytic_contactPoint);
   WALBERLA_CHECK_FLOAT_EQUAL(epa_normal, analytic_normal);
   WALBERLA_CHECK_FLOAT_EQUAL(epa_penetrationDepth, analytic_penetrationDepth);
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
