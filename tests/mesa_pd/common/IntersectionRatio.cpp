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
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//
//======================================================================================================================

#include "core/DataTypes.h"
#include "core/Environment.h"
#include "core/debug/TestSubsystem.h"
#include "core/math/all.h"

#include "mesa_pd/common/RayParticleIntersection.h"
#include "mesa_pd/data/ParticleAccessorWithShape.h"
#include "mesa_pd/data/ParticleStorage.h"
#include "mesa_pd/data/ShapeStorage.h"
#include "mesa_pd/data/DataTypes.h"
#include <mesa_pd/data/shape/Ellipsoid.h>
#include <mesa_pd/data/shape/HalfSpace.h>
#include <mesa_pd/data/shape/Sphere.h>
#include "mesa_pd/kernel/SingleCast.h"

namespace intersection_ratio_test
{
using namespace walberla;
using mesa_pd::Vec3;

/*!\brief Tests the ray-particle intersection ratio functionality implemented in mesa_pd/common/RayParticleIntersection.h
 *
 * Currently the following shapes are tested:
 *  - sphere
 *  - halfspace
 *  - ellipsoid (default and rotated)
 *
 * Additionally, the default variant with the bisection line search is tested with the help of a sphere.
 */

//////////
// MAIN //
//////////
int main( int argc, char **argv )
{
   debug::enterTestMode();

   mpi::Environment env( argc, argv );

   auto ps = std::make_shared<mesa_pd::data::ParticleStorage>(1);
   auto shapeStorage = std::make_shared<mesa_pd::data::ShapeStorage>();
   using ParticleAccessor = mesa_pd::data::ParticleAccessorWithShape;
   ParticleAccessor accessor(ps, shapeStorage);

   const real_t epsilon( 1e-4_r );

   mesa_pd::kernel::SingleCast singleCast;
   mesa_pd::RayParticleIntersectionRatioFunctor intersectionRatioFctr;

   ////////////
   // SPHERE //
   ////////////
   {
      real_t sphereRadius = 1_r;
      auto sphereShape = shapeStorage->create<mesa_pd::data::Sphere>( sphereRadius );

      Vec3 position(1_r, 0_r, 0_r);

      mesa_pd::data::Particle&& p = *ps->create();
      p.setPosition(position);
      p.setShapeID(sphereShape);
      auto idx = p.getIdx();

      Vec3 pos1(-0.5_r, 0_r, 0_r);
      Vec3 dir1(1_r, 0_r, 0_r);
      real_t delta1 = singleCast(idx, accessor, intersectionRatioFctr, accessor, pos1, dir1, epsilon );
      WALBERLA_CHECK_FLOAT_EQUAL(delta1, 0.5_r, "Intersection ratio 1 with sphere wrong!");

      Vec3 pos2(1_r, 1_r, 1_r);
      Vec3 dir2(0_r, -1_r, -1_r);
      real_t delta2 = singleCast(idx, accessor, intersectionRatioFctr, accessor, pos2, dir2, epsilon );
      WALBERLA_CHECK_FLOAT_EQUAL(delta2, (std::sqrt(2_r) - 1_r) / std::sqrt(2_r), "Intersection ratio 2 with sphere wrong!");
   }

   ///////////////
   // HALFSPACE //
   ///////////////
   {
      Vec3 position(1_r, 0_r, 0_r);
      Vec3 normal(0_r, 1_r, 1_r);

      auto planeShape = shapeStorage->create<mesa_pd::data::HalfSpace>( normal.getNormalized() );

      mesa_pd::data::Particle&& p = *ps->create(true);
      p.setPosition(position);
      p.setShapeID(planeShape);
      auto idx = p.getIdx();

      Vec3 pos1(1_r, 0.5_r, 0.5_r);
      Vec3 dir1(0_r, -1_r, -1_r);
      real_t delta1 = singleCast(idx, accessor, intersectionRatioFctr, accessor, pos1, dir1, epsilon );
      WALBERLA_CHECK_FLOAT_EQUAL(delta1, 0.5_r, "Intersection ratio 1 with half space wrong!");

      Vec3 dir2(0_r, 0_r, -2_r);
      real_t delta2 = singleCast(idx, accessor, intersectionRatioFctr, accessor, pos1, dir2, epsilon );
      WALBERLA_CHECK_FLOAT_EQUAL(delta2, 0.5_r, "Intersection ratio 2 with half space wrong!");

      Vec3 dir3(0_r, -3_r, 0_r);
      real_t delta3 = singleCast(idx, accessor, intersectionRatioFctr, accessor, pos1, dir3, epsilon );
      WALBERLA_CHECK_FLOAT_EQUAL(delta3, 1_r/3_r, "Intersection ratio 3 with half space wrong!");
   }

   ///////////////
   // ELLIPSOID //
   ///////////////
   {
      Vec3 semiAxes{1_r, 1_r, 2_r};
      auto ellipsoidShape = shapeStorage->create<mesa_pd::data::Ellipsoid>( semiAxes );

      Vec3 position(1_r, 0_r, 0_r);

      mesa_pd::data::Particle&& p = *ps->create();
      p.setPosition(position);
      p.setShapeID(ellipsoidShape);
      auto idx = p.getIdx();

      Vec3 pos1(-0.5_r, 0_r, 0_r);
      Vec3 dir1(1_r, 0_r, 0_r);
      real_t delta1 = singleCast(idx, accessor, intersectionRatioFctr, accessor, pos1, dir1, epsilon );
      WALBERLA_CHECK_FLOAT_EQUAL(delta1, 0.5_r, "Intersection ratio 1 with ellipsoid wrong!");

      Vec3 pos2(1_r, 1.5_r, 0_r);
      Vec3 dir2(0_r, -1_r, 0_r);
      real_t delta2 = singleCast(idx, accessor, intersectionRatioFctr, accessor, pos2, dir2, epsilon );
      WALBERLA_CHECK_FLOAT_EQUAL(delta2, 0.5_r, "Intersection ratio 2 with ellipsoid wrong!");

      auto rotation = p.getRotation();
      rotation.rotate( Vec3(1_r,0_r,0_r), math::pi / 2_r ); // rotate by 90° around x axis
      p.setRotation(rotation);

      real_t delta1Rot = singleCast(idx, accessor, intersectionRatioFctr, accessor, pos1, dir1, epsilon );
      WALBERLA_CHECK_FLOAT_EQUAL(delta1Rot, 0.5_r, "Intersection ratio 1 with rotated ellipsoid wrong!");

      Vec3 pos2Rot(1_r, 0_r, -1.5_r);
      Vec3 dir2Rot(0_r, 0_r, 1_r);
      real_t delta2Rot = singleCast(idx, accessor, intersectionRatioFctr, accessor, pos2Rot, dir2Rot, epsilon );
      WALBERLA_CHECK_FLOAT_EQUAL(delta2Rot, 0.5_r, "Intersection ratio 2 with rotated ellipsoid wrong!");

   }
   

   ///////////////////////////
   // Bisection Line Search //
   ///////////////////////////
   {
      // note: here tested with a sphere and therefore called explicitly because otherwise the sphere specialization would be selected

      real_t sphereRadius = 1_r;
      auto sphereShape = shapeStorage->create<mesa_pd::data::Sphere>( sphereRadius );

      Vec3 position(1_r, 0_r, 0_r);

      mesa_pd::data::Particle&& p = *ps->create();
      p.setPosition(position);
      p.setShapeID(sphereShape);
      auto idx = p.getIdx();

      Vec3 pos1(-0.5_r, 0_r, 0_r);
      Vec3 dir1(1_r, 0_r, 0_r);
      real_t delta1 = mesa_pd::intersectionRatioBisection(idx, accessor, pos1, dir1, epsilon);
      WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(delta1, 0.5_r, epsilon, "Intersection ratio 1 with bisection line search for sphere wrong!");

      Vec3 pos2(1_r, 1_r, 1_r);
      Vec3 dir2(0_r, -1_r, -1_r);
      real_t delta2 = mesa_pd::intersectionRatioBisection(idx, accessor, pos2, dir2, epsilon);
      WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(delta2, (std::sqrt(2_r) - 1_r) / std::sqrt(2_r), epsilon, "Intersection ratio 2 with bisection line search for sphere wrong!");
   }

   return 0;

}

} //namespace intersection_ratio_test

int main( int argc, char **argv ){
   intersection_ratio_test::main(argc, argv);
}
