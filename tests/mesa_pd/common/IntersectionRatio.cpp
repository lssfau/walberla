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
//! \file IntersectionRatio.cpp
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//
//======================================================================================================================

#include "core/DataTypes.h"
#include "core/Environment.h"
#include "core/debug/TestSubsystem.h"
#include "core/math/all.h"

#include "mesa_pd/common/RayParticleIntersection.h"
#include "mesa_pd/data/ParticleAccessor.h"
#include "mesa_pd/data/ParticleStorage.h"
#include "mesa_pd/data/ShapeStorage.h"
#include "mesa_pd/data/DataTypes.h"
#include "mesa_pd/kernel/SingleCast.h"

namespace intersection_ratio_test
{
using namespace walberla;


class ParticleAccessorWithShape : public mesa_pd::data::ParticleAccessor
{
public:
   ParticleAccessorWithShape(std::shared_ptr<mesa_pd::data::ParticleStorage>& ps, std::shared_ptr<mesa_pd::data::ShapeStorage>& ss)
         : ParticleAccessor(ps)
         , ss_(ss)
   {}

   mesa_pd::data::BaseShape* getShape(const size_t p_idx) const {return ss_->shapes[ps_->getShapeID(p_idx)].get();}
private:
   std::shared_ptr<mesa_pd::data::ShapeStorage> ss_;
};


/*!\brief Tests the ray-particle intersection ratio functionality of the RPDUtility.h in the lbm_rpd_coupling module
 *
 * Currently the following shapes are tested:
 *  - sphere
 *  - halfspace ( default and rotated )
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
   using ParticleAccessor = ParticleAccessorWithShape;
   ParticleAccessor accessor(ps, shapeStorage);

   const real_t epsilon( real_t(1e-4) );

   mesa_pd::kernel::SingleCast singleCast;
   mesa_pd::RayParticleIntersectionRatioFunctor intersectionRatioFctr;

   ////////////
   // SPHERE //
   ////////////
   {
      real_t sphereRadius = real_t(1);
      auto sphereShape = shapeStorage->create<mesa_pd::data::Sphere>( sphereRadius );

      Vector3<real_t> position(real_t(1), real_t(0), real_t(0));

      mesa_pd::data::Particle&& p = *ps->create();
      p.setPosition(position);
      p.setShapeID(sphereShape);
      auto idx = p.getIdx();

      Vector3<real_t> pos1(real_t(-0.5), real_t(0), real_t(0));
      Vector3<real_t> dir1(real_t(1), real_t(0), real_t(0));
      real_t delta1 = singleCast(idx, accessor, intersectionRatioFctr, accessor, pos1, dir1, epsilon );
      WALBERLA_CHECK_FLOAT_EQUAL(delta1, real_t(0.5), "Intersection ratio 1 with sphere wrong!");

      Vector3<real_t> pos2(real_t(1), real_t(1), real_t(1));
      Vector3<real_t> dir2(real_t(0), -real_t(1), -real_t(1));
      real_t delta2 = singleCast(idx, accessor, intersectionRatioFctr, accessor, pos2, dir2, epsilon );
      WALBERLA_CHECK_FLOAT_EQUAL(delta2, (std::sqrt(2) - real_t(1)) / std::sqrt(2), "Intersection ratio 2 with sphere wrong!");
   }

   ///////////////
   // HALFSPACE //
   ///////////////
   {
      Vector3<real_t> position(real_t(1), real_t(0), real_t(0));
      Vector3<real_t> normal(real_t(0), real_t(1), real_t(1));

      auto planeShape = shapeStorage->create<mesa_pd::data::HalfSpace>( normal.getNormalized() );

      mesa_pd::data::Particle&& p = *ps->create(true);
      p.setPosition(position);
      p.setShapeID(planeShape);
      auto idx = p.getIdx();

      Vector3<real_t> pos1(real_t(1), real_t(0.5), real_t(0.5));
      Vector3<real_t> dir1(real_t(0), -real_t(1), -real_t(1));
      real_t delta1 = singleCast(idx, accessor, intersectionRatioFctr, accessor, pos1, dir1, epsilon );
      WALBERLA_CHECK_FLOAT_EQUAL(delta1, real_t(0.5), "Intersection ratio 1 with half space wrong!");

      Vector3<real_t> dir2(real_t(0), real_t(0), -real_t(2));
      real_t delta2 = singleCast(idx, accessor, intersectionRatioFctr, accessor, pos1, dir2, epsilon );
      WALBERLA_CHECK_FLOAT_EQUAL(delta2, real_t(0.5), "Intersection ratio 2 with half space wrong!");

      Vector3<real_t> dir3(real_t(0), -real_t(3), real_t(0));
      real_t delta3 = singleCast(idx, accessor, intersectionRatioFctr, accessor, pos1, dir3, epsilon );
      WALBERLA_CHECK_FLOAT_EQUAL(delta3, real_t(1)/real_t(3), "Intersection ratio 3 with half space wrong!");
   }

   /////////////////////////
   // HALFSPACE (rotated) //
   /////////////////////////
   {
      Vector3<real_t> position(real_t(1), real_t(0), real_t(0));
      Vector3<real_t> normal(real_t(0), real_t(0), real_t(1));

      auto planeShape = shapeStorage->create<mesa_pd::data::HalfSpace>( normal.getNormalized() );

      // rotate to same position as half space before
      Vector3<real_t> rotationAngles( -math::M_PI / real_t(4), real_t(0), real_t(0));
      Quaternion<real_t> quat( rotationAngles );

      mesa_pd::data::Particle&& p = *ps->create(true);
      p.setPosition(position);
      p.setShapeID(planeShape);
      p.setRotation(quat);
      auto idx = p.getIdx();

      Vector3<real_t> pos1(real_t(1), real_t(0.5), real_t(0.5));
      Vector3<real_t> dir1(real_t(0), -real_t(1), -real_t(1));
      real_t delta1 = singleCast(idx, accessor, intersectionRatioFctr, accessor, pos1, dir1, epsilon );
      WALBERLA_CHECK_FLOAT_EQUAL(delta1, real_t(0.5), "Intersection ratio 1 with rotated half space wrong!");

      Vector3<real_t> dir2(real_t(0), real_t(0), -real_t(2));
      real_t delta2 = singleCast(idx, accessor, intersectionRatioFctr, accessor, pos1, dir2, epsilon );
      WALBERLA_CHECK_FLOAT_EQUAL(delta2, real_t(0.5), "Intersection ratio 2 with rotated half space wrong!");

      Vector3<real_t> dir3(real_t(0), -real_t(3), real_t(0));
      real_t delta3 = singleCast(idx, accessor, intersectionRatioFctr, accessor, pos1, dir3, epsilon );
      WALBERLA_CHECK_FLOAT_EQUAL(delta3, real_t(1)/real_t(3), "Intersection ratio 3 with rotated half space wrong!");
   }

   ///////////////////////////
   // Bisection Line Search //
   ///////////////////////////
   {
      // note: here tested with a sphere and therefore called explicitly because otherwise the sphere specialization would be selected

      real_t sphereRadius = real_t(1);
      auto sphereShape = shapeStorage->create<mesa_pd::data::Sphere>( sphereRadius );

      Vector3<real_t> position(real_t(1), real_t(0), real_t(0));

      mesa_pd::data::Particle&& p = *ps->create();
      p.setPosition(position);
      p.setShapeID(sphereShape);
      auto idx = p.getIdx();

      Vector3<real_t> pos1(real_t(-0.5), real_t(0), real_t(0));
      Vector3<real_t> dir1(real_t(1), real_t(0), real_t(0));
      real_t delta1 = mesa_pd::intersectionRatioBisection(idx, accessor, pos1, dir1, epsilon);
      WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(delta1, real_t(0.5), epsilon, "Intersection ratio 1 with bisection line search for sphere wrong!");

      Vector3<real_t> pos2(real_t(1), real_t(1), real_t(1));
      Vector3<real_t> dir2(real_t(0), -real_t(1), -real_t(1));
      real_t delta2 = mesa_pd::intersectionRatioBisection(idx, accessor, pos2, dir2, epsilon);
      WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(delta2, (std::sqrt(2) - real_t(1)) / std::sqrt(2), epsilon, "Intersection ratio 2 with bisection line search for sphere wrong!");
   }

   return 0;

}

} //namespace intersection_ratio_test

int main( int argc, char **argv ){
   intersection_ratio_test::main(argc, argv);
}
