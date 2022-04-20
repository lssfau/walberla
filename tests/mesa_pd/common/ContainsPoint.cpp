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

#include "mesa_pd/common/Contains.h"
#include "mesa_pd/data/ParticleAccessorWithShape.h"
#include "mesa_pd/data/ParticleStorage.h"
#include "mesa_pd/data/ShapeStorage.h"
#include "mesa_pd/data/DataTypes.h"
#include <mesa_pd/data/shape/Box.h>
#include <mesa_pd/data/shape/CylindricalBoundary.h>
#include <mesa_pd/data/shape/Ellipsoid.h>
#include <mesa_pd/data/shape/HalfSpace.h>
#include <mesa_pd/data/shape/Sphere.h>
#include "mesa_pd/kernel/SingleCast.h"

namespace contains_point_test
{
using namespace walberla;
using mesa_pd::Vec3;


/*!\brief Tests the contains point functionality implemented in mesa_pd/common/Contains.h
 *
 * Currently the following shapes are tested:
 *  - sphere
 *  - halfspace
 *  - box (default and rotated)
 *  - ellipsoid (default and rotated)
 *  - cylindrical boundary
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


   mesa_pd::kernel::SingleCast singleCast;
   mesa_pd::ContainsPointFunctor containsPointFctr;

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

      std::vector<Vec3> testPositions {position, position + Vec3(0,0,1.01_r*sphereRadius)};
      std::vector<bool> shouldBeContained {true, false};

      for(size_t i = 0; i < testPositions.size(); ++i)
      {
         bool isContained = singleCast(idx, accessor, containsPointFctr, accessor, testPositions[i] );
         WALBERLA_CHECK_EQUAL(isContained, shouldBeContained[i], "Sphere check, case " << i << ": Wrong containment info for position " << testPositions[i] );
      }
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

      std::vector<Vec3> testPositions {position - normal, position + normal};
      std::vector<bool> shouldBeContained {true, false};

      for(size_t i = 0; i < testPositions.size(); ++i)
      {
         bool isContained = singleCast(idx, accessor, containsPointFctr, accessor, testPositions[i] );
         WALBERLA_CHECK_EQUAL(isContained, shouldBeContained[i], "Halfspace check, case " << i << ": Wrong containment info for position " << testPositions[i] );
      }
   }

   /////////
   // BOX //
   /////////
   {
      Vec3 position(1_r, 0_r, 0_r);
      Vec3 edgeLength{1_r};

      auto boxShape = shapeStorage->create<mesa_pd::data::Box>( edgeLength );

      mesa_pd::data::Particle&& p = *ps->create();
      p.setPosition(position);
      p.setShapeID(boxShape);
      auto idx = p.getIdx();

      std::vector<Vec3> testPositions {position, position + 0.51_r*Vec3(0_r,0_r,edgeLength[2])};
      std::vector<bool> shouldBeContained {true, false};

      for(size_t i = 0; i < testPositions.size(); ++i)
      {
         bool isContained = singleCast(idx, accessor, containsPointFctr, accessor, testPositions[i] );
         WALBERLA_CHECK_EQUAL(isContained, shouldBeContained[i], "Box check, case " << i << ": Wrong containment info for position " << testPositions[i] );
      }

      auto rotation = p.getRotation();
      rotation.rotate( Vec3(1_r,0_r,0_r), math::pi / 4_r ); // rotate by 45° around x axis
      p.setRotation(rotation);

      shouldBeContained[1] = true;

      for(size_t i = 0; i < testPositions.size(); ++i)
      {
         bool isContained = singleCast(idx, accessor, containsPointFctr, accessor, testPositions[i] );
         WALBERLA_CHECK_EQUAL(isContained, shouldBeContained[i], "Box rotation check, case " << i << ": Wrong containment info for position " << testPositions[i] );
      }
   }

   ///////////////
   // ELLIPSOID //
   ///////////////
   {
      Vec3 position(1_r, 0_r, 0_r);
      Vec3 semiAxes{2_r,1_r, 1_r};

      auto ellipsoidShape = shapeStorage->create<mesa_pd::data::Ellipsoid>( semiAxes );

      mesa_pd::data::Particle&& p = *ps->create();
      p.setPosition(position);
      p.setShapeID(ellipsoidShape);
      auto idx = p.getIdx();

      std::vector<Vec3> testPositions {position,
                                       position + 0.9_r*Vec3(semiAxes[0],0_r,0_r),
                                       position + 1.1_r*Vec3(0_r,semiAxes[1], 0_r),
                                       position + 1.1_r*Vec3(0_r,0_r,semiAxes[2])};
      std::vector<bool> shouldBeContained {true, true, false, false};

      for(size_t i = 0; i < testPositions.size(); ++i)
      {
         bool isContained = singleCast(idx, accessor, containsPointFctr, accessor, testPositions[i] );
         WALBERLA_CHECK_EQUAL(isContained, shouldBeContained[i], "Ellipsoid check, case " << i << ": Wrong containment info for position " << testPositions[i] );
      }

      auto rotation = p.getRotation();
      rotation.rotate( Vec3(0_r,1_r,0_r), math::pi / 2_r ); // rotate by 90° around y axis
      p.setRotation(rotation);

      shouldBeContained[1] = false;
      shouldBeContained[3] = true;

      for(size_t i = 0; i < testPositions.size(); ++i)
      {
         bool isContained = singleCast(idx, accessor, containsPointFctr, accessor, testPositions[i] );
         WALBERLA_CHECK_EQUAL(isContained, shouldBeContained[i], "Ellipsoid rotation check, case " << i << ": Wrong containment info for position " << testPositions[i] );
      }
   }

   //////////////////////////
   // CYLINDRICAL BOUNDARY //
   //////////////////////////
   {
      Vec3 position(1_r, 0_r, 0_r);
      Vec3 axis(0_r, 1_r, 1_r);
      real_t radius = 1_r;

      auto cylindricalBoundaryShape = shapeStorage->create<mesa_pd::data::CylindricalBoundary>( radius, axis.getNormalized() );

      mesa_pd::data::Particle&& p = *ps->create(true);
      p.setPosition(position);
      p.setShapeID(cylindricalBoundaryShape);
      auto idx = p.getIdx();

      Vec3 orthogonalAxis(0_r, -1_r, 1_r);
      std::vector<Vec3> testPositions {position, position + axis, position + 1.1_r*radius*orthogonalAxis.getNormalized()};
      std::vector<bool> shouldBeContained {false, false, true};

      for(size_t i = 0; i < testPositions.size(); ++i)
      {
         bool isContained = singleCast(idx, accessor, containsPointFctr, accessor, testPositions[i] );
         WALBERLA_CHECK_EQUAL(isContained, shouldBeContained[i], "Cylindrical boundary check, case " << i << ": Wrong containment info for position " << testPositions[i] );
      }
   }

   return 0;

}

} //namespace contains_point_test

int main( int argc, char **argv ){
   contains_point_test::main(argc, argv);
}
