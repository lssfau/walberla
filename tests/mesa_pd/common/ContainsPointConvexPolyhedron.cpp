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
#include "mesa_pd/kernel/SingleCast.h"

#include <mesa_pd/data/shape/ConvexPolyhedron.h>
#include <mesh_common/QHull.h>

namespace contains_point_convex_polyhedron_test
{
using namespace walberla;
using mesa_pd::Vec3;


/*!\brief Tests the contains point functionality for ConvexPolyhedron shape, implemented in mesa_pd/common/Contains.h
 * See tests/mesa_pd/common/ContainsPoint.cpp for other shapes
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

   /////////////////////////////////
   // CONVEX POLYHEDRON ( = BOX ) //
   /////////////////////////////////
   {
      Vec3 position(1_r, 0_r, 0_r);
      Vec3 edgeLength{1_r};
      Vec3 startVertex = - 0.5_r * edgeLength;

      // manually construct triangle mesh of a cube
      std::vector<Vec3> cornerVertices{startVertex,
                                       startVertex + Vec3{edgeLength[0], 0, 0},
                                       startVertex + Vec3{0, edgeLength[1], 0},
                                       startVertex + Vec3{0, 0, edgeLength[2]},
                                       startVertex + Vec3{edgeLength[0], edgeLength[1], 0},
                                       startVertex + Vec3{edgeLength[0], 0, edgeLength[2]},
                                       startVertex + Vec3{0, edgeLength[1], edgeLength[2]},
                                       startVertex + Vec3{edgeLength[0], edgeLength[1], edgeLength[2]}};

      mesh::QHull<mesh::TriangleMesh> qHull(cornerVertices);
      qHull.run();

      auto cpShape = shapeStorage->create<mesa_pd::data::ConvexPolyhedron>( qHull.mesh() );

      mesa_pd::data::Particle&& p = *ps->create();
      p.setPosition(position);
      p.setShapeID(cpShape);
      auto idx = p.getIdx();

      std::vector<Vec3> testPositions {position, position + 0.51_r*Vec3(0_r,0_r,edgeLength[2])};
      std::vector<bool> shouldBeContained {true, false};

      for(size_t i = 0; i < testPositions.size(); ++i)
      {
         bool isContained = singleCast(idx, accessor, containsPointFctr, accessor, testPositions[i] );
         WALBERLA_CHECK_EQUAL(isContained, shouldBeContained[i], "Convex Polyhedron check, case " << i << ": Wrong containment info for position " << testPositions[i] );
      }

      auto rotation = p.getRotation();
      rotation.rotate( Vec3(1_r,0_r,0_r), math::pi / 4_r ); // rotate by 45° around x axis
      p.setRotation(rotation);

      shouldBeContained[1] = true;

      for(size_t i = 0; i < testPositions.size(); ++i)
      {
         bool isContained = singleCast(idx, accessor, containsPointFctr, accessor, testPositions[i] );
         WALBERLA_CHECK_EQUAL(isContained, shouldBeContained[i], "Convex Polyhedron rotation check, case " << i << ": Wrong containment info for position " << testPositions[i] );
      }
   }

   return 0;

}

} //namespace contains_point_convex_polyhedron_test

int main( int argc, char **argv ){
   contains_point_convex_polyhedron_test::main(argc, argv);
}
