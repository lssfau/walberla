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
//! \file MeshPeRaytracing.cpp
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include <core/debug/TestSubsystem.h>
#include <core/logging/Logging.h>
#include <core/mpi/Environment.h>

#include <core/DataTypes.h>
#include <core/math/Vector3.h>
#include <mesh/pe/raytracing/Intersects.h>
#include <mesh/pe/rigid_body/ConvexPolyhedron.h>
#include <mesh/QHull.h>
#include <pe/Materials.h>
#include <pe/raytracing/Intersects.h>
#include <pe/rigidbody/Box.h>

#include <vector>

namespace walberla{
namespace mesh {
namespace pe {

int CpRayIntersectionTest(const int resolution = 10)
{
   using namespace walberla::math;
   using namespace walberla::pe::raytracing;

   std::vector<Vector3<real_t>> points;
   points.emplace_back( real_t(-1), real_t(-1), real_t(-1) );
   points.emplace_back( real_t(-1), real_t(-1), real_t( 1) );
   points.emplace_back( real_t(-1), real_t( 1), real_t(-1) );
   points.emplace_back( real_t(-1), real_t( 1), real_t( 1) );
   points.emplace_back( real_t( 1), real_t(-1), real_t(-1) );
   points.emplace_back( real_t( 1), real_t(-1), real_t( 1) );
   points.emplace_back( real_t( 1), real_t( 1), real_t(-1) );
   points.emplace_back( real_t( 1), real_t( 1), real_t( 1) );

   shared_ptr< TriangleMesh > mesh = make_shared<TriangleMesh>();
   mesh::QHull<TriangleMesh> qhull( points, mesh );
   qhull.run();

   const Vec3 center(1,2,3);

   ConvexPolyhedron cp(0, 0, center,Quat(), *mesh, Material::find("iron"), false, true, true);
   cp.rotate(real_t(1), real_t(2), real_t(3));
   Box bx(0, 0, center, Quat(), Vec3(2,2,2), Material::find("iron"), false, true, true);
   bx.rotate(real_t(1), real_t(2), real_t(3));

   real_t dx = real_t(1.0) / static_cast<real_t>(resolution);
   //rays pointed at center of body
   for (int x = 0; x < resolution; ++x)
   {
      WALBERLA_LOG_INFO("[" << x+1 << " / " << resolution << "]" );
      const real_t rand1 = real_c(x) * dx;
      for (int y = 0; y < resolution; ++y)
      {
         const real_t rand2 = real_c(y) * dx;
         real_t theta = real_t(2) * real_t(math::pi) * rand1;
         real_t phi = std::acos(real_t(1) - real_t(2) * rand2);
         Vec3 dir(std::sin(phi) * std::cos(theta), std::sin(phi) * std::sin(theta), std::cos(phi));

         Ray ray( center + dir*real_t(5), -dir);
         real_t bx_t, cp_t;
         Vec3   bx_n, cp_n;
         WALBERLA_CHECK( intersects(&bx, ray, bx_t, bx_n) );
         WALBERLA_CHECK( intersects(&cp, ray, cp_t, cp_n) );
         WALBERLA_CHECK_FLOAT_EQUAL(bx_t, cp_t);
         WALBERLA_CHECK_FLOAT_EQUAL(bx_n, cp_n);
      }
   }

   //rays emitted form a point outside of the body
   for (int x = 0; x < resolution; ++x)
   {
      WALBERLA_LOG_INFO("[" << x+1 << " / " << resolution << "]" );
      const real_t rand1 = real_c(x) * dx;
      for (int y = 0; y < resolution; ++y)
      {
         const real_t rand2 = real_c(y) * dx;
         real_t theta = real_t(2) * math::pi * rand1;
         real_t phi = std::acos(real_t(1) - real_t(2) * rand2);
         Vec3 dir(std::sin(phi) * std::cos(theta), std::sin(phi) * std::sin(theta), std::cos(phi));

         Ray ray( Vec3(real_t(5),real_t(5),real_t(5)), -dir);
         real_t bx_t, cp_t;
         Vec3   bx_n, cp_n;
         const bool bx_intersects = intersects(&bx, ray, bx_t, bx_n);
         const bool cp_intersects = intersects(&cp, ray, cp_t, cp_n);
         WALBERLA_CHECK_EQUAL( bx_intersects, cp_intersects );
         if (bx_intersects)
         {
            WALBERLA_CHECK_FLOAT_EQUAL(bx_t, cp_t);
            WALBERLA_CHECK_FLOAT_EQUAL(bx_n, cp_n);
         }
      }
   }

   return EXIT_SUCCESS;
}

} //namespace pe
} //namespace mesh
} //namespace walberla

int main( int argc, char * argv[] )
{
   walberla::debug::enterTestMode();
   walberla::mpi::Environment mpiEnv( argc, argv );
   walberla::mpi::MPIManager::instance()->useWorldComm();

   return walberla::mesh::pe::CpRayIntersectionTest(10);
}
