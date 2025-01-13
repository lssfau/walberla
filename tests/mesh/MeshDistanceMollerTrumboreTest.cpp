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
//! \file MeshDistanceMollerTrumboreTest.cpp
//! \ingroup mesh
//! \author Brendan Waters <brendan.waters@sydney.edu.au>
//
//======================================================================================================================

#include "core/debug/TestSubsystem.h"
#include "core/logging/Logging.h"
#include "core/mpi/Environment.h"

#include "geometry/mesh/TriangleMesh.h"
#include "geometry/mesh/TriangleMeshIO.h"

#include "mesh_common/TriangleMeshes.h"
#include "mesh_common/MeshOperations.h"
#include "mesh_common/DistanceComputations.h"
#include "mesh_common/distance_octree/DistanceOctree.h"
#include "mesh_common/MeshIO.h"

#include <vector>
#include <string>

namespace walberla {
namespace mesh {

template<typename MeshType>
std::vector<std::tuple<typename MeshType::Point, typename MeshType::Point, float>> 
generateRayInfo( const Vector3<real_t> & min, const Vector3<real_t> & max, std::vector< real_t > & offsets) 
{
   // All combinations of min/max for each axis
   std::vector< real_t > xValues = {min[0], max[0]};
   std::vector< real_t > yValues = {min[1], max[1]};
   std::vector< real_t > zValues = {min[2], max[2]};

   std::vector<std::tuple<typename MeshType::Point, typename MeshType::Point, float>> result;

   for (const auto& offset : offsets) {
      for (const auto& x : xValues) {
         auto x_dir = floatIsEqual(x , max[0]) ? -1 : floatIsEqual(x , min[0]) ? 1 : 0;
        
         typename MeshType::Point origin = { x - offset * x_dir,  0 , 0 };
         typename MeshType::Point direction = { x_dir * (floatIsEqual(origin[0] , 0 ) ? 0 : 1), 0, 0};

         result.emplace_back(origin, direction, offset);
      }

      for (const auto& y : yValues) {
         auto y_dir = floatIsEqual(y , max[1]) ? -1 : floatIsEqual(y , min[1]) ? 1 : 0;

         typename MeshType::Point origin = { 0 ,  y - offset * y_dir , 0  };
         typename MeshType::Point direction = { 0, y_dir * (floatIsEqual(origin[1] , 0 ) ? 0 : 1), 0};

         result.emplace_back(origin, direction, offset);
      }

      for (const auto& z : zValues) {
         auto z_dir = floatIsEqual(z , max[2]) ? -1 : floatIsEqual(z , min[2]) ? 1 : 0;

         typename MeshType::Point origin = { 0 , 0 , z - offset * z_dir };
         typename MeshType::Point direction = { 0, 0, z_dir * (floatIsEqual(origin[2] , 0 ) ? 0 : 1)};

         result.emplace_back(origin, direction, offset);
      }

      for (const auto& y : yValues) {
         for (const auto& x : xValues) {
            
            auto x_dir = floatIsEqual(x , max[0]) ? -1 : floatIsEqual(x , min[0]) ? 1 : 0;
            auto y_dir = floatIsEqual(y , max[1]) ? -1 : floatIsEqual(y , min[1]) ? 1 : 0;

            typename MeshType::Point origin = { x - offset * x_dir, y - offset * y_dir, 0 };

            typename MeshType::Point direction = { x_dir * (floatIsEqual(origin[0] , 0 ) ? 0 : 1), 
                                                   y_dir * (floatIsEqual(origin[1] , 0 ) ? 0 : 1),  0};

            result.emplace_back(origin, direction, offset);
         }
      }

      for (const auto& z : zValues) {
         for (const auto& x : xValues) {
            
            auto x_dir = floatIsEqual(x , max[0]) ? -1 : floatIsEqual(x , min[0]) ? 1 : 0;
            auto z_dir = floatIsEqual(z , max[2]) ? -1 : floatIsEqual(z , min[2]) ? 1 : 0;

            typename MeshType::Point origin = { x - offset * x_dir, 0, z - offset * z_dir };

            typename MeshType::Point direction = { x_dir * (floatIsEqual(origin[0] , 0 ) ? 0 : 1), 0, 
                                                   z_dir * (floatIsEqual(origin[2] , 0 ) ? 0 : 1)};

            result.emplace_back(origin, direction, offset);
         }
      }
   
      for (const auto& z : zValues) {
         for (const auto& y : yValues) {   

            auto y_dir = floatIsEqual(y , max[1]) ? -1 : floatIsEqual(y , min[1]) ? 1 : 0;
            auto z_dir = floatIsEqual(z , max[2]) ? -1 : floatIsEqual(z , min[2]) ? 1 : 0;

            typename MeshType::Point origin = { 0, y - offset * y_dir, z - offset * z_dir };

            typename MeshType::Point direction = { 0, y_dir * (floatIsEqual(origin[1] , 0 ) ? 0 : 1), 
                                                      z_dir * (floatIsEqual(origin[2] , 0 ) ? 0 : 1)};

            result.emplace_back(origin, direction, offset);
         
         }
      }
   
      for (const auto& z : zValues) {
         for (const auto& y : yValues) {
            for (const auto& x : xValues) {
               
               auto x_dir = floatIsEqual(x , max[0]) ? -1 : floatIsEqual(x , min[0]) ? 1 : 0;
               auto y_dir = floatIsEqual(y , max[1]) ? -1 : floatIsEqual(y , min[1]) ? 1 : 0;
               auto z_dir = floatIsEqual(z , max[2]) ? -1 : floatIsEqual(z , min[2]) ? 1 : 0;

               typename MeshType::Point origin = { x - offset * x_dir, 
                                                   y - offset * y_dir, 
                                                   z - offset * z_dir };

               typename MeshType::Point direction = { x_dir * (floatIsEqual(origin[0] , 0 ) ? 0 : 1), 
                                                      y_dir * (floatIsEqual(origin[1] , 0 ) ? 0 : 1), 
                                                      z_dir * (floatIsEqual(origin[2] , 0 ) ? 0 : 1)};

               result.emplace_back(origin, direction, offset);
            }
         }
      }
   
   }

   return result;
}

template<typename MeshType>
void test( const std::string & meshFile )
{
   auto mesh = make_shared<MeshType>();
   mesh::readAndBroadcast( meshFile, *mesh);

   const AABB meshAABB = computeAABB( *mesh );
   const Vector3<real_t> meshMin  = meshAABB.min();
   const Vector3<real_t> meshMax  = meshAABB.max();
   
   auto translationVector = 0.5*( meshMin - meshMax );

   translate( *mesh, translationVector );

   const AABB aabb = computeAABB( *mesh );
   const Vector3<real_t> min  = aabb.min();
   const Vector3<real_t> max = aabb.max();

   auto triDist = make_shared< mesh::TriangleDistance<MeshType> >( mesh );
   DistanceOctree<MeshType> distanceOctree( triDist );

   std::vector< real_t > offsets = {0.0, 0.25, 0.5, 0.77, 1.0, 1.5};

   std::vector<std::tuple<typename MeshType::Point, typename MeshType::Point, float>> 
   offsetCombos = generateRayInfo<MeshType>( min, max, offsets);

   for (auto& [origin, direction, offset] : offsetCombos) 
   {  
      const real_t ray_length = direction.norm();

      // As `direction' has not been normalised yet, q = offset.
      real_t q = distanceOctree.getRayDistanceToMeshObject(origin, direction);

      WALBERLA_CHECK_FLOAT_EQUAL( q, offset); 

      // With normalisation this should produce the actual length
      q = distanceOctree.getRayDistanceToMeshObject(origin, direction.normalize());
      const real_t q_normalised = q/ray_length; // dividing by the ray_length so that q_normalised = offset.

      WALBERLA_CHECK_FLOAT_EQUAL( q_normalised, offset); 

      WALBERLA_LOG_DETAIL_ON_ROOT("Wall Distance: " << q/ray_length << " origin: " << origin << " dir: " << direction << " offset: " << offset )
   }
}

int main( int argc, char * argv[] )
{
   debug::enterTestMode();
   mpi::Environment mpiEnv( argc, argv );
   mpi::MPIManager::instance()->useWorldComm();

   std::vector<std::string> args( argv, argv + argc );
   if( args.size() != 2 )
      WALBERLA_ABORT_NO_DEBUG_INFO( "USAGE: " << args[0] << " cube.obj" );

   const std::string & meshFile = args[1];

   test<mesh::TriangleMesh>( meshFile );
   test<mesh::FloatTriangleMesh>( meshFile );
   test<mesh::PythonTriangleMesh>( meshFile );

   return EXIT_SUCCESS;
}


} // namespace mesh
} // namespace walberla

int main( int argc, char * argv[] )
{
   return walberla::mesh::main( argc, argv );
}