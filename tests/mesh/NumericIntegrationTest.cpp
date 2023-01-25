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
//! \file NumericItegrationTest.cpp
//! \ingroup mesh
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//
//======================================================================================================================

#include "core/debug/TestSubsystem.h"
#include "core/grid_generator/SCIterator.h"
#include "core/logging/Logging.h"
#include "core/math/AABB.h"
#include "core/math/KahanSummation.h"
#include "core/mpi/Environment.h"


#include "geometry/containment_octree/ContainmentOctree.h"
#include "geometry/mesh/TriangleMesh.h"
#include "geometry/mesh/TriangleMeshIO.h"

#include "mesh_common/MeshOperations.h"
#include "mesh_common/TriangleMeshes.h"
#include "mesh_common/DistanceComputations.h"
#include "mesh_common/MeshIO.h"
#include "mesh_common/distance_octree/DistanceOctree.h"


namespace walberla {
namespace mesh {

template< typename ContainmentT >
real_t volumeNumeric( const ContainmentT & body, const AABB & aabb, const real_t spacing )
{
   Vector3<real_t> pointOfReference = aabb.min() + Vector3<real_t>( real_t(0.5) * spacing );
   
   uint_t volume = 0;
   
   for(grid_generator::SCIterator it( aabb, pointOfReference, spacing ); it != grid_generator::SCIterator(); ++it)
   {
      if( body.contains( ContainmentT::toPoint( *it ) ) )
         ++volume;
   }

   return real_t(volume) * spacing * spacing * spacing;
}


template< typename ContainmentT >
Vector3<real_t> centroidNumeric( const ContainmentT & body, const AABB & aabb, const real_t spacing )
{
   Vector3<real_t> pointOfReference = aabb.min() + Vector3<real_t>( real_t(0.5) * spacing );

   math::KahanAccumulator<real_t> centroid[3];
   uint_t numPoints = 0;

   for(grid_generator::SCIterator it( aabb, pointOfReference, spacing ); it != grid_generator::SCIterator(); ++it)
   {
      if(body.contains( ContainmentT::toPoint( *it ) ))
      {
         centroid[0] += (*it)[0];
         centroid[1] += (*it)[1];
         centroid[2] += (*it)[2];
         ++numPoints;
      }
   }

   return Vector3<real_t>( centroid[0].get(), centroid[1].get(), centroid[2].get() ) / real_t(numPoints);
}


template< typename ContainmentT >
Matrix3<real_t> inertiaTensorNumeric( const ContainmentT & body, const AABB & aabb, const real_t spacing )
{
   Vector3<real_t> pointOfReference = aabb.min() + Vector3<real_t>( real_t(0.5) * spacing );

   math::KahanAccumulator<real_t> inertiaTensor[6];

   for(grid_generator::SCIterator it( aabb, pointOfReference, spacing ); it != grid_generator::SCIterator(); ++it)
   {
      if(body.contains( ContainmentT::toPoint( *it ) ))
      {
         const Vector3<real_t> p = *it;
         const real_t & x = p[0];
         const real_t & y = p[1];
         const real_t & z = p[2];

         inertiaTensor[0] += y*y + z*z;
         inertiaTensor[1] += -x*y;
         inertiaTensor[2] += -x*z;
         inertiaTensor[3] += x*x + z*z;
         inertiaTensor[4] += -y*z;
         inertiaTensor[5] += x*x + y*y;
      }
   }

   return Matrix3<real_t>( inertiaTensor[0].get(), inertiaTensor[1].get(), inertiaTensor[2].get(),
                           inertiaTensor[1].get(), inertiaTensor[3].get(), inertiaTensor[4].get(),
                           inertiaTensor[2].get(), inertiaTensor[4].get(), inertiaTensor[5].get() ) * (spacing * spacing * spacing);
}


template< typename MeshType >
void testNumeric( const shared_ptr<MeshType> & mesh )
{
   auto triDist = make_shared< mesh::TriangleDistance<mesh::TriangleMesh> >( mesh );
   auto distanceOctree = make_shared< mesh::DistanceOctree<mesh::TriangleMesh> >( triDist );
   auto containmentOctree = make_shared< geometry::ContainmentOctree< mesh::DistanceOctree<mesh::TriangleMesh> > >( distanceOctree );

   AABB aabb = computeAABB(*mesh);
   uint_t numPoints = 1000000;
   real_t spacing = std::pow( aabb.volume() / real_t(numPoints), real_t(1) / real_t(3) );

   WALBERLA_LOG_INFO("Computing numeric volume");
   real_t numericVolume = volumeNumeric(*containmentOctree, aabb, spacing );
   WALBERLA_LOG_INFO("Computing geometrical volume");
   real_t geometricalVolume = computeVolume(*mesh);
   WALBERLA_LOG_INFO("Numerical volume:   " << numericVolume << "\n" <<
                     "Geometrical volume: " << geometricalVolume << "\n" <<
                     "Difference:         " << numericVolume - geometricalVolume );
   WALBERLA_CHECK( std::fabs( numericVolume - geometricalVolume ) < real_t(0.001) || std::fabs( real_t(1) - numericVolume / geometricalVolume ) < real_t(0.001) );


   WALBERLA_LOG_INFO("Computing numeric centroid");
   Vector3<real_t> numericCentroid = centroidNumeric(*containmentOctree, aabb, spacing );
   WALBERLA_LOG_INFO("Computing geometrical centroid");
   Vector3<real_t> geometricalCentroid = toWalberla( computeCentroid(*mesh) );
   WALBERLA_LOG_INFO("Numerical centroid:   " << numericCentroid << "\n" <<
                     "Geometrical centroid: " << geometricalCentroid << "\n" <<
                     "Difference:           " << numericCentroid - geometricalCentroid );

   WALBERLA_CHECK( std::fabs( numericCentroid[0] - geometricalCentroid[0] ) < real_t(0.001) || std::fabs( real_t(1) - numericCentroid[0] / geometricalCentroid[0] ) < real_t(0.001) );
   WALBERLA_CHECK( std::fabs( numericCentroid[0] - geometricalCentroid[0] ) < real_t(0.001) || std::fabs( real_t(1) - numericCentroid[1] / geometricalCentroid[1] ) < real_t(0.001) );
   WALBERLA_CHECK( std::fabs( numericCentroid[0] - geometricalCentroid[0] ) < real_t(0.001) || std::fabs( real_t(1) - numericCentroid[2] / geometricalCentroid[2] ) < real_t(0.001) );

   WALBERLA_LOG_INFO("Computing numeric inertia tensor");
   Matrix3<real_t> numericTensor = inertiaTensorNumeric(*containmentOctree, aabb, spacing );
   WALBERLA_LOG_INFO("Computing geometrical inertia tensor");
   Matrix3<real_t> geometricalTensor = computeInertiaTensor(*mesh);
   WALBERLA_LOG_INFO("Numerical tensor:\n"   << numericTensor << "\n" <<
                     "Geometrical tensor:\n" << geometricalTensor << "\n" <<
                     "Difference:\n"         << numericTensor - geometricalTensor );

   WALBERLA_CHECK( std::fabs( numericTensor[0] - geometricalTensor[0] ) < real_t(0.001) || std::fabs( real_t(1) - numericTensor[0] / geometricalTensor[0] ) < real_t(0.001) );
   WALBERLA_CHECK( std::fabs( numericTensor[1] - geometricalTensor[1] ) < real_t(0.001) || std::fabs( real_t(1) - numericTensor[1] / geometricalTensor[1] ) < real_t(0.001) );
   WALBERLA_CHECK( std::fabs( numericTensor[2] - geometricalTensor[2] ) < real_t(0.001) || std::fabs( real_t(1) - numericTensor[2] / geometricalTensor[2] ) < real_t(0.001) );
   WALBERLA_CHECK( std::fabs( numericTensor[3] - geometricalTensor[3] ) < real_t(0.001) || std::fabs( real_t(1) - numericTensor[3] / geometricalTensor[3] ) < real_t(0.001) );
   WALBERLA_CHECK( std::fabs( numericTensor[4] - geometricalTensor[4] ) < real_t(0.001) || std::fabs( real_t(1) - numericTensor[4] / geometricalTensor[4] ) < real_t(0.001) );
   WALBERLA_CHECK( std::fabs( numericTensor[5] - geometricalTensor[5] ) < real_t(0.001) || std::fabs( real_t(1) - numericTensor[5] / geometricalTensor[5] ) < real_t(0.001) );
   WALBERLA_CHECK( std::fabs( numericTensor[6] - geometricalTensor[6] ) < real_t(0.001) || std::fabs( real_t(1) - numericTensor[6] / geometricalTensor[6] ) < real_t(0.001) );
   WALBERLA_CHECK( std::fabs( numericTensor[7] - geometricalTensor[7] ) < real_t(0.001) || std::fabs( real_t(1) - numericTensor[7] / geometricalTensor[7] ) < real_t(0.001) );
   WALBERLA_CHECK( std::fabs( numericTensor[8] - geometricalTensor[8] ) < real_t(0.001) || std::fabs( real_t(1) - numericTensor[8] / geometricalTensor[8] ) < real_t(0.001) );
}


int main( int argc, char * argv[] )
{
   debug::enterTestMode();
   mpi::Environment mpiEnv( argc, argv );
   mpi::MPIManager::instance()->useWorldComm();

   std::vector<std::string> args( argv, argv + argc );
   if(args.size() != size_t( 2 ))
   {
      WALBERLA_ABORT_NO_DEBUG_INFO( "USAGE:\n" << args[0] << "<MeshFileName>" );
   }

   auto mesh = make_shared< mesh::TriangleMesh >();
   mesh::readAndBroadcast(args[1], *mesh);

   testNumeric(mesh);

   return EXIT_SUCCESS;
}


} // namespace mesh
} // namespace walberla

int main( int argc, char * argv[] )
{
   return walberla::mesh::main( argc, argv );
}
