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
//! \file QHullTest.cpp
//! \ingroup mesh
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//
//======================================================================================================================

#include "core/debug/TestSubsystem.h"
#include "core/logging/Logging.h"
#include "core/mpi/Environment.h"
#include "core/timing/Timer.h"
#include "core/math/Utility.h"

#include "mesh/TriangleMeshes.h"
#include "mesh/QHull.h"
#include "mesh/vtk/VTKMeshWriter.h"
#include "mesh/vtk/CommonDataSources.h"

#include "vtk/VTKOutput.h"
#include "vtk/PointDataSource.h"

#include <random>

#include <vector>
#include <string>

namespace walberla {
namespace mesh {

class PointCloudDataSource : public vtk::PointDataSource
{
public:
   PointCloudDataSource( const std::vector<Vector3<real_t>> & pointCloud ) : pointCloud_( pointCloud ) {}

   std::vector< Attributes > getAttributes() const override { return std::vector< Attributes >(); }
   std::vector< Vector3< real_t > > getPoints() override { return pointCloud_; }
   void configure() override {};

   void push( std::ostream& /*os*/,  const uint_t /*data*/, const uint_t /*point*/, const uint_t /*component*/ ) override {};
   void push( vtk::Base64Writer& /*b64*/, const uint_t /*data*/, const uint_t /*point*/, const uint_t /*component*/ ) override {};

private:
   std::vector<Vector3<real_t>> pointCloud_;
};



template< typename MeshType >
void test( const std::string & testName, const std::vector<Vector3<real_t>> & pointCloud, const bool doVTKOutput )
{
   // Textfile output for qhull comparison
   //std::ofstream ofs(testName + ".txt");
   //ofs << "3\n"
   //    << pointCloud.size() << "\n";
   //for( const auto & p : pointCloud )
   //   ofs << std::setprecision(16) << p[0] << " " << p[1] << " " << p[2] << "\n";

   mesh::QHull<MeshType> qhull( pointCloud );
   
   if( doVTKOutput )
   {
      auto pcds = walberla::make_shared<PointCloudDataSource>( pointCloud );
      vtk::createVTKOutput_PointData(pcds, testName + "PointCloud", 1)->write();
      qhull.enableDebugVTKOutput( testName );   
   }

   WcTimer timer;
   uint_t iterations = qhull.run();
   timer.end();
   WALBERLA_LOG_INFO( "QHull on \"" << testName << "\":\n" << 
                      "   Point cloud size: " << pointCloud.size() << "\n" << 
                      "   Iterations:       " << iterations << "\n" <<                
                      "   Num hull points:  " << qhull.mesh().n_vertices() << "\n" << 
                      "   Num hull faces:   " << qhull.mesh().n_faces() << "\n" << 
                      "   Time:             " << timer.last() << "s\n" <<                   
                      "   VTK output:       " << (doVTKOutput ? "enabled" : "disabled") );

   const MeshType & mesh = qhull.mesh();

   // Test wether result mesh is watertight
   for( auto eh : mesh.edges())
   {
      WALBERLA_CHECK( !mesh.is_boundary( eh ) );
   }

   // Test wether alls Points are outside of the polyhedron
   for(auto fh : mesh.faces())
   {
      const typename MeshType::Normal & n = mesh.normal(fh);
      const typename MeshType::Point    pointOnFace = mesh.point( *mesh.cfv_begin(fh) );
      for(const auto & p : pointCloud)
      {
         using Scalar = typename MeshType::Scalar;
         WALBERLA_CHECK_LESS( (toOpenMeshNumericCast<Scalar>(p) - pointOnFace) | n, real_comparison::Epsilon<Scalar>::value,
                                 "Point: " << p << " Face normal: " << n << " v: " << toOpenMeshNumericCast<Scalar>(p) - pointOnFace );
      }

   }
}


std::vector<Vector3<real_t>> generatePointCloudCube()
{
   std::vector<Vector3<real_t>> points;
   points.emplace_back( real_t(-1), real_t(-1), real_t(-1) );
   points.emplace_back( real_t(-1), real_t(-1), real_t( 1) );
   points.emplace_back( real_t(-1), real_t( 1), real_t(-1) );
   points.emplace_back( real_t(-1), real_t( 1), real_t( 1) );
   points.emplace_back( real_t( 1), real_t(-1), real_t(-1) );
   points.emplace_back( real_t( 1), real_t(-1), real_t( 1) );
   points.emplace_back( real_t( 1), real_t( 1), real_t(-1) );
   points.emplace_back( real_t( 1), real_t( 1), real_t( 1) );

   return points;
}


std::vector<Vector3<real_t>> generatePointCloudTetrahedron()
{
   std::vector<Vector3<real_t>> points;
   points.emplace_back( real_t( 1), real_t( 1), real_t(-1) );
   points.emplace_back( real_t(-1), real_t(-1), real_t(-1) );
   points.emplace_back( real_t(-1), real_t( 1), real_t( 1) );
   points.emplace_back( real_t( 1), real_t(-1), real_t( 1) );

   return points;
}

std::vector<Vector3<real_t>> generatePointCloudOctahedron()
{
   std::vector<Vector3<real_t>> points;

   for( auto one : {real_t(-1), real_t(1)} )
   {
      points.emplace_back( one,   0,   0 );
      points.emplace_back(   0, one,   0 );
      points.emplace_back(   0,   0, one );
   }

   return points;
}

std::vector<Vector3<real_t>> generatePointCloudIcosahedron()
{
   std::vector<Vector3<real_t>> points;
   
   static const real_t PHI = ( real_t(1) + std::sqrt( real_t(5) ) ) / real_t(2);

   for( auto one : {real_t(-1), real_t(1)} )
      for( auto phi : {-PHI, PHI} )
      {
         points.emplace_back( real_t(  0), real_t(one), real_t(phi) );
         points.emplace_back( real_t(one), real_t(phi), real_t(  0) );
         points.emplace_back( real_t(phi), real_t(  0), real_t(one) );
      }

   return points;
}

std::vector<Vector3<real_t>> generatePointCloudDodecahedron()
{
   std::vector<Vector3<real_t>> points = generatePointCloudCube();

   static const real_t PHI = ( real_t(1) + std::sqrt( real_t(5) ) ) / real_t(2);
   static const real_t PHI_INV = real_t(1) / PHI;

   for( auto phi : {-PHI, PHI} )
      for( auto piv : {-PHI_INV, PHI_INV} )
      {
         points.emplace_back( real_t(  0), real_t(piv), real_t(phi) );
         points.emplace_back( real_t(piv), real_t(phi), real_t(  0) );
         points.emplace_back( real_t(phi), real_t(  0), real_t(piv) );
      }

   return points;
}

std::vector<Vector3<real_t>> generatePointCloudInAABB( const math::AABB & aabb, const uint_t numPoints )
{
   std::mt19937 rng(42);
   
   std::vector<Vector3<real_t>> pointCloud( numPoints );
   for( auto & p : pointCloud )
      p = aabb.randomPoint(rng);
   
   return pointCloud;
}


std::vector<Vector3<real_t>> generatPointCloudOnSphere( const real_t radius, const uint_t numPoints )
{
   std::mt19937 rng(42);
   std::uniform_real_distribution<real_t> distribution;

   std::vector<Vector3<real_t>> pointCloud( numPoints );
   for( auto & p : pointCloud )
   {
      real_t theta = 2 * real_t(M_PI) * distribution(rng);
      real_t phi = std::acos( real_t(1.0) - real_t(2.0) * distribution(rng) );
      p[0] = std::sin(phi) * std::cos(theta) * radius;
      p[1] = std::sin(phi) * std::sin(theta) * radius;
      p[2] = std::cos(phi) * radius;
   }

   return pointCloud;
}


template< typename MeshType >
void runTests( const uint_t numPoints, const bool doVTKOutput )
{
   test<MeshType>( "tetrahedron", generatePointCloudTetrahedron(), doVTKOutput );
   test<MeshType>( "cube", generatePointCloudCube(), doVTKOutput );
   test<MeshType>( "octahedron", generatePointCloudOctahedron(), doVTKOutput );
   test<MeshType>( "icosahedron", generatePointCloudIcosahedron(), doVTKOutput );
   test<MeshType>( "dodecahedron", generatePointCloudDodecahedron(), doVTKOutput );

   math::AABB aabb( real_t(-1), real_t(-1), real_t(-1), real_t(1), real_t(1), real_t(1) );
   test<MeshType>( "aabb", generatePointCloudInAABB(aabb, numPoints), doVTKOutput );

   test<MeshType>( "sphere", generatPointCloudOnSphere(real_t(1), numPoints), doVTKOutput );
}


int main( int argc, char * argv[] )
{
   debug::enterTestMode();
   mpi::Environment mpiEnv( argc, argv );
   mpi::MPIManager::instance()->useWorldComm();

   std::vector<std::string> args( argv, argv + argc );

   bool doVTKOutput = false;
   auto it = std::find( args.begin(), args.end(), std::string("--vtk") );
   if(it != args.end())
   {
      doVTKOutput = true;
      args.erase( it );
   }

   bool useFloatMesh = false;
   it = std::find( args.begin(), args.end(), std::string("--floatMesh") );
   if(it != args.end())
   {
      useFloatMesh = true;
      args.erase( it );
   }

   if(args.size() != size_t( 2 ))
   {
      WALBERLA_ABORT_NO_DEBUG_INFO( "USAGE:\n" << args[0] << " [--vtk] [--floatMesh] NUM_POINTS" );
   }
   
   uint_t numPoints;
   try {
      numPoints = std::stoul( args[1] );
   }
   catch(std::exception &)
   {
      WALBERLA_ABORT_NO_DEBUG_INFO( "USAGE:\n" << args[0] << " [--vtk] [--floatMesh] NUM_POINTS" );
   }

   if(useFloatMesh)
      runTests<FloatTriangleMesh>( numPoints, doVTKOutput );
   else
      runTests<TriangleMesh>( numPoints, doVTKOutput );

   return EXIT_SUCCESS;
}


} // namespace mesh
} // namespace walberla

int main( int argc, char * argv[] )
{
   return walberla::mesh::main( argc, argv );
}