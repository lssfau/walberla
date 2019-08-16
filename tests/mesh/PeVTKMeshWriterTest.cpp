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
//! \file   PeVTKMeshWriterTest.cpp
//! \ingroup mesh
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//
//======================================================================================================================

#include <pe/basic.h>
#include <pe/statistics/BodyStatistics.h>
#include <pe/vtk/SphereVtkOutput.h>
#include <pe/fcd/GJKEPACollideFunctor.h>

#include <mesh/pe/vtk/PeVTKMeshWriter.h>
#include <mesh/pe/vtk/CommonDataSources.h>
#include <mesh/pe/DefaultTesselation.h>
#include <mesh/pe/rigid_body/ConvexPolyhedronFactory.h>
#include <mesh/pe/communication/ConvexPolyhedron.h>
#include <mesh/PolyMeshes.h>

#include <blockforest/Initialization.h>
#include <core/Abort.h>
#include <core/Environment.h>
#include <core/grid_generator/HCPIterator.h>
#include <core/grid_generator/SCIterator.h>
#include <core/logging/Logging.h>
#include <core/timing/TimingTree.h>
#include <core/waLBerlaBuildInfo.h>
#include <core/math/Random.h>
#include <core/math/Utility.h>
#include <vtk/VTKOutput.h>

#include <functional>
#include <random>
#include <tuple>

using namespace walberla;
using namespace walberla::pe;
using namespace walberla::mesh::pe;

typedef std::tuple<ConvexPolyhedron, Plane> BodyTuple ;

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


template< typename RNG >
std::vector<Vector3<real_t>> generatPointCloudOnSphere( const real_t radius, const uint_t numPoints, RNG & rng )
{
   std::uniform_real_distribution<real_t> distribution;

   std::vector<Vector3<real_t>> pointCloud( numPoints );
   for( auto & p : pointCloud )
   {
      real_t theta = 2 * real_t(math::pi) * distribution(rng);
      real_t phi = std::acos( real_t(1.0) - real_t(2.0) * distribution(rng) );
      p[0] = std::sin(phi) * std::cos(theta) * radius;
      p[1] = std::sin(phi) * std::sin(theta) * radius;
      p[2] = std::cos(phi) * radius;
   }

   return pointCloud;
}

int main( int argc, char ** argv )
{
   Environment env(argc, argv);

   logging::Logging::instance()->setStreamLogLevel(logging::Logging::INFO);
   logging::Logging::instance()->setFileLogLevel(logging::Logging::INFO);

   math::seedRandomGenerator( static_cast<unsigned int>(1337 * mpi::MPIManager::instance()->worldRank()) );

   real_t spacing = real_c(1);
   real_t radius = real_c(0.45);
   real_t vMax = real_c(1);
   int simulationSteps = 10000;
   real_t dt = real_c(0.01);
   int visSpacing = 10;

   WALBERLA_DEBUG_SECTION() 
   {
      simulationSteps = 2;
      visSpacing = 1;
   }

   shared_ptr<BodyStorage> globalBodyStorage = make_shared<BodyStorage>();

   auto forest = blockforest::createBlockForest( AABB( real_t(0), real_t(0), real_t(0), real_t(6), real_t(6), real_t(6) ),
                                                 Vector3<uint_t>( uint_t(2) ), Vector3<bool>( false ) );

   SetBodyTypeIDs<BodyTuple>::execute();

   // add block data
   auto storageID           = forest->addBlockData(createStorageDataHandling<BodyTuple>(), "Storage");
   auto ccdID               = forest->addBlockData(ccd::createHashGridsDataHandling( globalBodyStorage, storageID ), "CCD");
   auto fcdID               = forest->addBlockData(fcd::createGenericFCDDataHandling<BodyTuple, fcd::GJKEPACollideFunctor>(), "FCD");

   cr::HCSITS cr(globalBodyStorage, forest, storageID, ccdID, fcdID);
   cr.setMaxIterations( 10 );
   cr.setRelaxationModel( cr::HardContactSemiImplicitTimesteppingSolvers::ApproximateInelasticCoulombContactByDecoupling );
   cr.setRelaxationParameter( real_t(0.7) );
   cr.setGlobalLinearAcceleration( Vec3(0,0,5) );

   std::function<void(void)> syncCall = std::bind( pe::syncNextNeighbors<BodyTuple>, std::ref(*forest), storageID, static_cast<WcTimingTree*>(nullptr), real_c(0.0), false );

   using OutputMeshType = mesh::FloatPolyMesh;
   using TesselationType = mesh::pe::DefaultTesselation<OutputMeshType>;
   TesselationType tesselation;
   mesh::pe::PeVTKMeshWriter<OutputMeshType, TesselationType> writer( forest, storageID, tesselation, "MeshOutput", uint_c(visSpacing) );

   shared_ptr<mesh::pe::PeVTKMeshWriter<OutputMeshType, TesselationType>::FaceDataSource<float>> linearVelocityFace = make_shared<mesh::pe::LinearVelocityFaceDataSource<OutputMeshType, TesselationType, float>>();
   shared_ptr<mesh::pe::PeVTKMeshWriter<OutputMeshType, TesselationType>::FaceDataSource<float>> angularVelocityFace = make_shared<mesh::pe::AngularVelocityFaceDataSource<OutputMeshType, TesselationType, float>>();
   shared_ptr<mesh::pe::PeVTKMeshWriter<OutputMeshType, TesselationType>::VertexDataSource<float>> surfaceVelocityVertex = make_shared<mesh::pe::SurfaceVelocityVertexDataSource<OutputMeshType, TesselationType, float>>();

   writer.addDataSource( linearVelocityFace );
   writer.addDataSource( angularVelocityFace );
   writer.addDataSource( surfaceVelocityVertex );

   writer.addFacePropertyRank();

   WALBERLA_LOG_INFO_ON_ROOT("*** SETUP - START ***");
   const real_t   static_cof  ( real_c(0.1) / 2 );   // Coefficient of static friction. Note: pe doubles the input coefficient of friction for material-material contacts.
   const real_t   dynamic_cof ( static_cof ); // Coefficient of dynamic friction. Similar to static friction for low speed friction.
   MaterialID     material = createMaterial( "granular", real_t( 1.0 ), 0, static_cof, dynamic_cof, real_t( 0.5 ), 1, 1, 0, 0 );

   auto simulationDomain = forest->getDomain();
   const auto& generationDomain = simulationDomain; // simulationDomain.getExtended(-real_c(0.5) * spacing);
   createPlane(*globalBodyStorage, 0, Vec3(1,0,0), simulationDomain.minCorner(), material );
   createPlane(*globalBodyStorage, 0, Vec3(-1,0,0), simulationDomain.maxCorner(), material );
   createPlane(*globalBodyStorage, 0, Vec3(0,1,0), simulationDomain.minCorner(), material );
   createPlane(*globalBodyStorage, 0, Vec3(0,-1,0), simulationDomain.maxCorner(), material );
   createPlane(*globalBodyStorage, 0, Vec3(0,0,1), simulationDomain.minCorner(), material );
   createPlane(*globalBodyStorage, 0, Vec3(0,0,-1), simulationDomain.maxCorner(), material );

   uint_t numParticles = uint_c(0);

   //std::vector<Vec3> pointCloud = generatePointCloudDodecahedron();
   //for( auto & p: pointCloud)
   //   p = p.getNormalized() * radius;

   std::mt19937 rng(42);
   for (auto blkIt = forest->begin(); blkIt != forest->end(); ++blkIt)
   {
      IBlock & currentBlock = *blkIt;
      for (auto it = grid_generator::SCIterator(currentBlock.getAABB().getIntersection(generationDomain), Vector3<real_t>(spacing, spacing, spacing) * real_c(0.5), spacing); it != grid_generator::SCIterator(); ++it)
      {
         std::vector<Vec3> pointCloud = generatPointCloudOnSphere(radius, 10, rng);
         //SphereID sp = pe::createSphere( *globalBodyStorage, *forest, storageID, 0, *it, radius, material);
         //BoxID sp = pe::createBox(  *globalBodyStorage, *forest, storageID, 0, *it, Vec3(radius), material );
         mesh::pe::ConvexPolyhedronID sp = mesh::pe::createConvexPolyhedron(  *globalBodyStorage, *forest, storageID, 0, *it, pointCloud, material );
         Vec3 rndVel(math::realRandom<real_t>(-vMax, vMax), math::realRandom<real_t>(-vMax, vMax), math::realRandom<real_t>(-vMax, vMax));
         if (sp != nullptr) sp->setLinearVel(rndVel);
         if (sp != nullptr) ++numParticles;
      }
   }
   mpi::reduceInplace(numParticles, mpi::SUM);
   WALBERLA_LOG_INFO_ON_ROOT("#particles created: " << numParticles);

   auto vtkDomainOutput = vtk::createVTKOutput_DomainDecomposition( forest, "domain_decomposition", 1, "vtk_out", "simulation_step" );

   writer();
   vtkDomainOutput->write();

   // synchronize particles
   syncCall();
   syncCall();


   WALBERLA_LOG_INFO_ON_ROOT("*** SIMULATION - START ***");
   for (int i=0; i < simulationSteps; ++i)
   {
      if( i % 10 == 0 )
      {
         WALBERLA_LOG_DEVEL_ON_ROOT( "Timestep " << i << " / " << simulationSteps );
      }

      cr.timestep( real_c(dt) );
      syncCall();
      writer();
   }
   WALBERLA_LOG_INFO_ON_ROOT("*** SIMULATION - END ***");

   return EXIT_SUCCESS;
}
