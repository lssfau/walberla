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
//! \file   02_ConfinedGasExtended.cpp
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include <pe/basic.h>
#include <pe/statistics/BodyStatistics.h>
#include <pe/vtk/SphereVtkOutput.h>
#include <pe/raytracing/Raytracer.h>

#include <blockforest/Initialization.h>
#include <core/Abort.h>
#include <core/Environment.h>
#include <core/grid_generator/HCPIterator.h>
#include <core/grid_generator/SCIterator.h>
#include <core/logging/Logging.h>
#include <core/math/Random.h>
#include <core/timing/TimingTree.h>
#include <core/waLBerlaBuildInfo.h>
#include <postprocessing/sqlite/SQLite.h>
#include <vtk/VTKOutput.h>

#include <functional>
#include <tuple>

namespace walberla {
using namespace walberla::pe;
using namespace walberla::timing;
using namespace walberla::pe::raytracing;

typedef std::tuple<Sphere, Plane> BodyTuple ;

int main( int argc, char ** argv )
{
   Environment env(argc, argv);

   logging::Logging::instance()->setStreamLogLevel(logging::Logging::INFO);
   logging::Logging::instance()->setFileLogLevel(logging::Logging::INFO);

   WALBERLA_LOG_INFO_ON_ROOT( "config file: " << argv[1] )
   WALBERLA_LOG_INFO_ON_ROOT( "waLBerla Revision: " << WALBERLA_GIT_SHA1 );

   math::seedRandomGenerator( static_cast<unsigned int>(1337 * mpi::MPIManager::instance()->worldRank()) );

   WcTimingTree tt;

   WALBERLA_LOG_INFO_ON_ROOT("*** READING COMMANDLINE ARGUMENTS ***");
   bool syncShadowOwners = false;

   for( int i = 1; i < argc; ++i )
   {
      if( std::strcmp( argv[i], "--syncShadowOwners" ) == 0 ) syncShadowOwners = true;
   }

   WALBERLA_LOG_INFO_ON_ROOT("*** READING CONFIG FILE ***");
   auto cfg = env.config();
   if (cfg == nullptr) WALBERLA_ABORT("No config specified!");
   const Config::BlockHandle mainConf  = cfg->getBlock( "ConfinedGasExtended" );

   const std::string sqlFile = mainConf.getParameter< std::string >( "sqlFile", "ConfinedGas.sqlite" );

   //! [SQLProperties]
   std::map< std::string, walberla::int64_t >     integerProperties;
   std::map< std::string, double >      realProperties;
   std::map< std::string, std::string > stringProperties;
   //! [SQLProperties]

   stringProperties["walberla_git"] = WALBERLA_GIT_SHA1;

   real_t spacing = mainConf.getParameter<real_t>("spacing", real_c(1.0) );
   WALBERLA_LOG_INFO_ON_ROOT("spacing: " << spacing);
   realProperties["spacing"] = double_c(spacing);

   real_t radius = mainConf.getParameter<real_t>("radius", real_c(0.4) );
   WALBERLA_LOG_INFO_ON_ROOT("radius: " << radius);
   realProperties["radius"] = double_c(radius);

   real_t vMax = mainConf.getParameter<real_t>("vMax", real_c(1.0) );
   WALBERLA_LOG_INFO_ON_ROOT("vMax: " << vMax);
   realProperties["vMax"] = vMax;

   int warmupSteps = mainConf.getParameter<int>("warmupSteps", 0 );
   WALBERLA_LOG_INFO_ON_ROOT("warmupSteps: " << warmupSteps);
   integerProperties["warmupSteps"] = warmupSteps;

   int simulationSteps = mainConf.getParameter<int>("simulationSteps", 10 );
   WALBERLA_LOG_INFO_ON_ROOT("simulationSteps: " << simulationSteps);
   integerProperties["simulationSteps"] = simulationSteps;

   real_t dt = mainConf.getParameter<real_t>("dt", real_c(0.01) );
   WALBERLA_LOG_INFO_ON_ROOT("dt: " << dt);
   realProperties["dt"] = dt;

   const int visSpacing = mainConf.getParameter<int>("visSpacing",  1000 );
   WALBERLA_LOG_INFO_ON_ROOT("visSpacing: " << visSpacing);
   const std::string path = mainConf.getParameter<std::string>("path",  "vtk_out" );
   WALBERLA_LOG_INFO_ON_ROOT("path: " << path);

   WALBERLA_LOG_INFO_ON_ROOT("syncShadowOwners: " << syncShadowOwners);
   integerProperties["syncShadowOwners"] = syncShadowOwners;

   WALBERLA_LOG_INFO_ON_ROOT("*** GLOBALBODYSTORAGE ***");
   shared_ptr<BodyStorage> globalBodyStorage = make_shared<BodyStorage>();

   WALBERLA_LOG_INFO_ON_ROOT("*** BLOCKFOREST ***");
   // create forest
   shared_ptr< BlockForest > forest = blockforest::createBlockForestFromConfig( mainConf );
   if (!forest)
   {
      WALBERLA_LOG_INFO_ON_ROOT( "No BlockForest created ... exiting!");
      return EXIT_SUCCESS;
   }

   WALBERLA_LOG_INFO_ON_ROOT("simulationDomain: " << forest->getDomain());
   integerProperties["sim_x"] = int64_c(forest->getDomain().maxCorner()[0]);
   integerProperties["sim_y"] = int64_c(forest->getDomain().maxCorner()[1]);
   integerProperties["sim_z"] = int64_c(forest->getDomain().maxCorner()[2]);

   WALBERLA_LOG_INFO_ON_ROOT("blocks: " << Vector3<uint_t>(forest->getXSize(), forest->getYSize(), forest->getZSize()) );
   integerProperties["blocks_x"] = int64_c(forest->getXSize());
   integerProperties["blocks_y"] = int64_c(forest->getYSize());
   integerProperties["blocks_z"] = int64_c(forest->getZSize());

   WALBERLA_LOG_INFO_ON_ROOT("*** BODYTUPLE ***");
   // initialize body type ids
   SetBodyTypeIDs<BodyTuple>::execute();

   WALBERLA_LOG_INFO_ON_ROOT("*** STORAGEDATAHANDLING ***");
   // add block data
   auto storageID           = forest->addBlockData(createStorageDataHandling<BodyTuple>(), "Storage");
   auto ccdID               = forest->addBlockData(ccd::createHashGridsDataHandling( globalBodyStorage, storageID ), "CCD");
   auto fcdID               = forest->addBlockData(fcd::createGenericFCDDataHandling<BodyTuple, fcd::AnalyticCollideFunctor>(), "FCD");

   WALBERLA_LOG_INFO_ON_ROOT("*** INTEGRATOR ***");
   cr::HCSITS cr(globalBodyStorage, forest, storageID, ccdID, fcdID);
   cr.setMaxIterations( 10 );
   cr.setRelaxationModel( cr::HardContactSemiImplicitTimesteppingSolvers::ApproximateInelasticCoulombContactByDecoupling );
   cr.setRelaxationParameter( real_t(0.7) );
   cr.setGlobalLinearAcceleration( Vec3(0,0,0) );

   WALBERLA_LOG_INFO_ON_ROOT("*** SYNCCALL ***");
   std::function<void(void)> syncCall;
   if (!syncShadowOwners)
   {
      syncCall = std::bind( pe::syncNextNeighbors<BodyTuple>, std::ref(*forest), storageID, &tt, real_c(0.0), false );
   } else
   {
      syncCall = std::bind( pe::syncShadowOwners<BodyTuple>, std::ref(*forest), storageID, &tt, real_c(0.0), false );
   }

   //! [Bind Sync Call]
   std::function<void(void)> syncCallWithoutTT;
   if (!syncShadowOwners)
   {
      syncCallWithoutTT = std::bind( pe::syncNextNeighbors<BodyTuple>, std::ref(*forest), storageID, static_cast<WcTimingTree*>(nullptr), real_c(0.0), false );
   } else
   {
      syncCallWithoutTT = std::bind( pe::syncShadowOwners<BodyTuple>, std::ref(*forest), storageID, static_cast<WcTimingTree*>(nullptr), real_c(0.0), false );
   }
   //! [Bind Sync Call]
   
   WALBERLA_LOG_INFO_ON_ROOT("*** RAYTRACER ***");
   //! [Raytracer Init]
   std::function<ShadingParameters (const BodyID body)> customShadingFunction = [](const BodyID body) {
      if (body->getTypeID() == Sphere::getStaticTypeID()) {
         return processRankDependentShadingParams(body).makeGlossy();
      }
      return defaultBodyTypeDependentShadingParams(body);
   };
   Raytracer raytracer(forest, storageID, globalBodyStorage, ccdID,
                       cfg->getBlock("Raytracing"),
                       customShadingFunction);
   //! [Raytracer Init]

   WALBERLA_LOG_INFO_ON_ROOT("*** VTK ***");
   //! [VTK Domain Output]
   auto vtkDomainOutput = vtk::createVTKOutput_DomainDecomposition( forest, "domain_decomposition", 1, "vtk_out", "simulation_step" );
   //! [VTK Domain Output]
   //! [VTK Sphere Output]
   auto vtkSphereHelper = make_shared<SphereVtkOutput>(storageID, *forest) ;
   auto vtkSphereOutput = vtk::createVTKOutput_PointData(vtkSphereHelper, "Bodies", 1, "vtk_out", "simulation_step", false, false);
   //! [VTK Sphere Output]

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

   for (auto blkIt = forest->begin(); blkIt != forest->end(); ++blkIt)
   {
      IBlock & currentBlock = *blkIt;
      for (auto it = grid_generator::SCIterator(currentBlock.getAABB().getIntersection(generationDomain), Vector3<real_t>(spacing, spacing, spacing) * real_c(0.5), spacing); it != grid_generator::SCIterator(); ++it)
      {
         SphereID sp = pe::createSphere( *globalBodyStorage, *forest, storageID, 0, *it, radius, material);
         Vec3 rndVel(math::realRandom<real_t>(-vMax, vMax), math::realRandom<real_t>(-vMax, vMax), math::realRandom<real_t>(-vMax, vMax));
         if (sp != nullptr) sp->setLinearVel(rndVel);
         if (sp != nullptr) ++numParticles;
      }
   }
   mpi::reduceInplace(numParticles, mpi::SUM);
   WALBERLA_LOG_INFO_ON_ROOT("#particles created: " << numParticles);

   WALBERLA_LOG_INFO_ON_ROOT("*** SETUP - END ***");

   // synchronize particles
   //! [TT Example]
   //WcTimingTree tt;

   tt.start("Initial Sync");
   syncCallWithoutTT();
   syncCallWithoutTT();
   tt.stop("Initial Sync");
   //! [TT Example]

   WALBERLA_LOG_INFO_ON_ROOT("*** SIMULATION - START ***");
   WcTimingPool tp;
   tt.start("Simulation Loop");
   tp["Total"].start();
   for (int i=0; i < simulationSteps; ++i)
   {
      if( i % 200 == 0 )
      {
         WALBERLA_LOG_DEVEL_ON_ROOT( "Timestep " << i << " / " << simulationSteps );
      }

      tp["Solver"].start();
      cr.timestep( real_c(dt) );
      tp["Solver"].end();
      tp["Sync"].start();
      syncCall();
      tp["Sync"].end();

      if( i % visSpacing == 0 )
      {
         //! [VTK Output]
         vtkDomainOutput->write( );
         vtkSphereOutput->write( );
         //! [VTK Output]
         //! [Image Output]
         raytracer.generateImage<BodyTuple>(size_t(i), &tt);
         //! [Image Output]
      }
   }
   tp["Total"].end();
   tt.stop("Simulation Loop");
   WALBERLA_LOG_INFO_ON_ROOT("*** SIMULATION - END ***");

   BodyStatistics bodyStats( forest, storageID );
   bodyStats();
   WALBERLA_LOG_INFO_ON_ROOT( bodyStats );
   integerProperties["numBodies"]        = int64_c(bodyStats.numBodies());
   integerProperties["numShadowCopies"]  = int64_c(bodyStats.numShadowCopies());

   //! [TT Log]
   auto temp = tt.getReduced( );
   WALBERLA_ROOT_SECTION()
   {
      std::cout << temp;
   }
   //! [TT Log]

   auto tpReduced = tp.getReduced();
   //! [SQL Save]
   WALBERLA_ROOT_SECTION()
   {
      auto runId = postprocessing::storeRunInSqliteDB( sqlFile, integerProperties, stringProperties, realProperties );
      postprocessing::storeTimingPoolInSqliteDB( sqlFile, runId, *tpReduced, "Timeloop" );
      postprocessing::storeTimingTreeInSqliteDB( sqlFile, runId, tt, "TimingTree" );
   }
   //! [SQL Save]

   return EXIT_SUCCESS;
}
} // namespace walberla

int main( int argc, char* argv[] )
{
  return walberla::main( argc, argv );
}
