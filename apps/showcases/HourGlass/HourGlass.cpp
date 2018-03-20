//======================================================================================================================
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
//! \file HourGlass.cpp
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include "blockforest/all.h"
#include "core/all.h"
#include "domain_decomposition/all.h"

#include <pe/basic.h>
#include <pe/rigidbody/BodyIterators.h>
#include <pe/statistics/BodyStatistics.h>
#include <pe/synchronization/ClearSynchronization.h>
#include <pe/synchronization/SyncNextNeighbors.h>
#include <pe/utility/CreateWorld.h>
#include <pe/vtk/BodyVtkOutput.h>
#include <pe/vtk/SphereVtkOutput.h>

#include <blockforest/loadbalancing/DynamicParMetis.h>
#include <blockforest/loadbalancing/PODPhantomData.h>
#include <core/Environment.h>
#include <core/grid_generator/SCIterator.h>
#include <core/logging/Logging.h>
#include <core/mpi/MPITextFile.h>
#include <core/timing/TimingTree.h>
#include <core/waLBerlaBuildInfo.h>
#include <postprocessing/sqlite/SQLite.h>
#include <vtk/VTKOutput.h>

#include <fstream>
#include <iomanip>

using namespace walberla;
using namespace walberla::pe;
using namespace walberla::timing;

typedef boost::tuple<Box, Capsule, Sphere, Plane> BodyTuple ;

int main( int argc, char ** argv )
{
   Environment env(argc, argv);
   WALBERLA_UNUSED(env);

   std::map< std::string, int64_t >     integerProperties;
   std::map< std::string, double >      realProperties;
   std::map< std::string, std::string > stringProperties;

   WALBERLA_LOG_INFO_ON_ROOT( "config file: " << argv[1] );
   stringProperties["config"]       = argv[1];

   WALBERLA_LOG_INFO_ON_ROOT( "waLBerla Revision: " << WALBERLA_GIT_SHA1 );
   stringProperties["WALBERLA_GIT"] = WALBERLA_GIT_SHA1;

   WcTimingTree tt;

   //   logging::Logging::instance()->setStreamLogLevel(logging::Logging::INFO);
   //   logging::Logging::instance()->includeLoggingToFile("HourGlass");
   //   logging::Logging::instance()->setFileLogLevel(logging::Logging::INFO);


   WALBERLA_LOG_INFO_ON_ROOT("*** CONFIG ***");

   // configuration file
   const Config::BlockHandle cfg  = env.config()->getBlock( "HourGlass" );
   if (!cfg)
   {
      WALBERLA_LOG_WARNING_ON_ROOT("Configuration file is broken... using default values");
   }

   const std::string sqlFile             = cfg.getParameter< std::string >( "sqlFile", "HourGlass.sqlite" );

   integerProperties["numMPIProcesses"]  = mpi::MPIManager::instance()->numProcesses();

   const uint_t numberOfRotations                = cfg.getParameter<uint_t>( "numberOfRotations", uint_c(12) );
   WALBERLA_LOG_INFO_ON_ROOT("numberOfRotations: " << numberOfRotations);
   integerProperties["openingRadius"] = int64_c(numberOfRotations);
   const real_t openingRadius                    = cfg.getParameter<real_t>( "openingRadius", real_t(0.8) );
   WALBERLA_LOG_INFO_ON_ROOT("openingRadius: " << openingRadius);
   realProperties["openingRadius"] = double_c(openingRadius);

   const real_t spacing                          = cfg.getParameter<real_t>( "spacing", real_t(0.45) );
   WALBERLA_LOG_INFO_ON_ROOT("spacing: " << spacing);
   realProperties["spacing"] = double_c(spacing);
   const real_t radius                           = cfg.getParameter<real_t>( "radius", real_t(0.2) );
   WALBERLA_LOG_INFO_ON_ROOT("radius: " << radius);
   realProperties["radius"] = double_c(radius);
   const real_t vMax                             = cfg.getParameter<real_t>( "vMax", real_t(0.05) );
   WALBERLA_LOG_INFO_ON_ROOT("vMax: " << vMax);
   realProperties["vMax"] = double_c(vMax);

   const real_t dt                               = cfg.getParameter<real_t>( "dt", real_t(0.01) );
   WALBERLA_LOG_INFO_ON_ROOT("dt: " << dt);
   realProperties["dt"] = double_c(dt);
   const uint_t simulationSteps                  = cfg.getParameter<uint_t>( "simulationSteps", uint_t(30000) );
   WALBERLA_LOG_INFO_ON_ROOT("simulationSteps: " << simulationSteps);
   integerProperties["simulationSteps"] = int64_c(simulationSteps);
   const uint_t vtkInterval                      = cfg.getParameter<uint_t>( "vtkInterval", uint_t(10000) );
   WALBERLA_LOG_INFO_ON_ROOT("vtkInterval: " << vtkInterval);
   integerProperties["vtkInterval"] = int64_c(vtkInterval);
   const std::string vtkPath                        = cfg.getParameter<std::string>("vtkPath",  "vtk_out" );
   WALBERLA_LOG_INFO_ON_ROOT("vtkPath: " << vtkPath);

   const std::string outputFilename              = cfg.getParameter<std::string>( "outputFilename", "balancing.txt" );
   WALBERLA_LOG_INFO_ON_ROOT("outputFilename: " << outputFilename);

   ///////////////////
   // Customization //
   ///////////////////

   uint_t     numParticles      = uint_c(0);

   WALBERLA_LOG_INFO_ON_ROOT("*** BODYTUPLE ***");
   // initialize body type ids
   SetBodyTypeIDs<BodyTuple>::execute();

   shared_ptr<BodyStorage> globalBodyStorage = make_shared<BodyStorage>();

   const real_t   static_cof  ( 0.4 / 2 );    // Coefficient of static friction. Roughly 0.85 with high variation depending on surface roughness for low stresses. Note: pe doubles the input coefficient of friction for material-material contacts.
   const real_t   dynamic_cof ( static_cof ); // Coefficient of dynamic friction. Similar to static friction for low speed friction.
   MaterialID     material = createMaterial( "granular", real_t( 1.0 ), 0, static_cof, dynamic_cof, real_t( 0.5 ), 1, 1, 0, 0 );

   // create forest
   shared_ptr< BlockForest > forest   = createBlockForestFromConfig( cfg );
   if (!forest)
   {
      WALBERLA_LOG_INFO_ON_ROOT( "No BlockForest created ... exiting!");
      return EXIT_SUCCESS;
   }
   BlockDataID storageID              = forest->addBlockData(createStorageDataHandling<BodyTuple>(), "Storage");
   BlockDataID ccdID                  = forest->addBlockData(ccd::createHashGridsDataHandling( globalBodyStorage, storageID ), "CCD");
   BlockDataID fcdID                  = forest->addBlockData(fcd::createGenericFCDDataHandling<BodyTuple, fcd::AnalyticCollideFunctor>(), "FCD");
   const math::AABB& simulationDomain = forest->getDomain();

   //***** SETUP SOLVER
   WALBERLA_LOG_INFO_ON_ROOT("*** INTEGRATOR ***");
   cr::HCSITS cr(globalBodyStorage, forest, storageID, ccdID, fcdID, &tt);
   configure(cfg, cr);

   WALBERLA_LOG_INFO_ON_ROOT("*** VTK ***");

   //write domain decomposition to file
   auto vtkDomainOutput = vtk::createVTKOutput_DomainDecomposition( forest, "domain_decomposition", 1, vtkPath, "simulation_step" );
   auto vtkOutput       = make_shared<SphereVtkOutput>(storageID, *forest) ;
   auto vtkWriter       = vtk::createVTKOutput_PointData(vtkOutput, "Bodies", 1, vtkPath, "simulation_step", false, false);

   WALBERLA_LOG_INFO_ON_ROOT("*** SETUP - START ***");


   //calculate center
   auto longestEdge = std::max(simulationDomain.size(0), std::max(simulationDomain.size(1), simulationDomain.size(2)));
   auto center      = real_t(0.5) * (simulationDomain.minCorner() + simulationDomain.maxCorner());
   //generate hour clock
   for (uint_t i = 0; i < numberOfRotations; ++i)
   {
      BoxID box;

      box = createBox(*globalBodyStorage,
                      *forest,
                      storageID,
                      0,
                      center,
                      Vec3(longestEdge),
                      material,
                      true,
                      false,
                      true);
      box->rotate( Vec3(0,1,0), math::PI * real_t(0.25) );
      box->setPosition( center + (longestEdge * sqrt(2) * real_t(0.5) + openingRadius ) * Vec3(1,0,0) );
      box->rotateAroundPoint(center, Vec3(0,0,1), math::PI * real_c(i) / real_c(numberOfRotations));

      box = createBox(*globalBodyStorage,
                      *forest,
                      storageID,
                      0,
                      center,
                      Vec3(longestEdge),
                      material,
                      true,
                      false,
                      true);
      box->rotate( Vec3(0,1,0), math::PI * real_t(0.25) );
      box->setPosition( center - (longestEdge * sqrt(2) * real_t(0.5) + openingRadius ) * Vec3(1,0,0) );
      box->rotateAroundPoint(center, Vec3(0,0,1), math::PI * real_c(i) / real_c(numberOfRotations));
   }
   createPlane( *globalBodyStorage, 0, Vec3(0,0,1), Vec3(0,0,0), material );

   WALBERLA_MPI_BARRIER();
   WALBERLA_LOG_INFO_ON_ROOT("*** PARTICLE GENERATION ***");
   tt.start("Particle Creation");

   for (auto blkIt = forest->begin(); blkIt != forest->end(); ++blkIt)
   {
      IBlock & currentBlock = *blkIt;
      for (auto it = grid_generator::SCIterator(currentBlock.getAABB(), Vector3<real_t>(spacing) * real_t(0.5), spacing); it != grid_generator::SCIterator(); ++it)
      {
         Vec3 tmp = *it - center;
         if (tmp[2] < real_t(0) ) continue;
         if ( (tmp[0] * tmp[0] + tmp[1] * tmp[1]) > (tmp[2] * tmp[2]) ) continue;
         SphereID sp = createSphere( *globalBodyStorage, *forest, storageID, 0, *it, radius, material );
         Vec3 rndVel(math::realRandom<real_t>(-vMax, vMax), math::realRandom<real_t>(-vMax, vMax), math::realRandom<real_t>(-vMax, 0));
         if (sp != NULL) sp->setLinearVel(rndVel);
         if (sp != NULL) ++numParticles;
      }
   }
   WALBERLA_MPI_BARRIER();
   WALBERLA_LOG_INFO_ON_ROOT("*** SYNCHRONIZATION ***");
   syncNextNeighbors<BodyTuple>(*forest, storageID, &tt);
   syncNextNeighbors<BodyTuple>(*forest, storageID, &tt);
   WALBERLA_MPI_BARRIER();
   WALBERLA_LOG_INFO_ON_ROOT("*** AFTER SYNCHRONIZATION ***");
   tt.stop("Particle Creation");
   mpi::reduceInplace(numParticles, mpi::SUM);
   WALBERLA_LOG_INFO_ON_ROOT("#particles created: " << numParticles);

   WALBERLA_LOG_INFO_ON_ROOT("*** SIMULATION LOOP***");

   std::stringstream ssTimestamp;
   std::stringstream ssNumBlocks;
   std::stringstream ssNumParticles;
   std::stringstream ssNumShadowParticles;
   std::stringstream ssNumContacts;
   auto currentTimestamp = timing::WcPolicy::getTimestamp();

   WALBERLA_ROOT_SECTION()
   {
      ssTimestamp << "# timestamps" << std::endl << "# three lines for every process: numBlocks, numParticles, numShadowParticles" << std::endl;
   }

   //   vtkWriter->write( true );
   //   vtkDomainOutput->write( true );

   WcTimingPool tp;
   tt.start("Simulation Loop");
   tp["Total"].start();
   for (uint_t i = 0; i < simulationSteps; ++i)
   {
      tp["VTKOutput"].start();
      if( i % vtkInterval == 0 )
      {
         vtkWriter->write( true );
         vtkDomainOutput->write( true );
      }
      tp["VTKOutput"].end();

      tp["Solver"].start();
      cr.timestep( real_c(dt) );
      tp["Solver"].end();
      tp["Sync"].start();
      syncNextNeighbors<BodyTuple>(*forest, storageID, &tt);
      tp["Sync"].end();
   }
   tp["Total"].end();
   tt.stop("Simulation Loop");
   WALBERLA_LOG_INFO_ON_ROOT("*** SIMULATION - END ***");

   ssTimestamp << std::endl << ssNumBlocks.str()
               << std::endl << ssNumParticles.str()
               << std::endl << ssNumShadowParticles.str()
               << std::endl << ssNumContacts.str()
               << std::endl;

   mpi::writeMPITextFile(outputFilename, ssTimestamp.str());

   auto ttReduced = tt.getReduced( );
   WALBERLA_ROOT_SECTION()
   {
      std::cout << ttReduced << std::endl;
   }
   auto tpReduced = tp.getReduced();
   auto tpLoadBalancing = forest->getRefreshTiming().getReduced();
   WALBERLA_ROOT_SECTION()
   {
      std::cout << *tpReduced << std::endl;
      std::cout << *tpLoadBalancing << std::endl;
   }

   WALBERLA_ROOT_SECTION()
   {
      // Log to SQL Database
      auto runId = postprocessing::storeRunInSqliteDB( sqlFile, integerProperties, stringProperties, realProperties );
      postprocessing::storeTimingPoolInSqliteDB( sqlFile, runId, *tpReduced, "Timeloop" );
      postprocessing::storeTimingPoolInSqliteDB( sqlFile, runId, *tpLoadBalancing, "LoadBalancing" );
      postprocessing::storeTimingTreeInSqliteDB( sqlFile, runId, ttReduced, "TimingTree" );
   }

   return EXIT_SUCCESS;
}
