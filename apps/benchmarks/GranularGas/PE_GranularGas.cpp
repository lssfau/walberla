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
//! \file   PE_GranularGas.cpp
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include "NodeTimings.h"
#include "Parameters.h"
#include "SQLProperties.h"

#include <pe/basic.h>
#include <pe/vtk/SphereVtkOutput.h>

#include <blockforest/Initialization.h>
#include <core/Abort.h>
#include <core/Environment.h>
#include <core/math/Random.h>
#include <core/grid_generator/SCIterator.h>
#include <core/logging/Logging.h>
#include <core/OpenMP.h>
#include <core/timing/TimingTree.h>
#include <core/waLBerlaBuildInfo.h>
#include <sqlite/SQLite.h>
#include <vtk/VTKOutput.h>

#include <functional>
#include <memory>
#include <tuple>

namespace walberla {
using namespace walberla::pe;
using namespace walberla::timing;

using BodyTuple = std::tuple<Sphere, Plane> ;

int main( int argc, char ** argv )
{
   WcTimingTree tt;
   Environment env(argc, argv);

   logging::Logging::instance()->setStreamLogLevel(logging::Logging::INFO);
   logging::Logging::instance()->setFileLogLevel(logging::Logging::INFO);

   WALBERLA_LOG_INFO_ON_ROOT( "config file: " << argv[1] )
   WALBERLA_LOG_INFO_ON_ROOT( "waLBerla Revision: " << WALBERLA_GIT_SHA1 );

   math::seedRandomGenerator( static_cast<unsigned int>(1337 * mpi::MPIManager::instance()->worldRank()) );

   WALBERLA_LOG_INFO_ON_ROOT("*** READING COMMANDLINE ARGUMENTS ***");
   bool bDEM = false;
   bool bHCSITS = false;

   bool bNN = false;
   bool bSO = false;

   bool bInelasticFrictionlessContact = false;
   bool bApproximateInelasticCoulombContactByDecoupling = false;
   bool bInelasticCoulombContactByDecoupling = false;
   bool bInelasticGeneralizedMaximumDissipationContact = false;

   for( int i = 1; i < argc; ++i )
   {
      if( std::strcmp( argv[i], "--DEM" ) == 0 ) bDEM = true;
      if( std::strcmp( argv[i], "--HCSITS" ) == 0 ) bHCSITS = true;

      if( std::strcmp( argv[i], "--syncNextNeighbor" ) == 0 ) bNN = true;
      if( std::strcmp( argv[i], "--syncShadowOwners" ) == 0 ) bSO = true;

      if( std::strcmp( argv[i], "--InelasticFrictionlessContact" ) == 0 ) bInelasticFrictionlessContact = true;
      if( std::strcmp( argv[i], "--ApproximateInelasticCoulombContactByDecoupling" ) == 0 ) bApproximateInelasticCoulombContactByDecoupling = true;
      if( std::strcmp( argv[i], "--InelasticCoulombContactByDecoupling" ) == 0 ) bInelasticCoulombContactByDecoupling = true;
      if( std::strcmp( argv[i], "--InelasticGeneralizedMaximumDissipationContact" ) == 0 ) bInelasticGeneralizedMaximumDissipationContact = true;
   }

   WALBERLA_LOG_INFO_ON_ROOT("*** READING CONFIG FILE ***");
   auto cfg = env.config();
   if (cfg == nullptr) WALBERLA_ABORT("No config specified!");
   const Config::BlockHandle mainConf  = cfg->getBlock( "GranularGas" );
   mesa_pd::Parameters params;
   loadFromConfig(params, mainConf);

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

   WALBERLA_LOG_INFO_ON_ROOT("blocks: " << Vector3<uint_t>(forest->getXSize(), forest->getYSize(), forest->getZSize()) );

   WALBERLA_LOG_INFO_ON_ROOT("*** BODYTUPLE ***");
   // initialize body type ids
   SetBodyTypeIDs<BodyTuple>::execute();

   WALBERLA_LOG_INFO_ON_ROOT("*** STORAGEDATAHANDLING ***");
   // add block data
   auto storageID           = forest->addBlockData(createStorageDataHandling<BodyTuple>(), "Storage");
   auto ccdID               = forest->addBlockData(ccd::createHashGridsDataHandling( globalBodyStorage, storageID ), "CCD");
   auto fcdID               = forest->addBlockData(fcd::createGenericFCDDataHandling<BodyTuple, fcd::AnalyticCollideFunctor>(), "FCD");

   WALBERLA_LOG_INFO_ON_ROOT("*** INTEGRATOR ***");
   std::unique_ptr<cr::ICR> cr;
   if (bDEM)
   {
      cr = std::make_unique<cr::DEM>(globalBodyStorage, forest, storageID, ccdID, fcdID, &tt);
      WALBERLA_LOG_INFO_ON_ROOT("Using DEM!");
   } else if (bHCSITS)
   {
      cr = std::make_unique<cr::HCSITS>(globalBodyStorage, forest, storageID, ccdID, fcdID, &tt);
      configure(mainConf, *static_cast<cr::HCSITS*>(cr.get()));
      WALBERLA_LOG_INFO_ON_ROOT("Using HCSITS!");

      cr::HCSITS* hcsits = static_cast<cr::HCSITS*>(cr.get());

      if (bInelasticFrictionlessContact)
      {
         hcsits->setRelaxationModel(cr::HCSITS::InelasticFrictionlessContact);
         WALBERLA_LOG_INFO_ON_ROOT("Using InelasticFrictionlessContact!");
      } else if (bApproximateInelasticCoulombContactByDecoupling)
      {
         hcsits->setRelaxationModel(cr::HCSITS::ApproximateInelasticCoulombContactByDecoupling);
         WALBERLA_LOG_INFO_ON_ROOT("Using ApproximateInelasticCoulombContactByDecoupling!");
      } else if (bInelasticCoulombContactByDecoupling)
      {
         hcsits->setRelaxationModel(cr::HCSITS::InelasticCoulombContactByDecoupling);
         WALBERLA_LOG_INFO_ON_ROOT("Using InelasticCoulombContactByDecoupling!");
      } else if (bInelasticGeneralizedMaximumDissipationContact)
      {
         hcsits->setRelaxationModel(cr::HCSITS::InelasticGeneralizedMaximumDissipationContact);
         WALBERLA_LOG_INFO_ON_ROOT("Using InelasticGeneralizedMaximumDissipationContact!");
      } else
      {
         WALBERLA_ABORT("Friction model could not be determined!");
      }
   } else
   {
      WALBERLA_ABORT("Model could not be determined!");
   }

   WALBERLA_LOG_INFO_ON_ROOT("*** SYNCCALL ***");
   std::function<void(void)> syncCallWithoutTT;
   if (bNN)
   {
      syncCallWithoutTT = std::bind( pe::syncNextNeighbors<BodyTuple>, std::ref(*forest), storageID, &tt, real_c(0.1), false );
      WALBERLA_LOG_INFO_ON_ROOT("Using NextNeighbor sync!");
   } else if (bSO)
   {
      syncCallWithoutTT = std::bind( pe::syncShadowOwners<BodyTuple>, std::ref(*forest), storageID, &tt, real_c(0.1), false );
      WALBERLA_LOG_INFO_ON_ROOT("Using ShadowOwner sync!");
   } else
   {
      WALBERLA_ABORT("Synchronization method could not be determined!");
   }

   WALBERLA_LOG_INFO_ON_ROOT("*** VTK ***");
   auto vtkDomainOutput = vtk::createVTKOutput_DomainDecomposition( forest, "domain_decomposition", 1, params.vtk_out, "simulation_step" );
   auto vtkSphereHelper = make_shared<SphereVtkOutput>(storageID, *forest) ;
   auto vtkSphereOutput = vtk::createVTKOutput_PointData(vtkSphereHelper, "Bodies", 1, params.vtk_out, "simulation_step", false, false);

   WALBERLA_LOG_INFO_ON_ROOT("*** SETUP - START ***");
   //const real_t   static_cof  ( real_c(0.1) / 2 );   // Coefficient of static friction. Note: pe doubles the input coefficient of friction for material-material contacts.
   //const real_t   dynamic_cof ( static_cof ); // Coefficient of dynamic friction. Similar to static friction for low speed friction.
   MaterialID     material = createMaterial( "granular", real_t( 1.0 ), 0, 0, 0, real_t( 0.5 ), 1, real_t(1e-6), 0, 0 );

   auto simulationDomain = forest->getDomain();
   const auto& generationDomain = simulationDomain; // simulationDomain.getExtended(-real_c(0.5) * spacing);
   int64_t numParticles = 0;

   for (auto& currentBlock : *forest)
   {
      for (auto it = grid_generator::SCIterator(currentBlock.getAABB().getIntersection(generationDomain),
                                                Vector3<real_t>(params.spacing) * real_c(0.5) + params.shift,
                                                params.spacing);
           it != grid_generator::SCIterator();
           ++it)
      {
         SphereID sp = pe::createSphere( *globalBodyStorage, *forest, storageID, 0, *it, params.radius, material);
         if (sp != nullptr) ++numParticles;
      }
   }
   mpi::reduceInplace(numParticles, mpi::SUM);
   WALBERLA_LOG_INFO_ON_ROOT("#particles created: " << numParticles);

   if (!forest->isPeriodic(0))
   {
      createPlane(*globalBodyStorage, 0, Vec3(+1,0,0), forest->getDomain().minCorner());
      createPlane(*globalBodyStorage, 0, Vec3(-1,0,0), forest->getDomain().maxCorner());
   }
   if (!forest->isPeriodic(1))
   {
      createPlane(*globalBodyStorage, 0, Vec3(0,+1,0), forest->getDomain().minCorner());
      createPlane(*globalBodyStorage, 0, Vec3(0,-1,0), forest->getDomain().maxCorner());
   }
   if (!forest->isPeriodic(2))
   {
      createPlane(*globalBodyStorage, 0, Vec3(0,0,+1), forest->getDomain().minCorner());
      createPlane(*globalBodyStorage, 0, Vec3(0,0,-1), forest->getDomain().maxCorner());
   }

   WALBERLA_LOG_INFO_ON_ROOT("*** SETUP - END ***");

   // synchronize particles
   syncCallWithoutTT();
   syncCallWithoutTT();

   WALBERLA_LOG_INFO_ON_ROOT("*** SIMULATION - START ***");
   for (int64_t outerIteration = 0; outerIteration < params.numOuterIterations; ++outerIteration)
   {
      WALBERLA_LOG_INFO_ON_ROOT("*** RUNNING OUTER ITERATION " << outerIteration << " ***");
      WcTimer timer;
      WcTimingPool tp;
      WALBERLA_MPI_BARRIER();
      timer.start();
      for (int64_t i=0; i < params.simulationSteps; ++i)
      {
         if( i % 200 == 0 )
         {
            WALBERLA_LOG_DEVEL_ON_ROOT( "Timestep " << i << " / " << params.simulationSteps );
         }

         tp["CR"].start();
         cr->timestep( real_c(params.dt) );
         tp["CR"].end();
         tp["Sync"].start();
         syncCallWithoutTT();
         tp["Sync"].end();

         //if( i % visSpacing == 0 )
         //{
         //   vtkDomainOutput->write( );
         //   vtkSphereOutput->write( );
         //}
      }
      timer.end();
      auto timer_reduced = walberla::timing::getReduced(timer, REDUCE_TOTAL, 0);
      double PUpS = 0.0;
      WALBERLA_ROOT_SECTION()
      {
         WALBERLA_LOG_INFO_ON_ROOT(*timer_reduced);
         WALBERLA_LOG_INFO_ON_ROOT("runtime: " << timer_reduced->max());
         PUpS = double_c(numParticles) * double_c(params.simulationSteps) / double_c(timer_reduced->max());
         WALBERLA_LOG_INFO_ON_ROOT("PUpS: " << PUpS);
      }

      auto tp_reduced = tp.getReduced();
      WALBERLA_LOG_INFO_ON_ROOT(*tp_reduced);
      WALBERLA_LOG_INFO_ON_ROOT("*** SIMULATION - END ***");

      auto temp = tt.getReduced( );
      WALBERLA_ROOT_SECTION()
      {
         std::cout << temp;
      }

      WALBERLA_LOG_INFO_ON_ROOT("*** CHECKING RESULT - START ***");
      numParticles = 0;
      int64_t numGhostParticles = 0;
      for (auto& currentBlock : *forest)
      {
         Storage * storage = currentBlock.getData< Storage >( storageID );
         BodyStorage& localStorage = (*storage)[0];
         BodyStorage& shadowStorage = (*storage)[1];
         numParticles += localStorage.size();
         numGhostParticles += shadowStorage.size();

         auto bodyIt = localStorage.begin();
         for (auto it = grid_generator::SCIterator(currentBlock.getAABB().getIntersection(generationDomain),
                                                   Vector3<real_t>(params.spacing) * real_c(0.5) + params.shift,
                                                   params.spacing);
              it != grid_generator::SCIterator();
              ++it, ++bodyIt)
         {
            WALBERLA_CHECK_UNEQUAL(bodyIt, localStorage.end());
            WALBERLA_CHECK_FLOAT_EQUAL(bodyIt->getPosition(), *it);
         }
      }
      mpi::reduceInplace(numParticles, mpi::SUM);
      mpi::reduceInplace(numGhostParticles, mpi::SUM);
      WALBERLA_LOG_INFO_ON_ROOT("*** CHECKING RESULT - END ***");

      uint_t runId = uint_c(-1);
      WALBERLA_ROOT_SECTION()
      {
         std::map< std::string, walberla::int64_t > integerProperties;
         std::map< std::string, double >            realProperties;
         std::map< std::string, std::string >       stringProperties;

         stringProperties["walberla_git"]         = WALBERLA_GIT_SHA1;
         stringProperties["tag"]                  = "pe";
         integerProperties["bDEM"]                = bDEM;
         integerProperties["bNN"]                 = bNN;
         integerProperties["mpi_num_processes"]   = mpi::MPIManager::instance()->numProcesses();
         integerProperties["omp_max_threads"]     = omp_get_max_threads();
         realProperties["PUpS"]                   = double_c(PUpS);
         realProperties["timer_min"]              = timer_reduced->min();
         realProperties["timer_max"]              = timer_reduced->max();
         realProperties["timer_average"]          = timer_reduced->average();
         realProperties["timer_total"]            = timer_reduced->total();
         integerProperties["num_particles"]       = numParticles;
         integerProperties["num_ghost_particles"] = numGhostParticles;

         mesa_pd::addBuildInfoToSQL( integerProperties, realProperties, stringProperties );
         saveToSQL(params, integerProperties, realProperties, stringProperties );
         mesa_pd::addDomainPropertiesToSQL(*forest, integerProperties, realProperties, stringProperties);
         mesa_pd::addSlurmPropertiesToSQL(integerProperties, realProperties, stringProperties);

         runId = sqlite::storeRunInSqliteDB( params.sqlFile, integerProperties, stringProperties, realProperties );
         sqlite::storeTimingPoolInSqliteDB( params.sqlFile, runId, *tp_reduced, "Timeloop" );
      }
      if (params.storeNodeTimings)
      {
         mesa_pd::storeNodeTimings(runId, params.sqlFile, "NodeTiming", tp);
      }
      WALBERLA_LOG_INFO_ON_ROOT("*** SQL OUTPUT - END ***");
   }

   return EXIT_SUCCESS;
}
} // namespace walberla

int main( int argc, char* argv[] )
{
   return walberla::main( argc, argv );
}
