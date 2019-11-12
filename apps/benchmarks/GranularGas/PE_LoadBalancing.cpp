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
//! \file   PE_LoadBalancing.cpp
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include "NodeTimings.h"
#include "Parameters.h"
#include "SQLProperties.h"

#include <pe/amr/InfoCollection.h>
#include <pe/amr/level_determination/MinMaxLevelDetermination.h>
#include <pe/amr/weight_assignment/MetisAssignmentFunctor.h>
#include <pe/amr/weight_assignment/WeightAssignmentFunctor.h>
#include <pe/basic.h>
#include <pe/synchronization/ClearSynchronization.h>
#include <pe/vtk/SphereVtkOutput.h>

#include <blockforest/Initialization.h>
#include <blockforest/loadbalancing/DynamicCurve.h>
#include <blockforest/loadbalancing/DynamicParMetis.h>
#include <blockforest/loadbalancing/PODPhantomData.h>
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

   std::map< std::string, walberla::int64_t > integerProperties;
   std::map< std::string, double >            realProperties;
   std::map< std::string, std::string >       stringProperties;

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

   forest->recalculateBlockLevelsInRefresh( params.recalculateBlockLevelsInRefresh );
   forest->alwaysRebalanceInRefresh( params.alwaysRebalanceInRefresh );
   forest->reevaluateMinTargetLevelsAfterForcedRefinement( params.reevaluateMinTargetLevelsAfterForcedRefinement );
   forest->allowRefreshChangingDepth( params.allowRefreshChangingDepth );

   forest->allowMultipleRefreshCycles( params.allowMultipleRefreshCycles );
   forest->checkForEarlyOutInRefresh( params.checkForEarlyOutInRefresh );
   forest->checkForLateOutInRefresh( params.checkForLateOutInRefresh );

   WALBERLA_LOG_INFO_ON_ROOT("simulationDomain: " << forest->getDomain());

   WALBERLA_LOG_INFO_ON_ROOT("blocks: " << Vector3<uint_t>(forest->getXSize(), forest->getYSize(), forest->getZSize()) );

   auto ic = make_shared<blockforest::InfoCollection>();

   pe::amr::MinMaxLevelDetermination regrid(ic, params.regridMin, params.regridMax);
   forest->setRefreshMinTargetLevelDeterminationFunction( regrid );

   bool bRebalance = true;
   if (params.LBAlgorithm == "None")
   {
      bRebalance = false;
   } else if (params.LBAlgorithm == "Morton")
   {
      forest->setRefreshPhantomBlockDataAssignmentFunction( pe::amr::WeightAssignmentFunctor( ic, params.baseWeight ) );
      forest->setRefreshPhantomBlockDataPackFunction( pe::amr::WeightAssignmentFunctor::PhantomBlockWeightPackUnpackFunctor() );
      forest->setRefreshPhantomBlockDataUnpackFunction( pe::amr::WeightAssignmentFunctor::PhantomBlockWeightPackUnpackFunctor() );

      auto prepFunc = blockforest::DynamicCurveBalance< pe::amr::WeightAssignmentFunctor::PhantomBlockWeight >( false, true, false );
      prepFunc.setMaxBlocksPerProcess( params.maxBlocksPerProcess );
      forest->setRefreshPhantomBlockMigrationPreparationFunction( prepFunc );
   } else if (params.LBAlgorithm == "Hilbert")
   {
      forest->setRefreshPhantomBlockDataAssignmentFunction( pe::amr::WeightAssignmentFunctor( ic, params.baseWeight ) );
      forest->setRefreshPhantomBlockDataPackFunction( pe::amr::WeightAssignmentFunctor::PhantomBlockWeightPackUnpackFunctor() );
      forest->setRefreshPhantomBlockDataUnpackFunction( pe::amr::WeightAssignmentFunctor::PhantomBlockWeightPackUnpackFunctor() );

      auto prepFunc = blockforest::DynamicCurveBalance< pe::amr::WeightAssignmentFunctor::PhantomBlockWeight >( true, true, false );
      prepFunc.setMaxBlocksPerProcess( params.maxBlocksPerProcess );
      forest->setRefreshPhantomBlockMigrationPreparationFunction( prepFunc );
   } else if (params.LBAlgorithm == "Metis")
   {
      auto assFunc = pe::amr::MetisAssignmentFunctor( ic, params.baseWeight );
      forest->setRefreshPhantomBlockDataAssignmentFunction( assFunc );
      forest->setRefreshPhantomBlockDataPackFunction( pe::amr::MetisAssignmentFunctor::PhantomBlockWeightPackUnpackFunctor() );
      forest->setRefreshPhantomBlockDataUnpackFunction( pe::amr::MetisAssignmentFunctor::PhantomBlockWeightPackUnpackFunctor() );

      auto alg     = blockforest::DynamicParMetis::stringToAlgorithm(    params.metisAlgorithm );
      auto vWeight = blockforest::DynamicParMetis::stringToWeightsToUse( params.metisWeightsToUse );
      auto eWeight = blockforest::DynamicParMetis::stringToEdgeSource(   params.metisEdgeSource );

      auto prepFunc = blockforest::DynamicParMetis( alg, vWeight, eWeight );
      prepFunc.setipc2redist(params.metisipc2redist);
      mesa_pd::addParMetisPropertiesToSQL(prepFunc, integerProperties, realProperties, stringProperties);
      forest->setRefreshPhantomBlockMigrationPreparationFunction( prepFunc );
   } else if (params.LBAlgorithm == "Diffusive")
   {
      forest->setRefreshPhantomBlockDataAssignmentFunction( pe::amr::WeightAssignmentFunctor( ic, params.baseWeight ) );
      forest->setRefreshPhantomBlockDataPackFunction( pe::amr::WeightAssignmentFunctor::PhantomBlockWeightPackUnpackFunctor() );
      forest->setRefreshPhantomBlockDataUnpackFunction( pe::amr::WeightAssignmentFunctor::PhantomBlockWeightPackUnpackFunctor() );
      auto prepFunc = blockforest::DynamicDiffusionBalance< pe::amr::WeightAssignmentFunctor::PhantomBlockWeight >( 1, 1, false );
      //configure(cfg, prepFunc);
      //addDynamicDiffusivePropertiesToSQL(prepFunc, integerProperties, realProperties, stringProperties);
      forest->setRefreshPhantomBlockMigrationPreparationFunction(prepFunc);
   } else
   {
      WALBERLA_ABORT("Unknown LBAlgorithm: " << params.LBAlgorithm);
   }

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

   auto center = forest->getDomain().center();
   for (auto& currentBlock : *forest)
   {
      for (auto it = grid_generator::SCIterator(currentBlock.getAABB().getIntersection(generationDomain),
                                                Vector3<real_t>(params.spacing) * real_c(0.5) + params.shift,
                                                params.spacing);
           it != grid_generator::SCIterator();
           ++it)
      {
         auto tmp = dot( (*it - center), params.normal );
         if (tmp < 0)
         {
            SphereID sp = pe::createSphere( *globalBodyStorage, *forest, storageID, 0, *it, params.radius, material);
            if (sp != nullptr) ++numParticles;
         }
      }
   }
   mpi::reduceInplace(numParticles, mpi::SUM);
   WALBERLA_LOG_INFO_ON_ROOT("#particles created: " << numParticles);

   WALBERLA_LOG_INFO_ON_ROOT("*** SETUP - END ***");

   // synchronize particles
   syncCallWithoutTT();
   syncCallWithoutTT();

   WALBERLA_LOG_INFO_ON_ROOT("*** SIMULATION - START ***");

   WcTimer      timerImbalanced;
   WcTimer      timerLoadBalancing;
   WcTimer      timerBalanced;
   WcTimingPool tpImbalanced;
   WcTimingPool tpBalanced;
   WALBERLA_MPI_BARRIER();
   timerImbalanced.start();
   for (int64_t i=0; i < params.simulationSteps; ++i)
   {
      if( i % 200 == 0 )
      {
         WALBERLA_LOG_DEVEL_ON_ROOT( "Timestep " << i << " / " << params.simulationSteps );
      }

      tpImbalanced["CR"].start();
      cr->timestep( real_c(params.dt) );
      tpImbalanced["CR"].end();
      tpImbalanced["Sync"].start();
      syncCallWithoutTT();
      tpImbalanced["Sync"].end();

      //      if( i % params.visSpacing == 0 )
      //      {
      //         vtkDomainOutput->write( );
      //         vtkSphereOutput->write( );
      //      }
   }
   timerImbalanced.end();

   if (bRebalance)
   {
      vtkDomainOutput->write( );
      vtkSphereOutput->write( );
      WALBERLA_MPI_BARRIER();
      timerLoadBalancing.start();
      WALBERLA_LOG_INFO_ON_ROOT("*** Rebalance ***");
      createWithNeighborhoodLocalShadow( *forest, storageID, *ic );
      clearSynchronization( *forest, storageID );
      forest->refresh();
      integerProperties["MigrationIterations1"] = int64_c(forest->phantomBlockMigrationIterations());
      syncNextNeighbors<BodyTuple>(*forest, storageID);
      for (auto blockIt = forest->begin(); blockIt != forest->end(); ++blockIt)
      {
         ccd::ICCD* ccd = blockIt->getData< ccd::ICCD >( ccdID );
         ccd->reloadBodies();
      }
      timerLoadBalancing.end();
      vtkDomainOutput->write( );
      vtkSphereOutput->write( );
   }

   WALBERLA_MPI_BARRIER();
   timerBalanced.start();
   for (int64_t i=0; i < params.simulationSteps; ++i)
   {
      if( i % 200 == 0 )
      {
         WALBERLA_LOG_DEVEL_ON_ROOT( "Timestep " << i << " / " << params.simulationSteps );
      }

      tpBalanced["CR"].start();
      cr->timestep( real_c(params.dt) );
      tpBalanced["CR"].end();
      tpBalanced["Sync"].start();
      syncCallWithoutTT();
      tpBalanced["Sync"].end();

      //      if( i % params.visSpacing == 0 )
      //      {
      //         vtkDomainOutput->write( );
      //         vtkSphereOutput->write( );
      //      }
   }
   timerBalanced.end();

   auto timerImbalancedReduced = walberla::timing::getReduced(timerImbalanced, REDUCE_TOTAL, 0);
   double PUpSImbalanced = 0.0;
   WALBERLA_ROOT_SECTION()
   {
      WALBERLA_LOG_INFO_ON_ROOT("IMBALANCED " << *timerImbalancedReduced);
      PUpSImbalanced = double_c(numParticles) * double_c(params.simulationSteps) / double_c(timerImbalancedReduced->max());
      WALBERLA_LOG_INFO_ON_ROOT("PUpS: " << PUpSImbalanced);
   }

   auto timerBalancedReduced = walberla::timing::getReduced(timerBalanced, REDUCE_TOTAL, 0);
   double PUpSBalanced = 0.0;
   WALBERLA_ROOT_SECTION()
   {
      WALBERLA_LOG_INFO_ON_ROOT("BALANCED " << *timerBalancedReduced);
      PUpSBalanced = double_c(numParticles) * double_c(params.simulationSteps) / double_c(timerBalancedReduced->max());
      WALBERLA_LOG_INFO_ON_ROOT("PUpS: " << PUpSBalanced);
   }

   auto timerLoadBalancingReduced = walberla::timing::getReduced(timerLoadBalancing, REDUCE_TOTAL, 0);

   auto tpImbalancedReduced = tpImbalanced.getReduced();
   WALBERLA_LOG_INFO_ON_ROOT(*tpImbalancedReduced);

   auto tpBalancedReduced = tpBalanced.getReduced();
   WALBERLA_LOG_INFO_ON_ROOT(*tpBalancedReduced);
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
   }
   auto minParticles = mpi::reduce(numParticles, mpi::MIN);
   auto maxParticles = mpi::reduce(numParticles, mpi::MAX);
   WALBERLA_LOG_DEVEL_ON_ROOT("particle ratio: " << minParticles << " / " << maxParticles);

   mpi::reduceInplace(numParticles, mpi::SUM);
   mpi::reduceInplace(numGhostParticles, mpi::SUM);
   WALBERLA_LOG_INFO_ON_ROOT("*** CHECKING RESULT - END ***");

   uint_t runId = uint_c(-1);
   WALBERLA_ROOT_SECTION()
   {
      stringProperties["walberla_git"]         = WALBERLA_GIT_SHA1;
      stringProperties["tag"]                  = "pe";
      integerProperties["bDEM"]                = bDEM;
      integerProperties["bNN"]                 = bNN;
      integerProperties["mpi_num_processes"]   = mpi::MPIManager::instance()->numProcesses();
      integerProperties["omp_max_threads"]     = omp_get_max_threads();
      realProperties["imbalanced_PUpS"]          = double_c(PUpSImbalanced);
      realProperties["imbalanced_timer_min"]     = timerImbalancedReduced->min();
      realProperties["imbalanced_timer_max"]     = timerImbalancedReduced->max();
      realProperties["imbalanced_timer_average"] = timerImbalancedReduced->average();
      realProperties["imbalanced_timer_total"]   = timerImbalancedReduced->total();
      realProperties["loadbalancing_timer_min"]     = timerLoadBalancingReduced->min();
      realProperties["loadbalancing_timer_max"]     = timerLoadBalancingReduced->max();
      realProperties["loadbalancing_timer_average"] = timerLoadBalancingReduced->average();
      realProperties["loadbalancing_timer_total"]   = timerLoadBalancingReduced->total();
      realProperties["balanced_PUpS"]          = double_c(PUpSBalanced);
      realProperties["balanced_timer_min"]     = timerBalancedReduced->min();
      realProperties["balanced_timer_max"]     = timerBalancedReduced->max();
      realProperties["balanced_timer_average"] = timerBalancedReduced->average();
      realProperties["balanced_timer_total"]   = timerBalancedReduced->total();
      integerProperties["num_particles"]       = numParticles;
      integerProperties["num_ghost_particles"] = numGhostParticles;
      integerProperties["minParticles"]        = minParticles;
      integerProperties["maxParticles"]        = maxParticles;

      mesa_pd::addBuildInfoToSQL( integerProperties, realProperties, stringProperties );
      saveToSQL(params, integerProperties, realProperties, stringProperties );
      mesa_pd::addDomainPropertiesToSQL(*forest, integerProperties, realProperties, stringProperties);
      mesa_pd::addSlurmPropertiesToSQL(integerProperties, realProperties, stringProperties);

      runId = sqlite::storeRunInSqliteDB( params.sqlFile, integerProperties, stringProperties, realProperties );
      sqlite::storeTimingPoolInSqliteDB( params.sqlFile, runId, *tpImbalancedReduced, "imbalanced" );
      sqlite::storeTimingPoolInSqliteDB( params.sqlFile, runId, *tpImbalancedReduced, "balanced" );
   }
   if (params.storeNodeTimings)
   {
      mesa_pd::storeNodeTimings(runId, params.sqlFile, "NodeTimingImbalanced", tpImbalanced);
      mesa_pd::storeNodeTimings(runId, params.sqlFile, "NodeTimingBalanced", tpBalanced);
   }
   WALBERLA_LOG_INFO_ON_ROOT("*** SQL OUTPUT - END ***");

   return EXIT_SUCCESS;
}
} // namespace walberla

int main( int argc, char* argv[] )
{
   return walberla::main( argc, argv );
}
