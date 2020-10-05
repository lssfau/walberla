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
//! \file   MESA_PD_LoadBalancing.cpp
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include "Accessor.h"
#include "check.h"
#include "Contact.h"
#include "CreateParticles.h"
#include "NodeTimings.h"
#include "Parameters.h"
#include "SelectProperty.h"
#include "sortParticleStorage.h"
#include "SQLProperties.h"

#include <mesa_pd/vtk/ParticleVtkOutput.h>

#include <mesa_pd/collision_detection/AnalyticContactDetection.h>
#include <mesa_pd/data/SparseLinkedCells.h>
#include <mesa_pd/data/ParticleAccessor.h>
#include <mesa_pd/data/ParticleStorage.h>
#include <mesa_pd/data/ShapeStorage.h>
#include <mesa_pd/data/STLOverloads.h>
#include <mesa_pd/domain/BlockForestDataHandling.h>
#include <mesa_pd/domain/BlockForestDomain.h>
#include <mesa_pd/domain/InfoCollection.h>
#include <mesa_pd/kernel/AssocToBlock.h>
#include <mesa_pd/kernel/DoubleCast.h>
#include <mesa_pd/kernel/ExplicitEuler.h>
#include <mesa_pd/kernel/InsertParticleIntoSparseLinkedCells.h>
#include <mesa_pd/kernel/ParticleSelector.h>
#include <mesa_pd/kernel/SpringDashpot.h>
#include <mesa_pd/mpi/ContactFilter.h>
#include <mesa_pd/mpi/ReduceProperty.h>
#include <mesa_pd/mpi/SyncNextNeighbors.h>
#include <mesa_pd/mpi/SyncNextNeighborsBlockForest.h>

#include <mesa_pd/mpi/notifications/ForceTorqueNotification.h>

#include <blockforest/BlockForest.h>
#include <blockforest/Initialization.h>
#include <blockforest/loadbalancing/DynamicCurve.h>
#include <blockforest/loadbalancing/DynamicParMetis.h>
#include <blockforest/loadbalancing/InfoCollection.h>
#include <blockforest/loadbalancing/PODPhantomData.h>
#include <core/Abort.h>
#include <core/Environment.h>
#include <core/math/Random.h>
#include <core/mpi/Reduce.h>
#include <core/mpi/MPITextFile.h>
#include <core/grid_generator/SCIterator.h>
#include <core/logging/Logging.h>
#include <core/OpenMP.h>
#include <core/timing/Timer.h>
#include <core/waLBerlaBuildInfo.h>
#include <pe/amr/level_determination/MinMaxLevelDetermination.h>
#include <pe/amr/weight_assignment/MetisAssignmentFunctor.h>
#include <pe/amr/weight_assignment/WeightAssignmentFunctor.h>
#include <sqlite/SQLite.h>
#include <vtk/VTKOutput.h>

#include <functional>
#include <memory>
#include <string>
#include <type_traits>

namespace walberla {
namespace mesa_pd {

template <typename PhantomBlockWeight_T>
void configure( const walberla::Config::BlockHandle& config, walberla::blockforest::DynamicDiffusionBalance<PhantomBlockWeight_T>& ddb )
{
   using namespace walberla;

   ddb.setMaxIterations  ( config.getParameter<uint_t>( "diffMaxIterations", uint_c(20)) );
   WALBERLA_LOG_INFO_ON_ROOT( "diffMaxIterations: " << ddb.getMaxIterations() );

   ddb.setFlowIterations ( config.getParameter<uint_t>( "diffFlowIterations", uint_c(12)) );
   WALBERLA_LOG_INFO_ON_ROOT( "diffFlowIterations: " << ddb.getFlowIterations() );

   ddb.checkForEarlyAbort( config.getParameter<bool>("bDiffAbortEarly", true ) );
   WALBERLA_LOG_INFO_ON_ROOT("bDiffAbortEarly: " << ddb.checkForEarlyAbort());

   ddb.adaptInflowWithGlobalInformation( config.getParameter<bool>("bDiffAdaptInflow", true ) );
   WALBERLA_LOG_INFO_ON_ROOT("bDiffAdaptInflow: " << ddb.adaptInflowWithGlobalInformation());

   ddb.adaptOutflowWithGlobalInformation( config.getParameter<bool>("bDiffAdaptOutflow", true ) );
   WALBERLA_LOG_INFO_ON_ROOT("bDiffAdaptOutflow: " << ddb.adaptOutflowWithGlobalInformation());

   std::string diffModeStr = config.getParameter<std::string>("diffMode", "push" );
   if (diffModeStr == "push")
   {
      ddb.setMode( blockforest::DynamicDiffusionBalance<PhantomBlockWeight_T>::DIFFUSION_PUSH );
   } else if (diffModeStr == "pull")
   {
      ddb.setMode( blockforest::DynamicDiffusionBalance<PhantomBlockWeight_T>::DIFFUSION_PULL );
   } else if (diffModeStr == "pushpull")
   {
      ddb.setMode( blockforest::DynamicDiffusionBalance<PhantomBlockWeight_T>::DIFFUSION_PUSHPULL );
   } else
   {
      WALBERLA_ABORT("Unknown Diffusion Mode: " << diffModeStr);
   }
   diffModeStr = "unknown";
   if (ddb.getMode() == blockforest::DynamicDiffusionBalance<PhantomBlockWeight_T>::DIFFUSION_PUSH) diffModeStr = "push";
   if (ddb.getMode() == blockforest::DynamicDiffusionBalance<PhantomBlockWeight_T>::DIFFUSION_PULL) diffModeStr = "pull";
   if (ddb.getMode() == blockforest::DynamicDiffusionBalance<PhantomBlockWeight_T>::DIFFUSION_PUSHPULL) diffModeStr = "pushpull";
   WALBERLA_LOG_INFO_ON_ROOT("diffMode: " << diffModeStr);
}


int main( int argc, char ** argv )
{
   using namespace walberla::timing;

   Environment env(argc, argv);
   auto mpiManager = walberla::mpi::MPIManager::instance();
   mpiManager->useWorldComm();

   WALBERLA_LOG_DEVEL_ON_ROOT("MESA_PD_LoadBalancing" );

   //   logging::Logging::instance()->setStreamLogLevel(logging::Logging::INFO);
   //   logging::Logging::instance()->includeLoggingToFile("LoadBalancing");
   //   logging::Logging::instance()->setFileLogLevel(logging::Logging::DETAIL);

   WALBERLA_LOG_INFO_ON_ROOT( "config file: " << argv[1] );
   WALBERLA_LOG_INFO_ON_ROOT( "waLBerla Revision: " << WALBERLA_GIT_SHA1 );

   math::seedRandomGenerator( static_cast<unsigned int>(1337 * mpiManager->worldRank()) );

   std::map< std::string, walberla::int64_t > integerProperties;
   std::map< std::string, double >            realProperties;
   std::map< std::string, std::string >       stringProperties;

   WALBERLA_LOG_INFO_ON_ROOT("*** READING CONFIG FILE ***");
   auto cfg = env.config();
   if (cfg == nullptr) WALBERLA_ABORT("No config specified!");
   const Config::BlockHandle mainConf  = cfg->getBlock( "GranularGas" );

   Parameters params;
   loadFromConfig(params, mainConf);

   WALBERLA_LOG_INFO_ON_ROOT("*** BLOCKFOREST ***");
   // create forest
   auto forest = blockforest::createBlockForestFromConfig( mainConf );
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
      addParMetisPropertiesToSQL(prepFunc, integerProperties, realProperties, stringProperties);
      forest->setRefreshPhantomBlockMigrationPreparationFunction( prepFunc );
   } else if (params.LBAlgorithm == "Diffusive")
   {
      forest->setRefreshPhantomBlockDataAssignmentFunction( pe::amr::WeightAssignmentFunctor( ic, params.baseWeight ) );
      forest->setRefreshPhantomBlockDataPackFunction( pe::amr::WeightAssignmentFunctor::PhantomBlockWeightPackUnpackFunctor() );
      forest->setRefreshPhantomBlockDataUnpackFunction( pe::amr::WeightAssignmentFunctor::PhantomBlockWeightPackUnpackFunctor() );
      auto prepFunc = blockforest::DynamicDiffusionBalance< pe::amr::WeightAssignmentFunctor::PhantomBlockWeight >( 1, 1, false );
      configure(mainConf, prepFunc);
      addDynamicDiffusivePropertiesToSQL(prepFunc, integerProperties, realProperties, stringProperties);
      forest->setRefreshPhantomBlockMigrationPreparationFunction(prepFunc);
   } else
   {
      WALBERLA_ABORT("Unknown LBAlgorithm: " << params.LBAlgorithm);
   }

   auto domain = std::make_shared<domain::BlockForestDomain>(forest);

   WALBERLA_LOG_INFO_ON_ROOT("*** SETUP - START ***");

   //init data structures
   auto ps = std::make_shared<data::ParticleStorage>(100);
   auto ss = std::make_shared<data::ShapeStorage>();
   ParticleAccessorWithShape accessor(ps, ss);
   auto lc = std::make_shared<data::SparseLinkedCells>(domain->getUnionOfLocalAABBs().getExtended(params.spacing), params.spacing );
   forest->addBlockData(domain::createBlockForestDataHandling(ps), "Storage");

   auto center = forest->getDomain().center();
   auto  smallSphere = ss->create<data::Sphere>( params.radius );
   ss->shapes[smallSphere]->updateMassAndInertia(real_t(2707));
   for (auto& iBlk : *forest)
   {
      for (auto pt : grid_generator::SCGrid(iBlk.getAABB(),
                                            Vector3<real_t>(params.spacing) * real_c(0.5) + params.shift,
                                            params.spacing))
      {
         WALBERLA_CHECK(iBlk.getAABB().contains(pt));
         auto tmp = dot( (pt - center), params.normal );
         if (tmp < 0)
         {
            createSphere(*ps, pt, params.radius, smallSphere);
         }
      }
   }
   int64_t numParticles = int64_c(ps->size());
   walberla::mpi::reduceInplace(numParticles, walberla::mpi::SUM);
   WALBERLA_LOG_INFO_ON_ROOT("#particles created: " << numParticles);

   WALBERLA_LOG_INFO_ON_ROOT("*** SETUP - END ***");

   WALBERLA_LOG_INFO_ON_ROOT("*** VTK ***");
   auto vtkDomainOutput = walberla::vtk::createVTKOutput_DomainDecomposition( forest, "domain_decomposition", 1, params.vtk_out, "simulation_step" );
   auto vtkOutput       = make_shared<mesa_pd::vtk::ParticleVtkOutput>(ps) ;
   auto vtkWriter       = walberla::vtk::createVTKOutput_PointData(vtkOutput, "Bodies", 1, params.vtk_out, "simulation_step", false, false);
   vtkOutput->addOutput<SelectRank>("rank");
   vtkOutput->addOutput<data::SelectParticleOwner>("owner");

   WALBERLA_LOG_INFO_ON_ROOT("*** SIMULATION - START ***");
   // Init kernels
   kernel::ExplicitEuler                 explicitEuler( params.dt );
   kernel::InsertParticleIntoSparseLinkedCells ipilc;
   kernel::SpringDashpot                 dem(1);
   dem.setStiffness(0, 0, real_t(0));
   dem.setDampingN (0, 0, real_t(0));
   dem.setDampingT (0, 0, real_t(0));
   dem.setFriction (0, 0, real_t(0));
   collision_detection::AnalyticContactDetection              acd;
   kernel::AssocToBlock                  assoc(forest);
   kernel::DoubleCast                    double_cast;
   mpi::ContactFilter                    contact_filter;
   mpi::ReduceProperty                   RP;
   mpi::SyncNextNeighborsBlockForest     SNN;

   // initial sync
   ps->forEachParticle(true, kernel::SelectLocal(), accessor, assoc, accessor);
   SNN(*ps, forest, domain);
   sortParticleStorage(*ps, params.sorting, lc->domain_, uint_c(lc->numCellsPerDim_[0]));

   WcTimer      timerImbalanced;
   WcTimer      timerLoadBalancing;
   WcTimer      timerBalanced;
   WcTimingPool tpImbalanced;
   WcTimingPool tpBalanced;

   auto    SNNBytesSent     = SNN.getBytesSent();
   auto    SNNBytesReceived = SNN.getBytesReceived();
   auto    SNNSends         = SNN.getNumberOfSends();
   auto    SNNReceives      = SNN.getNumberOfReceives();
   auto    RPBytesSent      = RP.getBytesSent();
   auto    RPBytesReceived  = RP.getBytesReceived();
   auto    RPSends          = RP.getNumberOfSends();
   auto    RPReceives       = RP.getNumberOfReceives();
   auto    SNNMinSends      = walberla::mpi::reduce(SNNSends, walberla::mpi::MIN);
   auto    SNNMaxSends      = walberla::mpi::reduce(SNNSends, walberla::mpi::MAX);
   WALBERLA_LOG_DEVEL_ON_ROOT( "SNN min/max communication partners: " << SNNMinSends << " / " << SNNMaxSends );
   walberla::mpi::reduceInplace(SNNBytesSent, walberla::mpi::SUM);
   walberla::mpi::reduceInplace(SNNBytesReceived, walberla::mpi::SUM);
   walberla::mpi::reduceInplace(SNNSends, walberla::mpi::SUM);
   walberla::mpi::reduceInplace(SNNReceives, walberla::mpi::SUM);
   walberla::mpi::reduceInplace(RPBytesSent, walberla::mpi::SUM);
   walberla::mpi::reduceInplace(RPBytesReceived, walberla::mpi::SUM);
   walberla::mpi::reduceInplace(RPSends, walberla::mpi::SUM);
   walberla::mpi::reduceInplace(RPReceives, walberla::mpi::SUM);
   WALBERLA_LOG_DEVEL_ON_ROOT( "SNN bytes communicated:   " << SNNBytesSent << " / " << SNNBytesReceived );
   WALBERLA_LOG_DEVEL_ON_ROOT( "SNN communication partners: " << SNNSends << " / " << SNNReceives );
   WALBERLA_LOG_DEVEL_ON_ROOT( "RP bytes communicated:  " << RPBytesSent << " / " << RPBytesReceived );
   WALBERLA_LOG_DEVEL_ON_ROOT( "RP communication partners: " << RPSends << " / " << RPReceives );

   int64_t imbalancedContactsChecked  = 0;
   int64_t imbalancedContactsDetected = 0;
   int64_t imbalancedContactsTreated  = 0;

   WALBERLA_MPI_BARRIER();
   WALBERLA_LOG_DEVEL_ON_ROOT("running imbalanced simulation");
   timerImbalanced.start();
   for (int64_t i=0; i < params.simulationSteps; ++i)
   {
      //WALBERLA_LOG_DEVEL_ON_ROOT("timestep: " << i << " / " << params.simulationSteps );
      //         if (i % params.visSpacing == 0)
      //         {
      //            vtkWriter->write();
      //         }

      tpImbalanced["AssocToBlock"].start();
      ps->forEachParticle(true, kernel::SelectLocal(), accessor, assoc, accessor);
      if (params.bBarrier) WALBERLA_MPI_BARRIER();
      tpImbalanced["AssocToBlock"].end();

      tpImbalanced["GenerateSparseLinkedCells"].start();
      lc->clear();
      ps->forEachParticle(true, kernel::SelectAll(), accessor, ipilc, accessor, *lc);
      if (params.bBarrier) WALBERLA_MPI_BARRIER();
      tpImbalanced["GenerateSparseLinkedCells"].end();

      tpImbalanced["DEM"].start();
      imbalancedContactsChecked  = 0;
      imbalancedContactsDetected = 0;
      imbalancedContactsTreated  = 0;
      lc->forEachParticlePairHalf(true,
                                  kernel::SelectAll(),
                                  accessor,
                                  [&](const size_t idx1, const size_t idx2, auto& ac)
      {
         ++imbalancedContactsChecked;
         if (double_cast(idx1, idx2, ac, acd, ac ))
         {
            ++imbalancedContactsDetected;
            if (contact_filter(acd.getIdx1(), acd.getIdx2(), ac, acd.getContactPoint(), *domain))
            {
               ++imbalancedContactsTreated;
               dem(acd.getIdx1(), acd.getIdx2(), ac, acd.getContactPoint(), acd.getContactNormal(), acd.getPenetrationDepth());
            }
         }
      },
      accessor );
      if (params.bBarrier) WALBERLA_MPI_BARRIER();
      tpImbalanced["DEM"].end();

      tpImbalanced["ReduceForce"].start();
      RP.operator()<ForceTorqueNotification>(*ps);
      if (params.bBarrier) WALBERLA_MPI_BARRIER();
      tpImbalanced["ReduceForce"].end();

      tpImbalanced["Euler"].start();
      //ps->forEachParticle(false, [&](const size_t idx){WALBERLA_CHECK_EQUAL(ps->getForce(idx), Vec3(0,0,0), *(*ps)[idx] << "\n" << idx);});
      ps->forEachParticle(true, kernel::SelectLocal(), accessor, explicitEuler, accessor);
      if (params.bBarrier) WALBERLA_MPI_BARRIER();
      tpImbalanced["Euler"].end();

      tpImbalanced["SNN"].start();
      SNN(*ps, forest, domain);
      if (params.bBarrier) WALBERLA_MPI_BARRIER();
      tpImbalanced["SNN"].end();
   }
   timerImbalanced.end();

   vtkDomainOutput->write( );
   vtkWriter->write();
   WALBERLA_MPI_BARRIER();
   timerLoadBalancing.start();
   if (bRebalance)
   {
      WALBERLA_LOG_DEVEL_ON_ROOT("running load balancing");
      domain::createWithNeighborhood( accessor, *forest, *ic );
      for (auto pIt = ps->begin(); pIt != ps->end(); )
      {
         using namespace walberla::mesa_pd::data::particle_flags;
         if (isSet(pIt->getFlags(), GHOST))
         {
            pIt = ps->erase(pIt);
         } else
         {
            pIt->getGhostOwnersRef().clear();
            ++pIt;
         }
      }
      forest->refresh();
      domain->refresh();
      lc = std::make_shared<data::SparseLinkedCells>(domain->getUnionOfLocalAABBs().getExtended(params.spacing), params.spacing );
      ps->forEachParticle(true, kernel::SelectLocal(), accessor, assoc, accessor);
      SNN(*ps, forest, domain);
      sortParticleStorage(*ps, params.sorting, lc->domain_, uint_c(lc->numCellsPerDim_[0]));
   }
   timerLoadBalancing.end();
   vtkDomainOutput->write( );
   vtkWriter->write();

   WALBERLA_MPI_BARRIER();
   WALBERLA_LOG_DEVEL_ON_ROOT("running balanced simulation");
   int64_t balancedContactsChecked  = 0;
   int64_t balancedContactsDetected = 0;
   int64_t balancedContactsTreated  = 0;
   timerBalanced.start();
   for (int64_t i=0; i < params.simulationSteps; ++i)
   {
      //WALBERLA_LOG_DEVEL_ON_ROOT("timestep: " << i << " / " << params.simulationSteps );
      //         if (i % params.visSpacing == 0)
      //         {
      //            vtkWriter->write();
      //         }

      tpBalanced["AssocToBlock"].start();
      ps->forEachParticle(true, kernel::SelectLocal(), accessor, assoc, accessor);
      if (params.bBarrier) WALBERLA_MPI_BARRIER();
      tpBalanced["AssocToBlock"].end();

      tpBalanced["GenerateSparseLinkedCells"].start();
      lc->clear();
      ps->forEachParticle(true, kernel::SelectAll(), accessor, ipilc, accessor, *lc);
      if (params.bBarrier) WALBERLA_MPI_BARRIER();
      tpBalanced["GenerateSparseLinkedCells"].end();

      tpBalanced["DEM"].start();
      balancedContactsChecked  = 0;
      balancedContactsDetected = 0;
      balancedContactsTreated  = 0;
      lc->forEachParticlePairHalf(true,
                                  kernel::SelectAll(),
                                  accessor,
                                  [&](const size_t idx1, const size_t idx2, auto& ac)
      {
         ++balancedContactsChecked;
         if (double_cast(idx1, idx2, ac, acd, ac ))
         {
            ++balancedContactsDetected;
            if (contact_filter(acd.getIdx1(), acd.getIdx2(), ac, acd.getContactPoint(), *domain))
            {
               ++balancedContactsTreated;
               dem(acd.getIdx1(), acd.getIdx2(), ac, acd.getContactPoint(), acd.getContactNormal(), acd.getPenetrationDepth());
            }
         }
      },
      accessor );
      if (params.bBarrier) WALBERLA_MPI_BARRIER();
      tpBalanced["DEM"].end();

      tpBalanced["ReduceForce"].start();
      RP.operator()<ForceTorqueNotification>(*ps);
      if (params.bBarrier) WALBERLA_MPI_BARRIER();
      tpBalanced["ReduceForce"].end();

      tpBalanced["Euler"].start();
      //ps->forEachParticle(false, [&](const size_t idx){WALBERLA_CHECK_EQUAL(ps->getForce(idx), Vec3(0,0,0), *(*ps)[idx] << "\n" << idx);});
      ps->forEachParticle(true, kernel::SelectLocal(), accessor, explicitEuler, accessor);
      if (params.bBarrier) WALBERLA_MPI_BARRIER();
      tpBalanced["Euler"].end();

      tpBalanced["SNN"].start();
      SNN(*ps, forest, domain);
      if (params.bBarrier) WALBERLA_MPI_BARRIER();
      tpBalanced["SNN"].end();
   }
   timerBalanced.end();

   SNNBytesSent     = SNN.getBytesSent();
   SNNBytesReceived = SNN.getBytesReceived();
   SNNSends         = SNN.getNumberOfSends();
   SNNReceives      = SNN.getNumberOfReceives();
   RPBytesSent      = RP.getBytesSent();
   RPBytesReceived  = RP.getBytesReceived();
   RPSends          = RP.getNumberOfSends();
   RPReceives       = RP.getNumberOfReceives();
   SNNMinSends      = walberla::mpi::reduce(SNNSends, walberla::mpi::MIN);
   SNNMaxSends      = walberla::mpi::reduce(SNNSends, walberla::mpi::MAX);
   WALBERLA_LOG_DEVEL_ON_ROOT( "SNN min/max communication partners: " << SNNMinSends << " / " << SNNMaxSends );
   walberla::mpi::reduceInplace(SNNBytesSent, walberla::mpi::SUM);
   walberla::mpi::reduceInplace(SNNBytesReceived, walberla::mpi::SUM);
   walberla::mpi::reduceInplace(SNNSends, walberla::mpi::SUM);
   walberla::mpi::reduceInplace(SNNReceives, walberla::mpi::SUM);
   walberla::mpi::reduceInplace(RPBytesSent, walberla::mpi::SUM);
   walberla::mpi::reduceInplace(RPBytesReceived, walberla::mpi::SUM);
   walberla::mpi::reduceInplace(RPSends, walberla::mpi::SUM);
   walberla::mpi::reduceInplace(RPReceives, walberla::mpi::SUM);
   auto cC = walberla::mpi::reduce(balancedContactsChecked, walberla::mpi::SUM);
   auto cD = walberla::mpi::reduce(balancedContactsDetected, walberla::mpi::SUM);
   auto cT = walberla::mpi::reduce(balancedContactsTreated, walberla::mpi::SUM);
   WALBERLA_LOG_DEVEL_ON_ROOT( "SNN bytes communicated:   " << SNNBytesSent << " / " << SNNBytesReceived );
   WALBERLA_LOG_DEVEL_ON_ROOT( "SNN communication partners: " << SNNSends << " / " << SNNReceives );
   WALBERLA_LOG_DEVEL_ON_ROOT( "RP bytes communicated:  " << RPBytesSent << " / " << RPBytesReceived );
   WALBERLA_LOG_DEVEL_ON_ROOT( "RP communication partners: " << RPSends << " / " << RPReceives );
   WALBERLA_LOG_DEVEL_ON_ROOT( "contacts checked/detected/treated: " << cC << " / " << cD << " / " << cT );

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

   if (params.checkSimulation)
   {
      check(*ps, *forest, params.spacing, params.shift);
   }

   WALBERLA_LOG_INFO_ON_ROOT("*** SQL OUTPUT - START ***");
   numParticles = 0;
   int64_t numGhostParticles = 0;
   ps->forEachParticle(false,
                       kernel::SelectAll(),
                       accessor,
                       [&numParticles, &numGhostParticles](const size_t idx, auto& ac)
   {
      if (data::particle_flags::isSet( ac.getFlagsRef(idx), data::particle_flags::GHOST))
      {
         ++numGhostParticles;
      } else
      {
         ++numParticles;
      }
   },
   accessor);
   auto minParticles = walberla::mpi::reduce(numParticles, walberla::mpi::MIN);
   auto maxParticles = walberla::mpi::reduce(numParticles, walberla::mpi::MAX);
   auto minGhostParticles = walberla::mpi::reduce(numGhostParticles, walberla::mpi::MIN);
   auto maxGhostParticles = walberla::mpi::reduce(numGhostParticles, walberla::mpi::MAX);
   auto minTotalParticles = walberla::mpi::reduce(numParticles + numGhostParticles, walberla::mpi::MIN);
   auto maxTotalParticles = walberla::mpi::reduce(numParticles + numGhostParticles, walberla::mpi::MAX);
   WALBERLA_LOG_DEVEL_ON_ROOT("particle ratio: " << minParticles << " / " << maxParticles);
   walberla::mpi::reduceInplace(numParticles, walberla::mpi::SUM);
   walberla::mpi::reduceInplace(numGhostParticles, walberla::mpi::SUM);
   walberla::mpi::reduceInplace(imbalancedContactsChecked, walberla::mpi::SUM);
   walberla::mpi::reduceInplace(imbalancedContactsDetected, walberla::mpi::SUM);
   walberla::mpi::reduceInplace(imbalancedContactsTreated, walberla::mpi::SUM);
   walberla::mpi::reduceInplace(balancedContactsChecked, walberla::mpi::SUM);
   walberla::mpi::reduceInplace(balancedContactsDetected, walberla::mpi::SUM);
   walberla::mpi::reduceInplace(balancedContactsTreated, walberla::mpi::SUM);
   double linkedCellsVolume = lc->domain_.volume();
   walberla::mpi::reduceInplace(linkedCellsVolume, walberla::mpi::SUM);
   size_t numLinkedCells = lc->cells_.size();
   walberla::mpi::reduceInplace(numLinkedCells, walberla::mpi::SUM);
   size_t local_aabbs         = domain->getNumLocalAABBs();
   size_t neighbor_subdomains = domain->getNumNeighborSubdomains();
   size_t neighbor_processes  = domain->getNumNeighborProcesses();
   walberla::mpi::reduceInplace(local_aabbs, walberla::mpi::SUM);
   walberla::mpi::reduceInplace(neighbor_subdomains, walberla::mpi::SUM);
   walberla::mpi::reduceInplace(neighbor_processes, walberla::mpi::SUM);

   uint_t runId = uint_c(-1);
   WALBERLA_ROOT_SECTION()
   {
      stringProperties["walberla_git"]                  = WALBERLA_GIT_SHA1;
      stringProperties["tag"]                           = "mesa_pd";
      integerProperties["mpi_num_processes"]            = mpiManager->numProcesses();
      integerProperties["omp_max_threads"]              = omp_get_max_threads();
      realProperties["imbalanced_PUpS"]                 = double_c(PUpSImbalanced);
      realProperties["imbalanced_timer_min"]            = timerImbalancedReduced->min();
      realProperties["imbalanced_timer_max"]            = timerImbalancedReduced->max();
      realProperties["imbalanced_timer_average"]        = timerImbalancedReduced->average();
      realProperties["imbalanced_timer_total"]          = timerImbalancedReduced->total();
      realProperties["loadbalancing_timer_min"]         = timerLoadBalancingReduced->min();
      realProperties["loadbalancing_timer_max"]         = timerLoadBalancingReduced->max();
      realProperties["loadbalancing_timer_average"]     = timerLoadBalancingReduced->average();
      realProperties["loadbalancing_timer_total"]       = timerLoadBalancingReduced->total();
      realProperties["balanced_PUpS"]                   = double_c(PUpSBalanced);
      realProperties["balanced_timer_min"]              = timerBalancedReduced->min();
      realProperties["balanced_timer_max"]              = timerBalancedReduced->max();
      realProperties["balanced_timer_average"]          = timerBalancedReduced->average();
      realProperties["balanced_timer_total"]            = timerBalancedReduced->total();
      integerProperties["num_particles"]                = numParticles;
      integerProperties["num_ghost_particles"]          = numGhostParticles;
      integerProperties["minParticles"]                 = minParticles;
      integerProperties["maxParticles"]                 = maxParticles;
      integerProperties["minGhostParticles"]            = minGhostParticles;
      integerProperties["maxGhostParticles"]            = maxGhostParticles;
      integerProperties["minTotalParticles"]            = minTotalParticles;
      integerProperties["maxTotalParticles"]            = maxTotalParticles;
      integerProperties["imbalancedContactsChecked"]    = imbalancedContactsChecked;
      integerProperties["imbalancedContactsDetected"]   = imbalancedContactsDetected;
      integerProperties["imbalancedContactsTreated"]    = imbalancedContactsTreated;
      integerProperties["balancedContactsChecked"]      = balancedContactsChecked;
      integerProperties["balancedContactsDetected"]     = balancedContactsDetected;
      integerProperties["balancedContactsTreated"]      = balancedContactsTreated;
      integerProperties["local_aabbs"]                  = int64_c(local_aabbs);
      integerProperties["neighbor_subdomains"]          = int64_c(neighbor_subdomains);
      integerProperties["neighbor_processes"]           = int64_c(neighbor_processes);
      integerProperties["SNNBytesSent"]                 = SNNBytesSent;
      integerProperties["SNNBytesReceived"]             = SNNBytesReceived;
      integerProperties["SNNSends"]                     = SNNSends;
      integerProperties["SNNReceives"]                  = SNNReceives;
      integerProperties["SNNMinSends"]                  = SNNMinSends;
      integerProperties["SNNMaxSends"]                  = SNNMaxSends;
      integerProperties["RPBytesSent"]                  = RPBytesSent;
      integerProperties["RPBytesReceived"]              = RPBytesReceived;
      integerProperties["RPSends"]                      = RPSends;
      integerProperties["RPReceives"]                   = RPReceives;
      realProperties["linkedCellsVolume"]               = linkedCellsVolume;
      integerProperties["numLinkedCells"]               = int64_c(numLinkedCells);

      addBuildInfoToSQL( integerProperties, realProperties, stringProperties );
      saveToSQL(params, integerProperties, realProperties, stringProperties );
      addDomainPropertiesToSQL(*forest, integerProperties, realProperties, stringProperties);
      addSlurmPropertiesToSQL(integerProperties, realProperties, stringProperties);

      runId = sqlite::storeRunInSqliteDB( params.sqlFile, integerProperties, stringProperties, realProperties );
      sqlite::storeTimingPoolInSqliteDB( params.sqlFile, runId, *tpImbalancedReduced, "imbalanced" );
      sqlite::storeTimingPoolInSqliteDB( params.sqlFile, runId, *tpBalancedReduced, "balanced" );
   }

   if (params.storeNodeTimings)
   {
      storeNodeTimings(runId, params.sqlFile, "NodeTimingImbalanced", tpImbalanced);
      storeNodeTimings(runId, params.sqlFile, "NodeTimingBalanced", tpBalanced);
   }
   WALBERLA_LOG_INFO_ON_ROOT("*** SQL OUTPUT - END ***");

   return EXIT_SUCCESS;
}

} // namespace mesa_pd
} // namespace walberla

int main( int argc, char* argv[] )
{
   return walberla::mesa_pd::main( argc, argv );
}
