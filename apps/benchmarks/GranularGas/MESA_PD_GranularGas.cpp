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
//! \file   MESA_PD_GranularGas.cpp
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
#include <mesa_pd/data/LinkedCells.h>
#include <mesa_pd/data/ParticleAccessor.h>
#include <mesa_pd/data/ParticleStorage.h>
#include <mesa_pd/data/ShapeStorage.h>
#include <mesa_pd/domain/BlockForestDomain.h>
#include <mesa_pd/kernel/AssocToBlock.h>
#include <mesa_pd/kernel/DoubleCast.h>
#include <mesa_pd/kernel/ExplicitEulerWithShape.h>
#include <mesa_pd/kernel/InsertParticleIntoLinkedCells.h>
#include <mesa_pd/kernel/ParticleSelector.h>
#include <mesa_pd/kernel/SpringDashpot.h>
#include <mesa_pd/mpi/ContactFilter.h>
#include <mesa_pd/mpi/ReduceProperty.h>
#include <mesa_pd/mpi/SyncNextNeighbors.h>
#include <mesa_pd/mpi/SyncNextNeighborsBlockForest.h>

#include <mesa_pd/mpi/notifications/ForceTorqueNotification.h>

#include <blockforest/BlockForest.h>
#include <blockforest/Initialization.h>
#include <core/Abort.h>
#include <core/Environment.h>
#include <core/math/Random.h>
#include <core/mpi/Reduce.h>
#include <core/grid_generator/SCIterator.h>
#include <core/logging/Logging.h>
#include <core/OpenMP.h>
#include <core/timing/Timer.h>
#include <core/waLBerlaBuildInfo.h>
#include <sqlite/SQLite.h>
#include <vtk/VTKOutput.h>

#include <functional>
#include <memory>
#include <string>
#include <type_traits>

namespace walberla {
namespace mesa_pd {

int main( int argc, char ** argv )
{
   using namespace walberla::timing;

   Environment env(argc, argv);
   auto mpiManager = walberla::mpi::MPIManager::instance();
   mpiManager->useWorldComm();

   logging::Logging::instance()->setStreamLogLevel(logging::Logging::INFO);
   logging::Logging::instance()->setFileLogLevel(logging::Logging::INFO);

   WALBERLA_LOG_INFO_ON_ROOT( "config file: " << argv[1] );
   WALBERLA_LOG_INFO_ON_ROOT( "waLBerla Revision: " << WALBERLA_GIT_SHA1 );

   math::seedRandomGenerator( static_cast<unsigned int>(1337 * mpiManager->worldRank()) );

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
   auto domain = std::make_shared<domain::BlockForestDomain> (forest);

   auto simulationDomain = forest->getDomain();
   auto localDomain = forest->begin()->getAABB();
   for (auto& blk : *forest)
   {
      localDomain.merge(blk.getAABB());
   }

   WALBERLA_LOG_INFO_ON_ROOT("*** SETUP - START ***");

   //init data structures
   auto ps = std::make_shared<data::ParticleStorage>(100);
   auto ss = std::make_shared<data::ShapeStorage>();
   ParticleAccessorWithShape accessor(ps, ss);
   data::LinkedCells     lc(localDomain.getExtended(params.spacing), params.spacing );

   auto  smallSphere = ss->create<data::Sphere>( params.radius );
   ss->shapes[smallSphere]->updateMassAndInertia(real_t(2707));
   for (auto& iBlk : *forest)
   {
      for (auto pt : grid_generator::SCGrid(iBlk.getAABB(),
                                            Vector3<real_t>(params.spacing) * real_c(0.5) + params.shift,
                                            params.spacing))
      {
         WALBERLA_CHECK(iBlk.getAABB().contains(pt));
         createSphere(*ps, pt, params.radius, smallSphere);
      }
   }
   int64_t numParticles = int64_c(ps->size());
   walberla::mpi::reduceInplace(numParticles, walberla::mpi::SUM);
   WALBERLA_LOG_INFO_ON_ROOT("#particles created: " << numParticles);

   auto shift = (params.spacing - params.radius - params.radius) * real_t(0.5);
   auto confiningDomain = simulationDomain.getExtended(shift);

   if (!forest->isPeriodic(0))
   {
      createPlane(*ps, *ss, confiningDomain.minCorner(), Vec3(+1,0,0));
      createPlane(*ps, *ss, confiningDomain.maxCorner(), Vec3(-1,0,0));
   }

   if (!forest->isPeriodic(1))
   {
      createPlane(*ps, *ss, confiningDomain.minCorner(), Vec3(0,+1,0));
      createPlane(*ps, *ss, confiningDomain.maxCorner(), Vec3(0,-1,0));
   }

   if (!forest->isPeriodic(2))
   {
      createPlane(*ps, *ss, confiningDomain.minCorner(), Vec3(0,0,+1));
      createPlane(*ps, *ss, confiningDomain.maxCorner(), Vec3(0,0,-1));
   }

   WALBERLA_LOG_INFO_ON_ROOT("*** SETUP - END ***");

   WALBERLA_LOG_INFO_ON_ROOT("*** VTK ***");
   auto vtkDomainOutput = walberla::vtk::createVTKOutput_DomainDecomposition( forest, "domain_decomposition", 1, params.vtk_out, "simulation_step" );
   auto vtkOutput       = make_shared<mesa_pd::vtk::ParticleVtkOutput>(ps) ;
   auto vtkWriter       = walberla::vtk::createVTKOutput_PointData(vtkOutput, "Bodies", 1, params.vtk_out, "simulation_step", false, false);
   vtkOutput->addOutput<SelectRank>("rank");
   vtkOutput->addOutput<data::SelectParticleOwner>("owner");
   //   vtkDomainOutput->write();

   WALBERLA_LOG_INFO_ON_ROOT("*** SIMULATION - START ***");
   // Init kernels
   kernel::ExplicitEulerWithShape        explicitEulerWithShape( params.dt );
   kernel::InsertParticleIntoLinkedCells ipilc;
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
   ps->forEachParticle(false, kernel::SelectLocal(), accessor, assoc, accessor);
   SNN(*ps, forest, domain);
   sortParticleStorage(*ps, params.sorting, lc.domain_, uint_c(lc.numCellsPerDim_[0]));
//   vtkWriter->write();

   for (int64_t outerIteration = 0; outerIteration < params.numOuterIterations; ++outerIteration)
   {
      WALBERLA_LOG_INFO_ON_ROOT("*** RUNNING OUTER ITERATION " << outerIteration << " ***");

      WcTimer      timer;
      WcTimingPool tp;
      auto    SNNBytesSent     = SNN.getBytesSent();
      auto    SNNBytesReceived = SNN.getBytesReceived();
      auto    SNNSends         = SNN.getNumberOfSends();
      auto    SNNReceives      = SNN.getNumberOfReceives();
      auto    RPBytesSent      = RP.getBytesSent();
      auto    RPBytesReceived  = RP.getBytesReceived();
      auto    RPSends          = RP.getNumberOfSends();
      auto    RPReceives       = RP.getNumberOfReceives();
      int64_t contactsChecked  = 0;
      int64_t contactsDetected = 0;
      int64_t contactsTreated  = 0;
      if (params.bBarrier) WALBERLA_MPI_BARRIER();
      timer.start();
      for (int64_t i=0; i < params.simulationSteps; ++i)
      {
         //      if (i % visSpacing == 0)
         //      {
         //         vtkWriter->write();
         //      }

         tp["AssocToBlock"].start();
         ps->forEachParticle(false, kernel::SelectLocal(), accessor, assoc, accessor);
         if (params.bBarrier) WALBERLA_MPI_BARRIER();
         tp["AssocToBlock"].end();

         tp["GenerateLinkedCells"].start();
         lc.clear();
         ps->forEachParticle(true, kernel::SelectAll(), accessor, ipilc, accessor, lc);
         if (params.bBarrier) WALBERLA_MPI_BARRIER();
         tp["GenerateLinkedCells"].end();

         tp["DEM"].start();
         contactsChecked  = 0;
         contactsDetected = 0;
         contactsTreated  = 0;
         lc.forEachParticlePairHalf(true,
                                    kernel::SelectAll(),
                                    accessor,
                                    [&](const size_t idx1, const size_t idx2, auto& ac)
         {
            ++contactsChecked;
            if (double_cast(idx1, idx2, ac, acd, ac ))
            {
               ++contactsDetected;
               if (contact_filter(acd.getIdx1(), acd.getIdx2(), ac, acd.getContactPoint(), *domain))
               {
                  ++contactsTreated;
                  dem(acd.getIdx1(), acd.getIdx2(), ac, acd.getContactPoint(), acd.getContactNormal(), acd.getPenetrationDepth());
               }
            }
         },
         accessor );
         if (params.bBarrier) WALBERLA_MPI_BARRIER();
         tp["DEM"].end();

         tp["ReduceForce"].start();
         RP.operator()<ForceTorqueNotification>(*ps);
         if (params.bBarrier) WALBERLA_MPI_BARRIER();
         tp["ReduceForce"].end();

         tp["Euler"].start();
         //ps->forEachParticle(false, [&](const size_t idx){WALBERLA_CHECK_EQUAL(ps->getForce(idx), Vec3(0,0,0), *(*ps)[idx] << "\n" << idx);});
         ps->forEachParticle(true, kernel::SelectLocal(), accessor, explicitEulerWithShape, accessor);
         if (params.bBarrier) WALBERLA_MPI_BARRIER();
         tp["Euler"].end();

         tp["SNN"].start();
         SNN(*ps, forest, domain);
         if (params.bBarrier) WALBERLA_MPI_BARRIER();
         tp["SNN"].end();
      }
      timer.end();

      SNNBytesSent     = SNN.getBytesSent();
      SNNBytesReceived = SNN.getBytesReceived();
      SNNSends         = SNN.getNumberOfSends();
      SNNReceives      = SNN.getNumberOfReceives();
      RPBytesSent      = RP.getBytesSent();
      RPBytesReceived  = RP.getBytesReceived();
      RPSends          = RP.getNumberOfSends();
      RPReceives       = RP.getNumberOfReceives();
      walberla::mpi::reduceInplace(SNNBytesSent, walberla::mpi::SUM);
      walberla::mpi::reduceInplace(SNNBytesReceived, walberla::mpi::SUM);
      walberla::mpi::reduceInplace(SNNSends, walberla::mpi::SUM);
      walberla::mpi::reduceInplace(SNNReceives, walberla::mpi::SUM);
      walberla::mpi::reduceInplace(RPBytesSent, walberla::mpi::SUM);
      walberla::mpi::reduceInplace(RPBytesReceived, walberla::mpi::SUM);
      walberla::mpi::reduceInplace(RPSends, walberla::mpi::SUM);
      walberla::mpi::reduceInplace(RPReceives, walberla::mpi::SUM);
      auto cC = walberla::mpi::reduce(contactsChecked, walberla::mpi::SUM);
      auto cD = walberla::mpi::reduce(contactsDetected, walberla::mpi::SUM);
      auto cT = walberla::mpi::reduce(contactsTreated, walberla::mpi::SUM);
      WALBERLA_LOG_DEVEL_ON_ROOT( "SNN bytes communicated:   " << SNNBytesSent << " / " << SNNBytesReceived );
      WALBERLA_LOG_DEVEL_ON_ROOT( "SNN communication partners: " << SNNSends << " / " << SNNReceives );
      WALBERLA_LOG_DEVEL_ON_ROOT( "RP bytes communicated:  " << RPBytesSent << " / " << RPBytesReceived );
      WALBERLA_LOG_DEVEL_ON_ROOT( "RP communication partners: " << RPSends << " / " << RPReceives );
      WALBERLA_LOG_DEVEL_ON_ROOT( "contacts checked/detected/treated: " << cC << " / " << cD << " / " << cT );

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

      if (params.checkSimulation)
      {
         check(*ps, *forest, params.spacing);
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
      walberla::mpi::reduceInplace(numParticles, walberla::mpi::SUM);
      walberla::mpi::reduceInplace(numGhostParticles, walberla::mpi::SUM);
      walberla::mpi::reduceInplace(contactsChecked, walberla::mpi::SUM);
      walberla::mpi::reduceInplace(contactsDetected, walberla::mpi::SUM);
      walberla::mpi::reduceInplace(contactsTreated, walberla::mpi::SUM);
      double linkedCellsVolume = lc.domain_.volume();
      walberla::mpi::reduceInplace(linkedCellsVolume, walberla::mpi::SUM);
      size_t numLinkedCells = lc.cells_.size();
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
         std::map< std::string, walberla::int64_t > integerProperties;
         std::map< std::string, double >            realProperties;
         std::map< std::string, std::string >       stringProperties;

         stringProperties["walberla_git"]         = WALBERLA_GIT_SHA1;
         stringProperties["tag"]                  = "mesa_pd";
         integerProperties["mpi_num_processes"]   = mpiManager->numProcesses();
         integerProperties["omp_max_threads"]     = omp_get_max_threads();
         realProperties["PUpS"]                   = double_c(PUpS);
         realProperties["timer_min"]              = timer_reduced->min();
         realProperties["timer_max"]              = timer_reduced->max();
         realProperties["timer_average"]          = timer_reduced->average();
         realProperties["timer_total"]            = timer_reduced->total();
         integerProperties["outerIteration"]      = int64_c(outerIteration);
         integerProperties["num_particles"]       = numParticles;
         integerProperties["num_ghost_particles"] = numGhostParticles;
         integerProperties["contacts_checked"]    = contactsChecked;
         integerProperties["contacts_detected"]   = contactsDetected;
         integerProperties["contacts_treated"]    = contactsTreated;
         integerProperties["local_aabbs"]         = int64_c(local_aabbs);
         integerProperties["neighbor_subdomains"] = int64_c(neighbor_subdomains);
         integerProperties["neighbor_processes"]  = int64_c(neighbor_processes);
         integerProperties["SNNBytesSent"]        = SNNBytesSent;
         integerProperties["SNNBytesReceived"]    = SNNBytesReceived;
         integerProperties["SNNSends"]            = SNNSends;
         integerProperties["SNNReceives"]         = SNNReceives;
         integerProperties["RPBytesSent"]         = RPBytesSent;
         integerProperties["RPBytesReceived"]     = RPBytesReceived;
         integerProperties["RPSends"]             = RPSends;
         integerProperties["RPReceives"]          = RPReceives;
         realProperties["linkedCellsVolume"]      = linkedCellsVolume;
         integerProperties["numLinkedCells"]      = int64_c(numLinkedCells);
         realProperties["PUpS"]                   = double_c(PUpS);
         realProperties["timer_min"]              = timer_reduced->min();
         realProperties["timer_max"]              = timer_reduced->max();
         realProperties["timer_average"]          = timer_reduced->average();
         realProperties["timer_total"]            = timer_reduced->total();

         addBuildInfoToSQL( integerProperties, realProperties, stringProperties );
         saveToSQL(params, integerProperties, realProperties, stringProperties );
         addDomainPropertiesToSQL(*forest, integerProperties, realProperties, stringProperties);
         addSlurmPropertiesToSQL(integerProperties, realProperties, stringProperties);

         runId = sqlite::storeRunInSqliteDB( params.sqlFile, integerProperties, stringProperties, realProperties );
         sqlite::storeTimingPoolInSqliteDB( params.sqlFile, runId, *tp_reduced, "Timeloop" );
      }

      if (params.storeNodeTimings)
      {
         storeNodeTimings(runId, params.sqlFile, "NodeTiming", tp);
      }
      WALBERLA_LOG_INFO_ON_ROOT("*** SQL OUTPUT - END ***");
   }

   return EXIT_SUCCESS;
}

} // namespace mesa_pd
} // namespace walberla

int main( int argc, char* argv[] )
{
   return walberla::mesa_pd::main( argc, argv );
}
