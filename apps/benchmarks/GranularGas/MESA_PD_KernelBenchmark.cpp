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
//! \file   MESA_PD_KernelBenchmark.cpp
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

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
#include <mesa_pd/data/ParticleAccessorWithShape.h>
#include <mesa_pd/data/ParticleStorage.h>
#include <mesa_pd/data/ShapeStorage.h>
#include <mesa_pd/domain/BlockForestDomain.h>
#include <mesa_pd/kernel/AssocToBlock.h>
#include <mesa_pd/kernel/DoubleCast.h>
#include <mesa_pd/kernel/ExplicitEuler.h>
#include <mesa_pd/kernel/InsertParticleIntoLinkedCells.h>
#include <mesa_pd/kernel/ParticleSelector.h>
#include <mesa_pd/kernel/SpringDashpot.h>
#include <mesa_pd/kernel/SpringDashpotSpring.h>
#include <mesa_pd/mpi/ContactFilter.h>
#include <mesa_pd/mpi/ReduceContactHistory.h>
#include <mesa_pd/mpi/ReduceProperty.h>
#include <mesa_pd/mpi/SyncNextNeighbors.h>
#include <mesa_pd/mpi/SyncNextNeighborsBlockForest.h>
#include <mesa_pd/mpi/notifications/ForceTorqueNotification.h>
#include <mesa_pd/sorting/HilbertCompareFunctor.h>
#include <mesa_pd/sorting/LinearizedCompareFunctor.h>

#include <blockforest/BlockForest.h>
#include <blockforest/Initialization.h>
#include <core/Abort.h>
#include <core/Environment.h>
#include <core/Hostname.h>
#include <core/math/Random.h>
#include <core/math/Sample.h>
#include <core/mpi/Gatherv.h>
#include <core/mpi/RecvBuffer.h>
#include <core/mpi/Reduce.h>
#include <core/mpi/SendBuffer.h>
#include <core/grid_generator/SCIterator.h>
#include <core/logging/Logging.h>
#include <core/OpenMP.h>
#include <core/perf_analysis/extern/likwid.h>
#include <core/timing/Timer.h>
#include <core/timing/TimingPool.h>
#include <core/waLBerlaBuildInfo.h>
#include <sqlite/SQLite.h>
#include <vtk/VTKOutput.h>

#include <functional>
#include <memory>
#include <string>
#include <type_traits>

namespace walberla {
namespace mesa_pd {

class ContactDetection
{
public:
   ContactDetection(const std::shared_ptr<domain::BlockForestDomain>& domain)
   : domain_(domain)
   {
      contacts_.reserve(4000000);
   }

   template <typename T>
   inline
   void operator()(const size_t idx1, const size_t idx2, T& ac)
   {
      ++contactsChecked_;
      if (double_cast_(idx1, idx2, ac, acd_, ac))
      {
         ++contactsDetected_;
         if (contact_filter_(acd_.getIdx1(), acd_.getIdx2(), ac, acd_.getContactPoint(),
                             *domain_))
         {
            ++contactsTreated_;
            contacts_.emplace_back(acd_.getIdx1(), acd_.getIdx2(), acd_.getContactPoint(),
                 acd_.getContactNormal(), acd_.getPenetrationDepth());
         }
      }
   }

   inline const auto& getContacts() const {return contacts_;}

   inline
   void resetCounters()
   {
      contactsChecked_ = 0;
      contactsDetected_ = 0;
      contactsTreated_ = 0;
      contacts_.clear();
   }

   inline
   int64_t getContactsChecked() const
   {
      return contactsChecked_;
   }

   inline
   int64_t getContactsDetected() const
   {
      return contactsDetected_;
   }

   inline
   int64_t getContactsTreated() const
   {
      return contactsTreated_;
   }

private:
   kernel::DoubleCast double_cast_;
   mpi::ContactFilter contact_filter_;
   std::shared_ptr<domain::BlockForestDomain> domain_;
   collision_detection::AnalyticContactDetection acd_;
   std::vector<Contact> contacts_;
   int64_t contactsChecked_ = 0;
   int64_t contactsDetected_ = 0;
   int64_t contactsTreated_ = 0;
};

template <typename T>
void reportOverRanks(const std::string& info, const T& value)
{
   math::Sample sample;
   sample.insert(value);
   sample.mpiGatherRoot();
   WALBERLA_LOG_INFO_ON_ROOT(info);
   WALBERLA_LOG_INFO_ON_ROOT(sample.format());
}

int main( int argc, char ** argv )
{
   LIKWID_MARKER_INIT;

   using namespace walberla::timing;

   Environment env(argc, argv);
   auto mpiManager = walberla::mpi::MPIManager::instance();
   mpiManager->useWorldComm();

   logging::Logging::instance()->setStreamLogLevel(logging::Logging::INFO);
   logging::Logging::instance()->setFileLogLevel(logging::Logging::INFO);
//   logging::Logging::instance()->includeLoggingToFile("KernelBenchmark");

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
   auto domain = std::make_shared<domain::BlockForestDomain>(forest);

   LIKWID_MARKER_THREADINIT;

   auto simulationDomain = forest->getDomain();
   auto localDomain = forest->begin()->getAABB();
   for (auto& blk : *forest)
   {
      localDomain.merge(blk.getAABB());
   }
//   WALBERLA_LOG_DEVEL_VAR(localDomain);

   WALBERLA_LOG_INFO_ON_ROOT("*** SETUP - START ***");

   //init data structures
   auto ss = std::make_shared<data::ShapeStorage>();
   auto ps = std::make_shared<data::ParticleStorage>(100);
   data::ParticleAccessorWithShape accessor(ps, ss);
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
//         WALBERLA_LOG_DEVEL_VAR_ON_ROOT(pt);
      }
   }
   int64_t numParticles = int64_c(ps->size());
   walberla::mpi::reduceInplace(numParticles, walberla::mpi::SUM);
   WALBERLA_LOG_INFO_ON_ROOT("#particles created: " << numParticles);

   auto shift = (params.spacing - params.radius - params.radius) * real_t(0.5);
   auto confiningDomain = simulationDomain.getExtended(shift);

   if (!forest->isPeriodic(0))
   {
      createPlane(*ps, *ss, confiningDomain.minCorner() + params.shift, Vec3(+1,0,0));
      createPlane(*ps, *ss, confiningDomain.maxCorner() + params.shift, Vec3(-1,0,0));
   }

   if (!forest->isPeriodic(1))
   {
      createPlane(*ps, *ss, confiningDomain.minCorner() + params.shift, Vec3(0,+1,0));
      createPlane(*ps, *ss, confiningDomain.maxCorner() + params.shift, Vec3(0,-1,0));
   }

   if (!forest->isPeriodic(2))
   {
      createPlane(*ps, *ss, confiningDomain.minCorner() + params.shift, Vec3(0,0,+1));
      createPlane(*ps, *ss, confiningDomain.maxCorner() + params.shift, Vec3(0,0,-1));
   }

   WALBERLA_LOG_INFO_ON_ROOT("*** SETUP - END ***");

   WALBERLA_LOG_INFO_ON_ROOT("*** VTK ***");
   auto vtkDomainOutput = walberla::vtk::createVTKOutput_DomainDecomposition( forest, "domain_decomposition", 1, params.vtk_out, "simulation_step" );
   auto vtkOutput       = make_shared<mesa_pd::vtk::ParticleVtkOutput>(ps) ;
   auto vtkWriter       = walberla::vtk::createVTKOutput_PointData(vtkOutput, "Bodies", 1, params.vtk_out, "simulation_step", false, false);
   vtkOutput->addOutput<SelectRank>("rank");
   vtkOutput->addOutput<data::SelectParticleOwner>("owner");
   vtkOutput->addOutput<SelectIdx>("idx");
//   vtkDomainOutput->write();

   WALBERLA_LOG_INFO_ON_ROOT("*** SIMULATION - START ***");
   // Init kernels
   kernel::ExplicitEuler                 explicitEuler( params.dt );
   kernel::AssocToBlock                  assoc(forest);
   kernel::InsertParticleIntoLinkedCells ipilc;
   kernel::SpringDashpot                 sd(1);
   sd.setDampingT (0, 0, real_t(0));
   sd.setFriction (0, 0, real_t(0));
   sd.setParametersFromCOR(0, 0, real_t(0.9), params.dt*real_t(20), ss->shapes[smallSphere]->getMass() * real_t(0.5));
   kernel::SpringDashpotSpring           sds(1);
   sds.setParametersFromCOR(0, 0, real_t(0.9), params.dt*real_t(20), ss->shapes[smallSphere]->getMass() * real_t(0.5));
   sds.setCoefficientOfFriction(0,0,real_t(0.4));
   sds.setStiffnessT(0,0,real_t(0.9) * sds.getStiffnessN(0,0));

   mpi::ReduceProperty                   RP;
   mpi::SyncNextNeighbors                SNN;
   mpi::ReduceContactHistory             RCH;
   ContactDetection                      CD(domain);

   // initial sync
   ps->forEachParticle(false, kernel::SelectLocal(), accessor, assoc, accessor);
   SNN(*ps, *domain);
   sortParticleStorage(*ps, params.sorting, lc.domain_, uint_c(lc.numCellsPerDim_[0]));
//   vtkWriter->write();

   for (int64_t outerIteration = 0; outerIteration < params.numOuterIterations; ++outerIteration)
   {
      WALBERLA_LOG_INFO_ON_ROOT("*** RUNNING OUTER ITERATION " << outerIteration << " ***");

      WcTimingPool tp;

      LIKWID_MARKER_REGISTER("SNN");
      WALBERLA_MPI_BARRIER();
      LIKWID_MARKER_START("SNN");
      tp["SNN"].start();
      for (int64_t i=0; i < params.simulationSteps; ++i)
      {
         SNN(*ps, *domain);
      }
      tp["SNN"].end();
      LIKWID_MARKER_STOP("SNN");

      LIKWID_MARKER_REGISTER("AssocToBlock");
      WALBERLA_MPI_BARRIER();
      LIKWID_MARKER_START("AssocToBlock");
      tp["AssocToBlock"].start();
      for (int64_t i=0; i < params.simulationSteps; ++i)
      {
         ps->forEachParticle(false, kernel::SelectLocal(), accessor, assoc, accessor);
      }
      tp["AssocToBlock"].end();
      LIKWID_MARKER_STOP("AssocToBlock");

      LIKWID_MARKER_REGISTER("GenerateLinkedCells");
      WALBERLA_MPI_BARRIER();
      LIKWID_MARKER_START("GenerateLinkedCells");
      tp["GenerateLinkedCells"].start();
      for (int64_t i=0; i < params.simulationSteps; ++i)
      {
         lc.clear();
         ps->forEachParticle(true, kernel::SelectAll(), accessor, ipilc, accessor, lc);
      }
      tp["GenerateLinkedCells"].end();
      LIKWID_MARKER_STOP("GenerateLinkedCells");

      LIKWID_MARKER_REGISTER("ContactDetection");
      WALBERLA_MPI_BARRIER();
      LIKWID_MARKER_START("ContactDetection");
      tp["ContactDetection"].start();
      for (int64_t i=0; i < params.simulationSteps; ++i)
      {
         CD.resetCounters();
         lc.forEachParticlePairHalf(true,
                                    kernel::SelectAll(),
                                    accessor,
                                    CD,
                                    accessor);
      }
      tp["ContactDetection"].end();
      LIKWID_MARKER_STOP("ContactDetection");

      LIKWID_MARKER_REGISTER("SpringDashpot");
      WALBERLA_MPI_BARRIER();
      LIKWID_MARKER_START("SpringDashpot");
      tp["SpringDashpot"].start();
      for (int64_t i=0; i < params.simulationSteps; ++i)
      {
         for (auto& c : CD.getContacts())
         {
            sd(c.idx1_, c.idx2_, accessor, c.contactPoint_, c.contactNormal_, c.penetrationDepth_);
         }
      }
      tp["SpringDashpot"].end();
      LIKWID_MARKER_STOP("SpringDashpot");

      LIKWID_MARKER_REGISTER("SpringDashpotSpring");
      WALBERLA_MPI_BARRIER();
      LIKWID_MARKER_START("SpringDashpotSpring");
      tp["SpringDashpotSpring"].start();
      for (int64_t i=0; i < params.simulationSteps; ++i)
      {
         for (auto& c : CD.getContacts())
         {
            sds(c.idx1_, c.idx2_, accessor, c.contactPoint_, c.contactNormal_, c.penetrationDepth_, params.dt);
         }
      }
      tp["SpringDashpotSpring"].end();
      LIKWID_MARKER_STOP("SpringDashpotSpring");

      LIKWID_MARKER_REGISTER("ReduceForce");
      WALBERLA_MPI_BARRIER();
      LIKWID_MARKER_START("ReduceForce");
      tp["ReduceForce"].start();
      for (int64_t i=0; i < params.simulationSteps; ++i)
      {
         RP.operator()<ForceTorqueNotification>(*ps);
      }
      tp["ReduceForce"].end();
      LIKWID_MARKER_STOP("ReduceForce");

      LIKWID_MARKER_REGISTER("ReduceContactHistory");
      WALBERLA_MPI_BARRIER();
      LIKWID_MARKER_START("ReduceContactHistory");
      tp["ReduceContactHistory"].start();
      for (int64_t i=0; i < params.simulationSteps; ++i)
      {
         RCH(*ps);
      }
      tp["ReduceContactHistory"].end();
      LIKWID_MARKER_STOP("ReduceContactHistory");

      LIKWID_MARKER_REGISTER("Euler");
      WALBERLA_MPI_BARRIER();
      LIKWID_MARKER_START("Euler");
      tp["Euler"].start();
      for (int64_t i=0; i < params.simulationSteps; ++i)
      {
         ps->forEachParticle(true, kernel::SelectLocal(), accessor, explicitEuler, accessor);
      }
      tp["Euler"].end();
      LIKWID_MARKER_STOP("Euler");

      WALBERLA_LOG_INFO_ON_ROOT("*** SIMULATION - END ***");

      if (params.checkSimulation)
      {
         //if you want to activate checking you have to deactivate sorting
         check(*ps, *forest, params.spacing, params.shift);
      }

      WALBERLA_LOG_INFO_ON_ROOT("*** SQL OUTPUT - START ***");
      int64_t contactsChecked = CD.getContactsChecked();
      int64_t contactsDetected = CD.getContactsDetected();
      int64_t contactsTreated = CD.getContactsTreated();

      reportOverRanks("Contacts Checked:", real_c(contactsChecked));
      reportOverRanks("Contacts Detected:", real_c(contactsDetected));
      reportOverRanks("Contacts Treated:", real_c(contactsTreated));

      auto SNNBytesSent     = SNN.getBytesSent();
      auto SNNBytesReceived = SNN.getBytesReceived();
      auto SNNSends         = SNN.getNumberOfSends();
      auto SNNReceives      = SNN.getNumberOfReceives();
      auto RPBytesSent      = RP.getBytesSent();
      auto RPBytesReceived  = RP.getBytesReceived();
      auto RPSends          = RP.getNumberOfSends();
      auto RPReceives       = RP.getNumberOfReceives();
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

      auto tp_reduced = tp.getReduced();
      WALBERLA_LOG_INFO_ON_ROOT(*tp_reduced);

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

      reportOverRanks("Total Particles:", real_c(ps->size()));
      reportOverRanks("Number of Particles:", real_c(numParticles));
      reportOverRanks("Number of Ghost Particles:", real_c(numGhostParticles));
      auto estGhostParticles = (localDomain.xSize()+2) * (localDomain.ySize()+2) * (localDomain.zSize()+2);
      estGhostParticles -= localDomain.xSize() * localDomain.ySize() * localDomain.zSize();
      estGhostParticles -= 4*localDomain.xSize() + 4*localDomain.ySize() + 4*localDomain.zSize();
      estGhostParticles -= 8;
      WALBERLA_LOG_DEVEL_ON_ROOT("Estimation: " << estGhostParticles);

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

   LIKWID_MARKER_CLOSE;

   return EXIT_SUCCESS;
}

} // namespace mesa_pd
} // namespace walberla

int main( int argc, char* argv[] )
{
   return walberla::mesa_pd::main( argc, argv );
}
