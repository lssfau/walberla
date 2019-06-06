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

#include <mesa_pd/vtk/ParticleVtkOutput.h>

#include <mesa_pd/collision_detection/AnalyticContactDetection.h>
#include <mesa_pd/data/LinkedCells.h>
#include <mesa_pd/data/ParticleAccessor.h>
#include <mesa_pd/data/ParticleStorage.h>
#include <mesa_pd/data/ShapeStorage.h>
#include <mesa_pd/domain/BlockForestDomain.h>
#include <mesa_pd/kernel/DoubleCast.h>
#include <mesa_pd/kernel/ExplicitEulerWithShape.h>
#include <mesa_pd/kernel/InsertParticleIntoLinkedCells.h>
#include <mesa_pd/kernel/ParticleSelector.h>
#include <mesa_pd/kernel/SpringDashpot.h>
#include <mesa_pd/mpi/ContactFilter.h>
#include <mesa_pd/mpi/ReduceProperty.h>
#include <mesa_pd/mpi/SyncNextNeighbors.h>

#include <mesa_pd/mpi/notifications/ForceTorqueNotification.h>

#include <blockforest/BlockForest.h>
#include <blockforest/Initialization.h>
#include <core/Abort.h>
#include <core/Environment.h>
#include <core/Hostname.h>
#include <core/math/Random.h>
#include <core/mpi/Gatherv.h>
#include <core/mpi/RecvBuffer.h>
#include <core/mpi/Reduce.h>
#include <core/mpi/SendBuffer.h>
#include <core/grid_generator/SCIterator.h>
#include <core/logging/Logging.h>
#include <core/OpenMP.h>
#include <core/timing/Timer.h>
#include <core/timing/TimingPool.h>
#include <core/waLBerlaBuildInfo.h>
#include <postprocessing/sqlite/SQLite.h>
#include <vtk/VTKOutput.h>

#include <functional>
#include <memory>
#include <string>
#include <type_traits>

namespace walberla {
namespace mesa_pd {

class SelectRank
{
public:
   using return_type = int;
   int operator()(const data::Particle& /*p*/) const { return rank_; }
   int operator()(const data::Particle&& /*p*/) const { return rank_; }
private:
   int rank_ = walberla::mpi::MPIManager::instance()->rank();
};

struct Contact
{
   Contact(const size_t idx1,
           const size_t idx2,
           const Vec3   contactPoint,
           const Vec3   contactNormal,
           const real_t penetrationDepth)
      : idx1_(idx1)
      , idx2_(idx2)
      , contactPoint_(contactPoint)
      , contactNormal_(contactNormal)
      , penetrationDepth_(penetrationDepth) {}

   size_t idx1_;
   size_t idx2_;
   Vec3   contactPoint_;
   Vec3   contactNormal_;
   real_t penetrationDepth_;
};


class ParticleAccessorWithShape : public data::ParticleAccessor
{
public:
   ParticleAccessorWithShape(std::shared_ptr<data::ParticleStorage>& ps, std::shared_ptr<data::ShapeStorage>& ss)
      : ParticleAccessor(ps)
      , ss_(ss)
   {}

   const walberla::real_t& getInvMass(const size_t p_idx) const {return ss_->shapes[ps_->getShapeIDRef(p_idx)]->getInvMass();}
   walberla::real_t& getInvMassRef(const size_t p_idx) {return ss_->shapes[ps_->getShapeIDRef(p_idx)]->getInvMass();}
   void setInvMass(const size_t p_idx, const walberla::real_t& v) { ss_->shapes[ps_->getShapeIDRef(p_idx)]->getInvMass() = v;}

   const auto& getInvInertiaBF(const size_t p_idx) const {return ss_->shapes[ps_->getShapeIDRef(p_idx)]->getInvInertiaBF();}
   auto& getInvInertiaBFRef(const size_t p_idx) {return ss_->shapes[ps_->getShapeIDRef(p_idx)]->getInvInertiaBF();}
   void setInvInertiaBF(const size_t p_idx, const Mat3& v) { ss_->shapes[ps_->getShapeIDRef(p_idx)]->getInvInertiaBF() = v;}

   data::BaseShape* getShape(const size_t p_idx) const {return ss_->shapes[ps_->getShapeIDRef(p_idx)].get();}
private:
   std::shared_ptr<data::ShapeStorage> ss_;
};

void createPlane( data::ParticleStorage& ps,
                  data::ShapeStorage& ss,
                  const Vec3& pos,
                  const Vec3& normal )
{
   auto p0              = ps.create(true);
   p0->getPositionRef() = pos;
   p0->getShapeIDRef()  = ss.create<data::HalfSpace>( normal );
   p0->getOwnerRef()    = walberla::mpi::MPIManager::instance()->rank();
   p0->getTypeRef()     = 0;
   data::particle_flags::set(p0->getFlagsRef(), data::particle_flags::INFINITE);
   data::particle_flags::set(p0->getFlagsRef(), data::particle_flags::FIXED);
   data::particle_flags::set(p0->getFlagsRef(), data::particle_flags::NON_COMMUNICATING);
}

std::string envToString(const char* env)
{
   return env != nullptr ? std::string(env) : "";
}

void storeNodeTimings( const uint_t                 runId,
                       const std::string          & dbFile,
                       const std::string          & tableName,
                       const WcTimingPool         & tp )
{
   std::map< std::string, walberla::int64_t > integerProperties;
   std::map< std::string, double >            realProperties;
   std::map< std::string, std::string >       stringProperties;

   walberla::mpi::SendBuffer sb;
   walberla::mpi::RecvBuffer rb;

   sb << walberla::getHostName();
   sb << int64_t(walberla::mpi::MPIManager::instance()->rank());
   sb << tp;

   walberla::mpi::gathervBuffer(sb, rb);

   WALBERLA_ROOT_SECTION()
   {
      while (!rb.isEmpty())
      {
         integerProperties.clear();
         realProperties.clear();
         stringProperties.clear();

         std::string  hostname;
         int64_t      rank;
         WcTimingPool cTP;
         rb >> hostname;
         rb >> rank;
         rb >> cTP;

         stringProperties["hostname"] = hostname;
         integerProperties["rank"]    = rank;
         for (auto& v : cTP)
         {
            realProperties[v.first] = v.second.average();
         }

         postprocessing::storeAdditionalRunInfoInSqliteDB( runId,
                                                           dbFile,
                                                           tableName,
                                                           integerProperties,
                                                           stringProperties,
                                                           realProperties );
      }
   }
}


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

   const std::string host = mainConf.getParameter<std::string>("host", "none" );
   WALBERLA_LOG_INFO_ON_ROOT("host: " << host);

   const int jobid = mainConf.getParameter<int>("jobid", 0 );
   WALBERLA_LOG_INFO_ON_ROOT("jobid: " << jobid);

   const real_t spacing = mainConf.getParameter<real_t>("spacing", real_t(1.0) );
   WALBERLA_LOG_INFO_ON_ROOT("spacing: " << spacing);

   const real_t radius = mainConf.getParameter<real_t>("radius", real_t(0.5) );
   WALBERLA_LOG_INFO_ON_ROOT("radius: " << radius);

   bool bBarrier = mainConf.getParameter<bool>("bBarrier", false );
   WALBERLA_LOG_INFO_ON_ROOT("bBarrier: " << bBarrier);

   int64_t numOuterIterations = mainConf.getParameter<int64_t>("numOuterIterations", 10 );
   WALBERLA_LOG_INFO_ON_ROOT("numOuterIterations: " << numOuterIterations);

   int64_t simulationSteps = mainConf.getParameter<int64_t>("simulationSteps", 10 );
   WALBERLA_LOG_INFO_ON_ROOT("simulationSteps: " << simulationSteps);

   real_t dt = mainConf.getParameter<real_t>("dt", real_c(0.01) );
   WALBERLA_LOG_INFO_ON_ROOT("dt: " << dt);

   const int visSpacing = mainConf.getParameter<int>("visSpacing",  1000 );
   WALBERLA_LOG_INFO_ON_ROOT("visSpacing: " << visSpacing);
   const std::string path = mainConf.getParameter<std::string>("path",  "vtk_out" );
   WALBERLA_LOG_INFO_ON_ROOT("path: " << path);

   const std::string sqlFile = mainConf.getParameter<std::string>("sqlFile",  "benchmark.sqlite" );
   WALBERLA_LOG_INFO_ON_ROOT("sqlFile: " << sqlFile);

   WALBERLA_LOG_INFO_ON_ROOT("*** BLOCKFOREST ***");
   // create forest
   auto forest = blockforest::createBlockForestFromConfig( mainConf );
   if (!forest)
   {
      WALBERLA_LOG_INFO_ON_ROOT( "No BlockForest created ... exiting!");
      return EXIT_SUCCESS;
   }
   domain::BlockForestDomain domain(forest);

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
   data::LinkedCells     lc(localDomain.getExtended(spacing), spacing );

   auto  smallSphere = ss->create<data::Sphere>( radius );
   ss->shapes[smallSphere]->updateMassAndInertia(real_t(2707));
   for (auto& iBlk : *forest)
   {
      for (auto pt : grid_generator::SCGrid(iBlk.getAABB(), Vector3<real_t>(spacing, spacing, spacing) * real_c(0.5), spacing))
      {
         WALBERLA_CHECK(iBlk.getAABB().contains(pt));

         auto p                       = ps->create();
         p->getPositionRef()          = pt;
         p->getInteractionRadiusRef() = radius;
         p->getShapeIDRef()           = smallSphere;
         p->getOwnerRef()             = mpiManager->rank();
         p->getTypeRef()              = 0;
      }
   }
   int64_t numParticles = int64_c(ps->size());
   walberla::mpi::reduceInplace(numParticles, walberla::mpi::SUM);
   WALBERLA_LOG_INFO_ON_ROOT("#particles created: " << numParticles);

   auto shift = (spacing - radius - radius) * real_t(0.5);
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
   auto vtkDomainOutput = walberla::vtk::createVTKOutput_DomainDecomposition( forest, "domain_decomposition", 1, "vtk_out", "simulation_step" );
   auto vtkOutput       = make_shared<mesa_pd::vtk::ParticleVtkOutput>(ps) ;
   auto vtkWriter       = walberla::vtk::createVTKOutput_PointData(vtkOutput, "Bodies", 1, "vtk", "simulation_step", false, false);
   vtkOutput->addOutput<SelectRank>("rank");
   vtkOutput->addOutput<data::SelectParticleOwner>("owner");
   //   vtkDomainOutput->write();

   WALBERLA_LOG_INFO_ON_ROOT("*** SIMULATION - START ***");
   // Init kernels
   kernel::ExplicitEulerWithShape        explicitEulerWithShape( dt );
   kernel::InsertParticleIntoLinkedCells ipilc;
   kernel::SpringDashpot                 dem(1);
   dem.setStiffness(0, 0, real_t(0));
   dem.setDampingN (0, 0, real_t(0));
   dem.setDampingT (0, 0, real_t(0));
   dem.setFriction (0, 0, real_t(0));
   collision_detection::AnalyticContactDetection              acd;
   kernel::DoubleCast                    double_cast;
   mpi::ContactFilter                    contact_filter;
   mpi::ReduceProperty                   RP;
   mpi::SyncNextNeighbors                SNN;
   std::vector<Contact>                  contacts;
   contacts.reserve(4000000);

   // initial sync
   SNN(*ps, domain);

   for (int64_t outerIteration = 0; outerIteration < numOuterIterations; ++outerIteration)
   {
      WALBERLA_LOG_INFO_ON_ROOT("*** RUNNING OUTER ITERATION " << outerIteration << " ***");

      WcTimingPool tp;

      WALBERLA_MPI_BARRIER();
      tp["GenerateLinkedCells"].start();
      for (int64_t i=0; i < simulationSteps; ++i)
      {
         lc.clear();
         ps->forEachParticle(true, kernel::SelectAll(), accessor, ipilc, accessor, lc);
      }
      tp["GenerateLinkedCells"].end();

      int64_t contactsChecked  = 0;
      int64_t contactsDetected = 0;
      int64_t contactsTreated  = 0;
      WALBERLA_MPI_BARRIER();
      tp["ContactDetection"].start();
      for (int64_t i=0; i < simulationSteps; ++i)
      {
         contacts.clear();
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
               if (contact_filter(acd.getIdx1(), acd.getIdx2(), ac, acd.getContactPoint(), domain))
               {
                  ++contactsTreated;
                  contacts.emplace_back(acd.getIdx1(), acd.getIdx2(), acd.getContactPoint(), acd.getContactNormal(), acd.getPenetrationDepth());
               }
            }
         },
         accessor );
      }
      tp["ContactDetection"].end();

      WALBERLA_MPI_BARRIER();
      tp["DEM"].start();
      for (int64_t i=0; i < simulationSteps; ++i)
      {
         for (auto& c : contacts)
         {
            dem(c.idx1_, c.idx2_, accessor, c.contactPoint_, c.contactNormal_, c.penetrationDepth_);
         }
      }
      tp["DEM"].end();

      WALBERLA_MPI_BARRIER();
      tp["ReduceForce"].start();
      for (int64_t i=0; i < simulationSteps; ++i)
      {
         RP.operator()<ForceTorqueNotification>(*ps);
      }
      tp["ReduceForce"].end();

      WALBERLA_MPI_BARRIER();
      tp["Euler"].start();
      for (int64_t i=0; i < simulationSteps; ++i)
      {
         ps->forEachParticle(true, kernel::SelectLocal(), accessor, explicitEulerWithShape, accessor);
      }
      tp["Euler"].end();

      WALBERLA_MPI_BARRIER();
      tp["SNN"].start();
      for (int64_t i=0; i < simulationSteps; ++i)
      {
         SNN(*ps, domain);
      }
      tp["SNN"].end();

      WALBERLA_LOG_INFO_ON_ROOT("*** SIMULATION - END ***");

      WALBERLA_LOG_INFO_ON_ROOT("*** CHECKING RESULT - START ***");
      auto pIt = ps->begin();
      for (auto& iBlk : *forest)
      {
         for (auto it = grid_generator::SCIterator(iBlk.getAABB(), Vector3<real_t>(spacing, spacing, spacing) * real_c(0.5), spacing);
              it != grid_generator::SCIterator();
              ++it, ++pIt)
         {
            WALBERLA_CHECK_UNEQUAL(pIt, ps->end());
            WALBERLA_CHECK_FLOAT_EQUAL((*pIt).getPositionRef(), *it);
         }
      }
      WALBERLA_LOG_INFO_ON_ROOT("*** CHECKING RESULT - END ***");

      WALBERLA_LOG_INFO_ON_ROOT("*** SQL OUTPUT - START ***");
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
      walberla::mpi::reduceInplace(numParticles, walberla::mpi::SUM);
      walberla::mpi::reduceInplace(numGhostParticles, walberla::mpi::SUM);
      walberla::mpi::reduceInplace(contactsChecked, walberla::mpi::SUM);
      walberla::mpi::reduceInplace(contactsDetected, walberla::mpi::SUM);
      walberla::mpi::reduceInplace(contactsTreated, walberla::mpi::SUM);
      double linkedCellsVolume = lc.domain_.volume();
      walberla::mpi::reduceInplace(linkedCellsVolume, walberla::mpi::SUM);
      size_t numLinkedCells = lc.cells_.size();
      walberla::mpi::reduceInplace(numLinkedCells, walberla::mpi::SUM);
      size_t local_aabbs         = domain.getNumLocalAABBs();
      size_t neighbor_subdomains = domain.getNumNeighborSubdomains();
      size_t neighbor_processes  = domain.getNumNeighborProcesses();
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
         stringProperties["host"]                 = host;
         integerProperties["jobid"]               = jobid;
         integerProperties["mpi_num_processes"]   = mpiManager->numProcesses();
         integerProperties["omp_max_threads"]     = omp_get_max_threads();
         integerProperties["outerIteration"]      = int64_c(outerIteration);
         integerProperties["numOuterIterations"]  = numOuterIterations;
         integerProperties["simulationSteps"]     = simulationSteps;
         integerProperties["bBarrier"]            = int64_c(bBarrier);
         integerProperties["num_particles"]       = numParticles;
         integerProperties["num_ghost_particles"] = numGhostParticles;
         integerProperties["contacts_checked"]    = contactsChecked;
         integerProperties["contacts_detected"]   = contactsDetected;
         integerProperties["contacts_treated"]    = contactsTreated;
         integerProperties["blocks_x"]            = int64_c(forest->getXSize());
         integerProperties["blocks_y"]            = int64_c(forest->getXSize());
         integerProperties["blocks_z"]            = int64_c(forest->getXSize());
         realProperties["domain_x"]               = double_c(forest->getDomain().xSize());
         realProperties["domain_y"]               = double_c(forest->getDomain().ySize());
         realProperties["domain_z"]               = double_c(forest->getDomain().zSize());
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
         stringProperties["SLURM_CLUSTER_NAME"]       = envToString(std::getenv( "SLURM_CLUSTER_NAME" ));
         stringProperties["SLURM_CPUS_ON_NODE"]       = envToString(std::getenv( "SLURM_CPUS_ON_NODE" ));
         stringProperties["SLURM_CPUS_PER_TASK"]      = envToString(std::getenv( "SLURM_CPUS_PER_TASK" ));
         stringProperties["SLURM_JOB_ACCOUNT"]        = envToString(std::getenv( "SLURM_JOB_ACCOUNT" ));
         stringProperties["SLURM_JOB_ID"]             = envToString(std::getenv( "SLURM_JOB_ID" ));
         stringProperties["SLURM_JOB_CPUS_PER_NODE"]  = envToString(std::getenv( "SLURM_JOB_CPUS_PER_NODE" ));
         stringProperties["SLURM_JOB_NAME"]           = envToString(std::getenv( "SLURM_JOB_NAME" ));
         stringProperties["SLURM_JOB_NUM_NODES"]      = envToString(std::getenv( "SLURM_JOB_NUM_NODES" ));
         stringProperties["SLURM_NTASKS"]             = envToString(std::getenv( "SLURM_NTASKS" ));
         stringProperties["SLURM_NTASKS_PER_CORE"]    = envToString(std::getenv( "SLURM_NTASKS_PER_CORE" ));
         stringProperties["SLURM_NTASKS_PER_NODE"]    = envToString(std::getenv( "SLURM_NTASKS_PER_NODE" ));
         stringProperties["SLURM_NTASKS_PER_SOCKET"]  = envToString(std::getenv( "SLURM_NTASKS_PER_SOCKET" ));
         stringProperties["SLURM_TASKS_PER_NODE"]     = envToString(std::getenv( "SLURM_TASKS_PER_NODE" ));


         runId = postprocessing::storeRunInSqliteDB( sqlFile, integerProperties, stringProperties, realProperties );
         postprocessing::storeTimingPoolInSqliteDB( sqlFile, runId, *tp_reduced, "Timeloop" );
      }
      storeNodeTimings(runId, sqlFile, "NodeTiming", tp);
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
