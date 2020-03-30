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
//! \file   Mixer.cpp
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================
#include <mesa_pd/vtk/ParticleVtkOutput.h>

#include <mesa_pd/collision_detection/BroadPhase.h>
#include <mesa_pd/collision_detection/AnalyticContactDetection.h>

#include <mesa_pd/data/LinkedCells.h>
#include <mesa_pd/data/ParticleAccessorWithShape.h>
#include <mesa_pd/data/ParticleStorage.h>

#include <mesa_pd/data/ShapeStorage.h>

#include <mesa_pd/domain/BlockForestDomain.h>

#include <mesa_pd/kernel/DoubleCast.h>
#include <mesa_pd/kernel/ExplicitEuler.h>
#include <mesa_pd/kernel/InsertParticleIntoLinkedCells.h>
#include <mesa_pd/kernel/ParticleSelector.h>
#include <mesa_pd/kernel/SpringDashpot.h>

#include <mesa_pd/mpi/ContactFilter.h>
#include <mesa_pd/mpi/SyncNextNeighbors.h>
#include <mesa_pd/mpi/ReduceProperty.h>

#include <mesa_pd/mpi/notifications/ForceTorqueNotification.h>

#include <mesa_pd/sorting/LinearizedCompareFunctor.h>

#include <blockforest/Initialization.h>
#include <core/Abort.h>
#include <core/Environment.h>
#include <core/math/Random.h>
#include <core/mpi/Reduce.h>
#include <core/grid_generator/SCIterator.h>
#include <core/logging/Logging.h>
#include <core/OpenMP.h>
#include <core/timing/Timer.h>
#include <core/timing/TimingPool.h>
#include <core/waLBerlaBuildInfo.h>
#include <sqlite/SQLite.h>
#include <vtk/VTKOutput.h>

#include <functional>
#include <memory>
#include <string>

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

class SelectIdx
{
public:
   using return_type = int;
   auto operator()(const data::Particle& p) const { return p.getIdx(); }
   auto operator()(const data::Particle&& p) const { return p.getIdx(); }
private:
};

class SelectGhost
{
public:
   using return_type = int;
   bool operator()(const data::Particle& p) const { return data::particle_flags::isSet(p.getFlags(), data::particle_flags::GHOST); }
   bool operator()(const data::Particle&& p) const { return data::particle_flags::isSet(p.getFlags(), data::particle_flags::GHOST); }
};

class SelectRotation
{
public:
   using return_type = Vec3;
   auto operator()(const data::Particle& p) const { return p.getRotation().getMatrix() * Vec3(1,0,0); }
   auto operator()(const data::Particle&& p) const { return p.getRotation().getMatrix() * Vec3(1,0,0); }
};

auto createBoundary( const std::shared_ptr<data::ParticleStorage>& ps,
                     const std::shared_ptr<data::ShapeStorage>& ss,
                     const Vec3&  pos,
                     const real_t radius,
                     const Vec3&  axis)
{
   auto p0              = ps->create(true);
   p0->setPosition( pos );
   auto cb = ss->create<data::CylindricalBoundary>( radius, axis );
   ss->shapes[cb]->updateMassAndInertia(std::numeric_limits<real_t>::infinity());
   p0->setShapeID( cb );
   p0->setOwner( walberla::mpi::MPIManager::instance()->rank() );
   p0->setType( 0 );
   data::particle_flags::set(p0->getFlagsRef(), data::particle_flags::GLOBAL);
   data::particle_flags::set(p0->getFlagsRef(), data::particle_flags::INFINITE);
   data::particle_flags::set(p0->getFlagsRef(), data::particle_flags::FIXED);
   data::particle_flags::set(p0->getFlagsRef(), data::particle_flags::NON_COMMUNICATING);
   return p0;
}

auto createPlane( const std::shared_ptr<data::ParticleStorage>& ps,
                  const std::shared_ptr<data::ShapeStorage>& ss,
                  const Vec3& pos,
                  const Vec3& normal)
{
   auto p0              = ps->create(true);
   p0->setPosition( pos );
   auto plane = ss->create<data::HalfSpace>( normal );
   ss->shapes[plane]->updateMassAndInertia(std::numeric_limits<real_t>::infinity());
   p0->setShapeID( plane );
   p0->setOwner( walberla::mpi::MPIManager::instance()->rank() );
   p0->setType( 0 );
   data::particle_flags::set(p0->getFlagsRef(), data::particle_flags::GLOBAL);
   data::particle_flags::set(p0->getFlagsRef(), data::particle_flags::INFINITE);
   data::particle_flags::set(p0->getFlagsRef(), data::particle_flags::FIXED);
   data::particle_flags::set(p0->getFlagsRef(), data::particle_flags::NON_COMMUNICATING);
   return p0;
}

int main( int argc, char ** argv )
{
   using namespace walberla::timing;

   Environment env(argc, argv);
   walberla::mpi::MPIManager::instance()->useWorldComm();

   logging::Logging::instance()->setStreamLogLevel(logging::Logging::INFO);
   logging::Logging::instance()->setFileLogLevel(logging::Logging::INFO);

   WALBERLA_LOG_INFO_ON_ROOT( "config file: " << argv[1] );
   WALBERLA_LOG_INFO_ON_ROOT( "waLBerla Revision: " << WALBERLA_GIT_SHA1 );

   math::seedRandomGenerator( static_cast<unsigned int>(1337 * walberla::mpi::MPIManager::instance()->worldRank()) );

   WALBERLA_LOG_INFO_ON_ROOT("*** READING CONFIG FILE ***");
   auto cfg = env.config();
   if (cfg == nullptr) WALBERLA_ABORT("No config specified!");
   const Config::BlockHandle mainConf  = cfg->getBlock( "MIXER" );

   const real_t spacing = mainConf.getParameter<real_t>("spacing", real_t(1.0) );
   WALBERLA_LOG_INFO_ON_ROOT("spacing: " << spacing);

   const real_t radius = mainConf.getParameter<real_t>("radius", real_t(0.5) );
   WALBERLA_LOG_INFO_ON_ROOT("radius: " << radius);

   const real_t vMax = mainConf.getParameter<real_t>("vMax", real_t(0.5) );
   WALBERLA_LOG_INFO_ON_ROOT("vMax: " << vMax);

   const real_t rotationSpeed = mainConf.getParameter<real_t>("rotationSpeed", real_t(0.5) );
   WALBERLA_LOG_INFO_ON_ROOT("rotationSpeed: " << rotationSpeed);

   const Vec3 gravity = mainConf.getParameter<Vec3>("gravity", Vec3(0, 0, real_t(9.81)));
   WALBERLA_LOG_INFO_ON_ROOT("rotationSpeed: " << rotationSpeed);

   int64_t simulationSteps = mainConf.getParameter<int64_t>("simulationSteps", 1000 );
   WALBERLA_LOG_INFO_ON_ROOT("simulationSteps: " << simulationSteps);

   real_t dt = mainConf.getParameter<real_t>("dt", real_c(0.01) );
   WALBERLA_LOG_INFO_ON_ROOT("dt: " << dt);

   real_t cor = mainConf.getParameter<real_t>("cor", real_c(0.95) );
   WALBERLA_LOG_INFO_ON_ROOT("cor: " << cor);

   real_t ct_in_dt = mainConf.getParameter<real_t>("ct_in_dt", real_c(20) );
   WALBERLA_LOG_INFO_ON_ROOT("ct_in_dt: " << ct_in_dt);

   real_t density = mainConf.getParameter<real_t>("density", real_c(2707) );
   WALBERLA_LOG_INFO_ON_ROOT("density: " << density);

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
   data::ParticleAccessorWithShape accessor(ps, ss);
   data::LinkedCells     lc(localDomain.getExtended(spacing), spacing+spacing );

   auto  smallSphere = ss->create<data::Sphere>( radius );
   ss->shapes[smallSphere]->updateMassAndInertia(density);
   for (auto& iBlk : *forest)
   {
      for (auto pt : grid_generator::SCGrid(iBlk.getAABB(), Vector3<real_t>(spacing, spacing, spacing) * real_c(0.5), spacing))
      {
         WALBERLA_CHECK(iBlk.getAABB().contains(pt));
         auto maxDist = real_t(simulationDomain.xSize()) * real_t(0.5) - spacing;
         auto dist = (pt - simulationDomain.center());
         dist[2] = real_t(0);
         if (dist.sqrLength() > maxDist*maxDist) continue;

         auto p                       = ps->create();
         p->getPositionRef()          = pt;
         p->getInteractionRadiusRef() = radius;
         p->setLinearVelocity( Vec3(math::realRandom(-vMax, vMax),math::realRandom(-vMax, vMax),math::realRandom(-vMax, vMax)) );
         p->getShapeIDRef()           = smallSphere;
         p->getOwnerRef()             = walberla::mpi::MPIManager::instance()->rank();
         p->getTypeRef()              = 0;
      }
   }
   int64_t numParticles = int64_c(ps->size());
   walberla::mpi::reduceInplace(numParticles, walberla::mpi::SUM);
   WALBERLA_LOG_INFO_ON_ROOT("#particles created: " << numParticles);

   auto boundary = createBoundary(ps,
                                  ss,
                                  simulationDomain.center(),
                                  real_t(simulationDomain.xSize()) * real_t(0.5),
                                  Vec3(0,0,1));
   WALBERLA_UNUSED(boundary);
   //boundary->setAngularVelocity(Vec3(real_t(0), real_t(0), rotationSpeed));
   createPlane(ps, ss, simulationDomain.minCorner(), Vec3(0,0,1));
   createPlane(ps, ss, simulationDomain.maxCorner(), Vec3(0,0,-1));

   Rot3 dRot;
   Vec3 origin = simulationDomain.center();
   Vec3 dp;

   auto mixingBlade = ss->create<data::Box>( Vec3(0.009,simulationDomain.ySize(),simulationDomain.zSize()) );
   ss->shapes[mixingBlade]->updateMassAndInertia(std::numeric_limits<real_t>::infinity());
   auto box0                       = ps->create();
   box0->getPositionRef()          = Vec3(simulationDomain.xSize() * real_t(0.5),0.0,0.0);
   box0->setShapeID( mixingBlade );
   box0->getOwnerRef()             = walberla::mpi::MPIManager::instance()->rank();
   box0->getTypeRef()              = 0;
   data::particle_flags::set(box0->getFlagsRef(), data::particle_flags::GLOBAL);
   data::particle_flags::set(box0->getFlagsRef(), data::particle_flags::INFINITE);
   data::particle_flags::set(box0->getFlagsRef(), data::particle_flags::FIXED);
   data::particle_flags::set(box0->getFlagsRef(), data::particle_flags::NON_COMMUNICATING);
   box0->getRotationRef().rotate(Vec3(0,1,0), math::pi * 0.15 );
//   box0->getRotationRef().rotate(Vec3(0,0,1), -math::pi * 0.15 );

   auto box1                       = ps->create();
   box1->getPositionRef()          = Vec3(simulationDomain.xSize() * real_t(0.5),0.0,0.0);
   box1->setShapeID( mixingBlade );
   box1->getOwnerRef()             = walberla::mpi::MPIManager::instance()->rank();
   box1->getTypeRef()              = 0;
   data::particle_flags::set(box1->getFlagsRef(), data::particle_flags::GLOBAL);
   data::particle_flags::set(box1->getFlagsRef(), data::particle_flags::INFINITE);
   data::particle_flags::set(box1->getFlagsRef(), data::particle_flags::FIXED);
   data::particle_flags::set(box1->getFlagsRef(), data::particle_flags::NON_COMMUNICATING);
   box1->getRotationRef().rotate(Vec3(0,1,0), math::pi * 0.15 );
//   box1->getRotationRef().rotate(Vec3(0,0,1), -math::pi * 0.15 );
   dp = ( box1->getPosition() - origin );
   dRot = Rot3(Vec3(real_t(0), real_t(0), math::pi * 0.5));
   box1->setPosition( origin + dRot.getMatrix() * dp );
   box1->getRotationRef().rotate(dRot);

   auto box2                       = ps->create();
   box2->getPositionRef()          = Vec3(simulationDomain.xSize() * real_t(0.5),0.0,0.0);
   box2->setShapeID( mixingBlade );
   box2->getOwnerRef()             = walberla::mpi::MPIManager::instance()->rank();
   box2->getTypeRef()              = 0;
   data::particle_flags::set(box2->getFlagsRef(), data::particle_flags::GLOBAL);
   data::particle_flags::set(box2->getFlagsRef(), data::particle_flags::INFINITE);
   data::particle_flags::set(box2->getFlagsRef(), data::particle_flags::FIXED);
   data::particle_flags::set(box2->getFlagsRef(), data::particle_flags::NON_COMMUNICATING);
   box2->getRotationRef().rotate(Vec3(0,1,0), math::pi * 0.15 );
//   box2->getRotationRef().rotate(Vec3(0,0,1), -math::pi * 0.15 );
   dp = ( box2->getPosition() - origin );
   dRot = Rot3(Vec3(real_t(0), real_t(0), math::pi));
   box2->setPosition( origin + dRot.getMatrix() * dp );
   box2->getRotationRef().rotate(dRot);

   auto box3                       = ps->create();
   box3->getPositionRef()          = Vec3(simulationDomain.xSize() * real_t(0.5),0.0,0.0);
   box3->setShapeID( mixingBlade );
   box3->getOwnerRef()             = walberla::mpi::MPIManager::instance()->rank();
   box3->getTypeRef()              = 0;
   data::particle_flags::set(box3->getFlagsRef(), data::particle_flags::GLOBAL);
   data::particle_flags::set(box3->getFlagsRef(), data::particle_flags::INFINITE);
   data::particle_flags::set(box3->getFlagsRef(), data::particle_flags::FIXED);
   data::particle_flags::set(box3->getFlagsRef(), data::particle_flags::NON_COMMUNICATING);
   box3->getRotationRef().rotate(Vec3(0,1,0), math::pi * 0.15 );
//   box3->getRotationRef().rotate(Vec3(0,0,1), -math::pi * 0.15 );
   dp = ( box3->getPosition() - origin );
   dRot = Rot3(Vec3(real_t(0), real_t(0), math::pi * 1.5));
   box3->setPosition( origin + dRot.getMatrix() * dp );
   box3->getRotationRef().rotate(dRot);

   WALBERLA_LOG_INFO_ON_ROOT("*** SETUP - END ***");

   WALBERLA_LOG_INFO_ON_ROOT("*** VTK ***");
   auto vtkDomainOutput = walberla::vtk::createVTKOutput_DomainDecomposition( forest, "domain_decomposition", 1, "vtk", "simulation_step" );
   auto vtkOutput       = make_shared<mesa_pd::vtk::ParticleVtkOutput>(ps) ;
   auto vtkWriter       = walberla::vtk::createVTKOutput_PointData(vtkOutput, "Bodies", 1, "vtk", "simulation_step", false, false);
   vtkOutput->addOutput<SelectGhost>("ghost");
   vtkOutput->addOutput<SelectRank>("rank");
   vtkOutput->addOutput<SelectIdx>("idx");
//   vtkOutput->addOutput<SelectRotation>("rot");
//   auto select_scaling = std::make_shared<vtk::OutputSelector<SelectScaling>>(SelectScaling(ss));
//   vtkOutput->addOutput("scale", select_scaling);
   vtkOutput->setParticleSelector([smallSphere](const data::ParticleStorage::iterator& pIt){ return pIt->getShapeID() == smallSphere;});
   vtkDomainOutput->write();

   WALBERLA_LOG_INFO_ON_ROOT("*** SIMULATION - START ***");
   // Init kernels
   kernel::ExplicitEuler                 explicitEuler( dt );
   kernel::InsertParticleIntoLinkedCells ipilc;
   kernel::SpringDashpot                 dem(1);
   const auto ct  = dt * ct_in_dt;
   const auto stiffness = (math::pi*math::pi + std::log(cor)*std::log(cor)) / (ct*ct);
   const auto damping   = (real_t(2) * std::log(cor)) / (ct);
   dem.setStiffness(0, 0, real_t(8.11e6));
   dem.setDampingN (0, 0, real_t(6.86e1));
   dem.setDampingT (0, 0, real_t(6.86e1));
   dem.setFriction (0, 0, real_t(0.4));

   collision_detection::AnalyticContactDetection acd;
   kernel::DoubleCast                 double_cast;
   mpi::ContactFilter                 contact_filter;
   mpi::ReduceProperty                RP;
   mpi::SyncNextNeighbors             SNN;

   // initial sync
   SNN(*ps, domain);

   std::ofstream fout;
   WALBERLA_ROOT_SECTION()
   {
      fout.open("timings.txt");
   }
   dRot = Rot3(Vec3(real_t(0), real_t(0), rotationSpeed) * dt);
   for (int64_t i=0; i < simulationSteps; ++i)
   {
      if (i % visSpacing == 0)
      {
//         vtkWriter->write();
         WALBERLA_ROOT_SECTION()
         {
            fout << std::setprecision( 20 ) << timing::WcPolicy::getTimestamp() << "\n";
         }
      }

      if (i == 80000)
      {
         WcTimer timer;
         sorting::LinearizedCompareFunctor linear(lc.domain_, lc.numCellsPerDim_);
         ps->sort(linear);
         timer.end();
         WALBERLA_LOG_DEVEL_ON_ROOT("time needed for sorting: " << timer.total());
      }

      ps->forEachParticle(false,
                          kernel::SelectLocal(),
                          accessor,
                          [&](const size_t idx, auto& ac){ac.setForce(idx, gravity );},
      accessor);

      lc.clear();
      ps->forEachParticle(true,
                          kernel::SelectAll(),
                          accessor,
                          ipilc,
                          accessor,
                          lc);

      lc.forEachParticlePairHalf(true,
                                 kernel::SelectAll(),
                                 accessor,
                                 [&](const size_t idx1, const size_t idx2, auto& ac)
      {
         if ((ac.getShapeID(idx1) != smallSphere) && (ac.getShapeID(idx2) != smallSphere))
         {

         } else
         {
            if (collision_detection::isInInteractionDistance(idx1, idx2, ac))
            {
               if (double_cast(idx1, idx2, ac, acd, ac ))
               {
                  if (contact_filter(acd.getIdx1(), acd.getIdx2(), ac, acd.getContactPoint(), domain))
                  {
                     auto meff = real_t(1) / (ac.getInvMass(idx1) + ac.getInvMass(idx2));
                     dem.setStiffness(0, 0, stiffness * meff);
                     dem.setDampingN(0, 0, damping * meff);
                     dem(acd.getIdx1(), acd.getIdx2(), ac, acd.getContactPoint(), acd.getContactNormal(), acd.getPenetrationDepth());
                  }
               }
            }
         }
      },
      accessor );

      RP.operator()<ForceTorqueNotification>(*ps);

      ps->forEachParticle(true,
                          kernel::SelectLocal(),
                          accessor,
                          explicitEuler,
                          accessor);

      ps->forEachParticle(true,
                          kernel::SelectAll(),
                          accessor,
                          [&](const size_t idx1, auto& ac)
      {
         if (ac.getShapeID(idx1) == mixingBlade)
         {
            //rotate
            dp = ( ac.getPosition(idx1) - origin );
            ac.setPosition(idx1, origin + dRot.getMatrix() * dp );
            ac.getRotationRef(idx1).rotate(dRot);
         }
      },
      accessor );

      //SNN(*ps, domain);

      if( i % 1000 == 0 )
      {
         WALBERLA_LOG_DEVEL_ON_ROOT( "Timestep " << i << " / " << simulationSteps );
      }
   }
   fout.close();
   WALBERLA_LOG_INFO_ON_ROOT("*** SIMULATION - END ***");

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
   WALBERLA_ROOT_SECTION()
   {
      std::map< std::string, walberla::int64_t > integerProperties;
      std::map< std::string, double >            realProperties;
      std::map< std::string, std::string >       stringProperties;

      stringProperties["walberla_git"]         = WALBERLA_GIT_SHA1;
      stringProperties["tag"]                  = "mesa_pd";
      integerProperties["mpi_num_processes"]   = walberla::mpi::MPIManager::instance()->numProcesses();
      integerProperties["omp_max_threads"]     = omp_get_max_threads();
      integerProperties["simulationSteps"]     = simulationSteps;
      integerProperties["num_particles"]       = numParticles;
      integerProperties["num_ghost_particles"] = numGhostParticles;
      integerProperties["blocks_x"]            = int64_c(forest->getXSize());
      integerProperties["blocks_y"]            = int64_c(forest->getYSize());
      integerProperties["blocks_z"]            = int64_c(forest->getZSize());
      realProperties["domain_x"]               = double_c(forest->getDomain().xSize());
      realProperties["domain_y"]               = double_c(forest->getDomain().ySize());
      realProperties["domain_z"]               = double_c(forest->getDomain().zSize());


      sqlite::storeRunInSqliteDB( sqlFile, integerProperties, stringProperties, realProperties );
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
