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
//! \file   HeatConduction.cpp
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

//======================================================================================================================
// This showcase combines a classic spring-dashpot interaction model with heat conduction and linear thermal expansion.
//
// First, a set of spheres is dropped onto a plate. The setup is periodic in horizontal directions. The plate as well
// as the spheres shares the same initial temperature. After the spheres have settled, the plate is heated and cooled
// down in a sinus wave. Through thermal conduction the heat is transported through the pile of spheres. At the same
// time a linear thermal expansion model makes the spheres grow and shrink. The parameters are chosen for a visual
// effect and are far to large for a realistic physical simulation. The evolution of the simulation is written to disk
// as vtk files which can be visualized with paraview.
//======================================================================================================================

#include <mesa_pd/vtk/ParticleVtkOutput.h>

#include <mesa_pd/collision_detection/BroadPhase.h>

#include <mesa_pd/data/LinkedCells.h>
#include <mesa_pd/data/ParticleAccessorWithShape.h>
#include <mesa_pd/data/ParticleStorage.h>
#include <mesa_pd/data/ShapeStorage.h>

#include <mesa_pd/domain/BlockForestDomain.h>

#include <mesa_pd/kernel/DoubleCast.h>
#include <mesa_pd/kernel/ExplicitEulerWithShape.h>
#include <mesa_pd/kernel/HeatConduction.h>
#include <mesa_pd/kernel/InsertParticleIntoLinkedCells.h>
#include <mesa_pd/kernel/SpringDashpot.h>
#include <mesa_pd/kernel/TemperatureIntegration.h>

#include <mesa_pd/mpi/ContactFilter.h>
#include <mesa_pd/mpi/SyncNextNeighbors.h>
#include <mesa_pd/mpi/ReduceProperty.h>

#include <mesa_pd/mpi/notifications/ForceTorqueNotification.h>
#include <mesa_pd/mpi/notifications/HeatFluxNotification.h>

#include "ContactDetection.h"
#include "ThermalExpansion.h"

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
#include <mesa_pd/kernel/ParticleSelector.h>

namespace walberla {
namespace mesa_pd {

class CustomParticleAccessor : public data::ParticleAccessorWithShape
{
public:
    CustomParticleAccessor(std::shared_ptr<data::ParticleStorage>& ps, std::shared_ptr<data::ShapeStorage>& ss, const real_t radius273K)
    : ParticleAccessorWithShape(ps, ss)
    , radius273K_(radius273K)
    {}

    real_t getRadius273K(const size_t /*p_idx*/) const {return radius273K_;}
private:
    real_t radius273K_;
};

auto createPlane( const std::shared_ptr<data::ParticleStorage>& ps,
                  const std::shared_ptr<data::ShapeStorage>& ss,
                  const Vec3& pos,
                  const Vec3& normal)
{
   auto p0              = ps->create(true);
   p0->setPosition( pos );
   p0->setShapeID( ss->create<data::HalfSpace>( normal ) );
   p0->setOwner( walberla::mpi::MPIManager::instance()->rank() );
   p0->setType( 0 );
   p0->setTemperature( 273 );
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

   WALBERLA_LOG_INFO_ON_ROOT( "config file: " << argv[1] );
   WALBERLA_LOG_INFO_ON_ROOT( "waLBerla Revision: " << WALBERLA_GIT_SHA1 );

   math::seedRandomGenerator( static_cast<unsigned int>(1337 * walberla::mpi::MPIManager::instance()->worldRank()) );

   WALBERLA_LOG_INFO_ON_ROOT("*** READING CONFIG FILE ***");
   auto cfg = env.config();
   if (cfg == nullptr) WALBERLA_ABORT("No config specified!");
   const Config::BlockHandle mainConf  = cfg->getBlock( "HeatConduction" );

   const real_t spacing = mainConf.getParameter<real_t>("spacing", real_t(1.0) );
   WALBERLA_LOG_INFO_ON_ROOT("spacing: " << spacing);

   const real_t radius = mainConf.getParameter<real_t>("radius", real_t(0.5) );
   WALBERLA_LOG_INFO_ON_ROOT("radius: " << radius);

   const real_t vMax = mainConf.getParameter<real_t>("vMax", real_t(0.5) );
   WALBERLA_LOG_INFO_ON_ROOT("vMax: " << vMax);

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
   CustomParticleAccessor accessor(ps, ss, radius);
   data::LinkedCells     lc(localDomain.getExtended(spacing), spacing+spacing );

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
         p->setLinearVelocity( Vec3(math::realRandom(-vMax, vMax),math::realRandom(-vMax, vMax),math::realRandom(-vMax, vMax)) );
         p->getShapeIDRef()           = smallSphere;
         p->getOwnerRef()             = walberla::mpi::MPIManager::instance()->rank();
         p->getTypeRef()              = 0;
         p->setTemperature( 273 );
         p->setRadiusAtTemperature(radius);
      }
   }
   int64_t numParticles = int64_c(ps->size());
   walberla::mpi::reduceInplace(numParticles, walberla::mpi::SUM);
   WALBERLA_LOG_INFO_ON_ROOT("#particles created: " << numParticles);

   auto plane = createPlane(ps, ss, simulationDomain.minCorner(), Vec3(0,0,1));
   createPlane(ps, ss, simulationDomain.maxCorner(), Vec3(0,0,-1));

   WALBERLA_LOG_INFO_ON_ROOT("*** SETUP - END ***");

   WALBERLA_LOG_INFO_ON_ROOT("*** VTK ***");
   auto vtkDomainOutput = walberla::vtk::createVTKOutput_DomainDecomposition( forest, "domain_partitioning", 1, "vtk", "simulation_step" );
   auto vtkOutput       = make_shared<mesa_pd::vtk::ParticleVtkOutput>(ps) ;
   auto vtkWriter       = walberla::vtk::createVTKOutput_PointData(vtkOutput, "particles", 1, "vtk", "simulation_step", false, false);
   vtkOutput->addOutput<data::SelectParticleOwner>("owner");
   vtkOutput->addOutput<data::SelectParticleTemperature>("temperature");
   vtkOutput->addOutput<data::SelectParticleRadiusAtTemperature>("radius");
   vtkOutput->addOutput<data::SelectParticleLinearVelocity>("linVel");
   vtkOutput->setParticleSelector([smallSphere](const data::ParticleStorage::iterator& pIt){ return pIt->getShapeID() == smallSphere;});
   vtkDomainOutput->write();

   WALBERLA_LOG_INFO_ON_ROOT("*** SIMULATION - START ***");
   // Init kernels
   kernel::ExplicitEuler                 explicitEuler( dt );
   kernel::InsertParticleIntoLinkedCells ipilc;
   kernel::HeatConduction                heatConduction(1);
   heatConduction.setConductance(0, 0, real_t(1));
   kernel::SpringDashpot                 dem(1);
   dem.setStiffness(0, 0, real_t(8.11e6));
   dem.setDampingN (0, 0, real_t(6.86e1));
   dem.setDampingT (0, 0, real_t(6.86e1));
   dem.setFriction (0, 0, real_t(1.2));
   kernel::TemperatureIntegration     temperatureIntegration(dt, 1);
   temperatureIntegration.setInvSpecificHeat(0, real_t(0.05) / dt / ss->shapes[smallSphere]->getInvMass());
   kernel::ThermalExpansion           thermalExpansion(1);
   thermalExpansion.setLinearExpansionCoefficient(0, real_t(0.0005));

   ContactDetection                   acd;
   kernel::DoubleCast                 double_cast;
   mpi::ContactFilter                 contact_filter;
   mpi::ReduceProperty                RP;
   mpi::SyncNextNeighbors             SNN;

   // initial sync
   SNN(*ps, domain);

   for (int64_t i=0; i < simulationSteps; ++i)
   {
      if (i % visSpacing == 0)
      {
         vtkWriter->write();
      }

      if (i > visSpacing * 100)
      {
         auto rad = (real_t(i) - real_t(10000)) * real_t(5e-5) * math::pi;
         plane->setTemperature(real_t(273) + real_t(300) * std::sin(rad));
      }

      ps->forEachParticle(false, kernel::SelectMaster(), accessor, [&](const size_t idx, auto& ac){ac.setForce(idx, Vec3(0,0,real_t(-9.81) ) );}, accessor);

      lc.clear();
      ps->forEachParticle(true, kernel::SelectAll(), accessor, ipilc, accessor, lc);

      lc.forEachParticlePairHalf(true, kernel::SelectAll(), accessor,
                                 [&](const size_t idx1, const size_t idx2, auto& ac)
      {
         if (collision_detection::isInInteractionDistance(idx1, idx2, ac))
         {
            if (double_cast(idx1, idx2, ac, acd, ac ))
            {
               if (contact_filter(acd.getIdx1(), acd.getIdx2(), ac, acd.getContactPoint(), domain))
               {
                  dem(acd.getIdx1(), acd.getIdx2(), ac, acd.getContactPoint(), acd.getContactNormal(), acd.getPenetrationDepth());
                  heatConduction(acd.getIdx1(), acd.getIdx2(), ac);
               }
            }
         }
      },
      accessor );

      RP.operator()<ForceTorqueNotification>(*ps);

      RP.operator()<HeatFluxNotification>(*ps);

      ps->forEachParticle(true, kernel::SelectMaster(), accessor, explicitEuler, accessor);

      ps->forEachParticle(true, kernel::SelectMaster(), accessor, temperatureIntegration, accessor);

      if( i % 10 == 0 )
      {
         ps->forEachParticle(true, kernel::SelectMaster(), accessor, thermalExpansion, accessor);
      }

      SNN(*ps, domain);
   }
   WALBERLA_LOG_INFO_ON_ROOT("*** SIMULATION - END ***");

   return EXIT_SUCCESS;
}

} // namespace mesa_pd
} // namespace walberla

int main( int argc, char* argv[] )
{
   return walberla::mesa_pd::main( argc, argv );
}
