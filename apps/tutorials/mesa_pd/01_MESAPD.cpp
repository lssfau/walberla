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
//! \file   01_ConfinedGas.cpp
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

//! [walberlaIncludes]
#include <blockforest/Initialization.h>

#include <core/Environment.h>
#include <core/grid_generator/HCPIterator.h>
#include <core/grid_generator/SCIterator.h>
#include <core/logging/Logging.h>
#include <core/math/Random.h>
//! [walberlaIncludes]

//! [mesapdIncludes]
#include <mesa_pd/collision_detection/AnalyticContactDetection.h>
#include <mesa_pd/data/LinkedCells.h>
#include <mesa_pd/data/ParticleAccessorWithShape.h>
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
#include <mesa_pd/vtk/ParticleVtkOutput.h>
//! [mesapdIncludes]

namespace walberla {
namespace mesa_pd {

//! [CreationHelper]
data::ParticleStorage::iterator createWall( data::ParticleStorage& ps,
                                            data::ShapeStorage& ss,
                                            const Vec3& pos,
                                            const Vec3& normal )
{
   using namespace walberla::mesa_pd::data::particle_flags;

   auto pl = ps.create(true);
   pl->setPosition          ( pos );
   pl->setInteractionRadius ( std::numeric_limits<real_t>::infinity() );
   pl->setShapeID           ( ss.create<data::HalfSpace>( normal ) );
   pl->setOwner             ( walberla::mpi::MPIManager::instance()->rank() );
   pl->setType              ( 1 );
   set(pl->getFlagsRef(), INFINITE);
   set(pl->getFlagsRef(), FIXED);
   set(pl->getFlagsRef(), NON_COMMUNICATING);
   return pl;
}

data::ParticleStorage::iterator createSphere( data::ParticleStorage& ps,
                                              const Vec3& pos,
                                              const real_t& radius,
                                              const size_t shapeID)
{
   auto sp = ps.create();
   sp->setPosition          ( pos );
   sp->setInteractionRadius ( radius );
   sp->setShapeID           ( shapeID );
   sp->setOwner             ( walberla::MPIManager::instance()->rank() );
   sp->setType              ( 0 );
   return sp;
}
//! [CreationHelper]

int main( int argc, char ** argv )
{
   using namespace walberla::mesa_pd;

   //! [Parameters]
   Environment env(argc, argv);
   WALBERLA_UNUSED(env);

   math::seedRandomGenerator( static_cast<unsigned int>(1337 * walberla::mpi::MPIManager::instance()->worldRank()) );

   real_t spacing          = real_c(1.0);
   real_t radius           = real_c(0.4);
   real_t density          = real_t(2707);
   real_t vMax             = real_c(1.0);
   int    simulationSteps  = 10;
   real_t dt               = real_c(0.01);
   //! [Parameters]

   WALBERLA_LOG_INFO_ON_ROOT("*** BLOCKFOREST ***");
   // create forest
   //! [BlockForest]
   auto forest = blockforest::createBlockForest( AABB(0,0,0,40,40,40), // simulation domain
                                                 Vector3<uint_t>(2,2,2), // blocks in each direction
                                                 Vector3<bool>(false, false, false) // periodicity
                                                 );
   if (!forest)
   {
      WALBERLA_LOG_INFO_ON_ROOT( "No BlockForest created ... exiting!");
      return EXIT_SUCCESS;
   }
   //! [BlockForest]

   auto simulationDomain = forest->getDomain();
   auto localDomain = forest->begin()->getAABB();
   for (auto& blk : *forest)
   {
      localDomain.merge(blk.getAABB());
   }

   auto domain = domain::BlockForestDomain(forest);

   WALBERLA_LOG_INFO_ON_ROOT("*** DATA STRUCTURES ***");
   //! [DataStructures]
   auto ps = std::make_shared<data::ParticleStorage>(100);
   auto ss = std::make_shared<data::ShapeStorage>();
   data::ParticleAccessorWithShape accessor(ps, ss);
   data::LinkedCells     lc(localDomain.getExtended(spacing), spacing );
   //! [DataStructures]

   const auto& generationDomain = simulationDomain; // simulationDomain.getExtended(-real_c(0.5) * spacing);
   //! [Walls]
   createWall(*ps, *ss, simulationDomain.minCorner(), Vec3(+1,0,0) );
   createWall(*ps, *ss, simulationDomain.maxCorner(), Vec3(-1,0,0) );
   createWall(*ps, *ss, simulationDomain.minCorner(), Vec3(0,+1,0) );
   createWall(*ps, *ss, simulationDomain.maxCorner(), Vec3(0,-1,0) );
   createWall(*ps, *ss, simulationDomain.minCorner(), Vec3(0,0,+1) );
   createWall(*ps, *ss, simulationDomain.maxCorner(), Vec3(0,0,-1) );
   //! [Walls]

   //! [Spheres]
   auto  sphereShape = ss->create<data::Sphere>( radius );
   ss->shapes[sphereShape]->updateMassAndInertia( density );
   uint_t numParticles = uint_c(0);
   for (auto blkIt = forest->begin(); blkIt != forest->end(); ++blkIt)
   {
      IBlock & currentBlock = *blkIt;
      for (auto it = grid_generator::SCIterator(currentBlock.getAABB().getIntersection(generationDomain),
                                                Vector3<real_t>(spacing, spacing, spacing) * real_c(0.5),
                                                spacing);
           it != grid_generator::SCIterator();
           ++it)
      {
         auto sp = createSphere( *ps, *it, radius, sphereShape);
         Vec3 rndVel(math::realRandom<real_t>(-vMax, vMax),
                     math::realRandom<real_t>(-vMax, vMax),
                     math::realRandom<real_t>(-vMax, vMax));
         sp->setLinearVelocity(rndVel);
         ++numParticles;
      }
   }
   walberla::mpi::reduceInplace(numParticles, walberla::mpi::SUM);
   WALBERLA_LOG_INFO_ON_ROOT("#particles created: " << numParticles);
   //! [Spheres]

   WALBERLA_LOG_INFO_ON_ROOT("*** SETUP - END ***");

   WALBERLA_LOG_INFO_ON_ROOT("*** KERNELS ***");
   //! [Kernels]
   kernel::ExplicitEulerWithShape        explicitEulerWithShape( dt );
   kernel::InsertParticleIntoLinkedCells ipilc;
   kernel::SpringDashpot                 dem( 2 );
   auto meffSphereSphere = real_t(0.5) * (real_c(4.0)/real_c(3.0) * math::pi) * radius * radius * radius * density;
   auto meffSphereWall   = real_t(1.0) * (real_c(4.0)/real_c(3.0) * math::pi) * radius * radius * radius * density;
   dem.setParametersFromCOR( 0, 0, real_t(0.9), real_t(20) * dt, meffSphereSphere );
   dem.setParametersFromCOR( 0, 1, real_t(0.9), real_t(20) * dt, meffSphereWall );
   collision_detection::AnalyticContactDetection              acd;
   kernel::DoubleCast                    double_cast;
   mpi::ContactFilter                    contact_filter;
   mpi::ReduceProperty                   RP;
   mpi::SyncNextNeighbors                SNN;
   //! [Kernels]

   SNN(*ps, domain);

   WALBERLA_LOG_INFO_ON_ROOT("*** SIMULATION - START ***");
   //! [Loop]
   for (int i=0; i < simulationSteps; ++i)
   {
      if( i % 10 == 0 )
      {
         WALBERLA_LOG_PROGRESS_ON_ROOT( "Timestep " << i << " / " << simulationSteps );
      }

      lc.clear();
      ps->forEachParticle(true, kernel::SelectAll(), accessor, ipilc, accessor, lc);

      lc.forEachParticlePairHalf(true,
                                 kernel::SelectAll(),
                                 accessor,
                                 [&](const size_t idx1, const size_t idx2, auto& ac)
      {
         if (double_cast(idx1, idx2, ac, acd, ac ))
         {
            if (contact_filter(acd.getIdx1(), acd.getIdx2(), ac, acd.getContactPoint(), domain))
            {
               dem(acd.getIdx1(), acd.getIdx2(), ac, acd.getContactPoint(), acd.getContactNormal(), acd.getPenetrationDepth());
            }
         }
      },
      accessor );

      RP.operator()<ForceTorqueNotification>(*ps);

      ps->forEachParticle(true, kernel::SelectLocal(), accessor, explicitEulerWithShape, accessor);

      SNN(*ps, domain);
   }
   //! [Loop]
   WALBERLA_LOG_INFO_ON_ROOT("*** SIMULATION - END ***");

   WALBERLA_LOG_INFO_ON_ROOT("*** GETTING STATISTICAL INFORMATION ***");
   //! [PostProcessing]
   auto meanVelocity = real_t(0);

   ps->forEachParticle(true,
                       kernel::SelectLocal(),
                       accessor,
                       [&meanVelocity](const size_t idx, auto& ac)
   {
      meanVelocity += length(ac.getLinearVelocity(idx));
   },
   accessor);

   walberla::mpi::reduceInplace(meanVelocity, walberla::mpi::SUM);
   meanVelocity /= real_c(numParticles);
   WALBERLA_LOG_INFO_ON_ROOT( "mean velocity: " << meanVelocity );
   //! [PostProcessing]

   return EXIT_SUCCESS;
}

} // namespace mesa_pd
} // namespace walberla

int main( int argc, char* argv[] )
{
  return walberla::mesa_pd::main( argc, argv );
}
