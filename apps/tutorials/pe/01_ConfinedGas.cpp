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

//! [Includes]
#include <pe/basic.h>

#include <blockforest/Initialization.h>
#include <core/Environment.h>
#include <core/grid_generator/HCPIterator.h>
#include <core/grid_generator/SCIterator.h>
#include <core/logging/Logging.h>
#include <core/math/Random.h>

#include <tuple>
//! [Includes]

namespace walberla {
using namespace walberla::pe;

//! [BodyTypeTuple]
using BodyTypeTuple = std::tuple<Sphere, Plane> ;
//! [BodyTypeTuple]

int main( int argc, char ** argv )
{
   //! [Parameters]
   Environment env(argc, argv);
   WALBERLA_UNUSED(env);

   math::seedRandomGenerator( static_cast<unsigned int>(1337 * mpi::MPIManager::instance()->worldRank()) );

   real_t spacing          = real_c(1.0);
   real_t radius           = real_c(0.4);
   real_t vMax             = real_c(1.0);
   int    simulationSteps  = 10;
   real_t dt               = real_c(0.01);
   //! [Parameters]

   WALBERLA_LOG_INFO_ON_ROOT("*** GLOBALBODYSTORAGE ***");
   //! [GlobalBodyStorage]
   shared_ptr<BodyStorage> globalBodyStorage = make_shared<BodyStorage>();
   //! [GlobalBodyStorage]

   WALBERLA_LOG_INFO_ON_ROOT("*** BLOCKFOREST ***");
   // create forest
   //! [BlockForest]
   auto forest = blockforest::createBlockForest( AABB(0,0,0,20,20,20), // simulation domain
                                                 Vector3<uint_t>(2,2,2), // blocks in each direction
                                                 Vector3<bool>(false, false, false) // periodicity
                                                 );
   //! [BlockForest]
   if (!forest)
   {
      WALBERLA_LOG_INFO_ON_ROOT( "No BlockForest created ... exiting!");
      return EXIT_SUCCESS;
   }

   WALBERLA_LOG_INFO_ON_ROOT("*** STORAGEDATAHANDLING ***");
   // add block data
   //! [StorageDataHandling]
   auto storageID           = forest->addBlockData(createStorageDataHandling<BodyTypeTuple>(), "Storage");
   //! [StorageDataHandling]
   //! [AdditionalBlockData]
   auto ccdID               = forest->addBlockData(ccd::createHashGridsDataHandling( globalBodyStorage, storageID ), "CCD");
   auto fcdID               = forest->addBlockData(fcd::createGenericFCDDataHandling<BodyTypeTuple, fcd::AnalyticCollideFunctor>(), "FCD");
   //! [AdditionalBlockData]

   WALBERLA_LOG_INFO_ON_ROOT("*** INTEGRATOR ***");
   //! [Integrator]
   cr::HCSITS cr(globalBodyStorage, forest, storageID, ccdID, fcdID);
   cr.setMaxIterations( 10 );
   cr.setRelaxationModel( cr::HardContactSemiImplicitTimesteppingSolvers::ApproximateInelasticCoulombContactByDecoupling );
   cr.setRelaxationParameter( real_t(0.7) );
   cr.setGlobalLinearAcceleration( Vec3(0,0,0) );
   //! [Integrator]

   WALBERLA_LOG_INFO_ON_ROOT("*** BodyTypeTuple ***");
   // initialize body type ids
   //! [SetBodyTypeIDs]
   SetBodyTypeIDs<BodyTypeTuple>::execute();
   //! [SetBodyTypeIDs]

   WALBERLA_LOG_INFO_ON_ROOT("*** SETUP - START ***");
   //! [Material]
   const real_t   static_cof  ( real_c(0.1) / 2 );   // Coefficient of static friction. Note: pe doubles the input coefficient of friction for material-material contacts.
   const real_t   dynamic_cof ( static_cof ); // Coefficient of dynamic friction. Similar to static friction for low speed friction.
   MaterialID     material = createMaterial( "granular", real_t( 1.0 ), 0, static_cof, dynamic_cof, real_t( 0.5 ), 1, 1, 0, 0 );
   //! [Material]

   auto simulationDomain = forest->getDomain();
   const auto& generationDomain = simulationDomain; // simulationDomain.getExtended(-real_c(0.5) * spacing);
   //! [Planes]
   createPlane(*globalBodyStorage, 0, Vec3(1,0,0), simulationDomain.minCorner(), material );
   createPlane(*globalBodyStorage, 0, Vec3(-1,0,0), simulationDomain.maxCorner(), material );
   createPlane(*globalBodyStorage, 0, Vec3(0,1,0), simulationDomain.minCorner(), material );
   createPlane(*globalBodyStorage, 0, Vec3(0,-1,0), simulationDomain.maxCorner(), material );
   createPlane(*globalBodyStorage, 0, Vec3(0,0,1), simulationDomain.minCorner(), material );
   createPlane(*globalBodyStorage, 0, Vec3(0,0,-1), simulationDomain.maxCorner(), material );
   //! [Planes]

   //! [Gas]
   uint_t numParticles = uint_c(0);
   for (auto blkIt = forest->begin(); blkIt != forest->end(); ++blkIt)
   {
      IBlock & currentBlock = *blkIt;
      for (auto it = grid_generator::SCIterator(currentBlock.getAABB().getIntersection(generationDomain), Vector3<real_t>(spacing, spacing, spacing) * real_c(0.5), spacing); it != grid_generator::SCIterator(); ++it)
      {
         SphereID sp = createSphere( *globalBodyStorage, *forest, storageID, 0, *it, radius, material);
         Vec3 rndVel(math::realRandom<real_t>(-vMax, vMax), math::realRandom<real_t>(-vMax, vMax), math::realRandom<real_t>(-vMax, vMax));
         if (sp != nullptr) sp->setLinearVel(rndVel);
         if (sp != nullptr) ++numParticles;
      }
   }
   WALBERLA_LOG_INFO_ON_ROOT("#particles created: " << numParticles);
   syncNextNeighbors<BodyTypeTuple>(*forest, storageID);
   //! [Gas]

   WALBERLA_LOG_INFO_ON_ROOT("*** SETUP - END ***");

   WALBERLA_LOG_INFO_ON_ROOT("*** SIMULATION - START ***");
   //! [GameLoop]
   for (int i=0; i < simulationSteps; ++i)
   {
      if( i % 10 == 0 )
      {
         WALBERLA_LOG_PROGRESS_ON_ROOT( "Timestep " << i << " / " << simulationSteps );
      }

      cr.timestep( real_c(dt) );
      syncNextNeighbors<BodyTypeTuple>(*forest, storageID);
   }
   //! [GameLoop]
   WALBERLA_LOG_INFO_ON_ROOT("*** SIMULATION - END ***");

   WALBERLA_LOG_INFO_ON_ROOT("*** GETTING STATISTICAL INFORMATION ***");
   //! [PostProcessing]
   Vec3 meanVelocity(0,0,0);
   for (auto blockIt = forest->begin(); blockIt != forest->end(); ++blockIt)
   {
      for (auto bodyIt = LocalBodyIterator::begin(*blockIt, storageID); bodyIt != LocalBodyIterator::end(); ++bodyIt)
      {
        meanVelocity += bodyIt->getLinearVel();
      }
   }
   meanVelocity /= numParticles;
   WALBERLA_LOG_INFO( "mean velocity: " << meanVelocity );
   //! [PostProcessing]

   return EXIT_SUCCESS;
}
} // namespace walberla

int main( int argc, char* argv[] )
{
  return walberla::main( argc, argv );
}
