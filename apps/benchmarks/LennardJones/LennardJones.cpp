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
//! \file   LennardJones.cpp
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include <pe/basic.h>
#include <pe/vtk/SphereVtkOutput.h>

#include <blockforest/Initialization.h>
#include <core/Abort.h>
#include <core/Environment.h>
#include <core/grid_generator/SCIterator.h>
#include <core/logging/Logging.h>
#include <core/math/Random.h>
#include <core/waLBerlaBuildInfo.h>
#include <vtk/VTKOutput.h>

#include <functional>
#include <memory>
#include <tuple>

namespace walberla {
using namespace walberla::pe;
using namespace walberla::timing;

using BodyTuple = std::tuple<Sphere> ;

inline
void LJ(RigidBody& bd1, RigidBody& bd2)
{
   if (bd1.getSystemID() != bd2.getSystemID())
   {
      Vec3 dir = bd1.getPosition() - bd2.getPosition();
      const real_t inv_d = real_t(1) / length(dir);
      dir *= inv_d;
      WALBERLA_ASSERT_FLOAT_EQUAL(length(dir),
                                  real_t(1),
                                  "direction not normalized: " << dir);
      const real_t sr = real_t(1.21) * inv_d;
      const real_t sr6 = (sr*sr*sr)*(sr*sr*sr);
      const real_t force = real_t(4) * real_t(2.12) * sr6 * ( sr6 - 1 );
      bd1.addForce(force * dir);
   }
}

inline
void integrate(RigidBody& bd)
{
      bd.setLinearVel( bd.getForce() * real_t(0.01) + bd.getLinearVel() );
      bd.setPosition( bd.getLinearVel() * real_t(0.01) + bd.getPosition() );
      bd.setForce( Vec3( real_t(0) ) );
}

int main( int argc, char ** argv )
{
   Environment env(argc, argv);

   logging::Logging::instance()->setStreamLogLevel(logging::Logging::INFO);
   logging::Logging::instance()->setFileLogLevel(logging::Logging::INFO);

   WALBERLA_LOG_INFO_ON_ROOT( "config file: " << argv[1] )
         WALBERLA_LOG_INFO_ON_ROOT( "waLBerla Revision: " << WALBERLA_GIT_SHA1 );

   math::seedRandomGenerator( static_cast<unsigned int>(1337 * mpi::MPIManager::instance()->worldRank()) );

   WALBERLA_LOG_INFO_ON_ROOT("*** GLOBALBODYSTORAGE ***");
   shared_ptr<BodyStorage> globalBodyStorage = make_shared<BodyStorage>();

   WALBERLA_LOG_INFO_ON_ROOT("*** BLOCKFOREST ***");
   //domain setup
   const real_t spacing = real_t(1.0);
   math::AABB domain( Vec3(real_t(-0.5), real_t(-0.5), real_t(-0.5)),
                      Vec3(real_t(+0.5), real_t(+0.5), real_t(+0.5)));
   domain.scale(real_t(20));

   // create forest
   auto forest = blockforest::createBlockForest(domain,
                                                Vector3<uint_t>(1,1,1),
                                                Vector3<bool>(false, false, false));
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
   WALBERLA_UNUSED(ccdID);
   WALBERLA_UNUSED(fcdID);

   WALBERLA_LOG_INFO_ON_ROOT("*** SYNCCALL ***");
   std::function<void(void)> syncCallWithoutTT;
   syncCallWithoutTT = std::bind( pe::syncNextNeighbors<BodyTuple>, std::ref(*forest), storageID, static_cast<WcTimingTree*>(nullptr), real_c(0.1), false );
   WALBERLA_LOG_INFO_ON_ROOT("Using NextNeighbor sync!");

   WALBERLA_LOG_INFO_ON_ROOT("*** VTK ***");
   auto vtkDomainOutput = vtk::createVTKOutput_DomainDecomposition( forest, "domain_decomposition", 1, "vtk_out", "simulation_step" );
   auto vtkSphereHelper = make_shared<SphereVtkOutput>(storageID, *forest) ;
   auto vtkSphereOutput = vtk::createVTKOutput_PointData(vtkSphereHelper, "Bodies", 1, "vtk_out", "simulation_step", false, false);

   WALBERLA_LOG_INFO_ON_ROOT("*** SETUP - START ***");
   auto simulationDomain = forest->getDomain();
   const auto& generationDomain = simulationDomain; // simulationDomain.getExtended(-real_c(0.5) * spacing);

   //initialize particles
   uint_t numParticles(0);
   for (auto it = grid_generator::SCIterator(generationDomain, Vec3(spacing, spacing, spacing) * real_c(0.5), spacing);
        it != grid_generator::SCIterator();
        ++it)
   {
      SphereID sp = createSphere(*globalBodyStorage, *forest, storageID, 0, *it, real_t(1));

      if (sp!=nullptr)
      {
         sp->setLinearVel( math::realRandom(real_t(-1), real_t(+1)),
                           math::realRandom(real_t(-1), real_t(+1)),
                           math::realRandom(real_t(-1), real_t(+1)));
         ++numParticles;
      }
   }
   mpi::reduceInplace(numParticles, mpi::SUM);
   WALBERLA_LOG_INFO_ON_ROOT("#particles created: " << numParticles);

   WALBERLA_LOG_INFO_ON_ROOT("*** SETUP - END ***");

   // synchronize particles
   //syncCallWithoutTT();
   //syncCallWithoutTT();

   WALBERLA_LOG_INFO_ON_ROOT("*** SIMULATION - START ***");
   WALBERLA_MPI_BARRIER();
   WcTimer timer;
   for (int i=0; i < 1000; ++i)
   {
      WALBERLA_LOG_DEVEL(i);

      for (auto& blk : *forest)
      {
         Storage* storage = blk.getData< Storage >( storageID );
         BodyStorage& localStorage = (*storage)[0];

         for (auto& bd1 : localStorage)
         {
            for (auto& bd2 : localStorage)
            {
               LJ(bd1, bd2);
            }
         }
         const real_t coeff = real_t(0.2);
         for (auto& bd : localStorage)
         {
            bd.addForce( -coeff*bd.getPosition() );
         }
      }

      for (auto& blk : *forest)
      {
         Storage* storage = blk.getData< Storage >( storageID );
         BodyStorage& localStorage = (*storage)[0];

         for (auto& bd : localStorage)
         {
            integrate(bd);
         }
      }

      vtkSphereOutput->write( );
      //syncCallWithoutTT();
   }
   WALBERLA_MPI_BARRIER();
   timer.end();
   WALBERLA_LOG_INFO_ON_ROOT("runtime: " << timer.average());
   WALBERLA_LOG_INFO_ON_ROOT("*** SIMULATION - END ***");
   return EXIT_SUCCESS;
}
} // namespace walberla

int main( int argc, char* argv[] )
{
   return walberla::main( argc, argv );
}
