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
//! \file BodyFlags.cpp
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================


#include "blockforest/all.h"
#include "core/all.h"
#include "domain_decomposition/all.h"

#include "pe/basic.h"
#include "pe/communication/ParseMessage.h"
#include "pe/synchronization/SyncNextNeighbors.h"
#include "pe/utility/GetBody.h"

#include "CheckVitalParameters.h"

#include "core/debug/TestSubsystem.h"

#include <tuple>

namespace walberla {
using namespace walberla::pe;

using BodyTuple = std::tuple<Sphere> ;

int main( int argc, char ** argv )
{
   walberla::debug::enterTestMode();

   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );

//   logging::Logging::instance()->setFileLogLevel( logging::Logging::DETAIL );
//   logging::Logging::instance()->includeLoggingToFile("BodyFlagsLog");

   shared_ptr<BodyStorage> globalStorage = make_shared<BodyStorage>();

   // create blocks
   shared_ptr< StructuredBlockForest > forest = blockforest::createUniformBlockGrid(
            uint_c( 2), uint_c( 2), uint_c( 2), // number of blocks in x,y,z direction
            uint_c( 1), uint_c( 1), uint_c( 1), // how many cells per block (x,y,z)
            real_c(10),                         // dx: length of one cell in physical coordinates
            0,                                  // max blocks per process
            false, false,                       // include metis / force metis
            false, false, false );                 // full periodicity

   SetBodyTypeIDs<BodyTuple>::execute();

   auto storageID           = forest->addBlockData(createStorageDataHandling<BodyTuple>(), "Storage");
   auto ccdID               = forest->addBlockData(ccd::createHashGridsDataHandling( globalStorage, storageID ), "CCD");
   auto fcdID               = forest->addBlockData(fcd::createGenericFCDDataHandling<BodyTuple, fcd::AnalyticCollideFunctor>(), "FCD");

   //cr::DEM    cr(globalStorage, forest->getBlockStorage(), storageID, ccdID, fcdID, NULL );
   cr::HCSITS cr(globalStorage, forest->getBlockStoragePointer(), storageID, ccdID, fcdID, nullptr );

   MaterialID iron = Material::find("iron");

   Sphere refGlobalSphere(1, 0, Vec3(9, 9, 9), Vec3(0,0,0), Quat(), 3, iron, true, false, true);
   refGlobalSphere.setLinearVel(Vec3(2,2,2));
   SphereID globalSphere = createSphere( *globalStorage, forest->getBlockStorage(), storageID, 0, Vec3(9,9,9), 3, iron, true, false, true);
   globalSphere->setLinearVel(Vec3(2,2,2));

   Sphere refFixedSphere(2, 0, Vec3(9,9,14), Vec3(0,0,0), Quat(), 3, iron, false, false, true);
   SphereID fixedSphere = createSphere( *globalStorage, forest->getBlockStorage(), storageID, 0, Vec3(9,9,14), 3, iron, false, false, true);
   walberla::id_t fixedSphereID = 0;
   if (fixedSphere != nullptr) fixedSphereID = fixedSphere->getSystemID();
   mpi::allReduceInplace(fixedSphereID, mpi::SUM);

   // synchronize particles
   syncShadowOwners<BodyTuple>( forest->getBlockForest(), storageID, nullptr, real_c(0.0), true);

   cr.setGlobalLinearAcceleration(Vec3(0, 0, real_c(-9.81)));

   checkVitalParameters(&refGlobalSphere, globalSphere);
   for (auto it = forest->begin(); it != forest->end(); ++it)
   {
      blockforest::Block& block = *(dynamic_cast<blockforest::Block*>(&(*it)));
      if (block.getAABB().intersects(refFixedSphere.getAABB()))
      {
         SphereID fixed = static_cast<SphereID> (getBody(*globalStorage, forest->getBlockStorage(), storageID, fixedSphereID));
         WALBERLA_ASSERT_NOT_NULLPTR(fixed);
         checkVitalParameters(&refFixedSphere,fixed);
         if (!block.getAABB().contains(refFixedSphere.getPosition()))
         {
            WALBERLA_ASSERT(fixed->isRemote());
         }
      }
   }

   WALBERLA_LOG_PROGRESS_ON_ROOT("*** SIMULATION - START ***");
   cr.timestep( real_c(1.0) );
   syncShadowOwners<BodyTuple>( forest->getBlockForest(), storageID, nullptr, real_c(0.0), false);
   WALBERLA_LOG_PROGRESS_ON_ROOT("*** SIMULATION - END ***");

   refGlobalSphere.setPosition(Vec3(11,11,11));
   checkVitalParameters(&refGlobalSphere, globalSphere);
   for (auto it = forest->begin(); it != forest->end(); ++it)
   {
      blockforest::Block& block = *(dynamic_cast<blockforest::Block*>(&(*it)));
      if (block.getAABB().intersects(refFixedSphere.getAABB()))
      {
         SphereID fixed = static_cast<SphereID> (getBody(*globalStorage, forest->getBlockStorage(), storageID, fixedSphereID));
         WALBERLA_ASSERT_NOT_NULLPTR(fixed);
         checkVitalParameters(&refFixedSphere, fixed);
         if (!block.getAABB().contains(refFixedSphere.getPosition()))
         {
            WALBERLA_ASSERT(fixed->isRemote());
         }
      }
   }

   return EXIT_SUCCESS;
}
} // namespace walberla

int main( int argc, char* argv[] )
{
  return walberla::main( argc, argv );
}
