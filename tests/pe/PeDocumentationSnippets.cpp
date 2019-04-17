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
//! \file DocumentationSnippets.cpp
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================


#include "blockforest/all.h"
#include "core/all.h"
#include "domain_decomposition/all.h"

#include "pe/basic.h"
#include "pe/rigidbody/BoxFactory.h"
#include "pe/rigidbody/CapsuleFactory.h"
#include "pe/rigidbody/SphereFactory.h"
#include "pe/rigidbody/UnionFactory.h"
#include "pe/ccd/SimpleCCDDataHandling.h"
#include "pe/synchronization/SyncNextNeighbors.h"
#include "pe/vtk/BodyVtkOutput.h"
#include "pe/vtk/SphereVtkOutput.h"

#include "CheckVitalParameters.h"

#include "core/debug/TestSubsystem.h"
#include "vtk/VTKOutput.h"

#include <tuple>

#include <algorithm>
#include <vector>

namespace walberla {
using namespace walberla::pe;

//! [Definition of Union Types]
using UnionT = Union<Box, Capsule, Sphere>;
using UnionID = UnionT *;
//! [Definition of Union Types]

//! [Definition BodyTypeTuple]
typedef std::tuple<Box, Capsule, Plane, Sphere, UnionT> BodyTypeTuple ;
//! [Definition BodyTypeTuple]

int main( int argc, char ** argv )
{
   walberla::debug::enterTestMode();

   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );

   shared_ptr<BodyStorage> globalBodyStorage = make_shared<BodyStorage>();

   // create blocks
   shared_ptr< StructuredBlockForest > forest = blockforest::createUniformBlockGrid(
            uint_c( 1), uint_c( 1), uint_c( 1), // number of blocks in x,y,z direction
            uint_c( 1), uint_c( 1), uint_c( 1), // how many cells per block (x,y,z)
            real_c(10),                         // dx: length of one cell in physical coordinates
            0,                                  // max blocks per process
            false, false,                       // include metis / force metis
            false, false, false );                 // full periodicity

   //! [Definition Setup TypeIds]
   SetBodyTypeIDs<BodyTypeTuple>::execute();
   //! [Definition Setup TypeIds]

   auto storageID           = forest->addBlockData(createStorageDataHandling<BodyTypeTuple>(), "Storage");
   auto ccdID               = forest->addBlockData(ccd::createHashGridsDataHandling( globalBodyStorage, storageID ), "CCD");
   auto fcdID               = forest->addBlockData(fcd::createGenericFCDDataHandling<BodyTypeTuple, fcd::AnalyticCollideFunctor>(), "FCD");

   //! [Setup HCSITS]
   cr::HCSITS hcsits(globalBodyStorage, forest->getBlockStoragePointer(), storageID, ccdID, fcdID);
   hcsits.setMaxIterations           ( uint_c(10) );
   hcsits.setRelaxationModel         ( cr::HCSITS::InelasticGeneralizedMaximumDissipationContact );
   hcsits.setRelaxationParameter     ( real_t(0.7) );
   hcsits.setErrorReductionParameter ( real_t(0.8) );
   hcsits.setGlobalLinearAcceleration( Vec3(0,0,-1) );
   //! [Setup HCSITS]

   //! [Create a Box]
   // Creating the iron box 1 with a side lengths of 2.5 at the global position (2,3,4).
   // Note that the box is
   // automatically added to the simulation world and is immediately part of the entire
   // simulation. The function returns a handle to the newly created box, which can
   // be used to for instance rotate the box around the global y-axis.
   BoxID box = createBox( *globalBodyStorage, forest->getBlockStorage(), storageID, 1, Vec3(2,3,4), Vec3(2.5,2.5,2.5) );
   if (box != nullptr)
      box->rotate( 0.0, real_c(math::M_PI/3.0), 0.0 );
   //! [Create a Box]

   //! [Create a Capsule]
   // Create a capsule and rotate it after successfull creation.
   CapsuleID capsule = createCapsule( *globalBodyStorage, forest->getBlockStorage(), storageID, 1, Vec3(2,3,4), real_t(1), real_t(1) );
   if (capsule != nullptr)
      capsule->rotate( 0.0, real_c(math::M_PI/3.0), 0.0 );
   //! [Create a Capsule]

   //! [Create a Plane]
   // Create a plane
   PlaneID plane = createPlane( *globalBodyStorage, 1, Vec3(2,3,4), Vec3(2,3,4) );
   //! [Create a Plane]
   WALBERLA_UNUSED( plane );

   //! [Create a Sphere]
   // Create a sphere and rotate it after successfull creation.
   SphereID sphere = createSphere( *globalBodyStorage, forest->getBlockStorage(), storageID, 1, Vec3(2,3,4), real_t(1) );
   if (sphere != nullptr)
      sphere->rotate( 0.0, real_c(math::M_PI/3.0), 0.0 );
   //! [Create a Sphere]

   //! [Create a Union]
   // Create a union and add a box, capsule and sphere.
   UnionID un = createUnion<Box, Capsule, Sphere>( *globalBodyStorage, forest->getBlockStorage(), storageID, 1, Vec3(2,3,4) );
   if (un != nullptr)
   {
      createBox    ( un, 1, Vec3(2,3,4), Vec3(2.5,2.5,2.5) );
      createCapsule( un, 1, Vec3(3,3,4), real_t(1), real_t(1) );
      createSphere ( un, 1, Vec3(4,3,4), real_t(1) );
   }
   //! [Create a Union]

   return EXIT_SUCCESS;
}
} // namespace walberla

int main( int argc, char* argv[] )
{
  return walberla::main( argc, argv );
}