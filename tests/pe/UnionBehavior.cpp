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
//! \file   UnionBehavior.cpp
//! \author Tobias Leemann <tobias.leemann@fau.de>
//! \brief Testcase checking whether a single body and a union representing the same geometry behave equally.
//======================================================================================================================

//! [Includes]
#include <pe/basic.h>

#include <core/Environment.h>
#include <core/grid_generator/HCPIterator.h>
#include <core/grid_generator/SCIterator.h>
#include <core/logging/Logging.h>
#include <core/math/Random.h>

#include "blockforest/StructuredBlockForest.h"
#include "blockforest/Initialization.h"
#include "pe/fcd/GenericFCD.h"

#include "pe/rigidbody/BoxFactory.h"
#include "vtk/VTKOutput.h"
#include "pe/vtk/BodyVtkOutput.h"


#include <iostream>
#include <sstream>
#include <fstream>
#include <cstring>
#include <tuple>


//! [Includes]

namespace walberla {
using namespace walberla::pe;


//! [BodyTypeTuple]
using BoxUnion = Union<Box>;
using BodyTypeTuple = std::tuple<Plane, Box, BoxUnion>;
//! [BodyTypeTuple]

// Override the collide function, so that the collisions between the two bodies themselves are not tracked,
// only collisions with the surrounding planes are handled
class PlaneOnlyFCD : public fcd::GenericFCD<BodyTypeTuple, fcd::AnalyticCollideFunctor>{
public:
   Contacts& generateContacts(PossibleContacts& possibleContacts) override
   {
      contacts_.clear();
      fcd::AnalyticCollideFunctor<decltype(contacts_)> func(contacts_);
      for (auto &c : possibleContacts)
      {
         if(c.first->getTypeID() == Plane::getStaticTypeID() || c.second->getTypeID() == Plane::getStaticTypeID()) {
            DoubleCast<BodyTypeTuple, BodyTypeTuple, fcd::AnalyticCollideFunctor<decltype(contacts_)>, bool>::execute(
                     c.first, c.second, func);
         }
      }
      return contacts_;
   }
};


int main( int argc, char ** argv )
{

   Environment env(argc, argv);
   WALBERLA_UNUSED(env);

   // Parameters
   auto dt = real_t(0.005);
   int simulationSteps = 4000;
   auto domainsize = real_t(100);
   // Size of the box
   auto boxlenght = real_t(6.0);
   auto boxwidth = real_t(2.0);

   Vec3 pos(real_t(domainsize/2.0),real_t(domainsize/2.0),real_t(0.8*domainsize));
   Vec3 v0(10,2,3); // Linear velocity of the objects
   Vec3 w0(4,15,1); // Some angular velocity of the objects (quite high)
   int visSpacing = 10;


   // Simulation part
   math::seedRandomGenerator( static_cast<unsigned int>(1337 * mpi::MPIManager::instance()->worldRank()) );



   //! [Parameters]

   WALBERLA_LOG_INFO_ON_ROOT("*** GLOBALBODYSTORAGE ***");
   //! [GlobalBodyStorage]
   shared_ptr<BodyStorage> globalBodyStorage = make_shared<BodyStorage>();
   //! [GlobalBodyStorage]

   WALBERLA_LOG_INFO_ON_ROOT("*** BLOCKFOREST ***");
   // create forest
   //! [BlockForest]
   auto forest = blockforest::createBlockForest( AABB(0,0,0,domainsize,domainsize,domainsize), // simulation domain
                                                 Vector3<uint_t>(3,1,1), // blocks in each direction
                                                 Vector3<bool>(true, false, false) // periodicity
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

   auto sharedFCDBDH        = make_shared<blockforest::AlwaysCreateBlockDataHandling<PlaneOnlyFCD> >( );

   auto fcdID               = forest->addBlockData(sharedFCDBDH, "FCD");
   //! [AdditionalBlockData]

   WALBERLA_LOG_INFO_ON_ROOT("*** INTEGRATOR ***");
   //! [Integrator]
   cr::DEM cr(globalBodyStorage, forest, storageID, ccdID, fcdID);
   cr.setGlobalLinearAcceleration( Vec3(0,0,-6) );
   //! [Integrator]

   WALBERLA_LOG_INFO_ON_ROOT("*** BodyTypeTuple ***");
   // initialize body type ids
   //! [SetBodyTypeIDs]
   SetBodyTypeIDs<BodyTypeTuple>::execute();
   //! [SetBodyTypeIDs]


   // UNCOMMENT THIS BLOCK FOR VTK OUTPUT
   /*WALBERLA_LOG_INFO_ON_ROOT("*** VTK OUTPUT ***");
   //! [VTK Domain Output]
   auto vtkDomainOutput = vtk::createVTKOutput_DomainDecomposition( forest, "domain_decomposition", 1, "VTK", "simulation_step" );
   vtkDomainOutput->write(true);
   auto vtkSphereHelper = make_shared<DefaultBodyVTKOutput>(storageID, *forest) ;
   auto vtkSphereOutput = vtk::createVTKOutput_PointData(vtkSphereHelper, "Bodies", 1, "VTK", "simulation_step", false, false);*/
   //! [VTK Domain Output]

   WALBERLA_LOG_INFO_ON_ROOT("*** SETUP - START ***");
   //! [Material]
   const real_t   static_cof  ( real_c(1.2) / 2 );   // Coefficient of static friction. Note: pe doubles the input coefficient of friction for material-material contacts.
   const real_t   dynamic_cof ( static_cof ); // Coefficient of dynamic friction. Similar to static friction for low speed friction.
   MaterialID     material = createMaterial( "granular", real_t( 1.0 ), 0, static_cof, dynamic_cof, real_t( 0.5 ), 1, real_t(8.11e5), real_t(6.86e1), real_t(6.86e1) );
   //! [Material]

   // Create surrounding planes in y,z directions, but not in x.
   auto simulationDomain = forest->getDomain();
   //! [Planes]
   createPlane(*globalBodyStorage, 100, Vec3(0,1,0), simulationDomain.minCorner(), material );
   createPlane(*globalBodyStorage, 101, Vec3(0,-1,0), simulationDomain.maxCorner(), material );
   createPlane(*globalBodyStorage, 102, Vec3(0,0,1), simulationDomain.minCorner(), material );
   createPlane(*globalBodyStorage, 103, Vec3(0,0,-1), simulationDomain.maxCorner(), material );
   //! [Planes]

   //! [Gas]
   uint_t numParticles = uint_c(0);

   // Create a union of a two half boxes with the center at pos
   BoxUnion* particleU = createUnion<Box>( *globalBodyStorage, *forest, storageID, 0, Vec3(10,10,10));
   // add 2 parts to the union
   createBox(particleU, 4, pos-Vec3(real_t(boxlenght/4.0),0,0), Vec3(real_t(boxlenght/2.0), boxwidth, boxwidth), material);
   createBox(particleU, 5, pos+Vec3(real_t(boxlenght/4.0),0,0), Vec3(real_t(boxlenght/2.0), boxwidth, boxwidth), material);

   if (particleU != nullptr) particleU->setLinearVel(v0);
   if (particleU != nullptr) particleU->setAngularVel(w0);
   if (particleU != nullptr) ++numParticles;

   // Generate the same box as one particle at pos
   BoxID particleS = createBox( *globalBodyStorage, *forest, storageID, 1, pos, Vec3(boxlenght, boxwidth, boxwidth), material );
   if (particleS != nullptr) particleS->setLinearVel(v0);
   if (particleS != nullptr) particleS->setAngularVel(w0);
   if (particleS != nullptr) ++numParticles;

   WALBERLA_LOG_INFO_ON_ROOT("#particles created per process: " << numParticles);
   syncNextNeighbors<BodyTypeTuple>(*forest, storageID);
   //! [Gas]

   WALBERLA_LOG_INFO_ON_ROOT("*** SETUP - END ***");

   WALBERLA_LOG_INFO_ON_ROOT("*** SIMULATION - START ***");
   //vtkDomainOutput->write(true);
   //! [GameLoop]
   for (int i=0; i < simulationSteps; ++i)
   {
      if( i % visSpacing == 0 )
      {
         WALBERLA_LOG_INFO_ON_ROOT( "Timestep " << i << " / " << simulationSteps );
         // UNCOMMENT FOR OUTPUT
         //vtkSphereOutput->write(true);
      }

      cr.timestep( real_c(dt) );

      int bdcount = 0;
      RigidBody* refU = nullptr;
      RigidBody* refS = nullptr;
      for (auto blockIt = forest->begin(); blockIt != forest->end(); ++blockIt)
      {
         for (auto bodyIt = LocalBodyIterator::begin(*blockIt, storageID); bodyIt != LocalBodyIterator::end(); ++bodyIt)
         {
            if(bodyIt->getID() == 0){
               refU = &*bodyIt;
            }
            if(bodyIt->getID() == 1){
               refS = &*bodyIt;
            }
            bdcount ++;
         }
      }

      // Compare with an quite large epsilon, as round of errors in collisions will amplify themselves over time
      auto eps = real_t(1e-4);
      WALBERLA_CHECK_GREATER_EQUAL(bdcount, 2, "Lost a body.");
      WALBERLA_CHECK_NOT_NULLPTR(refS, "Single body was not found.");
      WALBERLA_CHECK_NOT_NULLPTR(refU, "Union was not found.");
      // norm(pos1-pos2) < 1e-4 -> norm*norm < 1e-8
      WALBERLA_CHECK_LESS((refU->getPosition()-refS->getPosition()).sqrLength(), eps*eps  );
      WALBERLA_CHECK_LESS((refU->getLinearVel()-refS->getLinearVel()).sqrLength(), eps*eps );
      WALBERLA_CHECK_LESS((refU->getAngularVel()-refS->getAngularVel()).sqrLength(), eps*eps);

      // Sums of matrix differences...
      WALBERLA_CHECK_LESS(std::abs(Vec3(1,1,1) * ((refU->getRotation()-refS->getRotation())*Vec3(1,1,1))), real_t(9.0*eps));
      WALBERLA_CHECK_LESS(std::abs(Vec3(1,1,1) * ((refU->getInertia()-refS->getInertia())*Vec3(1,1,1))), real_t(9.0*eps));
      //WALBERLA_LOG_INFO_ON_ROOT("Timestep " << i << ": " <<refU->getPosition() << "/" << refS->getPosition() << (Vec3(1,1,1) * ((refU->getInertia()-refS->getInertia())*Vec3(1,1,1)))<< " " << abs(Vec3(1,1,1) * ((refU->getRotation()-refS->getRotation())*Vec3(1,1,1))));


      syncNextNeighbors<BodyTypeTuple>(*forest, storageID);
   }
   //! [GameLoop]
   WALBERLA_LOG_INFO_ON_ROOT("*** SIMULATION - END ***");

   return EXIT_SUCCESS;
}
} // namespace walberla

int main( int argc, char* argv[] )
{
   return walberla::main( argc, argv );
}
