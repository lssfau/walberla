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

#include <core/Environment.h>
#include <core/grid_generator/HCPIterator.h>
#include <core/grid_generator/SCIterator.h>
#include <core/logging/Logging.h>

#include "gui/Gui.h"
#include "timeloop/SweepTimeloop.h"

#include <pe/vtk/SphereVtkOutput.h>
#include "vtk/VTKOutput.h"

#include <pe/raytracing/Raytracer.h>

#include <core/mpi/all.h>

//! [Includes]

using namespace walberla;
using namespace walberla::pe;
using namespace walberla::pe::raytracing;

//! [BodyTypeTuple]
typedef boost::tuple<Sphere, Plane, Box> BodyTypeTuple ;
//! [BodyTypeTuple]

void testRayTracing () {
   shared_ptr<BodyStorage> globalBodyStorage = make_shared<BodyStorage>();
   shared_ptr<BlockForest> forest = createBlockForest(AABB(0,-5,-5,10,5,5), Vec3(1,1,1), Vec3(false, false, false));
   auto storageID = forest->addBlockData(createStorageDataHandling<BodyTypeTuple>(), "Storage");
   SetBodyTypeIDs<BodyTypeTuple>::execute();
   createSphere(*globalBodyStorage, *forest, storageID, 2, Vec3(6,4.5,4.5), real_t(0.5));
   createSphere(*globalBodyStorage, *forest, storageID, 3, Vec3(3.5,-2,0), real_t(1));
   //createSphere(*globalBodyStorage, *forest, storageID, 6, Vec3(3,2,0), real_t(1));
   BoxID box = createBox(*globalBodyStorage, *forest, storageID, 7, Vec3(5,0,0), Vec3(2,4,3));
   box->rotate(0,math::M_PI/4,math::M_PI/4);
   createBox(*globalBodyStorage, *forest, storageID, 7, Vec3(5,-4,3), Vec3(2,2,2));
   //createSphere(*globalBodyStorage, *forest, storageID, 5, Vec3(1,0,0), real_t(0.1));

   //Raytracer raytracer(forest, storageID, uint8_t(640), uint8_t(480), 49.13, Vec3(-5,0,0), Vec3(-1,0,0), Vec3(0,0,1));
   //raytracer.rayTrace<BodyTypeTuple>(0);
}

int main( int argc, char ** argv )
{
   //! [Parameters]
   Environment env(argc, argv);
   WALBERLA_UNUSED(env);
   
   math::seedRandomGenerator( static_cast<unsigned int>(1337 * mpi::MPIManager::instance()->worldRank()) );
   
   //testRayTracing();
   //return 0;
   
   real_t spacing          = real_c(1.0);
   real_t radius           = real_c(0.4);
   real_t vMax             = real_c(1.0);
   size_t    simulationSteps  = 10;
   size_t raytraceSkippedSteps = 10;
   real_t dt               = real_c(0.01);
   
   uint_t blocks_x = 2, blocks_y = 2, blocks_z = 2;
   
   auto cfg = env.config();
   if (cfg != NULL) {
      const Config::BlockHandle confBlock = cfg->getBlock("gas");
      
      spacing = confBlock.getParameter<real_t>("spacing", spacing);
      radius = confBlock.getParameter<real_t>("radius", radius);
      simulationSteps = confBlock.getParameter<size_t>("simulationSteps", simulationSteps);
      raytraceSkippedSteps = confBlock.getParameter<size_t>("raytraceSkippedSteps", raytraceSkippedSteps);
      blocks_x = confBlock.getParameter<uint_t>("blocks_x", blocks_x);
      blocks_y = confBlock.getParameter<uint_t>("blocks_y", blocks_y);
      blocks_z = confBlock.getParameter<uint_t>("blocks_z", blocks_z);
      
      WALBERLA_LOG_INFO_ON_ROOT("spacing: " << spacing);
      WALBERLA_LOG_INFO_ON_ROOT("radius: " << radius);
      WALBERLA_LOG_INFO_ON_ROOT("blocks_x: " << blocks_x);
      WALBERLA_LOG_INFO_ON_ROOT("blocks_y: " << blocks_y);
      WALBERLA_LOG_INFO_ON_ROOT("blocks_z: " << blocks_z);
   }
   
   //! [Parameters]

   WALBERLA_LOG_INFO_ON_ROOT("*** GLOBALBODYSTORAGE ***");
   //! [GlobalBodyStorage]
   shared_ptr<BodyStorage> globalBodyStorage = make_shared<BodyStorage>();
   //! [GlobalBodyStorage]

   WALBERLA_LOG_INFO_ON_ROOT("*** BLOCKFOREST ***");
   // create forest
   //! [BlockForest]
   shared_ptr< BlockForest > forest = createBlockForest( AABB(0,0,0,20,20,20), // simulation domain
                                                         Vector3<uint_t>(blocks_x,blocks_y,blocks_z), // blocks in each direction
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
   
   WALBERLA_LOG_INFO_ON_ROOT("*** RAYTRACER ***");
   if (cfg == NULL) {
      WALBERLA_ABORT("raytracer needs a working config");
   }
   Raytracer raytracer(forest, storageID, globalBodyStorage, cfg->getBlock("raytracing"));

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
   auto generationDomain = simulationDomain; // simulationDomain.getExtended(-real_c(0.5) * spacing);
   //! [Planes]
   PlaneID xPosPlane = createPlane(*globalBodyStorage, 0, Vec3(1,0,0), simulationDomain.minCorner(), material );
   createPlane(*globalBodyStorage, 0, Vec3(-1,0,0), simulationDomain.maxCorner(), material );
   createPlane(*globalBodyStorage, 0, Vec3(0,1,0), simulationDomain.minCorner(), material );
   createPlane(*globalBodyStorage, 0, Vec3(0,-1,0), simulationDomain.maxCorner(), material );
   createPlane(*globalBodyStorage, 0, Vec3(0,0,1), simulationDomain.minCorner(), material );
   createPlane(*globalBodyStorage, 0, Vec3(0,0,-1), simulationDomain.maxCorner(), material );
   //! [Planes]
   
   raytracer.setBodyInvisible(xPosPlane);

   //! [Gas]
   uint_t numParticles = uint_c(0);
   for (auto blkIt = forest->begin(); blkIt != forest->end(); ++blkIt)
   {
      IBlock & currentBlock = *blkIt;
      for (auto it = grid_generator::SCIterator(currentBlock.getAABB().getIntersection(generationDomain), Vector3<real_t>(spacing, spacing, spacing) * real_c(0.5), spacing); it != grid_generator::SCIterator(); ++it)
      {
         SphereID sp = createSphere( *globalBodyStorage, *forest, storageID, 0, *it, radius, material);
         Vec3 rndVel(math::realRandom<real_t>(-vMax, vMax), math::realRandom<real_t>(-vMax, vMax), math::realRandom<real_t>(-vMax, vMax));
         if (sp != NULL) sp->setLinearVel(rndVel);
         if (sp != NULL) ++numParticles;
      }
   }
   WALBERLA_LOG_INFO_ON_ROOT("#particles created: " << numParticles);
   syncNextNeighbors<BodyTypeTuple>(*forest, storageID);
   //! [Gas]

   WALBERLA_LOG_INFO_ON_ROOT("*** SETUP - END ***");

   WALBERLA_LOG_INFO_ON_ROOT("*** SIMULATION - START ***");
   //! [GameLoop]
   for (size_t i=0; i < simulationSteps; ++i)
   {
      if( i % 10 == 0 )
      {
         WALBERLA_LOG_PROGRESS_ON_ROOT( "Timestep " << i << " / " << simulationSteps );
      }

      cr.timestep( real_c(dt) );
      syncNextNeighbors<BodyTypeTuple>(*forest, storageID);
      
      if (i%raytraceSkippedSteps == 0) {
         raytracer.rayTrace<BodyTypeTuple>(i);
      }
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
   
   //WALBERLA_LOG_INFO_ON_ROOT("*** RAYTRACING - START ***");
   //raytracer.rayTrace<BodyTypeTuple>(0);
   //WALBERLA_LOG_INFO_ON_ROOT("*** RAYTRACING - END ***");

   // für einzelne sphere vtks: -> SphereVtkOutput.cpp
   
   auto vtkDomainOutput = vtk::createVTKOutput_DomainDecomposition( forest, "domain_decomposition", 1, "vtk_out", "simulation_step" );
   auto vtkSphereHelper = make_shared<SphereVtkOutput>(storageID, *forest) ;
   auto vtkSphereOutput = vtk::createVTKOutput_PointData(vtkSphereHelper, "Bodies", 1, "vtk_out", "simulation_step", false, false);
   vtkDomainOutput->write();
   vtkSphereOutput->write();
   
   return EXIT_SUCCESS;
}
