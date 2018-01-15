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

#include <pe/raytracing/Intersects.h>
#include <pe/raytracing/Ray.h>

//! [Includes]

using namespace walberla;
using namespace walberla::pe;
using namespace walberla::pe::raytracing;

//! [BodyTypeTuple]
typedef boost::tuple<Sphere, Plane, Box> BodyTypeTuple ;
//! [BodyTypeTuple]


real_t deg2rad(real_t deg) {
   return deg * math::M_PI / real_t(180.0);
}

void writeTBufferToFile(real_t buffer[], walberla::id_t idBuffer[], size_t width, size_t height, std::string fileName) {
   real_t t_max = 1;
   real_t t_min = INFINITY;
   walberla::id_t maxId = 0;
   for (size_t i = 0; i < width*height; i++) {
      if (buffer[i] > t_max && buffer[i] != INFINITY) {
         t_max = buffer[i];
      }
      if (buffer[i] < t_min) {
         //t_min = buffer[i];
      }
      if (idBuffer[i] > maxId) {
         maxId = idBuffer[i];
      }
   }
   if (t_min == INFINITY) t_min = 0;
   
   std::map<walberla::id_t, char> idToColors;
   
   std::ofstream ofs(fileName, std::ios::out | std::ios::binary);
   ofs << "P6\n" << width << " " << height << "\n255\n";
   for (size_t i = height*width-1; i > 0; i--) {
      char r = 0, g = 0, b = 0;
      real_t t = buffer[i];
      walberla::id_t id = idBuffer[i];
      if (t == INFINITY) {
         r = g = b = (char)255;
      } else {
         //r = g = b = (char)(255 * (t-t_min)/t_max);
         //b = (char)(255 * (real_t(id)/real_t(maxId)));
         //if (b > (char)250) b = (char)250;
         if (idToColors.count(id) == 0) {
            idToColors.insert(std::make_pair(id, math::intRandom(0, 240)));
         }
         b = (char)(idToColors.at(id) * ((t-t_min)/t_max));
         //b = (char)(200 * ((t-t_min)/t_max));
         r = (char)(200 * ((t-t_min)/t_max));
      }
      ofs << r << g << b;
   }
   
   ofs.close();
}

void rayTrace (shared_ptr<BlockForest> forest, BlockDataID storageID) {
   // - settings
   size_t pixelsHorizontal = 640;
   size_t pixelsVertical = 480;
   real_t fov_vertical = 49.13; // in degrees, in vertical direction
   // camera settings
   Vec3 cameraPosition(-5,0,0);
   Vec3 lookAtPoint(1,0,0);
   Vec3 upVector(0,1,0);
   
   //real_t scale = tan(deg2rad(fov * 0.5));
   
   // - viewing plane construction
   // eye cos setup
   Vec3 n = (cameraPosition - lookAtPoint).getNormalized(); // normal vector of viewing plane
   Vec3 u = (upVector % n).getNormalized(); // u and ...
   Vec3 v = n % u; // ... v span the viewing plane
   // image plane setup
   real_t d = (cameraPosition - lookAtPoint).length(); // distance from camera to viewing plane
   real_t aspectRatio = real_t(pixelsHorizontal) / real_t(pixelsVertical);
   real_t imagePlaneHeight = tan(deg2rad(fov_vertical)/real_t(2.)) * real_t(2.) * d;
   real_t imagePlaneWidth = imagePlaneHeight * aspectRatio;
   Vec3 imagePlaneOrigin = lookAtPoint - u*imagePlaneWidth/real_t(2.) - v*imagePlaneHeight/real_t(2.);
   real_t pixelWidth = imagePlaneWidth / pixelsHorizontal;
   real_t pixelHeight = imagePlaneHeight / pixelsVertical;
   
   // - raytracing
   real_t frameBuffer[pixelsHorizontal*pixelsVertical];
   walberla::id_t idBuffer[pixelsHorizontal*pixelsVertical];
   std::map<RigidBody*, bool> visibleBodies;
   
   real_t t, t_closest;
   walberla::id_t id_closest;
   for (size_t x = 0; x < pixelsHorizontal; x++) {
      for (size_t y = 0; y < pixelsVertical; y++) {
         //WALBERLA_LOG_INFO(x << "/" << y);
         Vec3 pixelLocation = imagePlaneOrigin + u*(x+real_t(0.5))*pixelWidth + v*(y+real_t(0.5))*pixelHeight;
         Vec3 direction = (pixelLocation - cameraPosition).getNormalized();
         Ray ray(cameraPosition, direction);
         
         t_closest = INFINITY;
         id_closest = 0;
         RigidBody* body_closest;
         for (auto blockIt = forest->begin(); blockIt != forest->end(); ++blockIt) {
            for (auto bodyIt = LocalBodyIterator::begin(*blockIt, storageID); bodyIt != LocalBodyIterator::end(); ++bodyIt) {
               IntersectsFunctor func(&ray, &t);
               bool intersects = SingleCast<BodyTypeTuple, IntersectsFunctor, bool>::execute(*bodyIt, func);
               
               if (intersects && t < t_closest) {
                  // body is shot by ray and currently closest to camera
                  //WALBERLA_LOG_INFO("object is closer than currently closest one");
                  
                  t_closest = t;
                  id_closest = bodyIt->getID();
                  body_closest = *bodyIt;
               }
            }
         }
         //std::cout << (t_closest != INFINITY ? size_t(t_closest) : 0) << " ";
         if (visibleBodies.find(body_closest) != visibleBodies.end()) {
            visibleBodies.insert(std::make_pair(body_closest, true));
         }
         frameBuffer[y*pixelsHorizontal + x] = t_closest;
         idBuffer[y*pixelsHorizontal + x] = id_closest;
      }
      //std::cout << std::endl;
   }
   
   writeTBufferToFile(frameBuffer, idBuffer, pixelsHorizontal, pixelsVertical, "/Users/ng/Desktop/tbuffer.ppm");
}

void testRayTracing () {
   shared_ptr<BodyStorage> globalBodyStorage = make_shared<BodyStorage>();
   shared_ptr<BlockForest> forest = createBlockForest(AABB(0,-5,-5,10,5,5), Vec3(1,1,1), Vec3(false, false, false));
   auto storageID = forest->addBlockData(createStorageDataHandling<BodyTypeTuple>(), "Storage");
   SetBodyTypeIDs<BodyTypeTuple>::execute();
   createSphere(*globalBodyStorage, *forest, storageID, 2, Vec3(6,0,0), real_t(3));
   createSphere(*globalBodyStorage, *forest, storageID, 3, Vec3(3,0,-3), real_t(1));
   createSphere(*globalBodyStorage, *forest, storageID, 6, Vec3(3,2,0), real_t(1));
   createBox(*globalBodyStorage, *forest, storageID, 7, Vec3(5,4,2), Vec3(2,3,4));
   //createSphere(*globalBodyStorage, *forest, storageID, 5, Vec3(1,0,0), real_t(0.1));

   rayTrace(forest, storageID);
}

int main( int argc, char ** argv )
{
   //! [Parameters]
   Environment env(argc, argv);
   WALBERLA_UNUSED(env);
   
   testRayTracing();
   return 0;
   
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
   shared_ptr< BlockForest > forest = createBlockForest( AABB(0,0,0,20,20,20), // simulation domain
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
   auto generationDomain = simulationDomain; // simulationDomain.getExtended(-real_c(0.5) * spacing);
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
   
   rayTrace(forest, storageID);
   
   // für einzelne sphere vtks: -> SphereVtkOutput.cpp
   
   auto vtkDomainOutput = vtk::createVTKOutput_DomainDecomposition( forest, "domain_decomposition", 1, "vtk_out", "simulation_step" );
   auto vtkSphereHelper = make_shared<SphereVtkOutput>(storageID, *forest) ;
   auto vtkSphereOutput = vtk::createVTKOutput_PointData(vtkSphereHelper, "Bodies", 1, "vtk_out", "simulation_step", false, false);
   vtkDomainOutput->write();
   vtkSphereOutput->write();
   
   return EXIT_SUCCESS;
}
