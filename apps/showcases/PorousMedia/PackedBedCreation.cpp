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
//! \file PackedBedCreation.cpp
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//
//======================================================================================================================

//======================================================================================================================
// This showcase creates a packed bed consisting of spherical particles by simulating a sedimentation process using
// rigid body dynamics. It must be given a parameter file as command-line argument, such as "PackedBedCreation.prm". The
// resulting porous medium is periodic in two directions and the location of the spherical particles are written to a
// file in which the particle positions are denoted relative to the particles' diameter.
// In the simulation, the spherical particles fall into a "box" due to a gravitational acceleration. The box is periodic
// in the side-directions (x- and y-directions) but bound by a solid plane at the bottom (in z-direction). The spheres
// are initially generated on top of the domain in a structured order with slight random disturbances. Once the
// particles have settled (according to their velocity), all particles with center point closer than three particle
// diameters to the bottom plane are removed. This is intended to reduce the influence of the bottom plane on the
// geometry creation of the packed bed. Finally, the packed bed is shifted upwards such that it is placed in the center
// of the domain.
//======================================================================================================================


#include "blockforest/all.h"

#include "core/Environment.h"
#include "core/grid_generator/all.h"
#include "core/mpi/MPITextFile.h"
#include "core/timing/RemainingTimeLogger.h"

#include "pe/basic.h"
#include "pe/synchronization/SyncShadowOwners.h"
#include "pe/utility/DestroyBody.h"
#include "pe/vtk/SphereVtkOutput.h"

#include "vtk/VTKOutput.h"

#include <string>
#include <tuple>

namespace walberla
{
///////////
// USING //
///////////

using BodyTypeTuple = std::tuple< pe::Plane, pe::Sphere >; // must contain all body types used in this simulation

///////////////
// UTILITIES //
///////////////

// write the particles' locations (non-dimensionalized with particle diameter) to file
void writePorousGeometry(BlockForest& blockForest, BlockDataID bodyStorageID, uint_t particleDiameter,
                         const std::string& filename);

// get the absolute maximum velocity of all particles
real_t getParticleMaxVel(BlockForest& blockForest, const BlockDataID& bodyStorageID);

////////////////
// PARAMETERS //
////////////////

// data structure class for storing basic parameters
struct Setup
{
   Vector3< uint_t > numBlocks;  // number of blocks (= processes) in x,y,z, direction
   Vector3< uint_t > domainSize; // domain size in x,y,z direction

   uint_t particleDiameter; // cells per diameter of the particles

   uint_t numInitialParticles; // number of initially created particles (some spheres at the bottom will be removed)
   uint_t maxTimesteps;        // maximum number of PE timesteps (if convergence criterion is not fulfilled)
   real_t convVelocity; // velocity at which the particles are considered stationary and the simulation is converged

   uint_t maxIterations;                       // number of maximum iterations of the physics engine (PE) contact solver
   real_t relaxationParameter;                 // relaxation parameter of the PE contact solver
   Vector3< real_t > globalLinearAcceleration; // global linear acceleration (= force) of the PE contact solver
   real_t dt;                                  // temporal resolution of the PE contact solver

   std::string porousGeometryFile; // filename of the output file that defines the packed bed
   uint_t vtkInterval;             // frequency at which vtk output is written

   void printSetup()
   {
      WALBERLA_LOG_DEVEL_VAR_ON_ROOT(numBlocks);
      WALBERLA_LOG_DEVEL_VAR_ON_ROOT(particleDiameter);
      WALBERLA_LOG_DEVEL_VAR_ON_ROOT(domainSize);
      WALBERLA_LOG_DEVEL_VAR_ON_ROOT(numInitialParticles);
      WALBERLA_LOG_DEVEL_VAR_ON_ROOT(maxTimesteps);
      WALBERLA_LOG_DEVEL_VAR_ON_ROOT(convVelocity);
      WALBERLA_LOG_DEVEL_VAR_ON_ROOT(maxIterations);
      WALBERLA_LOG_DEVEL_VAR_ON_ROOT(relaxationParameter);
      WALBERLA_LOG_DEVEL_VAR_ON_ROOT(globalLinearAcceleration);
      WALBERLA_LOG_DEVEL_VAR_ON_ROOT(dt);
      WALBERLA_LOG_DEVEL_VAR_ON_ROOT(porousGeometryFile);
      WALBERLA_LOG_DEVEL_VAR_ON_ROOT(vtkInterval);
   }
}; // struct Setup

//////////
// MAIN //
//////////

int main(int argc, char** argv)
{
   Environment walberlaEnv(argc, argv);

   if (argc < 2) { WALBERLA_ABORT("Please specify a parameter file as input argument."); }

   Setup setup;

   ///////////////////////////
   // SIMULATION PROPERTIES //
   ///////////////////////////

   // get parameters from config file
   auto domainParameters  = walberlaEnv.config()->getOneBlock("DomainParameters");
   setup.numBlocks        = domainParameters.getParameter< Vector3< uint_t > >("numBlocks");
   setup.particleDiameter = domainParameters.getParameter< uint_t >("particleDiameter");
   setup.domainSize       = domainParameters.getParameter< Vector3< uint_t > >("relDomainSize");
   setup.domainSize *= setup.particleDiameter;

   auto simulationParameters = walberlaEnv.config()->getOneBlock("SimulationParameters");
   setup.numInitialParticles = simulationParameters.getParameter< uint_t >("numInitialParticles");
   setup.maxTimesteps        = simulationParameters.getParameter< uint_t >("maxTimesteps");
   setup.convVelocity        = simulationParameters.getParameter< real_t >("convVelocity");

   auto contactSolverParameters = walberlaEnv.config()->getOneBlock("ContactSolverParameters");
   setup.maxIterations          = contactSolverParameters.getParameter< uint_t >("maxIterations");
   setup.relaxationParameter    = contactSolverParameters.getParameter< real_t >("relaxationParameter");
   setup.globalLinearAcceleration =
      contactSolverParameters.getParameter< Vector3< real_t > >("globalLinearAcceleration");
   setup.dt = contactSolverParameters.getParameter< real_t >("dt");

   auto outputParameters    = walberlaEnv.config()->getOneBlock("OutputParameters");
   setup.porousGeometryFile = outputParameters.getParameter< std::string >("porousGeometryFile");
   setup.vtkInterval        = outputParameters.getParameter< uint_t >("vtkInterval");

   if (setup.numBlocks[0] < uint_c(3) || setup.numBlocks[1] < uint_c(3))
   {
      WALBERLA_ABORT("The number of blocks in x- and y-direction must not be smaller than 3. This is a strict "
                     "requirement of the PE when using periodic boundaries in these directions.");
   }

   setup.printSetup();

   //////////////////////
   // SIMULATION SETUP //
   //////////////////////

   // axis-aligned bounding box (AABB) that represents the final domain
   const auto domainAABB = AABB(real_c(0), real_c(0), real_c(0), real_c(setup.domainSize[0]),
                                real_c(setup.domainSize[1]), real_c(2) * real_c(setup.domainSize[2]));

   // AABB in which the particles are initially positioned before falling down (above domainAABB)
   const auto genAABB =
      AABB(real_c(setup.domainSize[0]), real_c(setup.domainSize[1]), real_c(2) * real_c(setup.domainSize[2]), real_c(0),
           real_c(0), real_c(0.5) * real_c(setup.domainSize[2]));

   // create blockforest with periodicity in side-directions
   auto forest = blockforest::createBlockForest(domainAABB, setup.numBlocks, Vector3< bool >(true, true, false));

   // generate IDs of specified PE body types
   pe::SetBodyTypeIDs< BodyTypeTuple >::execute();

   // add global body storage for PE bodies
   shared_ptr< pe::BodyStorage > globalBodyStorage = make_shared< pe::BodyStorage >();

   // add block-local body storage
   const auto bodyStorageID = forest->addBlockData(pe::createStorageDataHandling< BodyTypeTuple >(), "bodyStorage");

   // add data-handling for coarse collision detection
   const auto ccdID =
      forest->addBlockData(pe::ccd::createHashGridsDataHandling(globalBodyStorage, bodyStorageID), "ccd");

   // add data handling for fine collision detection
   const auto fcdID = forest->addBlockData(
      pe::fcd::createGenericFCDDataHandling< BodyTypeTuple, pe::fcd::AnalyticCollideFunctor >(), "fcd");

   // add PE contact solver and set its parameters
   const auto cr = make_shared< pe::cr::HCSITS >(globalBodyStorage, forest, bodyStorageID, ccdID, fcdID, nullptr);
   cr->setMaxIterations(setup.maxIterations);
   cr->setRelaxationModel(
      pe::cr::HardContactSemiImplicitTimesteppingSolvers::ApproximateInelasticCoulombContactByDecoupling);
   cr->setRelaxationParameter(setup.relaxationParameter);
   cr->setGlobalLinearAcceleration(setup.globalLinearAcceleration);

   // create a structure of regularly aligned spherical particles with slight random offsets using a loop; the loop
   // starts from a certain position and iterates over a pre-defined spacing

   // point at which particle generation starts
   const Vector3< real_t > genPoint(real_c(0.55) * real_c(setup.particleDiameter),
                                    real_c(0.55) * real_c(setup.particleDiameter),
                                    real_c(setup.domainSize[2]) - real_c(0.55) * real_c(setup.particleDiameter));

   uint_t numParticles = uint_t(0);                                     // number of created particles
   const real_t spacing(real_c(1.25) * real_c(setup.particleDiameter)); // spacing between the particles

   for (auto it = grid_generator::SCIterator(genAABB, genPoint, spacing); it != grid_generator::SCIterator(); ++it)
   {
      // create slight random disturbances in the positions of the created particles
      Vector3< real_t > rndPosDist(real_c(0.1) * math::realRandom< real_t >(-spacing, spacing),
                                   real_c(0.1) * math::realRandom< real_t >(-spacing, spacing),
                                   real_c(0.1) * math::realRandom< real_t >(-spacing, spacing));

      // create a sphere particle with default material (iron, cf. src/pe/Materials.h)
      pe::createSphere(*globalBodyStorage, *forest, bodyStorageID, 0, *it + rndPosDist,
                       real_c(setup.particleDiameter) * real_c(0.5));
      ++numParticles;
      if (numParticles >= setup.numInitialParticles) { break; }
   }

   // create the bottom plane
   pe::createPlane(*globalBodyStorage, 1, Vector3< real_t >(0, 0, 1), domainAABB.minCorner());

   // communicate the created particles between processes
   pe::syncShadowOwners< BodyTypeTuple >(*forest, bodyStorageID);

   // add vtk output for the particles
   const auto vtkOutput = make_shared< pe::SphereVtkOutput >(bodyStorageID, *forest);
   const auto vtkWriter = vtk::createVTKOutput_PointData(vtkOutput, "spheres");

   // add vtk output for the domain decomposition
   vtk::writeDomainDecomposition(forest);

   ///////////////
   // TIME LOOP //
   ///////////////

   for (uint_t timestep = uint_c(0); timestep < setup.maxTimesteps; ++timestep)
   {
      // check for convergence every 500 time steps
      if (timestep % uint_c(500) == uint_c(0) && timestep > uint_c(0))
      {
         // get the absolute maximum velocity of all particles
         const auto particleMaxVel = getParticleMaxVel(*forest, bodyStorageID);
         WALBERLA_LOG_INFO_ON_ROOT("Time step: " << timestep
                                                 << "\t\tMaximum velocity of all particles: " << particleMaxVel);

         if (particleMaxVel < setup.convVelocity)
         {
            WALBERLA_LOG_INFO_ON_ROOT("Converged after " << timestep << " timesteps.")
            break;
         }
      }
      // perform a PE timestep
      cr->timestep(setup.dt);

      pe::syncShadowOwners< BodyTypeTuple >(*forest, bodyStorageID); // communicate

      if (timestep % setup.vtkInterval == uint_c(0)) { vtkWriter->write(); }
   }

   pe::syncShadowOwners< BodyTypeTuple >(*forest, bodyStorageID); // communicate

   ////////////////////
   // POSTPROCESSING //
   ////////////////////

   // write settled packed bed to vtk output
   vtkWriter->write();

   // delete particles the center point of which is below 3 * particleDiameter
   for (auto blockIt = forest->begin(); blockIt != forest->end(); ++blockIt)
   {
      for (auto bodyIt = pe::LocalBodyIterator::begin(*blockIt, bodyStorageID); bodyIt != pe::LocalBodyIterator::end();
           ++bodyIt)
      {
         if (bodyIt->getPosition()[2] < real_c(setup.particleDiameter) * real_c(3)) { bodyIt->markForDeletion(); }
      }
   }

   pe::syncShadowOwners< BodyTypeTuple >(*forest, bodyStorageID); // communicate

   // write packed bed after removal of bottom particles to vtk output
   vtkWriter->write();

   // find position of the lowest and highest particle (closest to bottom and top, respectively) to get the center of
   // the packed bed
   real_t lowestPosition  = real_c(domainAABB.maxCorner()[2]); // initialize with the opposite maximum value
   real_t highestPosition = real_c(domainAABB.minCorner()[2]);
   for (auto blockIt = forest->begin(); blockIt != forest->end(); ++blockIt)
   {
      for (auto bodyIt = pe::LocalBodyIterator::begin(*blockIt, bodyStorageID); bodyIt != pe::LocalBodyIterator::end();
           ++bodyIt)
      {
         real_t zPos = bodyIt->getPosition()[2];
         if (zPos < lowestPosition) { lowestPosition = zPos; }
         if (zPos > highestPosition) { highestPosition = zPos; }
      }
   }

   // communicate lowest and highest particle positions to all participating processes
   mpi::allReduceInplace< real_t >(lowestPosition, mpi::MIN);
   mpi::allReduceInplace< real_t >(highestPosition, mpi::MAX);

   // shift packed bed such that the packed bed's center in z-direction lays in the domain's center in z-direction
   for (auto blockIt = forest->begin(); blockIt != forest->end(); ++blockIt)
   {
      for (auto bodyIt = pe::LocalBodyIterator::begin(*blockIt, bodyStorageID); bodyIt != pe::LocalBodyIterator::end();
           ++bodyIt)
      {
         bodyIt->translate(real_c(0), real_c(0),
                           real_c(setup.domainSize[2]) * real_c(0.5) -
                              (lowestPosition + (highestPosition - lowestPosition) * real_c(0.5)));
      }
   }
   pe::syncShadowOwners< BodyTypeTuple >(*forest, bodyStorageID); // communicate

   // write final packed bed to vtk output
   vtkWriter->write();

   // write particle positions (relative to particle diameter) to file
   writePorousGeometry(*forest, bodyStorageID, setup.particleDiameter, setup.porousGeometryFile);

   return EXIT_SUCCESS;
}

void writePorousGeometry(BlockForest& blockForest, BlockDataID bodyStorageID, uint_t particleDiameter,
                         const std::string& filename)
{
   std::string bodyStats;

   // create file header
   WALBERLA_ROOT_SECTION() { bodyStats += std::string("x/Dp\t") + std::string("y/Dp\t") + std::string("z/Dp\n"); }

   // write to file from each process
   for (auto blockIt = blockForest.begin(); blockIt != blockForest.end(); ++blockIt)
   {
      for (auto bodyIt = pe::LocalBodyIterator::begin(*blockIt, bodyStorageID); bodyIt != pe::LocalBodyIterator::end();
           ++bodyIt)
      {
         Vector3< real_t > bodyPos = bodyIt->getPosition();
         bodyPos /= real_c(particleDiameter);
         bodyStats +=
            std::to_string(bodyPos[0]) + "\t" + std::to_string(bodyPos[1]) + "\t" + std::to_string(bodyPos[2]) + "\n";
      }
   }

   // write text file
   mpi::writeMPITextFile(filename, bodyStats);
}

real_t getParticleMaxVel(BlockForest& blockForest, const BlockDataID& bodyStorageID)
{
   real_t maxVel = real_c(0);

   for (auto blockIt = blockForest.begin(); blockIt != blockForest.end(); ++blockIt)
   {
      for (auto bodyIt = pe::LocalBodyIterator::begin(*blockIt, bodyStorageID); bodyIt != pe::LocalBodyIterator::end();
           ++bodyIt)
      {
         const Vector3< real_t > bodyVel = bodyIt->getLinearVel();
         real_t bodyVelAbsMax;
         if (math::abs(bodyVel.max()) > math::abs(bodyVel.min())) { bodyVelAbsMax = math::abs(bodyVel.max()); }
         else
         {
            bodyVelAbsMax = math::abs(bodyVel.min());
         }
         if (bodyVelAbsMax > maxVel) { maxVel = bodyVelAbsMax; }
      }
   }

   mpi::allReduceInplace< real_t >(maxVel, mpi::MAX);
   return maxVel;
}
} // namespace walberla

int main(int argc, char** argv) { return walberla::main(argc, argv); }