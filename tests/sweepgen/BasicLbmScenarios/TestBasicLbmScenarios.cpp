
#include "blockforest/all.h"

#include "core/all.h"

#include "geometry/all.h"

#include "vtk/all.h"

#include <iostream>
#include <limits>
#include <random>
#include <map>

#include "SimDomain.hpp"

#include "walberla/experimental/Sweep.hpp"

namespace BasicLbmScenarios
{
using namespace walberla;
using namespace walberla::experimental;

using TestFunction = std::function< void(mpi::Environment&) >;

std::string vtk_folder{"vtk_out"};

/**
 * Fully periodic force-driven flow.
 * The velocity in each cell should be steadily increasing.
 */
void fullyPeriodic(mpi::Environment&)
{
   uint_t numProcesses = uint_c(MPIManager::instance()->numProcesses());
   Vector3< uint_t > numBlocks{ math::getFactors3D(numProcesses) };

   SimDomain dom{ SimDomainBuilder{ .blocks        = { numBlocks[0], numBlocks[1], numBlocks[2] },
                                    .cellsPerBlock = { 16, 16, 16 },
                                    .periodic      = { true, true, true } }
                     .build() };

   const Vector3< real_t > force{ 0.005, 0., 0. };

   dom.initConstant(1.0, Vector3< real_t >{ 0.0 }, force);

   auto streamCollide = dom.streamCollideSweep(1.0, force);

   sweep::SerialSweeper sweeper{ dom.blocks };

   for (uint_t t = 0; t < 10; ++t)
   {
      for(auto& b: *dom.blocks) { streamCollide(&b); }
      dom.syncGhostLayers();

      dom.fields2host();

      for(auto& block: *dom.blocks) {
         const VectorField_T& velField = *block.template getData< VectorField_T >(dom.cpuFields.uId);

         sweeper.forAllCells([&](Cell c) {
            real_t expected{ real_c(t) * force[0] };
            real_t actual{ velField.get(c, 0) };
            WALBERLA_CHECK_FLOAT_EQUAL(expected, actual);
         });
      }
   }
}

/**
 * Periodic channel flow with a no-slip boundary at the top
 * and a symmetry plane at the bottom implemented using the free-slip boundary.
 * Flow is governed by the Hagen-Poiseuille-law,
 * with the maximum at the bottom.
 */
void mirroredHalfChannel(mpi::Environment&)
{
   size_t zCells{ 64 };
   uint_t numProcesses = uint_c(MPIManager::instance()->numProcesses());
   std::vector< uint_t > numBlocksXY{ math::getFactors(numProcesses, 2u) };

   SimDomain dom{ SimDomainBuilder{ .blocks        = { numBlocksXY[0], numBlocksXY[1], 1 },
                                    .cellsPerBlock = { 4, 4, zCells },
                                    .periodic      = { true, true, false } }
                     .build() };

   /* Hagen-Poiseuille-law in lattice units */
   const real_t u_max{ 0.025 };
   const real_t reynolds{ 10.0 };
   const real_t L_z{ real_c(2 * zCells) };
   const real_t radius{ L_z / 2.0 };
   const real_t r_squared{ radius * radius };
   const real_t lattice_viscosity{ L_z * u_max / reynolds };
   const real_t omega{ 2. / (6. * lattice_viscosity + 1.) };
   const real_t acceleration{ (u_max * 2.0 * lattice_viscosity) / r_squared };

   const Vector3< real_t > force{ acceleration, 0., 0. };

   auto velocityProfile = [&](real_t z) -> real_t {
      return acceleration / (2.0 * lattice_viscosity) * (r_squared - z * z);
   };

   sweep::SerialSweeper sweeper{ dom.blocks };

   for(auto& block: *dom.blocks) {
      ScalarField_T& densityField  = *block.getData< ScalarField_T >(dom.cpuFields.rhoId);
      VectorField_T& velocityField = *block.getData< VectorField_T >(dom.cpuFields.uId);

      sweeper.forAllCells([&](Cell c) {
         Cell globalCell{ c };
         dom.blocks->transformBlockLocalToGlobalCell(globalCell, block);
         Vector3< real_t > cellCenter{ dom.blocks->getCellCenter(globalCell) };

         densityField.get(c)     = 1.0;
         velocityField.get(c, 0) = velocityProfile(cellCenter[2]);
         velocityField.get(c, 1) = velocityField.get(c, 2) = 0.;
      });
   }

   dom.fields2device();
   dom.initFromFields(force);

   auto streamCollide  = dom.streamCollideSweep(omega, force);
   auto noSlipTop      = dom.noSlipTop();
   auto freeSlipBottom = dom.freeSlipBottom();
   auto velOutput      = field::createVTKOutput< VectorField_T >(dom.cpuFields.uId, *dom.blocks, "vel", uint_t{1u}, uint_t{0u}, false, vtk_folder);

   for (uint_t t = 0; t < 50; ++t)
   {
      for(auto& b: *dom.blocks) { streamCollide(&b); }
      dom.syncGhostLayers();

      dom.fields2host();

      for(auto& block: *dom.blocks) {
         const VectorField_T& velField = *block.template getData< VectorField_T >(dom.cpuFields.uId);

         sweeper.forAllCells([&](Cell c) {
            Cell globalCell{ c };
            dom.blocks->transformBlockLocalToGlobalCell(globalCell, block);
            Vector3< real_t > cellCenter{ dom.blocks->getCellCenter(globalCell) };

            real_t expected{ velocityProfile(cellCenter[2]) };
            real_t actual{ velField.get(c, 0) };
            WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(expected, actual, 1e-5);
         });
      }

      velOutput();

      for(auto& b: *dom.blocks) {
         noSlipTop(&b);
         freeSlipBottom(&b);
      }
   }
}

/**
 * Pipe flow with circular cross-section and free-slip walls.
 * The pipe flow is initialized with a constant velocity in x-direction.
 * As free-slip walls do not impose any friction, velocity should remain constant in time.
 */
void freeSlipPipe(mpi::Environment&)
{
   uint_t numProcesses = uint_c(MPIManager::instance()->numProcesses());

   SimDomain dom{ SimDomainBuilder{
      .blocks = { numProcesses, 1, 1 }, .cellsPerBlock = { 4, 32, 32 }, .periodic = { true, false, false } }
                     .build() };

   const FlagUID fluidFlagUid{ "Fluid" };
   const FlagUID freeSlipFlagUID{ "FreeSlip" };

   const CellInterval allCells{ { 0, 0, 0 },
                                { dom.blocks->getNumberOfXCellsPerBlock() - 1,
                                  dom.blocks->getNumberOfYCellsPerBlock() - 1,
                                  dom.blocks->getNumberOfZCellsPerBlock() - 1 } };
   const CellInterval allCellsWithGl{ { -1, 0, 0 },
                                      { dom.blocks->getNumberOfXCellsPerBlock(),
                                        dom.blocks->getNumberOfYCellsPerBlock() - 1,
                                        dom.blocks->getNumberOfZCellsPerBlock() - 1 } };

   const real_t pipeRadius{ 14.0 };
   const Vector3< real_t > pipeAnchor{ 0.0, 16.0, 16.0 };
   const real_t maxVelocity{ 0.02 };
   const Vector3< real_t > force{ 0., 0., 0. };

   for(auto& block: *dom.blocks) {
      FlagField_T& flagField = *block.getData< FlagField_T >(dom.cpuFields.flagFieldId);
      const uint8_t freeSlipFlag{ flagField.getOrRegisterFlag(freeSlipFlagUID) };

      ScalarField_T& densField = *block.getData< ScalarField_T >(dom.cpuFields.rhoId);
      VectorField_T& velField  = *block.getData< VectorField_T >(dom.cpuFields.uId);

      for (Cell c : allCellsWithGl)
      {
         Cell globalCell{ c };
         dom.blocks->transformBlockLocalToGlobalCell(globalCell, block);
         Vector3< real_t > cellCenter{ dom.blocks->getCellCenter(globalCell) };
         cellCenter[0] = 0.0;

         Vector3< real_t > initVelocity;
         real_t radialDistance = (cellCenter - pipeAnchor).length();
         if (radialDistance > pipeRadius)
         {
            flagField.addFlag(c, freeSlipFlag);
            initVelocity = Vector3< real_t >{ NAN };
         }
         else { initVelocity = Vector3< real_t >{ maxVelocity, 0., 0. }; }

         densField.get(c, 0) = 1.0;
         velField.get(c, 0)  = initVelocity[0];
         velField.get(c, 1)  = initVelocity[1];
         velField.get(c, 2)  = initVelocity[2];
      }
   }

   geometry::setNonBoundaryCellsToDomain< FlagField_T >(*dom.blocks, dom.cpuFields.flagFieldId, fluidFlagUid);

   auto flagsOutput = field::createVTKOutput< FlagField_T >(dom.cpuFields.flagFieldId, *dom.blocks, "flags", uint_t{1u}, uint_t{0u}, false, vtk_folder);
   flagsOutput();

   dom.fields2device();
   dom.initFromFields(force);

   gen::bc_sparse::FreeSlipIrregular freeSlip{ dom.irregularFreeSlipFactory().fromFlagField< FlagField_T >(
      dom.cpuFields.flagFieldId, freeSlipFlagUID, fluidFlagUid) };

   auto streamCollide = dom.streamCollideSweep(1.0, force);

   auto velOutput = field::createVTKOutput< VectorField_T >(dom.cpuFields.uId, *dom.blocks, "vel", uint_t{1u}, uint_t{0u}, false, vtk_folder);

   sweep::SerialSweeper sweeper{ dom.blocks };

   for (uint_t i = 0; i < 10; ++i)
   {
      for(auto& block: *dom.blocks) { streamCollide(&block); }
      dom.syncGhostLayers();

      dom.fields2host();

      for(auto& block: *dom.blocks) {
         VectorField_T& velField = *block.getData< VectorField_T >(dom.cpuFields.uId);
         FlagField_T& flagField  = *block.getData< FlagField_T >(dom.cpuFields.flagFieldId);
         uint8_t fluidFlag       = flagField.getFlag(fluidFlagUid);

         sweeper.forAllCells([&](Cell c) {
            if (flagField.isFlagSet(c, fluidFlag))
            {
               WALBERLA_CHECK_FLOAT_EQUAL(velField.get(c, 0), maxVelocity);
               WALBERLA_CHECK_FLOAT_EQUAL(velField.get(c, 1), 0.);
               WALBERLA_CHECK_FLOAT_EQUAL(velField.get(c, 2), 0.);
            }
         });
      }

      velOutput();

      for (auto& block : *dom.blocks)
      {
         freeSlip(&block);
      }
   }

   velOutput();
}

std::map< std::string, TestFunction > TESTS{
   { "FullyPeriodic", fullyPeriodic },
   { "MirroredHalfChannel", mirroredHalfChannel },
   { "FreeSlipPipe", freeSlipPipe },
};

} // namespace BasicLbmScenarios

int main(int argc, char** argv)
{
   walberla::mpi::Environment env{ argc, argv };

   if (argc < 1)
   {
      std::cerr << "No Test ID was specified." << std::endl;
      return -1;
   }

   std::string testId{ argv[1] };
   if (BasicLbmScenarios::TESTS.contains(testId))
   {
      std::random_device rd;
      std::uniform_int_distribution<std::size_t> generator(std::size_t{0u}, std::numeric_limits<std::size_t>::max() / 1024ul);
      auto unique_id = generator(rd) * 1024ul;
      unique_id += walberla::uint_c(walberla::mpi::MPIManager::instance()->numProcesses());
      walberla::mpi::broadcastObject(unique_id);
      BasicLbmScenarios::vtk_folder += "_" + testId + "_" + std::to_string(unique_id);
      BasicLbmScenarios::TESTS.at(testId)(env);
      return EXIT_SUCCESS;
   }

   WALBERLA_ABORT("Invalid test ID: " << testId);
}
