//========================================================================================================================
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
//! \file GPUBlockSelectorCommunicationTest.cpp
//! \ingroup gpu
//! \author Helen Schottenhamml <helen.schottenhamml@fau.de>
//! \brief Short communication test for the usage of block selectors in UniformGPUScheme.
//
//========================================================================================================================

#include <blockforest/GlobalLoadBalancing.h>
#include <blockforest/Initialization.h>
#include <blockforest/SetupBlockForest.h>

#include <core/DataTypes.h>
#include <core/Environment.h>
#include <core/debug/TestSubsystem.h>
#include <core/math/Random.h>

#include <domain_decomposition/BlockDataID.h>

#include <field/AddToStorage.h>
#include <field/GhostLayerField.h>

#include "gpu/GPUWrapper.h"
#include <gpu/AddGPUFieldToStorage.h>
#include <gpu/FieldCopy.h>
#include <gpu/GPUField.h>
#include <gpu/communication/MemcpyPackInfo.h>
#include <gpu/communication/UniformGPUScheme.h>

#include <stencil/D3Q27.h>
#include <stencil/Directions.h>

namespace walberla
{
using Type_T = int;

using Stencil_T        = stencil::D3Q27;
using ScalarField_T    = field::GhostLayerField< Type_T, 1 >;
using GPUScalarField_T = gpu::GPUField< Type_T >;

const Set< SUID > requiredBlockSelector("communication");
const Set< SUID > incompatibleBlockSelector("no communication");

void suidAssignmentFunction(blockforest::SetupBlockForest& forest)
{
   for (auto& sblock : forest)
   {
      if (forest.atDomainXMinBorder(sblock)) { sblock.addState(incompatibleBlockSelector); }
      else { sblock.addState(requiredBlockSelector); }
      sblock.setWorkload(walberla::numeric_cast< walberla::workload_t >(1));
   }
}

void initScalarField(std::shared_ptr< StructuredBlockForest >& blocks, const BlockDataID& fieldID)
{
   for (auto& block : *blocks)
   {
      Type_T val;
      if (blocks->atDomainXMinBorder(block)) { val = Type_T(-1); }
      else if (blocks->atDomainXMaxBorder(block)) { val = Type_T(1); }
      else { val = Type_T(0); }

      auto* field = block.getData< ScalarField_T >(fieldID);
      WALBERLA_ASSERT_NOT_NULLPTR(field)

      const auto cells = field->xyzSizeWithGhostLayer();

      for (auto cell : cells)
      {
         field->get(cell) = val;
      }
   }
}

std::shared_ptr< StructuredBlockForest >
   createSelectorBlockGrid(const uint_t numberOfXBlocks, const uint_t numberOfYBlocks, const uint_t numberOfZBlocks,
                           const uint_t numberOfXCellsPerBlock, const uint_t numberOfYCellsPerBlock,
                           const uint_t numberOfZCellsPerBlock, const real_t dx, const bool xPeriodic,
                           const bool yPeriodic, const bool zPeriodic, const bool keepGlobalBlockInformation)
{
   // initialize SetupBlockForest = determine domain decomposition

   SetupBlockForest sforest;

   sforest.addWorkloadMemorySUIDAssignmentFunction(suidAssignmentFunction);

   AABB const domainAABB{ real_c(0),
                          real_c(0),
                          real_c(0),
                          dx * real_c(numberOfXBlocks * numberOfXCellsPerBlock),
                          dx * real_c(numberOfYBlocks * numberOfYCellsPerBlock),
                          dx * real_c(numberOfZBlocks * numberOfZCellsPerBlock) };
   sforest.init(domainAABB, numberOfXBlocks, numberOfYBlocks, numberOfZBlocks, xPeriodic, yPeriodic, zPeriodic);

   // calculate process distribution

   const memory_t memoryLimit = numeric_cast< memory_t >(sforest.getNumberOfBlocks());

   blockforest::GlobalLoadBalancing::MetisConfiguration< SetupBlock > const metisConfig(
      true, false,
      std::bind(blockforest::cellWeightedCommunicationCost, std::placeholders::_1, std::placeholders::_2,
                numberOfXCellsPerBlock, numberOfYCellsPerBlock, numberOfZCellsPerBlock));

   sforest.calculateProcessDistribution_Default(uint_c(MPIManager::instance()->numProcesses()), memoryLimit, "hilbert",
                                                10, false, metisConfig);

   if (!MPIManager::instance()->rankValid()) MPIManager::instance()->useWorldComm();

   // create StructuredBlockForest (encapsulates a newly created BlockForest)

   auto bf =
      std::make_shared< BlockForest >(uint_c(MPIManager::instance()->rank()), sforest, keepGlobalBlockInformation);

   auto sbf = std::make_shared< StructuredBlockForest >(bf, numberOfXCellsPerBlock, numberOfYCellsPerBlock,
                                                        numberOfZCellsPerBlock);
   sbf->createCellBoundingBoxes();

   return sbf;
}

int main(int argc, char** argv)
{
   debug::enterTestMode();
   walberla::Environment const walberlaEnv(argc, argv);

   const Vector3< uint_t > nBlocks{ 3, 1, 1 };
   const Vector3< uint_t > cells{ 2, 2, 1 };
   Vector3< real_t > domainSize;
   for (uint_t d = 0; d < 3; ++d)
   {
      domainSize[d] = real_c(cells[d] * nBlocks[d]);
   }

   auto blocks = createSelectorBlockGrid(nBlocks[0], nBlocks[1], nBlocks[2], cells[0], cells[1], cells[2], 1, false,
                                         true, true, true);

   BlockDataID const fieldID = field::addToStorage< ScalarField_T >(blocks, "scalar", Type_T(0), field::fzyx, uint_t(1));
   initScalarField(blocks, fieldID);

   BlockDataID const gpuFieldID = gpu::addGPUFieldToStorage< ScalarField_T >(blocks, fieldID, "GPU scalar");

   // Setup communication schemes for GPUPackInfo
   gpu::communication::UniformGPUScheme< Stencil_T > communication(blocks, requiredBlockSelector, incompatibleBlockSelector);
   communication.addPackInfo(std::make_shared< gpu::communication::MemcpyPackInfo< GPUScalarField_T > >(gpuFieldID));

   // Perform one communication step
   communication();
   WALBERLA_GPU_CHECK(gpuDeviceSynchronize())


   // Copy to CPU
   gpu::fieldCpy< ScalarField_T, GPUScalarField_T >( blocks, fieldID, gpuFieldID );

   // Check for correct data in ghost layers of middle block
   auto middleBlock = blocks->getBlock(domainSize[0] / real_c(2), domainSize[1] / real_c(2), domainSize[2] / real_c(2));
   auto cpuField    = middleBlock->getData< ScalarField_T >(fieldID);
   WALBERLA_ASSERT_NOT_NULLPTR(cpuField)

   // avoid unused variable warning in release mode
   (void) cpuField;
   // WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ(cpuField, WALBERLA_LOG_DEVEL_VAR(cpuField->get(x, y, z)))

   // check for missing communication with left neighbour (first block, incompatible selector)
   WALBERLA_ASSERT_EQUAL(cpuField->get(-1, 0, 0), 0, "Communication with left neighbor detected.")
   WALBERLA_ASSERT_EQUAL(cpuField->get(-1, 1, 0), 0, "Communication with left neighbor detected.")

   // check for correct communication with right neighbor (third block, required selector)
   WALBERLA_ASSERT_EQUAL(cpuField->get(cell_idx_t(cells[0]), 0, 0), 1, "No communication with right neighbor detected.")
   WALBERLA_ASSERT_EQUAL(cpuField->get(cell_idx_t(cells[0]), 1, 0), 1, "No communication with right neighbor detected.")

   return EXIT_SUCCESS;
}

} // namespace walberla

int main(int argc, char** argv) { return walberla::main(argc, argv); }
