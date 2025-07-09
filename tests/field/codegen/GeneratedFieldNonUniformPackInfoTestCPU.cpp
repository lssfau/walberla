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
//! \file GeneratedFieldNonUniformPackInfoTestCPU.cpp
//! \ingroup field
//! \author Philipp Suffa <philipp.suffa@fau.de>
//! \brief Tests if a CPU Field is correctly communicated using generated nonuniform pack info
//
//======================================================================================================================

#include "field/AddToStorage.h"
#include "field/GhostLayerField.h"
#include "field/refinement/PackInfo.h"

#include "blockforest/Initialization.h"
#include "blockforest/SetupBlockForest.h"
#include "blockforest/loadbalancing/StaticCurve.h"
#include "blockforest/communication/NonUniformBufferedScheme.h"
#include "blockforest/communication/GeneratedNonUniformFieldPackInfo.h"

#include "core/debug/TestSubsystem.h"
#include "core/Environment.h"

#include "stencil/D3Q27.h"

#include "ScalarFieldNonUniformCommunicationCPU.h"

namespace walberla {

using Stencil_T = stencil::D3Q27;
using d_type = double;
using Field_T = field::GhostLayerField< d_type, 1 >;


static void refinementSelectionFunction(SetupBlockForest& forest)
{
   const AABB & domain = forest.getDomain();
   const AABB rightCorner( domain.xMax() - domain.xSize() * real_c(0.1), real_c(0.0), real_c(0.0), domain.xMax(), domain.yMax() , domain.zMax() );
   for(auto & block : forest)
   {
      auto & aabb = block.getAABB();
      if( rightCorner.intersects( aabb ) )
      {
         if( block.getLevel() < 1)
            block.setMarker( true );
      }
   }
}

void testValuesInGL(shared_ptr<StructuredBlockForest> blocks, BlockDataID scalarFieldId) {
   for( auto & iBlock : *blocks )
   {
      if (blocks->getLevel(iBlock) == 0)
      {
         auto cpuField = iBlock.getData< Field_T >(scalarFieldId);

         // check periodic communication in z-dir
         WALBERLA_CHECK_FLOAT_EQUAL(cpuField->get(0, 0, 8), 0.0)
         WALBERLA_CHECK_FLOAT_EQUAL(cpuField->get(0, 1, 8), 1.0)
         WALBERLA_CHECK_FLOAT_EQUAL(cpuField->get(1, 0, 8), 1.0)
         WALBERLA_CHECK_FLOAT_EQUAL(cpuField->get(1, 1, 8), 2.0)

         // check fine to coarse communication
         WALBERLA_CHECK_FLOAT_EQUAL(cpuField->get(8, 0, 0), 1.5)
         WALBERLA_CHECK_FLOAT_EQUAL(cpuField->get(8, 1, 0), 3.5)
         WALBERLA_CHECK_FLOAT_EQUAL(cpuField->get(8, 2, 0), 5.5)
         WALBERLA_CHECK_FLOAT_EQUAL(cpuField->get(8, 3, 0), 7.5)
         WALBERLA_CHECK_FLOAT_EQUAL(cpuField->get(-1, 4, 0), 7.5)
         WALBERLA_CHECK_FLOAT_EQUAL(cpuField->get(-1, 5, 0), 9.5)
         WALBERLA_CHECK_FLOAT_EQUAL(cpuField->get(-1, 6, 0), 11.5)
         WALBERLA_CHECK_FLOAT_EQUAL(cpuField->get(-1, 7, 0), 13.5)

         // find fine block in direction E
         auto block             = dynamic_cast< Block* >(&iBlock);
         auto dir               = stencil::Direction::E;
         const auto neighborIdx = blockforest::getBlockNeighborhoodSectionIndex(dir);

         WALBERLA_CHECK_EQUAL(block->getNeighborhoodSectionSize(neighborIdx), uint_t(4))
         WALBERLA_ASSERT(block->neighborhoodSectionHasSmallerBlocks(neighborIdx))

         const BlockID& fineReceiverId = block->getNeighborId(neighborIdx, 0); // bottom left block
         auto fineReceiverBlock        = blocks->getBlock(fineReceiverId);
         auto fineCpuField             = fineReceiverBlock->getData< Field_T >(scalarFieldId);

         // check coarse to fine communication
         WALBERLA_CHECK_FLOAT_EQUAL(fineCpuField->get(-2, 0, 0), 7.0)
         WALBERLA_CHECK_FLOAT_EQUAL(fineCpuField->get(-2, 0, 1), 7.0)
         WALBERLA_CHECK_FLOAT_EQUAL(fineCpuField->get(-2, 1, 0), 7.0)
         WALBERLA_CHECK_FLOAT_EQUAL(fineCpuField->get(-2, 1, 1), 7.0)
         WALBERLA_CHECK_FLOAT_EQUAL(fineCpuField->get(-1, 0, 0), 7.0)
         WALBERLA_CHECK_FLOAT_EQUAL(fineCpuField->get(-1, 0, 1), 7.0)
         WALBERLA_CHECK_FLOAT_EQUAL(fineCpuField->get(-1, 1, 0), 7.0)
         WALBERLA_CHECK_FLOAT_EQUAL(fineCpuField->get(-1, 1, 1), 7.0)
      }
   }
}

int main(int argc, char **argv) {
   debug::enterTestMode();
   const mpi::Environment env(argc, argv);
   mpi::MPIManager::instance()->useWorldComm();

   // Create a BlockForest with 2x2x2 cells per block
   SetupBlockForest setupBfs;
   const AABB domain(real_t(0.0), real_t(0.0), real_t(0.0), real_t(1.0), real_t(1.0), real_t(1.0));
   const Vector3<uint_t>rootBlocks(2,1,1);
   const Vector3<bool>periodic(true,true,true);
   const Vector3<uint_t>cellsPerBlock(8,8,8);

   setupBfs.addRefinementSelectionFunction(refinementSelectionFunction);
   setupBfs.addWorkloadMemorySUIDAssignmentFunction(blockforest::uniformWorkloadAndMemoryAssignment);
   setupBfs.init(domain, rootBlocks[0], rootBlocks[1], rootBlocks[2], periodic[0], periodic[1], periodic[2]);
   setupBfs.balanceLoad(blockforest::StaticLevelwiseCurveBalanceWeighted(), 1);

   auto bfs = std::make_shared< BlockForest >(uint_c(MPIManager::instance()->worldRank()), setupBfs);
   auto blocks = std::make_shared< StructuredBlockForest >(bfs, cellsPerBlock[0], cellsPerBlock[1], cellsPerBlock[2]);

   // Create a Field with the same number of cells as the block
   const BlockDataID scalarFieldForGeneratedPackinfoId = field::addToStorage< Field_T >(blocks, "CPUField", 0.0, field::fzyx, uint_c(2));
   const BlockDataID scalarFieldId = field::addToStorage< Field_T >(blocks, "CPUField", 0.0, field::fzyx, uint_c(2));

   blockforest::communication::NonUniformBufferedScheme< Stencil_T > us{ blocks, false };

   //initialization with value=x+y+z on all blocks
   for( auto & block : *blocks )
   {
      auto cpuField = block.getData< Field_T >(scalarFieldId);
      auto cpuFieldGenerated = block.getData< Field_T >(scalarFieldForGeneratedPackinfoId);
      WALBERLA_FOR_ALL_CELLS_XYZ(cpuField,
                                 cpuField->get(x, y, z) = real_t(x + y + z);
                                 cpuFieldGenerated->get(x, y, z) = real_t(x + y + z);)
   }

   auto packInfoGenerated = make_shared< GeneratedNonUniformFieldPackInfo< Field_T, pystencils::ScalarFieldNonUniformCommunicationCPU > >(scalarFieldForGeneratedPackinfoId);
   us.addPackInfo(packInfoGenerated);

   auto packInfo = make_shared< field::refinement::PackInfo< Field_T, Stencil_T > >(scalarFieldId);
   us.addPackInfo(packInfo);

   // communicate
   us.communicateEqualLevel(0);
   us.communicateEqualLevel(1);
   us.communicateCoarseToFine(1);
   us.communicateFineToCoarse(1);

   //check received values in GL cells
   testValuesInGL(blocks, scalarFieldForGeneratedPackinfoId);
   testValuesInGL(blocks, scalarFieldId);

   //Check if generated packing kernels and handwritten packing kernels match
   for( auto & block : *blocks )
   {
      auto cpuField = block.getData< Field_T >(scalarFieldId);
      auto cpuFieldGenerated = block.getData< Field_T >(scalarFieldForGeneratedPackinfoId);
      WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ(cpuField,
                                 WALBERLA_CHECK_FLOAT_EQUAL(cpuField->get(x, y, z), cpuFieldGenerated->get(x, y, z)))
   }

   return EXIT_SUCCESS;
}
} // namespace walberla

int main( int argc, char* argv[] ) {
   return walberla::main( argc, argv );
}