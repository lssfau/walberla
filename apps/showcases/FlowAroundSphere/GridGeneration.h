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
//! \file GridGeneration.h
//! \author Markus Holzer <markus.holzer@fau.de>
//
//======================================================================================================================
#pragma once

#include "blockforest/Initialization.h"
#include "blockforest/SetupBlock.h"
#include "blockforest/SetupBlockForest.h"
#include "blockforest/loadbalancing/StaticCurve.h"

#include "core/Environment.h"
#include "core/logging/Initialization.h"
#include "core/timing/RemainingTimeLogger.h"
#include "core/timing/TimingPool.h"

#include <string>
#include <fstream>

#include "Types.h"
#include "Setup.h"
#include "Sphere.h"

using namespace walberla;

uint_t blockCells(const Setup& setup, const bool withGhostLayers){
   uint_t cells;
   if(!withGhostLayers){
      cells = setup.cellsPerBlock[0] * setup.cellsPerBlock[1] * setup.cellsPerBlock[2];
   }
   else{
      cells = (setup.cellsPerBlock[0] + 2 * setup.numGhostLayers) *
              (setup.cellsPerBlock[1] + 2 * setup.numGhostLayers) *
              (setup.cellsPerBlock[2] + 2 * setup.numGhostLayers);
   }
   return cells;
}


void createSetupBlockForest(SetupBlockForest& setupBfs, const Setup& setup, const bool useMPIManager=false)
{
   WALBERLA_LOG_INFO_ON_ROOT("Generating SetupBlockForest...")

   uint_t numProcesses = setup.numProcesses;
   const std::string blockForestFilestem = setup.blockForestFilestem;

   if(useMPIManager) {numProcesses = uint_c(mpi::MPIManager::instance()->numProcesses());}

   Sphere Sphere( setup );
   SphereRefinementSelection SphereRefinementSelection( Sphere, setup.refinementLevels );
   SphereBlockExclusion SphereBlockExclusion( Sphere );

   setupBfs.addRefinementSelectionFunction(std::function<void(SetupBlockForest &)>(SphereRefinementSelection));
   setupBfs.addBlockExclusionFunction(SphereBlockExclusion);

   const AABB domain(real_t(0.0), real_t(0.0), real_t(0.0), setup.domainSize[0], setup.domainSize[1], setup.domainSize[2]);
   setupBfs.addWorkloadMemorySUIDAssignmentFunction(blockforest::uniformWorkloadAndMemoryAssignment);
   setupBfs.init(domain, setup.rootBlocks[0], setup.rootBlocks[1], setup.rootBlocks[2],
                         setup.periodic[0], setup.periodic[1], setup.periodic[2]);
   setupBfs.balanceLoad(blockforest::StaticLevelwiseCurveBalanceWeighted(), numProcesses);

   WALBERLA_LOG_INFO_ON_ROOT("===========================  BLOCK FOREST STATISTICS ============================");
   WALBERLA_LOG_INFO_ON_ROOT("Blocks created: " << setupBfs.getNumberOfBlocks())
   for (uint_t level = 0; level <= setup.refinementLevels; level++){
      const uint_t numberOfBlocks = setupBfs.getNumberOfBlocks(level);
      WALBERLA_LOG_INFO_ON_ROOT("Level " << level << " Blocks: " << numberOfBlocks)
   }

   const real_t avgBlocksPerProc = real_c(setupBfs.getNumberOfBlocks()) / real_c(setupBfs.getNumberOfProcesses());
   WALBERLA_LOG_INFO_ON_ROOT("Average blocks per process: " << avgBlocksPerProc);

   const uint_t cellsPerBlock        = blockCells(setup, false);
   const uint_t totalNumberCells     = setupBfs.getNumberOfBlocks() * cellsPerBlock;
   const real_t averageCellsPerGPU   = avgBlocksPerProc * real_c(cellsPerBlock);
   const uint_t cellsPerBlockGL      = blockCells(setup, true);
   const uint_t totalNumberCellsGL   = setupBfs.getNumberOfBlocks() * cellsPerBlockGL;
   const real_t averageCellsPerGPUGL = avgBlocksPerProc * real_c(cellsPerBlockGL);

   const uint_t PDFsPerCell            = StorageSpecification_T::inplace ? Stencil_T::Q : 2 * Stencil_T::Q;
   const uint_t valuesPerCell          = (PDFsPerCell + VelocityField_T::F_SIZE + ScalarField_T::F_SIZE);
   const uint_t sizePerValue           = sizeof(StorageSpecification_T::value_type);
   const double expectedMemory         = double_c(totalNumberCells * valuesPerCell * sizePerValue) * 1e-9;
   const double expectedMemoryPerGPU   = double_c(averageCellsPerGPU * valuesPerCell * sizePerValue) * 1e-9;
   const double expectedMemoryGL       = double_c(totalNumberCellsGL * valuesPerCell * sizePerValue) * 1e-9;
   const double expectedMemoryPerGPUGL = double_c(averageCellsPerGPUGL * valuesPerCell * sizePerValue) * 1e-9;

   WALBERLA_LOG_INFO_ON_ROOT( "Total number of cells will be " << totalNumberCells << " fluid cells (in total on all levels)")
   WALBERLA_LOG_INFO_ON_ROOT( "Expected total memory demand will be " << expectedMemory << " GB")
   WALBERLA_LOG_INFO_ON_ROOT( "Average memory demand per GPU will be " << expectedMemoryPerGPU << " GB")
   WALBERLA_LOG_INFO_ON_ROOT( "Expected total memory demand (with GL) will be " << expectedMemoryGL << " GB")
   WALBERLA_LOG_INFO_ON_ROOT( "Average memory demand per GPU (with GL) will be " << expectedMemoryPerGPUGL << " GB")

   WALBERLA_LOG_INFO_ON_ROOT("=================================================================================")

   if(mpi::MPIManager::instance()->numProcesses() > 1)
      return;

   if(setup.writeSetupForestAndReturn){
      std::ostringstream oss;
      oss << blockForestFilestem << ".bfs";
      setupBfs.saveToFile(oss.str().c_str());

      setupBfs.writeVTKOutput(blockForestFilestem);
   }
}

void createBlockForest(shared_ptr< BlockForest >& bfs, const Setup& setup){
   if (mpi::MPIManager::instance()->numProcesses() > 1){
      // Load structured block forest from file
      std::ostringstream oss;
      oss << setup.blockForestFilestem << ".bfs";
      const std::string setupBlockForestFilepath = oss.str();
      std::ifstream infile(setupBlockForestFilepath.c_str());
      if(!infile.good()){
         WALBERLA_LOG_WARNING_ON_ROOT("Blockforest was not created beforehand and thus needs to be created on the fly. For large simulation runs this can be a severe problem!")
         SetupBlockForest setupBfs;
         createSetupBlockForest(setupBfs, setup, true);
         bfs = std::make_shared< BlockForest >(uint_c(MPIManager::instance()->worldRank()), setupBfs);
      }
      else{
         bfs = std::make_shared< BlockForest >(uint_c(MPIManager::instance()->worldRank()),
                                               setupBlockForestFilepath.c_str(), false);
      }
   }
   else{
      SetupBlockForest setupBfs;
      createSetupBlockForest(setupBfs, setup);
      bfs = std::make_shared< BlockForest >(uint_c(MPIManager::instance()->worldRank()), setupBfs);
   }
}