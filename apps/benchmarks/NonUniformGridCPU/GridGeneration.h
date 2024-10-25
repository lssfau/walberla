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

#include "LdcSetup.h"
#include "NonUniformGridCPUInfoHeader.h"

using StorageSpecification_T = lbm::NonUniformGridCPUStorageSpecification;
using Stencil_T              = StorageSpecification_T::Stencil;

using namespace walberla;

void createSetupBlockForest(SetupBlockForest& setupBfs,
                            const Config::BlockHandle& domainSetup, const Config::BlockHandle& blockForestSetup,
                            const bool useMPIManager=false)
{
   WALBERLA_LOG_INFO_ON_ROOT("Generating SetupBlockForest...")

   Vector3<real_t> domainSize = domainSetup.getParameter<Vector3<real_t> >("domainSize");
   Vector3< uint_t > cellsPerBlock = domainSetup.getParameter< Vector3< uint_t > >("cellsPerBlock");
   Vector3<uint_t> rootBlocks = domainSetup.getParameter<Vector3<uint_t> >("rootBlocks");
   Vector3<bool> periodic = domainSetup.getParameter<Vector3<bool> >("periodic");

   const uint_t refinementDepth = blockForestSetup.getParameter< uint_t >("refinementDepth", uint_c(1));
   uint_t numProcesses = blockForestSetup.getParameter< uint_t >( "numProcesses");
   const std::string blockForestFilestem = blockForestSetup.getParameter< std::string > ("blockForestFilestem", "blockforest");
   const bool writeVtk = blockForestSetup.getParameter< bool >("writeVtk", false);
   const bool outputStatistics = blockForestSetup.getParameter< bool >("outputStatistics", false);

   if(useMPIManager)
      numProcesses = uint_c(mpi::MPIManager::instance()->numProcesses());

   const LDC ldc(refinementDepth);

   auto refSelection = ldc.refinementSelector();
   setupBfs.addRefinementSelectionFunction(std::function<void(SetupBlockForest &)>(refSelection));
   const AABB domain(real_t(0.0), real_t(0.0), real_t(0.0), domainSize[0], domainSize[1], domainSize[2]);
   setupBfs.addWorkloadMemorySUIDAssignmentFunction(blockforest::uniformWorkloadAndMemoryAssignment);
   setupBfs.init(domain, rootBlocks[0], rootBlocks[1], rootBlocks[2], periodic[0], periodic[1], periodic[2]);
   setupBfs.balanceLoad(blockforest::StaticLevelwiseCurveBalanceWeighted(), numProcesses);

   if(mpi::MPIManager::instance()->numProcesses() > 1)
      return;

   {
      std::ostringstream oss;
      oss << blockForestFilestem << ".bfs";
      setupBfs.saveToFile(oss.str().c_str());
   }

   if(writeVtk){
      setupBfs.writeVTKOutput(blockForestFilestem);
   }

   if(outputStatistics){
      WALBERLA_LOG_INFO_ON_ROOT("===========================  BLOCK FOREST STATISTICS ============================");
      WALBERLA_LOG_INFO_ON_ROOT("Blocks created: " << setupBfs.getNumberOfBlocks())
      for (uint_t level = 0; level <= refinementDepth; level++)
      {
         const uint_t numberOfBlocks = setupBfs.getNumberOfBlocks(level);
         WALBERLA_LOG_INFO_ON_ROOT("Level " << level << " Blocks: " << numberOfBlocks)
      }

      const real_t avgBlocksPerProc = real_c(setupBfs.getNumberOfBlocks()) / real_c(setupBfs.getNumberOfProcesses());
      WALBERLA_LOG_INFO_ON_ROOT("Average blocks per process: " << avgBlocksPerProc);

      const uint_t totalNumberCells = setupBfs.getNumberOfBlocks() * cellsPerBlock[0] * cellsPerBlock[1] * cellsPerBlock[2];
      const real_t averageCellsPerGPU = avgBlocksPerProc * real_c(cellsPerBlock[0] * cellsPerBlock[1] * cellsPerBlock[2]);

      const uint_t PDFsPerCell = StorageSpecification_T::inplace ? Stencil_T::Q : 2 * Stencil_T::Q;
      const uint_t valuesPerCell = (PDFsPerCell + VelocityField_T::F_SIZE + ScalarField_T::F_SIZE);
      const uint_t sizePerValue = sizeof(StorageSpecification_T::value_type);
      const double expectedMemory = double_c(totalNumberCells * valuesPerCell * sizePerValue) * 1e-9;
      const double expectedMemoryPerGPU = double_c(averageCellsPerGPU * valuesPerCell * sizePerValue) * 1e-9;

      WALBERLA_LOG_INFO_ON_ROOT( "Total number of cells will be " << totalNumberCells << " fluid cells (in total on all levels)")
      WALBERLA_LOG_INFO_ON_ROOT( "Expected total memory demand will be " << expectedMemory << " GB")
      WALBERLA_LOG_INFO_ON_ROOT( "Average memory demand per GPU will be " << expectedMemoryPerGPU << " GB")

      WALBERLA_LOG_INFO_ON_ROOT("=================================================================================");
   }
}

void createBlockForest(shared_ptr< BlockForest >& bfs,
                       const Config::BlockHandle& domainSetup, const Config::BlockHandle& blockForestSetup)
{
   if (mpi::MPIManager::instance()->numProcesses() > 1)
   {
      const std::string blockForestFilestem =
         blockForestSetup.getParameter< std::string >("blockForestFilestem", "blockforest");
      // Load structured block forest from file
      std::ostringstream oss;
      oss << blockForestFilestem << ".bfs";
      const std::string setupBlockForestFilepath = oss.str();
      std::ifstream infile(setupBlockForestFilepath.c_str());
      if(!infile.good())
      {
         WALBERLA_LOG_WARNING_ON_ROOT("Blockforest was not created beforehand and thus needs to be created on the fly. For large simulation runs this can be a severe problem!")
         SetupBlockForest setupBfs;
         createSetupBlockForest(setupBfs, domainSetup, blockForestSetup, true);
         bfs = std::make_shared< BlockForest >(uint_c(MPIManager::instance()->worldRank()), setupBfs);
      }
      else
      {
         bfs = std::make_shared< BlockForest >(uint_c(MPIManager::instance()->worldRank()),
                                               setupBlockForestFilepath.c_str(), false);
      }
   }
   else
   {
      SetupBlockForest setupBfs;
      createSetupBlockForest(setupBfs, domainSetup, blockForestSetup);
      bfs = std::make_shared< BlockForest >(uint_c(MPIManager::instance()->worldRank()), setupBfs);
   }
}