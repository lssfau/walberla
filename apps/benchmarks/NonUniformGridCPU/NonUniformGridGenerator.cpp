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
//! \file NonUniformGridGenerator.cpp
//! \author Frederik Hennig <frederik.hennig@fau.de>
//
//======================================================================================================================

#include "blockforest/Initialization.h"
#include "blockforest/SetupBlock.h"
#include "blockforest/SetupBlockForest.h"
#include "blockforest/loadbalancing/StaticCurve.h"

#include "core/all.h"

#include "python_coupling/CreateConfig.h"

#include <string>

#include "LdcSetup.h"

using namespace walberla;


int main(int argc, char ** argv){
   const mpi::Environment env(argc, argv);
   mpi::MPIManager::instance()->useWorldComm();

   if(mpi::MPIManager::instance()->numProcesses() > 1){
      WALBERLA_ABORT("Commandment: Thou shalt not run thy grid generator with more than one process.");
   }

   for (auto cfg = python_coupling::configBegin(argc, argv); cfg != python_coupling::configEnd(); ++cfg)
   {
      auto config = *cfg;
      auto domainSetup = config->getOneBlock("DomainSetup");

      Vector3<real_t> domainSize = domainSetup.getParameter<Vector3<real_t> >("domainSize");
      Vector3<uint_t> rootBlocks = domainSetup.getParameter<Vector3<uint_t> >("rootBlocks");
      Vector3<bool> periodic = domainSetup.getParameter<Vector3<bool> >("periodic");

      auto blockForestSetup = config->getOneBlock("SetupBlockForest");
      const uint_t refinementDepth = blockForestSetup.getParameter< uint_t >("refinementDepth", uint_c(1));
      const uint_t numProcesses = blockForestSetup.getParameter< uint_t >( "numProcesses");
      const std::string blockForestFilestem = blockForestSetup.getParameter< std::string > ("blockForestFilestem", "blockforest");
      const bool writeVtk = blockForestSetup.getParameter< bool >("writeVtk", false);
      const bool outputStatistics = blockForestSetup.getParameter< bool >("outputStatistics", false);

      const LDC ldc(refinementDepth);
      SetupBlockForest setupBfs;

      auto refSelection = ldc.refinementSelector();
      setupBfs.addRefinementSelectionFunction(std::function<void(SetupBlockForest &)>(refSelection));
      const AABB domain(real_t(0.0), real_t(0.0), real_t(0.0), domainSize[0], domainSize[1], domainSize[2]);
      setupBfs.addWorkloadMemorySUIDAssignmentFunction(blockforest::uniformWorkloadAndMemoryAssignment);
      setupBfs.init(domain, rootBlocks[0], rootBlocks[1], rootBlocks[2], periodic[0], periodic[1], periodic[2]);
      setupBfs.balanceLoad(blockforest::StaticLevelwiseCurveBalanceWeighted(), numProcesses);

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

         const uint_t avgBlocksPerProc = setupBfs.getNumberOfBlocks() / setupBfs.getNumberOfProcesses();
         WALBERLA_LOG_INFO_ON_ROOT("Average blocks per process: " << avgBlocksPerProc);
         WALBERLA_LOG_INFO_ON_ROOT("=================================================================================");
      }


      WALBERLA_LOG_INFO_ON_ROOT("Ending program")
   }
}
