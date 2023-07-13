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
//! \file NonUniformGridCPU.cpp
//! \author Markus Holzer <markus.holzer@fau.de>
//
//======================================================================================================================

#include "blockforest/Initialization.h"
#include "blockforest/SetupBlockForest.h"
#include "blockforest/loadbalancing/StaticCurve.h"

#include "core/Environment.h"
#include "core/logging/Initialization.h"
#include "core/timing/RemainingTimeLogger.h"
#include "core/timing/TimingPool.h"

#include "field/AddToStorage.h"
#include "field/FlagField.h"
#include "field/vtk/VTKWriter.h"

#include "geometry/InitBoundaryHandling.h"

#include "lbm_generated/communication/NonuniformGeneratedPdfPackInfo.h"
#include "lbm_generated/field/AddToStorage.h"
#include "lbm_generated/field/PdfField.h"
#include "lbm_generated/refinement/BasicRecursiveTimeStep.h"
#include "lbm_generated/evaluation/PerformanceEvaluation.h"

#include "python_coupling/CreateConfig.h"
#include "python_coupling/PythonCallback.h"

#include "timeloop/SweepTimeloop.h"

#include <cmath>

#include "NonUniformGridCPUInfoHeader.h"

using namespace walberla;

using StorageSpecification_T = lbm::NonUniformGridCPUStorageSpecification;
using Stencil_T              = StorageSpecification_T::Stencil;
using CommunicationStencil_T = StorageSpecification_T::CommunicationStencil;

using PdfField_T           = lbm_generated::PdfField< StorageSpecification_T >;
using FlagField_T          = FlagField< uint8_t >;
using BoundaryCollection_T = lbm::NonUniformGridCPUBoundaryCollection< FlagField_T >;

using SweepCollection_T = lbm::NonUniformGridCPUSweepCollection;

using blockforest::communication::NonUniformBufferedScheme;
using RefinementSelectionFunctor = SetupBlockForest::RefinementSelectionFunction;


class LDCRefinement
{
private:
   const uint_t refinementDepth_;

public:
   explicit LDCRefinement(const uint_t depth) : refinementDepth_(depth){};

   void operator()(SetupBlockForest& forest) const
   {
      const AABB & domain = forest.getDomain();

      real_t xSize = ( domain.xSize() / real_t(12) ) * real_c( 0.99 );
      real_t ySize = ( domain.ySize() / real_t(12) ) * real_c( 0.99 );

      AABB leftCorner( domain.xMin(), domain.yMin(), domain.zMin(),
                       domain.xMin() + xSize, domain.yMin() + ySize, domain.zMax() );

      AABB rightCorner( domain.xMax() - xSize, domain.yMin(), domain.zMin(),
                        domain.xMax(), domain.yMin() + ySize, domain.zMax() );

      for(auto & block : forest)
      {
         auto & aabb = block.getAABB();
         if( leftCorner.intersects( aabb ) || rightCorner.intersects( aabb ) )
         {
            if( block.getLevel() < refinementDepth_)
               block.setMarker( true );
         }
      }
   }
};

class LDC
{
private:
   const std::string refinementProfile_;
   const uint_t refinementDepth_;

   const FlagUID noSlipFlagUID_;
   const FlagUID ubbFlagUID_;

public:
   explicit LDC(const uint_t depth) : refinementDepth_(depth), noSlipFlagUID_("NoSlip"), ubbFlagUID_("UBB"){};

   RefinementSelectionFunctor refinementSelector() const
   {
      return LDCRefinement(refinementDepth_);
   }

   void setupBoundaryFlagField(StructuredBlockForest& sbfs, const BlockDataID flagFieldID)
   {
      for (auto bIt = sbfs.begin(); bIt != sbfs.end(); ++bIt)
      {
         auto& b           = dynamic_cast< Block& >(*bIt);
         const uint_t level       = b.getLevel();
         auto flagField     = b.getData< FlagField_T >(flagFieldID);
         const uint8_t noslipFlag = flagField->registerFlag(noSlipFlagUID_);
         const uint8_t ubbFlag    = flagField->registerFlag(ubbFlagUID_);
         for (auto cIt = flagField->beginWithGhostLayerXYZ(2); cIt != flagField->end(); ++cIt)
         {
            const Cell localCell = cIt.cell();
            Cell globalCell(localCell);
            sbfs.transformBlockLocalToGlobalCell(globalCell, b);
            if (globalCell.y() >= cell_idx_c(sbfs.getNumberOfYCells(level))) { flagField->addFlag(localCell, ubbFlag); }
            else if (globalCell.z() < 0 || globalCell.y() < 0 || globalCell.x() < 0 ||
                     globalCell.x() >= cell_idx_c(sbfs.getNumberOfXCells(level)) || globalCell.z() >= cell_idx_c(sbfs.getNumberOfZCells(level)))
            {
               flagField->addFlag(localCell, noslipFlag);
            }
         }
      }
   }
};

static void createSetupBlockForest(SetupBlockForest& setupBfs, const Config::BlockHandle& domainSetup, LDC& ldcSetup, const uint_t numProcesses=uint_c(MPIManager::instance()->numProcesses()))
{
   Vector3<real_t> domainSize = domainSetup.getParameter<Vector3<real_t> >("domainSize");
   Vector3<uint_t> rootBlocks = domainSetup.getParameter<Vector3<uint_t> >("rootBlocks");
   Vector3<bool> periodic = domainSetup.getParameter<Vector3<bool> >("periodic");

   auto refSelection = ldcSetup.refinementSelector();
   setupBfs.addRefinementSelectionFunction(std::function<void(SetupBlockForest &)>(refSelection));
   const AABB domain(real_t(0.0), real_t(0.0), real_t(0.0), domainSize[0], domainSize[1], domainSize[2]);
   setupBfs.addWorkloadMemorySUIDAssignmentFunction(blockforest::uniformWorkloadAndMemoryAssignment);
   setupBfs.init(domain, rootBlocks[0], rootBlocks[1], rootBlocks[2], periodic[0], periodic[1], periodic[2]);
   setupBfs.balanceLoad(blockforest::StaticLevelwiseCurveBalanceWeighted(), numProcesses);
}

int main(int argc, char** argv)
{
   const mpi::Environment env(argc, argv);
   mpi::MPIManager::instance()->useWorldComm();

   for (auto cfg = python_coupling::configBegin(argc, argv); cfg != python_coupling::configEnd(); ++cfg)
   {
      WALBERLA_MPI_BARRIER()
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      ///                                        SETUP AND CONFIGURATION                                             ///
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      auto config = *cfg;
      logging::configureLogging(config);
      auto domainSetup = config->getOneBlock("DomainSetup");

      // Reading parameters
      auto parameters              = config->getOneBlock("Parameters");
      const real_t omega       = parameters.getParameter< real_t >("omega", real_c(1.4));
      const uint_t refinementDepth = parameters.getParameter< uint_t >("refinementDepth", uint_c(1));
      const uint_t timesteps   = parameters.getParameter< uint_t >("timesteps", uint_c(50));
      const bool writeSetupForestAndReturn = parameters.getParameter< bool >("writeSetupForestAndReturn", false);
      const bool benchmarkKernelOnly = parameters.getParameter< bool >("benchmarkKernelOnly", false);
      const uint_t numProcesses = parameters.getParameter< uint_t >( "numProcesses");

      auto ldc = std::make_shared< LDC >(refinementDepth);
      SetupBlockForest setupBfs;
      if (writeSetupForestAndReturn)
      {
         WALBERLA_LOG_INFO_ON_ROOT("Creating SetupBlockForest for " << numProcesses << " processes")
         WALBERLA_LOG_INFO_ON_ROOT("Generating SetupBlockForest...")
         createSetupBlockForest(setupBfs, domainSetup, *ldc, numProcesses);

         WALBERLA_ROOT_SECTION() { setupBfs.writeVTKOutput("SetupBlockForest"); }

         WALBERLA_LOG_INFO_ON_ROOT("Blocks created: " << setupBfs.getNumberOfBlocks())
         for (uint_t level = 0; level <= refinementDepth; level++)
         {
            const uint_t numberOfBlocks = setupBfs.getNumberOfBlocks(level);
            WALBERLA_LOG_INFO_ON_ROOT("Level " << level << " Blocks: " << numberOfBlocks)
         }

         WALBERLA_LOG_INFO_ON_ROOT("Ending program")
         return EXIT_SUCCESS;
      }

      WALBERLA_LOG_INFO_ON_ROOT("Generating SetupBlockForest...")
      createSetupBlockForest(setupBfs, domainSetup, *ldc);

      // Create structured block forest
      Vector3< uint_t > cellsPerBlock = domainSetup.getParameter< Vector3< uint_t > >("cellsPerBlock");
      WALBERLA_LOG_INFO_ON_ROOT("Creating structured block forest...")
      auto bfs = std::make_shared< BlockForest >(uint_c(MPIManager::instance()->worldRank()), setupBfs);
      auto blocks =
         std::make_shared< StructuredBlockForest >(bfs, cellsPerBlock[0], cellsPerBlock[1], cellsPerBlock[2]);
      blocks->createCellBoundingBoxes();

      WALBERLA_LOG_INFO_ON_ROOT("Blocks created: " << blocks->getNumberOfBlocks())
      for (uint_t level = 0; level <= refinementDepth; level++)
      {
         WALBERLA_LOG_INFO_ON_ROOT("Level " << level << " Blocks: " << blocks->getNumberOfBlocks(level))
      }

      // Creating fields
      const StorageSpecification_T StorageSpec = StorageSpecification_T();
      const BlockDataID pdfFieldID =
         lbm_generated::addPdfFieldToStorage(blocks, "pdfs", StorageSpec, uint_c(2), field::fzyx);
      const BlockDataID velFieldID =
         field::addToStorage< VelocityField_T >(blocks, "vel", real_c(0.0), field::fzyx, uint_c(2));
      const BlockDataID densityFieldID =
         field::addToStorage< ScalarField_T >(blocks, "density", real_c(1.0), field::fzyx, uint_c(2));
      const BlockDataID flagFieldID =
         field::addFlagFieldToStorage< FlagField_T >(blocks, "Boundary Flag Field", uint_c(3));

      const Cell innerOuterSplit =
         Cell(parameters.getParameter< Vector3< cell_idx_t > >("innerOuterSplit", Vector3< cell_idx_t >(1, 1, 1)));
      SweepCollection_T sweepCollection(blocks, pdfFieldID, densityFieldID, velFieldID, omega, innerOuterSplit);
      for (auto& block : *blocks)
      {
         sweepCollection.initialise(&block, 2);
      }
      WALBERLA_MPI_BARRIER()
      WALBERLA_LOG_INFO_ON_ROOT("Initialisation done")
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      ///                                      LB SWEEPS AND BOUNDARY HANDLING                                       ///
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      const FlagUID fluidFlagUID("Fluid");
      ldc->setupBoundaryFlagField(*blocks, flagFieldID);
      geometry::setNonBoundaryCellsToDomain< FlagField_T >(*blocks, flagFieldID, fluidFlagUID, 2);
      BoundaryCollection_T boundaryCollection(blocks, flagFieldID, pdfFieldID, fluidFlagUID);

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      ///                                           COMMUNICATION SCHEME                                             ///
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      WALBERLA_LOG_INFO_ON_ROOT("Setting up communication...")
      auto communication = std::make_shared< NonUniformBufferedScheme< CommunicationStencil_T > >(blocks);
      auto packInfo      = lbm_generated::setupNonuniformPdfCommunication< PdfField_T >(blocks, pdfFieldID);
      communication->addPackInfo(packInfo);

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      ///                                          TIME STEP DEFINITIONS                                             ///
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      lbm_generated::BasicRecursiveTimeStep< PdfField_T, SweepCollection_T, BoundaryCollection_T > LBMMeshRefinement(
         blocks, pdfFieldID, sweepCollection, boundaryCollection, communication, packInfo);

      SweepTimeloop timeLoop(blocks->getBlockStorage(), timesteps);

      if(benchmarkKernelOnly){
         timeLoop.add() << Sweep(sweepCollection.streamCollide(SweepCollection_T::ALL), "LBM StreamCollide");
      }
      else{
         LBMMeshRefinement.addRefinementToTimeLoop(timeLoop);
      }

      // VTK
      const uint_t vtkWriteFrequency = parameters.getParameter< uint_t >("vtkWriteFrequency", 0);
      if (vtkWriteFrequency > 0)
      {
         auto vtkOutput = vtk::createVTKOutput_BlockData(*blocks, "vtk", vtkWriteFrequency, 0, false, "vtk_out",
                                                         "simulation_step", false, true, true, false, 0);
         auto velWriter = make_shared< field::VTKWriter< VelocityField_T > >(velFieldID, "vel");
         vtkOutput->addCellDataWriter(velWriter);

         vtkOutput->addBeforeFunction([&]() {
            for (auto& block : *blocks)
               sweepCollection.calculateMacroscopicParameters(&block);
         });
         timeLoop.addFuncAfterTimeStep(vtk::writeFiles(vtkOutput), "VTK Output");
      }

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      ///                                               BENCHMARK                                                    ///
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      auto remainingTimeLoggerFrequency =
         parameters.getParameter< real_t >("remainingTimeLoggerFrequency", real_c(-1.0)); // in seconds
      if (remainingTimeLoggerFrequency > 0)
      {
         auto logger = timing::RemainingTimeLogger(timeLoop.getNrOfTimeSteps(), remainingTimeLoggerFrequency);
         timeLoop.addFuncAfterTimeStep(logger, "remaining time logger");
      }

      lbm_generated::PerformanceEvaluation<FlagField_T> const performance(blocks, flagFieldID, fluidFlagUID);
      field::CellCounter< FlagField_T > fluidCells( blocks, flagFieldID, fluidFlagUID );
      fluidCells();

      WALBERLA_LOG_INFO_ON_ROOT( "Non uniform Grid benchmark with " << fluidCells.numberOfCells() << " fluid cells (in total on all levels)")

      WcTimingPool timeloopTiming;
      WcTimer simTimer;

      WALBERLA_MPI_BARRIER()
      WALBERLA_LOG_INFO_ON_ROOT("Starting benchmark with " << timesteps << " time steps")
      WALBERLA_MPI_BARRIER()

      simTimer.start();
      timeLoop.run(timeloopTiming);
      WALBERLA_MPI_BARRIER()
      simTimer.end();

      WALBERLA_LOG_INFO_ON_ROOT("Benchmark finished")
      double time = simTimer.max();
      WALBERLA_MPI_SECTION() { walberla::mpi::reduceInplace(time, walberla::mpi::MAX); }
      performance.logResultOnRoot(timesteps, time);

      const auto reducedTimeloopTiming = timeloopTiming.getReduced();
      WALBERLA_LOG_RESULT_ON_ROOT("Time loop timing:\n" << *reducedTimeloopTiming)
   }
   return EXIT_SUCCESS;
}