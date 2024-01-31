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
//! \file FreeSlipRefinement.cpp
//! \author Markus Holzer <markus.holzer@fau.de>
//! \brief Channel flow with inlet and outlet on West and East. The rest of the Boundaries consist of FreeSlip. The
//!        Channel flow should reach the inlet velocity in the whole domain because the FreeSlip BC will not provide a
//!        restriction on it. With the D3Q27 stencil this only works if the BC is also set on the first fluid node
//
//======================================================================================================================

#include "blockforest/Initialization.h"
#include "blockforest/communication/NonUniformBufferedScheme.h"

#include "core/DataTypes.h"
#include "core/Environment.h"
#include "core/debug/TestSubsystem.h"
#include "core/logging/Initialization.h"
#include "core/math/Vector3.h"
#include "core/timing/RemainingTimeLogger.h"

#include "field/AddToStorage.h"
#include "field/FlagField.h"
#include "field/GhostLayerField.h"
#include "field/vtk/VTKWriter.h"

#include "geometry/InitBoundaryHandling.h"

#include "timeloop/SweepTimeloop.h"

#include "lbm_generated/communication/NonuniformGeneratedPdfPackInfo.h"
#include "lbm_generated/field/AddToStorage.h"
#include "lbm_generated/field/PdfField.h"
#include "lbm_generated/refinement/BasicRecursiveTimeStep.h"

// include the generated header file. It includes all generated classes
#include "FreeSlipRefinementInfoHeader.h"

using namespace walberla;
using namespace std::placeholders;

using StorageSpecification_T = lbm::FreeSlipRefinementStorageSpecification;
using Stencil_T              = StorageSpecification_T::Stencil;
using CommunicationStencil_T = StorageSpecification_T::CommunicationStencil;
using PdfField_T             = lbm_generated::PdfField< StorageSpecification_T >;

using SweepCollection_T = lbm::FreeSlipRefinementSweepCollection;

using VectorField_T = GhostLayerField< real_t, StorageSpecification_T::Stencil::D >;
using ScalarField_T = GhostLayerField< real_t, 1 >;

using flag_t               = walberla::uint8_t;
using FlagField_T          = FlagField< flag_t >;
using BoundaryCollection_T = lbm::FreeSlipRefinementBoundaryCollection< FlagField_T >;

using RefinementSelectionFunctor = SetupBlockForest::RefinementSelectionFunction;

class ChannelRefinement
{
 public:
   ChannelRefinement(const uint_t depth) : refinementDepth_(depth){};

   void operator()(SetupBlockForest& forest)
   {
      std::vector< SetupBlock* > blocks;
      forest.getBlocks(blocks);
      AABB refinementAABB = AABB(forest.getDomain().xSize() / 2 - 1, forest.getDomain().yMin(), forest.getDomain().zSize() / 2 - 1,
                                 forest.getDomain().xSize() / 2 + 1, forest.getDomain().yMin() + 1, forest.getDomain().zSize() / 2 + 1 );

      for (auto b : blocks)
      {
         if (refinementAABB.intersects(b->getAABB()) && b->getLevel() < refinementDepth_)
         {
            b->setMarker(true);
         }
      }
   }

 private:
   const uint_t refinementDepth_;
};

class Channel
{
 public:
   Channel(const uint_t depth) : refinementDepth_(depth), freeSlipFlagUID_("FreeSlip"), ubbFlagUID_("UBB"), outflowFlagUID_("Outflow"){};

   RefinementSelectionFunctor refinementSelector() { return ChannelRefinement(refinementDepth_); }
   void setupBoundaryFlagField(StructuredBlockForest& sbfs, const BlockDataID flagFieldID)
   {
      for (auto bIt = sbfs.begin(); bIt != sbfs.end(); ++bIt)
      {
         Block& b           = dynamic_cast< Block& >(*bIt);
         uint_t level       = b.getLevel();
         auto flagField     = b.getData< FlagField_T >(flagFieldID);
         uint8_t freeSlipFlag = flagField->registerFlag(freeSlipFlagUID_);
         uint8_t ubbFlag    = flagField->registerFlag(ubbFlagUID_);
         uint8_t outflowFlag    = flagField->registerFlag(outflowFlagUID_);
         for (auto cIt = flagField->beginWithGhostLayerXYZ(2); cIt != flagField->end(); ++cIt)
         {
            Cell localCell = cIt.cell();
            Cell globalCell(localCell);
            sbfs.transformBlockLocalToGlobalCell(globalCell, b);
            if (globalCell.x() < 0) { flagField->addFlag(localCell, ubbFlag); }
            else if (globalCell.x() >= cell_idx_c(sbfs.getNumberOfXCells(level))) { flagField->addFlag(localCell, outflowFlag); }
            else if (globalCell.y() >= cell_idx_c(sbfs.getNumberOfYCells(level))) { flagField->addFlag(localCell, freeSlipFlag); }
            else if (globalCell.y() < cell_idx_c(1 << level)) {flagField->addFlag(localCell, freeSlipFlag);}
         }
      }
   }

 private:
   const std::string refinementProfile_;
   const uint_t refinementDepth_;

   const FlagUID freeSlipFlagUID_;
   const FlagUID ubbFlagUID_;
   const FlagUID outflowFlagUID_;
};

static void createSetupBlockForest(SetupBlockForest& setupBfs, const Config::BlockHandle& domainSetup, Channel& setup)
{
   Vector3< real_t > domainSize = domainSetup.getParameter< Vector3< real_t > >("domainSize");
   Vector3< uint_t > rootBlocks = domainSetup.getParameter< Vector3< uint_t > >("rootBlocks");
   Vector3< bool > periodic     = domainSetup.getParameter< Vector3< bool > >("periodic");

   auto refSelection = setup.refinementSelector();
   setupBfs.addRefinementSelectionFunction(std::function< void(SetupBlockForest&) >(refSelection));
   AABB domain(real_t(0.0), real_t(0.0), real_t(0.0), domainSize[0], domainSize[1], domainSize[2]);
   setupBfs.init(domain, rootBlocks[0], rootBlocks[1], rootBlocks[2], periodic[0], periodic[1], periodic[2]);
   setupBfs.balanceLoad(blockforest::StaticLevelwiseCurveBalance(true), uint_c(MPIManager::instance()->numProcesses()));
}

int main(int argc, char** argv)
{
   walberla::Environment walberlaEnv(argc, argv);
   mpi::MPIManager::instance()->useWorldComm();

   logging::configureLogging(walberlaEnv.config());

   // read parameters
   auto domainSetup = walberlaEnv.config()->getOneBlock("DomainSetup");
   auto parameters  = walberlaEnv.config()->getOneBlock("Parameters");

   const real_t omega           = parameters.getParameter< real_t >("omega");
   const real_t inletVelocity   = parameters.getParameter< real_t >("inletVelocity");
   const uint_t timesteps       = parameters.getParameter< uint_t >("timesteps", uint_c(10)) + uint_c(1);
   const uint_t refinementDepth = parameters.getParameter< uint_t >("refinementDepth", uint_c(1));

   auto loggingParameters = walberlaEnv.config()->getOneBlock("Logging");
   bool writeSetupForestAndReturn = loggingParameters.getParameter<bool>("writeSetupForestAndReturn", false);

   auto remainingTimeLoggerFrequency =
      parameters.getParameter< real_t >("remainingTimeLoggerFrequency", real_c(3.0)); // in seconds

   auto flowSetup = std::make_shared< Channel >(refinementDepth);

   SetupBlockForest setupBfs;
   WALBERLA_LOG_INFO_ON_ROOT("Generating SetupBlockForest...")
   createSetupBlockForest(setupBfs, domainSetup, *flowSetup);

   // Create structured block forest
   Vector3< uint_t > cellsPerBlock = domainSetup.getParameter< Vector3< uint_t > >("cellsPerBlock");

   WALBERLA_LOG_INFO_ON_ROOT("Creating structured block forest...")
   auto bfs    = std::make_shared< BlockForest >(uint_c(MPIManager::instance()->worldRank()), setupBfs);
   auto blocks = std::make_shared< StructuredBlockForest >(bfs, cellsPerBlock[0], cellsPerBlock[1], cellsPerBlock[2]);
   blocks->createCellBoundingBoxes();

   if (writeSetupForestAndReturn)
   {
      WALBERLA_LOG_INFO_ON_ROOT("Writing SetupBlockForest to VTK file")
      WALBERLA_ROOT_SECTION() { vtk::writeDomainDecomposition(blocks, "FreeSlipRefinementDomainDecomposition", "vtk_out"); }

      WALBERLA_LOG_INFO_ON_ROOT("Blocks created: " << setupBfs.getNumberOfBlocks())
      for (uint_t level = 0; level <= refinementDepth; level++)
      {
         WALBERLA_LOG_INFO_ON_ROOT("Level " << level << " Blocks: " << setupBfs.getNumberOfBlocks(level))
      }
      WALBERLA_LOG_INFO_ON_ROOT("Ending program")
      return EXIT_SUCCESS;
   }

   WALBERLA_LOG_INFO_ON_ROOT("Blocks created: " << setupBfs.getNumberOfBlocks())
   for (uint_t level = 0; level <= refinementDepth; level++)
   {
      WALBERLA_LOG_INFO_ON_ROOT("Level " << level << " Blocks: " << setupBfs.getNumberOfBlocks(level))
   }

   StorageSpecification_T StorageSpec = StorageSpecification_T();
   BlockDataID pdfFieldId = lbm_generated::addPdfFieldToStorage(blocks, "pdf field", StorageSpec, uint_c(2));
   BlockDataID velFieldId = field::addToStorage< VectorField_T >(blocks, "Velocity", real_c(0.0), field::fzyx);

   BlockDataID flagFieldId = field::addFlagFieldToStorage< FlagField_T >(blocks, "flag field", uint_c(3));

   SweepCollection_T sweepCollection(blocks, pdfFieldId, velFieldId, omega);
   for (auto& block : *blocks)
   {
      sweepCollection.initialise(&block);
   }

   const FlagUID fluidFlagUID("Fluid");
   flowSetup->setupBoundaryFlagField(*blocks, flagFieldId);
   geometry::setNonBoundaryCellsToDomain< FlagField_T >(*blocks, flagFieldId, fluidFlagUID, 2);
   BoundaryCollection_T boundaryCollection(blocks, flagFieldId, pdfFieldId, fluidFlagUID, inletVelocity);

   WALBERLA_LOG_INFO_ON_ROOT("Setting up communication...")
   auto comm =
      std::make_shared< blockforest::communication::NonUniformBufferedScheme< CommunicationStencil_T > >(blocks);
   auto packInfo = lbm_generated::setupNonuniformPdfCommunication< PdfField_T >(blocks, pdfFieldId);
   comm->addPackInfo(packInfo);

   lbm_generated::BasicRecursiveTimeStep< PdfField_T, SweepCollection_T, BoundaryCollection_T > timestep(
      blocks, pdfFieldId, sweepCollection, boundaryCollection, comm, packInfo);

   SweepTimeloop timeloop(blocks->getBlockStorage(), timesteps);
   uint_t vtkWriteFrequency = parameters.getParameter< uint_t >("vtkWriteFrequency", 0);
   if (vtkWriteFrequency > 0)
   {
      auto vtkOutput = vtk::createVTKOutput_BlockData(*blocks, "FreeSlipRefinementVTK", vtkWriteFrequency, 0, false, "vtk_out",
                                                      "simulation_step", false, true, true, false, 0);

      auto velWriter = make_shared< field::VTKWriter< VectorField_T > >(velFieldId, "velocity");
      vtkOutput->addBeforeFunction([&]() {
         for (auto& block : *blocks)
         {
            sweepCollection.calculateMacroscopicParameters(&block);
         }
      });

      vtkOutput->addCellDataWriter(velWriter);
      timeloop.addFuncBeforeTimeStep(vtk::writeFiles(vtkOutput), "VTK Output");
   }
   timeloop.addFuncAfterTimeStep(timestep);

   // log remaining time
   if (remainingTimeLoggerFrequency > 1.0)
   {
      timeloop.addFuncAfterTimeStep(timing::RemainingTimeLogger(timeloop.getNrOfTimeSteps(), remainingTimeLoggerFrequency),
                                    "remaining time logger");
   }

   WALBERLA_LOG_INFO_ON_ROOT("Starting Simulation with " << timesteps << " timesteps")

   timeloop.run();


   for (auto& block : *blocks)
   {
      sweepCollection.calculateMacroscopicParameters(&block);
      Block& b           = dynamic_cast< Block& >(block);
      uint_t level       = b.getLevel();

      auto velField  = b.getData< VectorField_T >(velFieldId);
      for( auto it = velField->beginXYZ(); it != velField->end(); ++it )
      {
         Cell localCell = it.cell();
         Cell globalCell(localCell);
         blocks->transformBlockLocalToGlobalCell(globalCell, b);

         if (globalCell.y() >= (cell_idx_c(1 << level)))
         {
            WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(it.getF(0), inletVelocity, real_c(1e-5));
         }
      }
   }
   return EXIT_SUCCESS;
}
