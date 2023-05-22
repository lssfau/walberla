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
//! \file Example.cpp
//! \author Markus Holzer <markus.holzer@fau.de>
//
//======================================================================================================================

#include "blockforest/all.h"

#include "core/all.h"

#include "domain_decomposition/all.h"

#include "field/all.h"

#include "geometry/all.h"

#include "timeloop/all.h"

#include "lbm_generated/communication/NonuniformGeneratedPdfPackInfo.h"
#include "lbm_generated/communication/UniformGeneratedPdfPackInfo.h"
#include "lbm_generated/field/AddToStorage.h"
#include "lbm_generated/field/PdfField.h"
#include "lbm_generated/refinement/BasicRecursiveTimeStep.h"

// include the generated header file. It includes all generated classes
#include "Example_InfoHeader.h"

using namespace walberla;
using namespace std::placeholders;

using StorageSpecification_T = lbm::LBMStorageSpecification;
using Stencil_T              = StorageSpecification_T::Stencil;
using CommunicationStencil_T = StorageSpecification_T::CommunicationStencil;
using PdfField_T             = lbm_generated::PdfField< StorageSpecification_T >;
using PackInfo_T             = lbm_generated::UniformGeneratedPdfPackInfo< PdfField_T >;

using SweepCollection_T = lbm::LBMSweepCollection;

using VectorField_T = GhostLayerField< real_t, StorageSpecification_T::Stencil::D >;
using ScalarField_T = GhostLayerField< real_t, 1 >;

using flag_t               = walberla::uint8_t;
using FlagField_T          = FlagField< flag_t >;
using BoundaryCollection_T = lbm::LBMBoundaryCollection< FlagField_T >;

using RefinementSelectionFunctor = SetupBlockForest::RefinementSelectionFunction;

class LDCRefinement
{
 public:
   LDCRefinement(const uint_t depth) : refinementDepth_(depth){};

   void operator()(SetupBlockForest& forest)
   {
      std::vector< SetupBlock* > blocks;
      forest.getBlocks(blocks);

      for (auto b : blocks)
      {
         if (forest.atDomainZMaxBorder(*b))
         {
            if (b->getLevel() < refinementDepth_) { b->setMarker(true); }
         }
      }
   }

 private:
   const uint_t refinementDepth_;
};

class LDC
{
 public:
   LDC(const uint_t depth) : refinementDepth_(depth), noSlipFlagUID_("NoSlip"), ubbFlagUID_("UBB"){};

   Vector3< real_t > acceleration() const { return Vector3< real_t >(0.0); }

   RefinementSelectionFunctor refinementSelector() { return LDCRefinement(refinementDepth_); }

   void setupBoundaryFlagField(StructuredBlockForest& sbfs, const BlockDataID flagFieldID)
   {
      for (auto bIt = sbfs.begin(); bIt != sbfs.end(); ++bIt)
      {
         Block& b           = dynamic_cast< Block& >(*bIt);
         uint_t level       = b.getLevel();
         auto flagField     = b.getData< FlagField_T >(flagFieldID);
         uint8_t noslipFlag = flagField->registerFlag(noSlipFlagUID_);
         uint8_t ubbFlag    = flagField->registerFlag(ubbFlagUID_);
         for (auto cIt = flagField->beginWithGhostLayerXYZ(2); cIt != flagField->end(); ++cIt)
         {
            Cell localCell = cIt.cell();
            Cell globalCell(localCell);
            sbfs.transformBlockLocalToGlobalCell(globalCell, b);
            if (globalCell.z() >= cell_idx_c(sbfs.getNumberOfZCells(level))) { flagField->addFlag(localCell, ubbFlag); }
            else if (globalCell.z() < 0 || globalCell.x() < 0 ||
                     globalCell.x() >= cell_idx_c(sbfs.getNumberOfXCells(level)))
            {
               flagField->addFlag(localCell, noslipFlag);
            }
         }
      }
   }

 private:
   const std::string refinementProfile_;
   const uint_t refinementDepth_;

   const FlagUID noSlipFlagUID_;
   const FlagUID ubbFlagUID_;
};

static void createSetupBlockForest(SetupBlockForest& setupBfs, const Config::BlockHandle& domainSetup, LDC& setup)
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

   // read parameters
   auto domainSetup = walberlaEnv.config()->getOneBlock("DomainSetup");
   auto parameters  = walberlaEnv.config()->getOneBlock("Parameters");

   auto omega           = parameters.getParameter< real_t >("omega", real_c(1.4));
   auto timesteps       = parameters.getParameter< uint_t >("timesteps", uint_c(10)) + uint_c(1);
   auto refinementDepth = parameters.getParameter< uint_t >("refinementDepth", uint_c(1));

   auto remainingTimeLoggerFrequency =
      parameters.getParameter< real_t >("remainingTimeLoggerFrequency", real_c(3.0)); // in seconds

   auto flowSetup = std::make_shared< LDC >(refinementDepth);

   SetupBlockForest setupBfs;
   WALBERLA_LOG_INFO_ON_ROOT("Generating SetupBlockForest...")
   createSetupBlockForest(setupBfs, domainSetup, *flowSetup);
   // domainSetup

   // Create structured block forest
   Vector3< uint_t > cellsPerBlock = domainSetup.getParameter< Vector3< uint_t > >("cellsPerBlock");

   WALBERLA_LOG_INFO_ON_ROOT("Creating structured block forest...")
   auto bfs    = std::make_shared< BlockForest >(uint_c(MPIManager::instance()->worldRank()), setupBfs);
   auto blocks = std::make_shared< StructuredBlockForest >(bfs, cellsPerBlock[0], cellsPerBlock[1], cellsPerBlock[2]);
   blocks->createCellBoundingBoxes();

   WALBERLA_ROOT_SECTION() { vtk::writeDomainDecomposition(blocks, "domainDecomposition", "vtk_out"); }

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
   BoundaryCollection_T boundaryCollection(blocks, flagFieldId, pdfFieldId, fluidFlagUID);

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
      auto vtkOutput = vtk::createVTKOutput_BlockData(*blocks, "ExampleVTK", vtkWriteFrequency, 0, false, "vtk_out",
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
   timeloop.addFuncAfterTimeStep(timing::RemainingTimeLogger(timeloop.getNrOfTimeSteps(), remainingTimeLoggerFrequency),
                                 "remaining time logger");

   WALBERLA_LOG_INFO_ON_ROOT("Starting Simulation with " << timesteps << " timesteps")

   timeloop.run();

   return EXIT_SUCCESS;
}
