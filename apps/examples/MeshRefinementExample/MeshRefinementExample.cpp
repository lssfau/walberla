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
//! \file MeshRefinementExample.cpp
//! \author Philipp Suffa   philipp.suffa@fau.de
//
//======================================================================================================================

#include "blockforest/all.h"
#include "core/all.h"
#include "domain_decomposition/all.h"
#include "field/all.h"
#include "geometry/all.h"
#include "lbm/all.h"
#include "lbm_generated/field/PdfField.h"
#include "lbm_generated/field/AddToStorage.h"
#include "lbm_generated/evaluation/PerformanceEvaluation.h"
#include "lbm_generated/communication/NonuniformGeneratedPdfPackInfo.h"
#include "lbm_generated/refinement/BasicRecursiveTimeStep.h"
#include "timeloop/all.h"

#include "InfoHeader.h"

namespace walberla
{

constexpr uint_t FieldGhostLayer{2};

using StorageSpecification_T     = lbm::MeshRefinementExampleStorageSpecification;
using LBMCommunicationStencil_T  = StorageSpecification_T::CommunicationStencil;
using PdfField_T                 = lbm_generated::PdfField< StorageSpecification_T >;

using SweepCollection_T = lbm::MeshRefinementExampleSweepCollection;

using flag_t      = walberla::uint32_t;
using FlagField_T = FlagField< flag_t >;

using BoundaryCollection_T = lbm::MeshRefinementExampleBoundaryCollection< FlagField_T >;
const FlagUID FluidFlagUID("Fluid");


class AABBRefinement
{
 public:
   AABBRefinement(AABB objectAABB, const uint_t depth) : objectAABB_(objectAABB), refinementDepth_(depth){};

   void operator()(SetupBlockForest& forest) const
   {
      auto extendedAABB = AABB(objectAABB_.minCorner()[0] - objectAABB_.xSize() * 0.00, objectAABB_.minCorner()[1] - objectAABB_.ySize() * 0.0, objectAABB_.minCorner()[2] - objectAABB_.zSize() * 0.0,
                              objectAABB_.maxCorner()[0] + 0.3 * objectAABB_.xSize(), objectAABB_.maxCorner()[1] + objectAABB_.ySize() * 0.0, objectAABB_.maxCorner()[2] + objectAABB_.zSize() * 0.0);
      for(auto & block : forest) {
         auto blockAABB = block.getAABB();
         if (extendedAABB.intersects(blockAABB)) {
            if( block.getLevel() < refinementDepth_)
               block.setMarker( true );
         }
      }
   }
 private:
   const AABB objectAABB_;
   const uint_t refinementDepth_;
};


int main(int argc, char** argv)
{
   Environment env(argc, argv);
   if (!env.config()) { WALBERLA_ABORT("No configuration file specified!"); }
   mpi::MPIManager::instance()->useWorldComm();
   WALBERLA_MPI_WORLD_BARRIER()

   // Get parameters from config file
   auto parameters = env.config()->getOneBlock("Parameters");
   const uint_t timesteps = parameters.getParameter< uint_t >("timesteps");
   const Vector3<uint_t> rootBlocks = parameters.getParameter< Vector3<uint_t> >("rootBlocks");
   Vector3<uint_t> cellsPerBlock = parameters.getParameter< Vector3< uint_t > >("cellsPerBlock");
   Vector3<real_t> sphereCenter = parameters.getParameter< Vector3< real_t > >("sphereCenter");

   const real_t sphereRadius = parameters.getParameter< real_t >("sphereRadius");
   const uint_t refinementLevels = parameters.getParameter< uint_t >("refinementLevels");
   const uint_t vtkWriteFrequency = parameters.getParameter< uint_t >("vtkWriteFrequency");

   const real_t omega = parameters.getParameter< real_t >("omega");


   auto aabb = AABB(0,0,0,1.0, 1.0, 1.0);
   geometry::Sphere sphere(sphereCenter, sphereRadius);


   // Create block forest
   const uint_t numProcs = uint_t(mpi::MPIManager::instance()->numProcesses());
   SetupBlockForest setupBfs;
   AABBRefinement refinement( sphere.boundingBox(), refinementLevels );
   setupBfs.addRefinementSelectionFunction(std::function<void(SetupBlockForest &)>(refinement));

   setupBfs.addWorkloadMemorySUIDAssignmentFunction(blockforest::uniformWorkloadAndMemoryAssignment);
   setupBfs.init(aabb, rootBlocks[0], rootBlocks[1], rootBlocks[2], false, false, false);
   setupBfs.balanceLoad(blockforest::StaticLevelwiseCurveBalanceWeighted(), numProcs);

   shared_ptr< BlockForest > forest = std::make_shared< BlockForest >(uint_c(MPIManager::instance()->worldRank()), setupBfs);
   auto blocks = std::make_shared< StructuredBlockForest >(forest, cellsPerBlock[0], cellsPerBlock[1], cellsPerBlock[2]);
   blocks->createCellBoundingBoxes();

   // Data fields
   const StorageSpecification_T StorageSpec = StorageSpecification_T();
   const BlockDataID pdfFieldCpuId      = lbm_generated::addPdfFieldToStorage(blocks, "pdfs", StorageSpec, FieldGhostLayer, field::fzyx);
   const BlockDataID velocityFieldCpuId = field::addToStorage< VelocityField_T >(blocks, "velocity", real_c(0.0), field::fzyx, FieldGhostLayer);
   const BlockDataID densityFieldCpuId  = field::addToStorage< ScalarField_T   >(blocks, "density", real_c(1.0), field::fzyx, FieldGhostLayer);
   const BlockDataID flagFieldId        = field::addFlagFieldToStorage< FlagField_T >(blocks, "flag field", FieldGhostLayer);

   SweepCollection_T sweepCollection( blocks, densityFieldCpuId, pdfFieldCpuId, velocityFieldCpuId, omega);

   // Initialize Velocity and PDF field
   Vector3<real_t> initialVel(0.1, 0, 0);
   for (auto& block : *blocks)
   {
      auto velField = block.getData<VelocityField_T>(velocityFieldCpuId);
      WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ(velField,
                                                       for (uint_t d = 0; d < LBMCommunicationStencil_T::D; ++d) {
                                                         velField->get(x,y,z,d) =  initialVel[d];
                                                      }
      )
      sweepCollection.initialise(&block, 0);
   }

   // Boundary handling
   auto boundariesConfig = env.config()->getOneBlock("Boundaries");
   geometry::initBoundaryHandling< FlagField_T >(*blocks, flagFieldId, boundariesConfig);
   FlagUID objBoundaryUID("NoSlip");
   for (auto& block : *blocks)
   {
      auto * flagField = block.getData<FlagField_T>(flagFieldId);
      if ( !flagField->flagExists(objBoundaryUID))
         flagField->registerFlag(objBoundaryUID);
      auto flag = flagField->getFlag(objBoundaryUID);

      WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ(flagField,
                                                      Cell cell(x,y,z);
                                                      auto midPoint = blocks->getGlobalCellCenterFromBlockLocalCell(cell, block);
                                                      if (contains(sphere, midPoint)) {
                                                         flagField->addFlag(cell, flag);
                                                      }
      )
   }

   geometry::setNonBoundaryCellsToDomain< FlagField_T >(*blocks, flagFieldId, FluidFlagUID);
   WALBERLA_MPI_WORLD_BARRIER();
   BoundaryCollection_T boundaryCollection( blocks, flagFieldId, pdfFieldCpuId, FluidFlagUID, initialVel[0] );

   lbm::BlockForestEvaluation<FlagField_T>( blocks, flagFieldId, FluidFlagUID ).logInfoOnRoot();

   // Communication
   auto communication = std::make_shared< blockforest::communication::NonUniformBufferedScheme< LBMCommunicationStencil_T > >(blocks);
   auto packInfo      = lbm_generated::setupNonuniformPdfCommunication< PdfField_T >(blocks, pdfFieldCpuId);
   communication->addPackInfo(packInfo);

   // Timeloop
   SweepTimeloop timeloop(blocks->getBlockStorage(), timesteps);
   lbm_generated::BasicRecursiveTimeStep< PdfField_T, SweepCollection_T, BoundaryCollection_T > LBMMeshRefinement(
      blocks, pdfFieldCpuId, sweepCollection, boundaryCollection, communication, packInfo);
   LBMMeshRefinement.addRefinementToTimeLoop(timeloop, 0);

   timeloop.addFuncAfterTimeStep(timing::RemainingTimeLogger(timeloop.getNrOfTimeSteps(), 10), "remaining time logger");

   // VTK output
   if (vtkWriteFrequency > 0)
   {
      auto vtkOutput = vtk::createVTKOutput_BlockData(*blocks, "vtk", vtkWriteFrequency, 0);
      auto velocityWriter = make_shared< field::VTKWriter< VelocityField_T > >(velocityFieldCpuId, "velocity");
      auto densityWriter  = make_shared< field::VTKWriter< ScalarField_T > >(densityFieldCpuId, "density");

      vtkOutput->addCellDataWriter(velocityWriter);
      vtkOutput->addCellDataWriter(densityWriter);

      field::FlagFieldCellFilter< FlagField_T > fluidFilter(flagFieldId);
      fluidFilter.addFlag(FluidFlagUID);
      vtkOutput->addCellInclusionFilter(fluidFilter);

      timeloop.addFuncAfterTimeStep(vtk::writeFiles(vtkOutput), "VTK Output");
      vtk::writeDomainDecomposition(blocks, "domain_decomposition", "vtk_out", "write_call", true, true, 0);
   }

   // Run Simulation
   WcTimer simTimer;
   WcTimingPool timingPool;

   simTimer.start();
   timeloop.run(timingPool);
   simTimer.end();
   double time = simTimer.max();

   // Performance metrics
   lbm_generated::PerformanceEvaluation<FlagField_T> const performance(blocks, flagFieldId, FluidFlagUID);
   performance.logResultOnRoot(timesteps, time);
   timingPool.unifyRegisteredTimersAcrossProcesses();
   timingPool.logResultOnRoot( timing::REDUCE_TOTAL, true );

   return EXIT_SUCCESS;
}
} // namespace walberla

int main(int argc, char** argv) { walberla::main(argc, argv); }