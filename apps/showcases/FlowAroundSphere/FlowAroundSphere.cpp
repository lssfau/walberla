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
//! \file FlowAroundSphere.cpp
//! \author Markus Holzer <markus.holzer@fau.de>
//
//======================================================================================================================

#include "blockforest/StructuredBlockForest.h"
#include "blockforest/communication/NonUniformBufferedScheme.h"

#include "core/Abort.h"
#include "core/DataTypes.h"
#include "core/SharedFunctor.h"
#include "core/debug/CheckFunctions.h"
#include "core/logging/Logging.h"
#include "core/math/Vector3.h"
#include "core/Environment.h"
#include "core/mpi/MPIManager.h"
#include "core/MemoryUsage.h"
#include "core/mpi/Reduce.h"
#include "core/timing/RemainingTimeLogger.h"
#include "core/timing/TimingPool.h"

#include "field/AddToStorage.h"
#include "field/CellCounter.h"
#include "field/FlagField.h"
#include "field/StabilityChecker.h"
#include "field/adaptors/AdaptorCreators.h"
#include "field/iterators/FieldIterator.h"
#include "field/vtk/VTKWriter.h"
#include "field/vtk/FlagFieldCellFilter.h"

#include "geometry/InitBoundaryHandling.h"

#if defined(WALBERLA_BUILD_WITH_GPU_SUPPORT)
#   include "gpu/AddGPUFieldToStorage.h"
#   include "gpu/DeviceSelectMPI.h"
#   include "gpu/ErrorChecking.h"
#   include "gpu/FieldCopy.h"
#   include "gpu/HostFieldAllocator.h"
#   include "gpu/ParallelStreams.h"
#   include "gpu/communication/NonUniformGPUScheme.h"
#   include "gpu/communication/UniformGPUScheme.h"
#endif

#include "lbm_generated/communication/NonuniformGeneratedPdfPackInfo.h"
#include "lbm_generated/communication/UniformGeneratedPdfPackInfo.h"
#include "lbm_generated/evaluation/PerformanceEvaluation.h"
#include "lbm_generated/field/AddToStorage.h"
#include "lbm_generated/field/PdfField.h"
#include "lbm_generated/refinement/BasicRecursiveTimeStep.h"
#include "lbm_generated/refinement/RefinementScaling.h"

#if defined(WALBERLA_BUILD_WITH_GPU_SUPPORT)
#   include "lbm_generated/gpu/AddToStorage.h"
#   include "lbm_generated/gpu/BasicRecursiveTimeStepGPU.h"
#   include "lbm_generated/gpu/GPUPdfField.h"
#   include "lbm_generated/gpu/NonuniformGeneratedGPUPdfPackInfo.h"
#   include "lbm_generated/gpu/UniformGeneratedGPUPdfPackInfo.h"
#endif

#include "timeloop/SweepTimeloop.h"
#include "vtk/VTKOutput.h"

#include <cstdlib>
#include <iostream>
#include <memory>

#include "Evaluation.h"
#include "GridGeneration.h"
#include "Setup.h"
#include "Sphere.h"
#include "Types.h"

#include "FlowAroundSphereStaticDefines.h"

using namespace walberla;

using BoundaryCollection_T = lbm::FlowAroundSphereBoundaryCollection< FlagField_T >;
using SweepCollection_T = lbm::FlowAroundSphereSweepCollection;

#if defined(WALBERLA_BUILD_WITH_GPU_SUPPORT)
using GPUPdfField_T = lbm_generated::GPUPdfField< StorageSpecification_T >;
using gpu::communication::NonUniformGPUScheme;
using gpu::communication::UniformGPUScheme;

using lbm_generated::NonuniformGeneratedGPUPdfPackInfo;
using lbm_generated::UniformGeneratedGPUPdfPackInfo;
#else
using PdfField_T = lbm_generated::PdfField< StorageSpecification_T >;
using blockforest::communication::NonUniformBufferedScheme;
using blockforest::communication::UniformBufferedScheme;

using lbm_generated::NonuniformGeneratedPdfPackInfo;
using lbm_generated::UniformGeneratedPdfPackInfo;
#endif

int main(int argc, char** argv)
{
   Environment env( argc, argv );
   mpi::MPIManager::instance()->useWorldComm();
#if defined(WALBERLA_BUILD_WITH_GPU_SUPPORT)
   gpu::selectDeviceBasedOnMpiRank();
#endif

   auto config = env.config();

   ///////////////////////
   /// PARAMETER INPUT ///
   ///////////////////////

   // read general simulation parameters
   auto parameters       = config->getOneBlock("Parameters");
   auto domainParameters = config->getOneBlock("DomainSetup");
   auto boundariesConfig = config->getBlock("Boundaries");
   auto loggingParameters= config->getOneBlock("Logging");
   Setup setup(parameters, domainParameters, loggingParameters, boundariesConfig, infoMap);

   bool writeSetupForestAndReturn              = loggingParameters.getParameter< bool >("writeSetupForestAndReturn", false);
   if (uint_c(MPIManager::instance()->numProcesses()) > 1) writeSetupForestAndReturn = false;

   WALBERLA_LOG_INFO_ON_ROOT("Diameter of the Sphere is resolved with " << setup.resolutionSphere << " lattice cells.")
   Sphere Sphere( setup );

   ///////////////////////////
   /// CREATE BLOCK FOREST ///
   ///////////////////////////

   shared_ptr< BlockForest > bfs;
   createBlockForest(bfs, setup);

   if (writeSetupForestAndReturn && mpi::MPIManager::instance()->numProcesses() == 1)
   {
      WALBERLA_LOG_INFO_ON_ROOT("BlockForest has been created and writen to file. Returning program")

      const uint_t nBlocks = bfs->getNumberOfBlocks();
      const uint_t nCells = nBlocks * (setup.cellsPerBlock[0] * setup.cellsPerBlock[1] * setup.cellsPerBlock[2]);
      const memory_t totalMemory = memory_t(nCells) * setup.memoryPerCell;

      WALBERLA_LOG_INFO_ON_ROOT( "Benchmark run data:"
                                "\n- simulation parameters:"   << std::setprecision(16) <<
                                "\n   + collision model:     " << infoMap["collisionOperator"] <<
                                "\n   + stencil:             " << infoMap["stencil"] <<
                                "\n   + streaming:           " << infoMap["streamingPattern"] <<
                                "\n   + compressible:        " << ( StorageSpecification_T::compressible ? "yes" : "no" ) <<
                                "\n   + mesh levels:         " << setup.refinementLevels + uint_c(1) <<
                                "\n   + resolution:          " << setup.dxC << " - on the coarsest grid" <<
                                "\n   + resolution:          " << setup.dxF << " - on the finest grid" <<
                                "\n- domain properties:      "
                                "\n   + # blocks:            " << nBlocks <<
                                "\n   + # cells:             " << uint_c(real_c(nCells) / real_c(1e6)) << " [1e6]"
                                "\n   + # cells sphere D:    " << setup.resolutionSphere <<
                                "\n   + total memory:        " << totalMemory / 1e9 << " [GB]" <<
                                "\n- simulation properties:  "
                                "\n   + sphere pos.(x):      " << setup.sphereXPosition * setup.dxC << " [m]" <<
                                "\n   + sphere pos.(y):      " << setup.sphereYPosition * setup.dxC << " [m]" <<
                                "\n   + sphere pos.(z):      " << setup.sphereZPosition * setup.dxC << " [m]" <<
                                "\n   + sphere radius:       " << setup.sphereRadius * setup.dxC << " [m]" <<
                                "\n   + kin. viscosity:      " << setup.viscosity * setup.dxC * setup.dxC / setup.dt << " [m^2/s] (on the coarsest grid)" <<
                                "\n   + viscosity LB:        " << setup.viscosity  << " [dx*dx/dt] (on the coarsest grid)" <<
                                "\n   + omega:               " << setup.omega << " (on the coarsest grid)" <<
                                "\n   + rho:                 " << setup.rho << " [kg/m^3]" <<
                                "\n   + inflow velocity:     " << setup.referenceVelocity << " [m/s]" <<
                                "\n   + lattice velocity:    " << setup.inflowVelocity << " [dx/dt]" <<
                                "\n   + Reynolds number:     " << setup.reynoldsNumber <<
                                "\n   + Mach number:         " << setup.machNumber <<
                                "\n   + dx (coarsest grid):  " << setup.dxC << " [m]" <<
                                "\n   + dt (coarsest grid):  " << setup.dt << " [s]"
                                "\n   + #time steps:         " << setup.timesteps << " (on the coarsest grid, " << ( real_c(1.0) / setup.dt ) << " for 1s of real time)" <<
                                "\n   + simulation time:     " << ( real_c( setup.timesteps ) * setup.dt ) << " [s]" )

      logging::Logging::printFooterOnStream();
      return EXIT_SUCCESS;
   }

   auto blocks = std::make_shared< StructuredBlockForest >(bfs, setup.cellsPerBlock[0], setup.cellsPerBlock[1], setup.cellsPerBlock[2]);
   blocks->createCellBoundingBoxes();

   ////////////////////////////////////
   /// CREATE AND INITIALIZE FIELDS ///
   ////////////////////////////////////

   // create fields
   const StorageSpecification_T StorageSpec = StorageSpecification_T();
   IDs ids;

#if defined(WALBERLA_BUILD_WITH_GPU_SUPPORT)
   auto PDFAllocator = make_shared< gpu::HostFieldAllocator< PdfField_T::value_type > >();
   auto allocator = make_shared< gpu::HostFieldAllocator< VelocityField_T::value_type > >();
   ids.pdfField = lbm_generated::addPdfFieldToStorage(blocks, "pdfs", StorageSpec, uint_c(2), field::fzyx, PDFAllocator);
   ids.velocityField = field::addToStorage< VelocityField_T >(blocks, "velocity", real_c(0.0), field::fzyx, uint_c(2), allocator);
   ids.densityField = field::addToStorage< ScalarField_T >(blocks, "density", real_c(1.0), field::fzyx, uint_c(2), allocator);
   ids.omegaField = field::addToStorage< ScalarField_T >(blocks, "omega", real_c(0.0), field::fzyx, uint_c(2), allocator);
   ids.flagField = field::addFlagFieldToStorage< FlagField_T >(blocks, "Boundary Flag Field", uint_c(3));

   ids.pdfFieldGPU = lbm_generated::addGPUPdfFieldToStorage< PdfField_T >(blocks, ids.pdfField, StorageSpec, "pdfs on GPU", true);
   ids.velocityFieldGPU = gpu::addGPUFieldToStorage< VelocityField_T >(blocks, ids.velocityField, "velocity on GPU", true);
   ids.densityFieldGPU = gpu::addGPUFieldToStorage< ScalarField_T >(blocks, ids.densityField, "density on GPU", true);
   ids.omegaFieldGPU = gpu::addGPUFieldToStorage< ScalarField_T >(blocks, ids.omegaField, "omega on GPU", true);

   WALBERLA_GPU_CHECK(gpuDeviceSynchronize())
   WALBERLA_GPU_CHECK(gpuPeekAtLastError())
#else
   ids.pdfField = lbm_generated::addPdfFieldToStorage(blocks, "pdfs", StorageSpec, uint_c(2), field::fzyx);
   ids.velocityField = field::addToStorage< VelocityField_T >(blocks, "vel", VelocityField_T::value_type(0.0), field::fzyx, uint_c(2));
   ids.densityField = field::addToStorage< ScalarField_T >(blocks, "density", ScalarField_T::value_type(1.0), field::fzyx, uint_c(2));
   ids.omegaField = field::addToStorage< ScalarField_T >(blocks, "omega", ScalarField_T::value_type(0.0), field::fzyx, uint_c(2));
   ids.flagField = field::addFlagFieldToStorage< FlagField_T >(blocks, "Boundary Flag Field", uint_c(3));
#endif

   auto spongeLayer     = config->getOneBlock("SpongeLayer");
   const bool deactivateSpongeLayer = spongeLayer.getParameter< bool >("deactivateSpongeLayer");
   const real_t targetOmega         = spongeLayer.getParameter< real_t >("targetOmega");
   const real_t spongeStart         = spongeLayer.getParameter< real_t >("spongeStart");

   const real_t omegaDiff   = setup.omega - targetOmega;
   const real_t spongeWidth = real_c(blocks->getDomain().xMax()) - spongeStart;
   const real_t spongeMid   = spongeStart + (spongeWidth / real_c(2.0));

   if(!deactivateSpongeLayer)
      WALBERLA_LOG_WARNING_ON_ROOT("Using a sponge layer at the Outlet. The sponge layer starts at " << spongeStart << " [m] and targets a relaxation rate of " << targetOmega << " at the outlet" )

   for (auto& block : *blocks)
   {
      Block& b = dynamic_cast< Block& >(block);
      const uint_t level = blocks->getLevel(b);

      auto * omegaField = b.getData< ScalarField_T > ( ids.omegaField );
      for( auto it = omegaField->beginWithGhostLayer(0); it != omegaField->end(); ++it ){
         if(deactivateSpongeLayer){
            omegaField->get(it.cell()) = ScalarField_T::value_type(lbm_generated::relaxationRateScaling(setup.omega, level));
         }
         else{
            Cell globalCell;
            blocks->transformBlockLocalToGlobalCell(globalCell, block, it.cell());
            Vector3<real_t> cellCenter = blocks->getCellCenter(globalCell, level);
            const real_t factor = real_c(0.5) + real_c(0.5) * std::tanh(real_c(2.0) * (cellCenter[0] - spongeMid) / spongeWidth);
            omegaField->get(it.cell()) = ScalarField_T::value_type(lbm_generated::relaxationRateScaling(setup.omega - (factor * omegaDiff), level));
         }
      }
   }

#if defined(WALBERLA_BUILD_WITH_GPU_SUPPORT)
   gpu::fieldCpy< gpu::GPUField< ScalarField_T::value_type >, ScalarField_T>(blocks, ids.omegaFieldGPU, ids.omegaField);
#endif

   WALBERLA_MPI_BARRIER()

   const Cell innerOuterSplit =
      Cell(parameters.getParameter< Vector3< cell_idx_t > >("innerOuterSplit", Vector3< cell_idx_t >(1, 1, 1)));
#if defined(WALBERLA_BUILD_WITH_GPU_SUPPORT)
   SweepCollection_T sweepCollection(blocks, ids.pdfFieldGPU, ids.omegaFieldGPU, ids.densityFieldGPU, ids.velocityFieldGPU, innerOuterSplit);
#else
   SweepCollection_T sweepCollection(blocks, ids.pdfField, ids.omegaField, ids.densityField, ids.velocityField, innerOuterSplit);
#endif

#if defined(WALBERLA_BUILD_WITH_GPU_SUPPORT)
   WALBERLA_LOG_INFO_ON_ROOT("Setting up communication...")
   const bool gpuEnabledMPI       = parameters.getParameter< bool >("gpuEnabledMPI", false);

   auto nonUniformCommunication = std::make_shared< NonUniformGPUScheme< CommunicationStencil_T > >(blocks, gpuEnabledMPI);
   auto nonUniformPackInfo      = lbm_generated::setupNonuniformGPUPdfCommunication< GPUPdfField_T >(blocks, ids.pdfFieldGPU);
   nonUniformCommunication->addPackInfo(nonUniformPackInfo);
#else
   WALBERLA_LOG_INFO_ON_ROOT("Setting up communication...")
   auto nonUniformCommunication = std::make_shared< NonUniformBufferedScheme< CommunicationStencil_T > >(blocks);
   auto nonUniformPackInfo      = lbm_generated::setupNonuniformPdfCommunication< PdfField_T >(blocks, ids.pdfField);
   nonUniformCommunication->addPackInfo(nonUniformPackInfo);

#endif
   WALBERLA_MPI_BARRIER()
   WALBERLA_LOG_INFO_ON_ROOT("Setting up communication done")

   /////////////////////////
   /// BOUNDARY HANDLING ///
   /////////////////////////
   WALBERLA_LOG_INFO_ON_ROOT("Start BOUNDARY HANDLING")
   // create and initialize boundary handling
   Sphere.setupBoundary(blocks, ids.flagField);
   geometry::setNonBoundaryCellsToDomain< FlagField_T >(*blocks, ids.flagField, setup.fluidUID, cell_idx_c(0));

   if(parameters.getParameter< bool >("initialiseWithInletVelocity"))
   {
      for (auto& block : *blocks)
      {
         auto * flagField = block.getData< FlagField_T > ( ids.flagField );
         auto * velField = block.getData< VelocityField_T > ( ids.velocityField );
         // auto domainFlag = flagField->getFlag(fluidFlagUID);

         for( auto it = flagField->beginWithGhostLayer(2); it != flagField->end(); ++it )
         {
            // if (!isFlagSet(it, domainFlag))
            //  continue;

            velField->get(it.cell(), 0) = VelocityField_T::value_type(setup.inflowVelocity);
         }
      }
   }

#if defined(WALBERLA_BUILD_WITH_GPU_SUPPORT)
   gpu::fieldCpy< gpu::GPUField< VelocityField_T::value_type >, VelocityField_T>(blocks, ids.velocityFieldGPU, ids.velocityField);
#endif

   for (auto& block : *blocks)
   {
      sweepCollection.initialise(&block, cell_idx_c(1));
   }

#if defined(WALBERLA_BUILD_WITH_GPU_SUPPORT)
   gpu::fieldCpy< PdfField_T, gpu::GPUField< PdfField_T::value_type > >(blocks, ids.pdfField, ids.pdfFieldGPU);
   sweepCollection.initialiseBlockPointer();
#endif

   std::function< VelocityField_T::value_type (const Cell&, const Cell&, const shared_ptr< StructuredBlockForest >&, IBlock&) >
      wallDistanceFunctor = wallDistance(Sphere);

   const real_t omegaFinestLevel = lbm_generated::relaxationRateScaling(setup.omega, setup.refinementLevels);

#if defined(WALBERLA_BUILD_WITH_GPU_SUPPORT)
   BoundaryCollection_T boundaryCollection(blocks, ids.flagField, ids.pdfFieldGPU, setup.fluidUID, omegaFinestLevel, setup.inflowVelocity, wallDistanceFunctor, ids.pdfField);
   WALBERLA_GPU_CHECK(gpuDeviceSynchronize())
   WALBERLA_GPU_CHECK(gpuPeekAtLastError())
#else
   BoundaryCollection_T boundaryCollection(blocks, ids.flagField, ids.pdfField, setup.fluidUID, omegaFinestLevel, setup.inflowVelocity, wallDistanceFunctor);
#endif
   WALBERLA_MPI_BARRIER()
   WALBERLA_LOG_INFO_ON_ROOT("BOUNDARY HANDLING done")

   //////////////////////////////////
   /// SET UP SWEEPS AND TIMELOOP ///
   //////////////////////////////////
   WALBERLA_LOG_INFO_ON_ROOT("Start SWEEPS AND TIMELOOP")
   // flow evaluation
   auto EvaluationParameters      = config->getOneBlock("Evaluation");
   const uint_t evaluationCheckFrequency       = EvaluationParameters.getParameter< uint_t >("evaluationCheckFrequency");
   const uint_t rampUpTime                     = EvaluationParameters.getParameter< uint_t >("rampUpTime");
   const bool evaluationLogToStream            = EvaluationParameters.getParameter< bool >("logToStream");
   const bool evaluationLogToFile              = EvaluationParameters.getParameter< bool >("logToFile");

   shared_ptr< Evaluation > evaluation( new Evaluation( blocks, evaluationCheckFrequency, rampUpTime, boundaryCollection,
                                                        ids, setup, evaluationLogToStream, evaluationLogToFile));

   // create time loop
   SweepTimeloop timeloop(blocks->getBlockStorage(), setup.timesteps);
#if defined(WALBERLA_BUILD_WITH_GPU_SUPPORT)
   int streamHighPriority = 0;
   int streamLowPriority  = 0;
   WALBERLA_GPU_CHECK(gpuDeviceGetStreamPriorityRange(&streamLowPriority, &streamHighPriority))
   sweepCollection.setOuterPriority(streamHighPriority);
   auto defaultStream = gpu::StreamRAII::newPriorityStream(streamLowPriority);

   lbm_generated::BasicRecursiveTimeStepGPU< GPUPdfField_T, SweepCollection_T, BoundaryCollection_T >
      LBMMeshRefinement(blocks, ids.pdfFieldGPU, sweepCollection, boundaryCollection, nonUniformCommunication, nonUniformPackInfo);

   LBMMeshRefinement.addRefinementToTimeLoop(timeloop, uint_c(0));
   LBMMeshRefinement.addPostBoundaryHandlingFunction(evaluation->forceCalculationFunctor());
#else
   std::shared_ptr< lbm_generated::BasicRecursiveTimeStep< PdfField_T, SweepCollection_T, BoundaryCollection_T > >
      LBMRefinement;

   LBMRefinement = std::make_shared<
   lbm_generated::BasicRecursiveTimeStep< PdfField_T, SweepCollection_T, BoundaryCollection_T > >(
   blocks, ids.pdfField, sweepCollection, boundaryCollection, nonUniformCommunication, nonUniformPackInfo);
   LBMRefinement->addPostBoundaryHandlingFunction(evaluation->forceCalculationFunctor());
   LBMRefinement->addRefinementToTimeLoop(timeloop);
#endif
   //////////////////
   /// VTK OUTPUT ///
   //////////////////
   WALBERLA_LOG_INFO_ON_ROOT("SWEEPS AND TIMELOOP done")

   auto VTKWriter           = config->getOneBlock("VTKWriter");
   const uint_t vtkWriteFrequency       = VTKWriter.getParameter< uint_t >("vtkWriteFrequency", 0);
   const uint_t vtkGhostLayers          = VTKWriter.getParameter< uint_t >("vtkGhostLayers", 0);
   const bool amrFileFormat             = VTKWriter.getParameter< bool >("amrFileFormat", false);
   const bool oneFilePerProcess         = VTKWriter.getParameter< bool >("oneFilePerProcess", false);
   const real_t samplingResolution      = VTKWriter.getParameter< real_t >("samplingResolution", real_c(-1.0));
   const uint_t initialWriteCallsToSkip = VTKWriter.getParameter< uint_t >("initialWriteCallsToSkip", uint_c(0.0));
   const real_t sliceThickness          = VTKWriter.getParameter< real_t >("sliceThickness", real_c(1.0));

   auto finalDomain = blocks->getDomain();
   if (vtkWriteFrequency > 0)
   {
      const std::string vtkFolder = "VTKSphereRE_" + std::to_string(uint64_c(setup.reynoldsNumber));
      auto vtkOutput =
         vtk::createVTKOutput_BlockData(*blocks, "vtk", vtkWriteFrequency, vtkGhostLayers, false, vtkFolder,
                                        "simulation_step", false, true, true, false, 0, amrFileFormat, oneFilePerProcess);

      vtkOutput->setInitialWriteCallsToSkip(initialWriteCallsToSkip);

      vtkOutput->addBeforeFunction([&]() {
         for (auto& block : *blocks)
         {
            sweepCollection.calculateMacroscopicParameters(&block);
         }

#if defined(WALBERLA_BUILD_WITH_GPU_SUPPORT)
         gpu::fieldCpy< VelocityField_T, gpu::GPUField< VelocityField_T::value_type > >(blocks, ids.velocityField, ids.velocityFieldGPU);
         gpu::fieldCpy< ScalarField_T, gpu::GPUField< ScalarField_T::value_type > >(blocks, ids.densityField, ids.densityFieldGPU);
#endif
      });

      vtkOutput->setSamplingResolution(samplingResolution);

      field::FlagFieldCellFilter<FlagField_T> fluidFilter( ids.flagField );
      fluidFilter.addFlag( setup.obstacleUID );
      vtkOutput->addCellExclusionFilter(fluidFilter);


      if (VTKWriter.getParameter< bool >("writeOnlySlice", true)){
         const AABB sliceXY(finalDomain.xMin(), finalDomain.yMin(), finalDomain.center()[2] - sliceThickness * setup.dxF,
                            finalDomain.xMax(), finalDomain.yMax(), finalDomain.center()[2] + sliceThickness * setup.dxF);
         vtkOutput->addCellInclusionFilter(vtk::AABBCellFilter(sliceXY));

         if (VTKWriter.getParameter< bool >("writeXZSlice", false)){
            const AABB sliceXZ(finalDomain.xMin(), finalDomain.center()[1] - setup.dxF, finalDomain.zMin(),
                               finalDomain.xMax(), finalDomain.center()[1] + setup.dxF, finalDomain.zMax());
            vtkOutput->addCellInclusionFilter(vtk::AABBCellFilter(sliceXZ));
         }
      }

      if (VTKWriter.getParameter< bool >("velocity"))
      {
         auto velWriter = make_shared< field::VTKWriter< VelocityField_T, float32 > >(ids.velocityField, "velocity");
         vtkOutput->addCellDataWriter(velWriter);
      }
      if (VTKWriter.getParameter< bool >("density"))
      {
         auto densityWriter = make_shared< field::VTKWriter< ScalarField_T, float32 > >(ids.densityField, "density");
         vtkOutput->addCellDataWriter(densityWriter);
      }
      if (VTKWriter.getParameter< bool >("omega"))
      {
         auto omegaWriter = make_shared< field::VTKWriter< ScalarField_T, float32 > >(ids.omegaField, "omega");
         vtkOutput->addCellDataWriter(omegaWriter);
      }
      if (VTKWriter.getParameter< bool >("flag"))
      {
         auto flagWriter = make_shared< field::VTKWriter< FlagField_T > >(ids.flagField, "flag");
         vtkOutput->addCellDataWriter(flagWriter);
      }
      timeloop.addFuncAfterTimeStep(vtk::writeFiles(vtkOutput), "VTK Output");
   }

   // log remaining time
   const real_t remainingTimeLoggerFrequency =
      loggingParameters.getParameter< real_t >("remainingTimeLoggerFrequency", 3.0); // in seconds
   if (uint_c(remainingTimeLoggerFrequency) > 0)
   {
      timeloop.addFuncAfterTimeStep(
         timing::RemainingTimeLogger(timeloop.getNrOfTimeSteps(), remainingTimeLoggerFrequency),
         "remaining time logger");
   }

   // LBM stability check
   auto CheckerParameters      = config->getOneBlock("StabilityChecker");
   const uint_t checkFrequency = CheckerParameters.getParameter< uint_t >("checkFrequency", uint_t(0));
   if (checkFrequency > 0)
   {
      auto checkFunction = [](PdfField_T::value_type value) {  return value < math::abs(PdfField_T::value_type(10)); };
      timeloop.addFuncAfterTimeStep(
         makeSharedFunctor(field::makeStabilityChecker< PdfField_T, FlagField_T >(
            config, blocks, ids.pdfField, ids.flagField, setup.fluidUID, checkFunction)),
         "Stability check");
   }

   timeloop.addFuncBeforeTimeStep( SharedFunctor< Evaluation >(evaluation), "evaluation" );

   WALBERLA_MPI_BARRIER()
#if defined(WALBERLA_BUILD_WITH_GPU_SUPPORT)
   WALBERLA_GPU_CHECK(gpuDeviceSynchronize())
   WALBERLA_GPU_CHECK(gpuPeekAtLastError())
#endif
   // WALBERLA_LOG_INFO_ON_ROOT("Execute single timestep to fully complete the preprocessing done")

   //////////////////////
   /// RUN SIMULATION ///
   //////////////////////
   const lbm_generated::PerformanceEvaluation< FlagField_T > performance(blocks, ids.flagField, setup.fluidUID);
   field::CellCounter< FlagField_T > fluidCells(blocks, ids.flagField, setup.fluidUID);
   fluidCells();

   WALBERLA_LOG_INFO_ON_ROOT("Blocks created: " << blocks->getNumberOfBlocks())
   for (uint_t level = 0; level <= setup.refinementLevels; level++)
   {
      WALBERLA_LOG_INFO_ON_ROOT("Level " << level << " Blocks: " << blocks->getNumberOfBlocks(level))
   }

   WALBERLA_LOG_INFO_ON_ROOT( "Benchmark run data:"
                             "\n- simulation parameters:"   << std::setprecision(16) <<
                             "\n   + collision model:     " << infoMap["collisionOperator"] <<
                             "\n   + stencil:             " << infoMap["stencil"] <<
                             "\n   + streaming:           " << infoMap["streamingPattern"] <<
                             "\n   + compressible:        " << ( StorageSpecification_T::compressible ? "yes" : "no" ) <<
                             "\n   + mesh levels:         " << setup.refinementLevels + uint_c(1) <<
                             "\n   + resolution:          " << setup.dxC << " - on the coarsest grid" <<
                             "\n   + resolution:          " << setup.dxF << " - on the finest grid" <<
                             "\n- simulation properties:  "
                             "\n   + sphere pos.(x):      " << setup.sphereXPosition * setup.dxC << " [m]" <<
                             "\n   + sphere pos.(y):      " << setup.sphereYPosition * setup.dxC << " [m]" <<
                             "\n   + sphere pos.(z):      " << setup.sphereZPosition * setup.dxC << " [m]" <<
                             "\n   + sphere radius:       " << setup.sphereRadius * setup.dxC << " [m]" <<
                             "\n   + kin. viscosity:      " << setup.viscosity * setup.dxC * setup.dxC / setup.dt << " [m^2/s] (on the coarsest grid)" <<
                             "\n   + viscosity LB:        " << setup.viscosity  << " [dx*dx/dt] (on the coarsest grid)" <<
                             "\n   + omega:               " << setup.omega << " (on the coarsest grid)" <<
                             "\n   + rho:                 " << setup.rho << " [kg/m^3]" <<
                             "\n   + inflow velocity:     " << setup.referenceVelocity << " [m/s]" <<
                             "\n   + lattice velocity:    " << setup.inflowVelocity << " [dx/dt]" <<
                             "\n   + Reynolds number:     " << setup.reynoldsNumber <<
                             "\n   + Mach number:         " << setup.machNumber <<
                             "\n   + dt (coarsest grid):  " << setup.dt << " [s]"
                             "\n   + #time steps:         " << timeloop.getNrOfTimeSteps() << " (on the coarsest grid, " << ( real_c(1.0) / setup.dt ) << " for 1s of real time)" <<
                             "\n   + simulation time:     " << ( real_c( timeloop.getNrOfTimeSteps() ) * setup.dt ) << " [s]" )

   WALBERLA_LOG_INFO_ON_ROOT("Starting Simulation")
   WcTimingPool timeloopTiming;
   WcTimer simTimer;

   WALBERLA_MPI_BARRIER()
#if defined(WALBERLA_BUILD_WITH_GPU_SUPPORT)
   WALBERLA_GPU_CHECK(gpuDeviceSynchronize())
   WALBERLA_GPU_CHECK(gpuPeekAtLastError())
#endif

   simTimer.start();
   timeloop.run(timeloopTiming);
#if defined(WALBERLA_BUILD_WITH_GPU_SUPPORT)
   WALBERLA_GPU_CHECK(gpuDeviceSynchronize())
   WALBERLA_GPU_CHECK(gpuPeekAtLastError())
#endif
   WALBERLA_MPI_BARRIER()
   simTimer.end();
   WALBERLA_LOG_INFO_ON_ROOT("Simulation finished")
   double time = double_c(simTimer.max());
   WALBERLA_MPI_SECTION() { walberla::mpi::reduceInplace(time, walberla::mpi::MAX); }
   performance.logResultOnRoot(setup.timesteps, time);

   const auto reducedTimeloopTiming = timeloopTiming.getReduced();
   WALBERLA_LOG_RESULT_ON_ROOT("Time loop timing:\n" << *reducedTimeloopTiming)

   printResidentMemoryStatistics();
   return EXIT_SUCCESS;
}