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
//! \file 03_AdvancedLBMCodegen.cpp
//! \author Frederik Hennig <frederik.hennig@fau.de>
//
//======================================================================================================================

#include "blockforest/all.h"

#include "core/all.h"

#if defined(WALBERLA_BUILD_WITH_CUDA)
#   include "cuda/AddGPUFieldToStorage.h"
#   include "cuda/DeviceSelectMPI.h"
#   include "cuda/HostFieldAllocator.h"
#   include "cuda/NVTX.h"
#   include "cuda/ParallelStreams.h"
#   include "cuda/communication/GPUPackInfo.h"
#   include "cuda/communication/UniformGPUScheme.h"
#endif

#include "domain_decomposition/all.h"

#include "field/all.h"
#include "field/vtk/VTKWriter.h"

#include "geometry/all.h"

#include "stencil/D2Q9.h"

#include "timeloop/all.h"

//    Codegen Includes
#include "CumulantMRTNoSlip.h"
#include "CumulantMRTPackInfo.h"
#include "CumulantMRTSweep.h"
#include "InitialPDFsSetter.h"

namespace walberla
{
///////////////////////
/// Typedef Aliases ///
///////////////////////

// Communication Pack Info
typedef pystencils::CumulantMRTPackInfo PackInfo_T;

// LB Method Stencil
typedef stencil::D2Q9 Stencil_T;

// PDF field type
typedef field::GhostLayerField< real_t, Stencil_T::Size > PdfField_T;

// Velocity Field Type
typedef field::GhostLayerField< real_t, Stencil_T::D > VectorField_T;

// Boundary Handling
typedef walberla::uint8_t flag_t;
typedef FlagField< flag_t > FlagField_T;
typedef lbm::CumulantMRTNoSlip NoSlip_T;

#if defined(WALBERLA_BUILD_WITH_CUDA)
typedef cuda::GPUField< double > GPUField;
#endif

//////////////////////////////////////////
/// Shear Flow Velocity Initialization ///
//////////////////////////////////////////

void initShearFlowVelocityField(const shared_ptr< StructuredBlockForest >& blocks, const BlockDataID& velocityFieldId,
                                const Config::BlockHandle& config)
{
   math::RealRandom< real_t > rng(config.getParameter< std::mt19937::result_type >("noiseSeed", 42));

   real_t velocityMagnitude = config.getParameter< real_t >("velocityMagnitude", real_c(0.08));
   real_t noiseMagnitude    = config.getParameter< real_t >("noiseMagnitude", real_c(0.1) * velocityMagnitude);

   real_t n_y = real_c(blocks->getNumberOfYCells());

   for (auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt)
   {
      auto u = (*blockIt).getData< VectorField_T >(velocityFieldId);

      for (auto cellIt = u->beginWithGhostLayerXYZ(); cellIt != u->end(); ++cellIt)
      {
         Cell globalCell(cellIt.cell());
         blocks->transformBlockLocalToGlobalCell(globalCell, *blockIt);

         real_t relative_y = real_c(globalCell.y()) / n_y;

         u->get(cellIt.cell(), 0) = relative_y < 0.3 || relative_y > 0.7 ? velocityMagnitude : -velocityMagnitude;

         u->get(cellIt.cell(), 1) = noiseMagnitude * rng();
      }
   }
}

/////////////////////
/// Main Function ///
/////////////////////

int main(int argc, char** argv)
{
   walberla::Environment walberlaEnv(argc, argv);

   if (!walberlaEnv.config()) { WALBERLA_ABORT("No configuration file specified!"); }

   ///////////////////////////////////////////////////////
   /// Block Storage Creation and Simulation Parameter ///
   ///////////////////////////////////////////////////////

   auto blocks = blockforest::createUniformBlockGridFromConfig(walberlaEnv.config());

   // read parameters
   auto parameters = walberlaEnv.config()->getOneBlock("Parameters");

   const uint_t timesteps = parameters.getParameter< uint_t >("timesteps", uint_c(10));
   const real_t omega     = parameters.getParameter< real_t >("omega", real_c(1.8));
   const double remainingTimeLoggerFrequency =
      parameters.getParameter< double >("remainingTimeLoggerFrequency", 3.0); // in seconds
   const uint_t VTKwriteFrequency = parameters.getParameter< uint_t >("VTKwriteFrequency", 1000);

   ////////////////////////////////////
   /// PDF Field and Velocity Setup ///
   ////////////////////////////////////

   // Common Fields
   BlockDataID velocityFieldId = field::addToStorage< VectorField_T >(blocks, "velocity", real_c(0.0), field::fzyx);
   BlockDataID flagFieldId     = field::addFlagFieldToStorage< FlagField_T >(blocks, "flag field");

#if defined(WALBERLA_BUILD_WITH_CUDA)
   // GPU Field for PDFs
   BlockDataID pdfFieldId = cuda::addGPUFieldToStorage< cuda::GPUField< real_t > >(
      blocks, "pdf field on GPU", Stencil_T::Size, field::fzyx, uint_t(1));

   // GPU Velocity Field
   BlockDataID velocityFieldIdGPU =
      cuda::addGPUFieldToStorage< VectorField_T >(blocks, velocityFieldId, "velocity on GPU", true);
#else
   // CPU Field for PDFs
   BlockDataID pdfFieldId = field::addToStorage< PdfField_T >(blocks, "pdf field", real_c(0.0), field::fzyx);
#endif

   // Velocity field setup
   auto shearFlowSetup = walberlaEnv.config()->getOneBlock("ShearFlowSetup");
   initShearFlowVelocityField(blocks, velocityFieldId, shearFlowSetup);

   real_t rho = shearFlowSetup.getParameter("rho", real_c(1.0));

   // pdfs setup
#if defined(WALBERLA_BUILD_WITH_CUDA)
   cuda::fieldCpy< GPUField, VectorField_T >(blocks, velocityFieldIdGPU, velocityFieldId);
   pystencils::InitialPDFsSetter pdfSetter(pdfFieldId, velocityFieldIdGPU, rho);
#else
   pystencils::InitialPDFsSetter pdfSetter(pdfFieldId, velocityFieldId, rho);
#endif

   for (auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt)
   {
      pdfSetter(&(*blockIt));
   }

   /////////////
   /// Sweep ///
   /////////////

#if defined(WALBERLA_BUILD_WITH_CUDA)
   pystencils::CumulantMRTSweep CumulantMRTSweep(pdfFieldId, velocityFieldIdGPU, omega);
#else
   pystencils::CumulantMRTSweep CumulantMRTSweep(pdfFieldId, velocityFieldId, omega);
#endif

   /////////////////////////
   /// Boundary Handling ///
   /////////////////////////

   const FlagUID fluidFlagUID("Fluid");

   auto boundariesConfig = walberlaEnv.config()->getOneBlock("Boundaries");

   NoSlip_T noSlip(blocks, pdfFieldId);

   geometry::initBoundaryHandling< FlagField_T >(*blocks, flagFieldId, boundariesConfig);
   geometry::setNonBoundaryCellsToDomain< FlagField_T >(*blocks, flagFieldId, fluidFlagUID);

   noSlip.fillFromFlagField< FlagField_T >(blocks, flagFieldId, FlagUID("NoSlip"), fluidFlagUID);

   /////////////////
   /// Time Loop ///
   /////////////////

   SweepTimeloop timeloop(blocks->getBlockStorage(), timesteps);

   // Communication
#if defined(WALBERLA_BUILD_WITH_CUDA)
   cuda::communication::UniformGPUScheme< Stencil_T > com(blocks, 0);
   com.addPackInfo(make_shared< PackInfo_T >(pdfFieldId));
   auto communication = std::function< void() >([&]() { com.communicate(nullptr); });
#else
   blockforest::communication::UniformBufferedScheme< Stencil_T > communication(blocks);
   communication.addPackInfo(make_shared< PackInfo_T >(pdfFieldId));
#endif

   // Timeloop
   timeloop.add() << BeforeFunction(communication, "communication") << Sweep(noSlip);
   timeloop.add() << Sweep(CumulantMRTSweep);

   // Time logger
   timeloop.addFuncAfterTimeStep(timing::RemainingTimeLogger(timeloop.getNrOfTimeSteps(), remainingTimeLoggerFrequency),
                                 "remaining time logger");

   if (VTKwriteFrequency > 0)
   {
      const std::string path = "vtk_out/tut03";
      auto vtkOutput = vtk::createVTKOutput_BlockData(*blocks, "cumulant_mrt_velocity_field", VTKwriteFrequency, 0,
                                                      false, path, "simulation_step", false, true, true, false, 0);

#if defined(WALBERLA_BUILD_WITH_CUDA)
      // Copy velocity data to CPU before output
      vtkOutput->addBeforeFunction(
         [&]() { cuda::fieldCpy< VectorField_T, GPUField >(blocks, velocityFieldId, velocityFieldIdGPU); });
#endif

      auto velWriter = make_shared< field::VTKWriter< VectorField_T > >(velocityFieldId, "Velocity");
      vtkOutput->addCellDataWriter(velWriter);

      timeloop.addFuncBeforeTimeStep(vtk::writeFiles(vtkOutput), "VTK Output");
   }

   timeloop.run();

   return EXIT_SUCCESS;
}

} // namespace walberla

int main(int argc, char** argv) { return walberla::main(argc, argv); }