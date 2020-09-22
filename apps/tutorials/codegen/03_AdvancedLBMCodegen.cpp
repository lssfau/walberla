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

#include "domain_decomposition/all.h"

#include "field/all.h"
#include "field/vtk/VTKWriter.h"

#include "geometry/all.h"

#include "lbm/vtk/VTKOutput.h"

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

   ////////////////////////////////////
   /// PDF Field and Velocity Setup ///
   ////////////////////////////////////

   BlockDataID velocityFieldId = field::addToStorage< VectorField_T >(blocks, "velocity", real_c(0.0), field::fzyx);

   BlockDataID pdfFieldId  = field::addToStorage< PdfField_T >(blocks, "pdf field", real_c(0.0), field::fzyx);
   BlockDataID flagFieldId = field::addFlagFieldToStorage< FlagField_T >(blocks, "flag field");

   // Velocity field setup
   auto shearFlowSetup = walberlaEnv.config()->getOneBlock("ShearFlowSetup");
   initShearFlowVelocityField(blocks, velocityFieldId, shearFlowSetup);

   real_t rho = shearFlowSetup.getParameter("rho", real_c(1.0));

   // pdfs setup
   pystencils::InitialPDFsSetter pdfSetter(pdfFieldId, velocityFieldId, rho);

   for (auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt)
   {
      pdfSetter(&(*blockIt));
   }

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
   blockforest::communication::UniformBufferedScheme< Stencil_T > communication(blocks);
   communication.addPackInfo(make_shared< PackInfo_T >(pdfFieldId));

   // Timeloop
   timeloop.add() << BeforeFunction(communication, "communication") << Sweep(noSlip);
   timeloop.add() << Sweep(pystencils::CumulantMRTSweep(pdfFieldId, velocityFieldId, omega));

   // Stability Checker
   timeloop.addFuncAfterTimeStep(makeSharedFunctor(field::makeStabilityChecker< PdfField_T, FlagField_T >(
                                    walberlaEnv.config(), blocks, pdfFieldId, flagFieldId, fluidFlagUID)),
                                 "LBM stability check");

   // Time logger
   timeloop.addFuncAfterTimeStep(timing::RemainingTimeLogger(timeloop.getNrOfTimeSteps(), remainingTimeLoggerFrequency),
                                 "remaining time logger");

   uint_t vtkWriteFrequency = uint_c(100);
   if (vtkWriteFrequency > 0)
   {
      const std::string path = "vtk_out/tut03";
      auto vtkOutput = vtk::createVTKOutput_BlockData(*blocks, "cumulant_mrt_velocity_field", vtkWriteFrequency, 0,
                                                      false, path, "simulation_step", false, true, true, false, 0);

      auto velWriter = make_shared< field::VTKWriter< VectorField_T > >(velocityFieldId, "Velocity");
      vtkOutput->addCellDataWriter(velWriter);

      timeloop.addFuncBeforeTimeStep(vtk::writeFiles(vtkOutput), "VTK Output");
   }

   timeloop.run();

   return EXIT_SUCCESS;
}

} // namespace walberla

int main(int argc, char** argv) { return walberla::main(argc, argv); }