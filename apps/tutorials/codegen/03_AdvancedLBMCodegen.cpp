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

#include "geometry/all.h"

#include "lbm/boundary/factories/DefaultBoundaryHandling.h"
#include "lbm/communication/PdfFieldPackInfo.h"
#include "lbm/field/AddToStorage.h"
#include "lbm/field/PdfField.h"
#include "lbm/field/initializer/all.h"
#include "lbm/vtk/VTKOutput.h"

#include "stencil/D2Q9.h"

#include "timeloop/all.h"

//    Codegen Includes
#include "CumulantMRTNoSlip.h"
#include "CumulantMRTPackInfo.h"
#include "CumulantMRTSweep.h"
#include "DensityAndVelocityFieldSetter.h"

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

////////////////////////
/// VTK Output Setup ///
////////////////////////

struct VTKOutputSetup
{
 private:
   const ConstBlockDataID velocityFieldId_;
   const ConstBlockDataID flagFieldId_;
   const FlagUID& domainFlag_;

 public:
   VTKOutputSetup(const BlockDataID& velocityFieldId, const BlockDataID& flagFieldId, const FlagUID& domainFlag)
      : velocityFieldId_(velocityFieldId), flagFieldId_(flagFieldId), domainFlag_(domainFlag)
   {}

   void operator()(std::vector< shared_ptr< vtk::BlockCellDataWriterInterface > >& writers,
                   std::map< std::string, vtk::VTKOutput::CellFilter >& filters,
                   std::map< std::string, vtk::VTKOutput::BeforeFunction >& /*beforeFunctions*/) const
   {
      // Add VTK writers for velocity and velocity magnitude
      writers.push_back(make_shared< field::VTKWriter< VectorField_T, float > >(velocityFieldId_, "Velocity"));

      // Add domain cell filter
      filters["DomainFilter"] = field::FlagFieldCellFilter< FlagField_T >(flagFieldId_, domainFlag_);
   }

   void initializeAndAdd(SweepTimeloop& timeloop, const shared_ptr< StructuredBlockForest >& blocks,
                         const shared_ptr<Config> & config)
   {
      std::map< std::string, vtk::SelectableOutputFunction > vtkOutputFunctions;
      vtk::initializeVTKOutput(vtkOutputFunctions, *this, blocks, config, "VTK");

      for (auto funcIt = vtkOutputFunctions.begin(); funcIt != vtkOutputFunctions.end(); ++funcIt)
      {
         timeloop.addFuncBeforeTimeStep(funcIt->second.outputFunction, funcIt->first,
                                        funcIt->second.requiredGlobalStates, funcIt->second.incompatibleGlobalStates);
      }
   }
};

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

   // First, set up the initial shear flow in the velocity field...

   auto shearFlowSetup = walberlaEnv.config()->getOneBlock("ShearFlowSetup");
   initShearFlowVelocityField(blocks, velocityFieldId, shearFlowSetup);

   real_t rho = shearFlowSetup.getParameter("rho", real_c(1.8));

   // ... and then use the generated PDF setter to initialize the PDF field.
   pystencils::DensityAndVelocityFieldSetter pdfSetter(pdfFieldId, velocityFieldId, rho);

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

   // VTK Output
   VTKOutputSetup vtkOutput(velocityFieldId, flagFieldId, fluidFlagUID);
   vtkOutput.initializeAndAdd(timeloop, blocks, walberlaEnv.config());

   timeloop.run();

   return EXIT_SUCCESS;
}

} // namespace walberla

int main(int argc, char** argv) { return walberla::main(argc, argv); }