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
//! \file 02_LBMLatticeModelGeneration.cpp
//! \author Frederik Hennig <frederik.hennig@fau.de>
//
//======================================================================================================================

#include "blockforest/all.h"

#include "core/all.h"

#include "domain_decomposition/all.h"

#include "field/all.h"

#include "geometry/all.h"

#include "lbm/communication/PdfFieldPackInfo.h"
#include "lbm/field/AddToStorage.h"
#include "lbm/field/PdfField.h"
#include "lbm/field/initializer/all.h"
#include "lbm/vtk/VTKOutput.h"

#include "timeloop/all.h"

//    Codegen Includes
#include "CumulantMRTLatticeModel.h"
#include "EntropicMRTLatticeModel.h"
#include "SRTLatticeModel.h"

//    Generated Communication Pack Infos
#include "CumulantMRTPackInfo.h"
#include "EntropicMRTPackInfo.h"
#include "SRTPackInfo.h"

//    Generated NoSlip Boundaries
#include "CumulantMRTNoSlip.h"
#include "EntropicMRTNoSlip.h"
#include "SRTNoSlip.h"

namespace walberla
{
///////////////////////
/// Typedef Aliases ///
///////////////////////

// Set a typedef alias for our generated lattice model and the boundary.
// Changing the Method only requires changing these aliases.

typedef lbm::SRTLatticeModel LatticeModel_T;
typedef lbm::SRTNoSlip NoSlip_T;

// Entropic Model
// typedef lbm::EntropicMRTLatticeModel LatticeModel_T;
// typedef lbm::EntropicMRTNoSlip NoSlip_T;

// Cumulant Model
// typedef lbm::CumulantMRTLatticeModel LatticeModel_T;
// typedef lbm::CumulantMRTNoSlip NoSlip_T;

// Communication Pack Info
typedef lbm::PdfFieldPackInfo< LatticeModel_T > PackInfo_T;

// Also set aliases for the stencils involved...
typedef LatticeModel_T::Stencil Stencil_T;
typedef LatticeModel_T::CommunicationStencil CommunicationStencil_T;

// ... and for the required field types. These include the PdfField for the simulation
typedef lbm::PdfField< LatticeModel_T > PdfField_T;

// and scalar- and vector-valued fields for output of macroscopic quantities.
typedef GhostLayerField< real_t, LatticeModel_T::Stencil::D > VectorField_T;
typedef GhostLayerField< real_t, 1 > ScalarField_T;

// Also, for boundary handling, a flag data type and flag field are required.
typedef walberla::uint8_t flag_t;
typedef FlagField< flag_t > FlagField_T;

/////////////////////////////////////////
/// Shear Flow Initialization Functor ///
/////////////////////////////////////////

struct ShearFlowInit
{
 private:
   lbm::initializer::ExprSystemInitFunction exprInitFunc_;
   const real_t noiseMagnitude_;
   math::RealRandom< real_t > rng_;

 public:
   ShearFlowInit(const shared_ptr< StructuredBlockForest >& blocks, const Config::BlockHandle& setup)
      : exprInitFunc_(blocks), noiseMagnitude_(setup.getParameter< real_t >("u_y_noise_magnitude")),
        rng_(setup.getParameter< std::mt19937::result_type >("noise_seed"))
   {
      if (!exprInitFunc_.parse(setup)) { WALBERLA_ABORT("Shear Flow Setup was incomplete."); }
   }

   std::vector< real_t > operator()(const Cell& globalCell)
   {
      std::vector< real_t > densityAndVelocity = exprInitFunc_(globalCell);
      real_t yPerturbation                     = noiseMagnitude_ * rng_();
      densityAndVelocity[2] += yPerturbation;

      return densityAndVelocity;
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

   ///////////////////
   /// Field Setup ///
   ///////////////////

   // TODO: Do we need a force field?
   // BlockDataID forceFieldId = field::addToStorage<VectorField_T>( blocks, "Force", real_t( 0.0 ));

   // Velocity output field
   // We don't need that either, do we?
   BlockDataID velFieldId = field::addToStorage< VectorField_T >(blocks, "Velocity", real_t(0.0));

   LatticeModel_T latticeModel = LatticeModel_T(velFieldId, omega);
   BlockDataID pdfFieldId      = lbm::addPdfFieldToStorage(blocks, "pdf field", latticeModel, field::zyxf);
   BlockDataID flagFieldId     = field::addFlagFieldToStorage< FlagField_T >(blocks, "flag field");

   ////////////////////////
   /// Shear Flow Setup ///
   ////////////////////////

   auto shearFlowSetup = walberlaEnv.config()->getOneBlock("ShearFlowSetup");
   ShearFlowInit shearFlowInitFunc(blocks, shearFlowSetup);
   lbm::initializer::PdfFieldInitializer< LatticeModel_T > fieldInit(pdfFieldId, blocks);
   fieldInit.initDensityAndVelocity(shearFlowInitFunc);

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
   blockforest::communication::UniformBufferedScheme< CommunicationStencil_T > communication(blocks);
   communication.addPackInfo(make_shared< PackInfo_T >(pdfFieldId));
   timeloop.add() << BeforeFunction(communication, "communication");

   // Boundaries and LBM Sweep
   timeloop.add() << Sweep(noSlip, "NoSlip Boundaries")
                  << Sweep(LatticeModel_T::Sweep(pdfFieldId), "LB stream & collide");

   // Stability Checker
   timeloop.addFuncAfterTimeStep(makeSharedFunctor(field::makeStabilityChecker< PdfField_T, FlagField_T >(
                                    walberlaEnv.config(), blocks, pdfFieldId, flagFieldId, fluidFlagUID)),
                                 "LBM stability check");

   // Time logger
   timeloop.addFuncAfterTimeStep(timing::RemainingTimeLogger(timeloop.getNrOfTimeSteps(), remainingTimeLoggerFrequency),
                                 "remaining time logger");

   // LBM VTK Output
   lbm::VTKOutput< LatticeModel_T, FlagField_T >::addToTimeloop(timeloop, blocks, walberlaEnv.config(), pdfFieldId,
                                                                flagFieldId, fluidFlagUID);

   timeloop.run();

   return EXIT_SUCCESS;
}

} // namespace walberla

int main(int argc, char** argv) { return walberla::main(argc, argv); }