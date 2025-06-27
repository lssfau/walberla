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

#include "lbm/boundary/factories/DefaultBoundaryHandling.h"
#include "lbm/field/AddToStorage.h"
#include "lbm/field/PdfField.h"
#include "lbm/field/initializer/all.h"
#include "lbm/vtk/VTKOutput.h"

#include "timeloop/all.h"

//    Codegen Includes
#include "SRTLatticeModel.h"
#include "SRTPackInfo.h"

namespace walberla
{
///////////////////////
/// Typedef Aliases ///
///////////////////////

// Typedef alias for the lattice model
using LatticeModel_T = lbm::SRTLatticeModel;

// Communication Pack Info
using PackInfo_T = pystencils::SRTPackInfo;

// Also set aliases for the stencils involved...
using Stencil_T = LatticeModel_T::Stencil;
using CommunicationStencil_T = LatticeModel_T::CommunicationStencil;

// ... and for the required field types.
using PdfField_T = lbm::PdfField< LatticeModel_T >;

// Also, for boundary handling, a flag data type and flag field are required.
using flag_t = walberla::uint8_t;
using FlagField_T = FlagField< flag_t >;
using BHFactory = lbm::DefaultBoundaryHandlingFactory< LatticeModel_T, FlagField_T >;

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
      : exprInitFunc_(blocks), noiseMagnitude_(setup.getParameter< real_t >("noise_magnitude")),
        rng_(setup.getParameter< std::mt19937::result_type >("noise_seed", 42))
   {
      if (!exprInitFunc_.parse(setup)) { WALBERLA_ABORT("Shear Flow Setup was incomplete."); }
   }

   std::vector< real_t > operator()(const Cell& globalCell)
   {
      std::vector< real_t > densityAndVelocity = exprInitFunc_(globalCell);
      auto yPerturbation                     = real_c(noiseMagnitude_ * rng_());
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
   const real_t remainingTimeLoggerFrequency =
      parameters.getParameter< real_t >("remainingTimeLoggerFrequency", real_c(3.0)); // in seconds

   ///////////////////
   /// Field Setup ///
   ///////////////////

   LatticeModel_T const latticeModel = LatticeModel_T(omega);
   BlockDataID const pdfFieldId      = lbm::addPdfFieldToStorage(blocks, "pdf field", latticeModel, field::fzyx);
   BlockDataID const flagFieldId     = field::addFlagFieldToStorage< FlagField_T >(blocks, "flag field");

   ////////////////////////
   /// Shear Flow Setup ///
   ////////////////////////

   auto shearFlowSetup = walberlaEnv.config()->getOneBlock("ShearFlowSetup");
   ShearFlowInit shearFlowInitFunc(blocks, shearFlowSetup);
   lbm::initializer::PdfFieldInitializer< LatticeModel_T > const fieldInit(pdfFieldId, blocks);
   fieldInit.initDensityAndVelocity(shearFlowInitFunc);

   /////////////////////////
   /// Boundary Handling ///
   /////////////////////////

   const FlagUID fluidFlagUID("Fluid");

   auto boundariesConfig = walberlaEnv.config()->getOneBlock("Boundaries");

   BlockDataID boundaryHandlingId =
      BHFactory::addBoundaryHandlingToStorage(blocks, "boundary handling", flagFieldId, pdfFieldId, fluidFlagUID,
                                              Vector3< real_t >(), Vector3< real_t >(), real_c(0.0), real_c(0.0));

   geometry::initBoundaryHandling< BHFactory::BoundaryHandling >(*blocks, boundaryHandlingId, boundariesConfig);
   geometry::setNonBoundaryCellsToDomain< BHFactory::BoundaryHandling >(*blocks, boundaryHandlingId);

   /////////////////
   /// Time Loop ///
   /////////////////

   SweepTimeloop timeloop(blocks->getBlockStorage(), timesteps);

   // Communication
   blockforest::communication::UniformBufferedScheme< CommunicationStencil_T > communication(blocks);
   communication.addPackInfo(make_shared< PackInfo_T >(pdfFieldId));

   // Timeloop
   timeloop.add() << BeforeFunction(communication, "communication")
                  << Sweep(BHFactory::BoundaryHandling::getBlockSweep(boundaryHandlingId), "Boundary Handling");
   timeloop.add() << Sweep(LatticeModel_T::Sweep(pdfFieldId), "LB stream & collide");

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