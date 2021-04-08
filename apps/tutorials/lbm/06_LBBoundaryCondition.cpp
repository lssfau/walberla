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
//! \file 06_LBBoundaryConditions.cpp
//! \author Helen Schottenhamml <helen.schottenhamml@fau.de>
//
//======================================================================================================================

#include "blockforest/all.h"

#include "core/all.h"

#include "domain_decomposition/all.h"

#include "field/all.h"

#include "geometry/all.h"

#include "gui/all.h"

#include "lbm/all.h"

#include "timeloop/all.h"

namespace walberla
{
using LatticeModel_T         = lbm::D2Q9< lbm::collision_model::SRT >;
using Stencil_T              = LatticeModel_T::Stencil;
using CommunicationStencil_T = LatticeModel_T::CommunicationStencil;

using PdfField_T = lbm::PdfField< LatticeModel_T >;

using flag_t      = walberla::uint16_t;
using FlagField_T = FlagField< flag_t >;

//! [variableDefines]
// number of ghost layers
const uint_t FieldGhostLayers = uint_t(1);

// unique identifiers for flags
const FlagUID FluidFlagUID("Fluid Flag");
const FlagUID NoSlipFlagUID("NoSlip Flag");
const FlagUID FreeSlipFlagUID("FreeSlip Flag");
const FlagUID SimpleUBBFlagUID("SimpleUBB Flag");
const FlagUID UBBFlagUID("UBB Flag");
const FlagUID DynamicUBBFlagUID("DynamicUBB Flag");
const FlagUID ParserUBBFlagUID("ParserUBB Flag");
const FlagUID SimplePressureFlagUID("SimplePressure Flag");
const FlagUID PressureFlagUID("Pressure Flag");
const FlagUID OutletFlagUID("Outlet Flag");
const FlagUID SimplePABFlagUID("SimplePAB Flag");
//! [variableDefines]

///////////////////////////
/// PARAMETER STRUCTURE ///
///////////////////////////

//! [setupStruct]
struct Setup
{
   std::string wallType;
   std::string inflowType;
   std::string outflowType;

   Vector3< real_t > inflowVelocity;
   real_t outflowPressure;

   // DynamicUBB
   real_t period;

   // ParserUBB
   Config::BlockHandle parser;

   // SimplePAB
   real_t omega;
};
//! [setupStruct]

////////////////////////////////
/// CUSTOM BOUNDARY HANDLING ///
////////////////////////////////

//! [VelocityFunctor]
class VelocityFunctor
{
 public:
   VelocityFunctor(const Vector3< real_t >& velocity, real_t period, real_t H)
      : velocity_(velocity), period_(period), H_(H)
   {
      constantTerm_ = real_t(4) * velocity_[0] / (H_ * H_);
   }

   void operator()(const real_t time)
   {
      amplitude_ = constantTerm_ * real_t(0.5) * (real_t(1) - std::cos(real_t(2) * math::pi * time / period_));
   }

   Vector3< real_t > operator()(const Vector3< real_t >& pos, const real_t)
   {
      return Vector3< real_t >(amplitude_ * pos[1] * (H_ - pos[1]), real_t(0), real_t(0));
   }

 private:
   const Vector3< real_t > velocity_;
   const real_t period_;
   const real_t H_;

   real_t constantTerm_; // part of the velocity that is constant in both time and space
   real_t amplitude_;
};
//! [VelocityFunctor]

//! [boundaryTypedefs]
/// boundary handling
using NoSlip_T    = lbm::NoSlip< LatticeModel_T, flag_t >;
using FreeSlip_T  = lbm::FreeSlip< LatticeModel_T, FlagField_T >;

using SimpleUBB_T    = lbm::SimpleUBB< LatticeModel_T, flag_t >;
using UBB_T          = lbm::UBB< LatticeModel_T, flag_t >;
using DynamicUBB_T   = lbm::DynamicUBB< LatticeModel_T, flag_t, VelocityFunctor >;
using ParserUBB_T    = lbm::ParserUBB< LatticeModel_T, flag_t >;

using SimplePressure_T  = lbm::SimplePressure< LatticeModel_T, flag_t >;
using Pressure_T        = lbm::Pressure< LatticeModel_T, flag_t >;
using Outlet_T          = lbm::Outlet< LatticeModel_T, FlagField_T >;
using SimplePAB_T       = lbm::SimplePAB< LatticeModel_T, FlagField_T >;

using BoundaryHandling_T = BoundaryHandling< FlagField_T, Stencil_T, NoSlip_T, FreeSlip_T,
                                             SimpleUBB_T, UBB_T, DynamicUBB_T, ParserUBB_T,
                                             SimplePressure_T, Pressure_T, Outlet_T, SimplePAB_T >;
//! [boundaryTypedefs]

//! [myBoundaryHandlingDeclaration]
class MyBoundaryHandling
{
 public:
   MyBoundaryHandling(const BlockDataID& flagFieldID, const BlockDataID& pdfFieldID, const Setup & setup,
                      const std::shared_ptr< lbm::TimeTracker >& timeTracker)
      : flagFieldID_(flagFieldID), pdfFieldID_(pdfFieldID), setup_(setup), timeTracker_(timeTracker)
   {}

   BoundaryHandling_T* operator()(IBlock* const block, const StructuredBlockStorage* const storage) const;

 private:
   const BlockDataID flagFieldID_;
   const BlockDataID pdfFieldID_;

   Setup setup_;
   std::shared_ptr< lbm::TimeTracker > timeTracker_;

}; // class MyBoundaryHandling
//! [myBoundaryHandlingDeclaration]

BoundaryHandling_T* MyBoundaryHandling::operator()(IBlock* const block,
                                                   const StructuredBlockStorage* const storage) const
{
   Vector3< real_t > domainSize(real_c(storage->getNumberOfXCells()), real_c(storage->getNumberOfYCells()),
                                real_c(storage->getNumberOfZCells()));

   real_t H = domainSize[1];

   VelocityFunctor velocity(setup_.inflowVelocity, setup_.period, H);

   WALBERLA_ASSERT_NOT_NULLPTR(block)

   //! [boundaryHandling_T fields]
   FlagField_T* flagField = block->getData< FlagField_T >(flagFieldID_);
   PdfField_T* pdfField   = block->getData< PdfField_T >(pdfFieldID_);

   const auto fluidFlag = flagField->getOrRegisterFlag(FluidFlagUID);
   //! [boundaryHandling_T fields]

   BoundaryHandling_T* handling = new BoundaryHandling_T(
      "Boundary Handling", flagField, fluidFlag,
      //! [handling_NoSlip]
      NoSlip_T("NoSlip", NoSlipFlagUID, pdfField),
      //! [handling_NoSlip]
      //! [handling_FreeSlip]
      FreeSlip_T("FreeSlip", FreeSlipFlagUID, pdfField, flagField, fluidFlag),
      //! [handling_FreeSlip]
      //! [handling_SimpleUBB]
      SimpleUBB_T("SimpleUBB", SimpleUBBFlagUID, pdfField, setup_.inflowVelocity),
      //! [handling_SimpleUBB]
      //! [handling_UBB]
      UBB_T("UBB", UBBFlagUID, pdfField),
      //! [handling_UBB]
      //! [handling_DynamicUBB]
      DynamicUBB_T("DynamicUBB", DynamicUBBFlagUID, pdfField, timeTracker_, storage->getLevel(*block), velocity,
                   block->getAABB()),
      //! [handling_DynamicUBB]
      //! [handling_ParserUBB]
      ParserUBB_T("ParserUBB", ParserUBBFlagUID, pdfField, flagField, timeTracker_, storage->getLevel(*block),
                  block->getAABB()),
      //! [handling_ParserUBB]
      //! [handling_SimplePressure]
      SimplePressure_T("SimplePressure", SimplePressureFlagUID, pdfField, setup_.outflowPressure),
      //! [handling_SimplePressure]
      //! [handling_Pressure]
      Pressure_T("Pressure", PressureFlagUID, pdfField),
      //! [handling_Pressure]
      //! [handling_Outlet]
      Outlet_T("Outlet", OutletFlagUID, pdfField, flagField, fluidFlag),
      //! [handling_Outlet]
      //! [handling_SimplePAB]
      SimplePAB_T("SimplePAB", SimplePABFlagUID, pdfField, flagField, setup_.outflowPressure, setup_.omega,
                  FluidFlagUID, NoSlipFlagUID));
      //! [handling_SimplePAB]

   //! [domainBB]
   CellInterval domainBB = storage->getDomainCellBB();
   storage->transformGlobalToBlockLocalCellInterval(domainBB, *block);
   //! [domainBB]

   //! [westBoundary]
   cell_idx_t ghost = cell_idx_t(FieldGhostLayers);

   domainBB.xMin() -= ghost;
   domainBB.xMax() += ghost;

   // WEST - Inflow
   CellInterval west(domainBB.xMin(), domainBB.yMin(),
                     domainBB.zMin(), domainBB.xMin(),
                     domainBB.yMax(), domainBB.zMin());

   //! [westBoundary]

   if (setup_.inflowType == "SimpleUBB") {
      //! [forceBoundary_SimpleUBB]
      handling->forceBoundary(SimpleUBBFlagUID, west);
      //! [forceBoundary_SimpleUBB]
   }
   else if (setup_.inflowType == "UBB")
   {
      //! [forceBoundary_UBB]
      Cell offset(0, 0, 0);
      storage->transformBlockLocalToGlobalCell(offset, *block);

      for (auto cellIt = west.begin(); cellIt != west.end(); ++cellIt)
      {
         Cell globalCell = *cellIt + offset;
         const real_t y = real_c(globalCell[1]);

         Vector3< real_t > ubbVel(0);
         ubbVel[0] = -real_t(4) * y * (y - H) / (H * H) * setup_.inflowVelocity[0];

         handling->forceBoundary(UBBFlagUID, cellIt->x(), cellIt->y(), cellIt->z(), UBB_T::Velocity(ubbVel));
      }
      //! [forceBoundary_UBB]
   }
   else if (setup_.inflowType == "DynamicUBB")
   {
      //! [forceBoundary_DynamicUBB]
      handling->forceBoundary(DynamicUBBFlagUID, west);
      //! [forceBoundary_DynamicUBB]
   }
   else if (setup_.inflowType == "ParserUBB")
   {
      //! [forceBoundary_ParserUBB_eqs]
      char x_eq[150];
      sprintf(x_eq, "0.1*4/%f/%f * y * (%f - y) * 0.5 * (1 - cos(2 * 3.1415926538 * t / %f));", H, H, H, setup_.period);

      std::array< std::string, 3 > eqs = { x_eq, "0", "0" };
      handling->forceBoundary(ParserUBBFlagUID, west, ParserUBB_T::Parser(eqs));
      //! [forceBoundary_ParserUBB_eqs]

      //! [forceBoundary_ParserUBB_config]
      handling->forceBoundary(ParserUBBFlagUID, west, ParserUBB_T::Parser(setup_.parser));
      //! [forceBoundary_ParserUBB_config]
   }
   else
   {
      WALBERLA_ABORT("Please specify a valid inflow type.")
   }

   // EAST - outflow
   CellInterval east(domainBB.xMax(), domainBB.yMin(),
                     domainBB.zMin(), domainBB.xMax(),
                     domainBB.yMax(), domainBB.zMin());

   if (setup_.outflowType == "SimplePressure") {
      //! [forceBoundary_SimplePressure]
      handling->forceBoundary(SimplePressureFlagUID, east);
      //! [forceBoundary_SimplePressure]
   }
   else if (setup_.outflowType == "Pressure")
   {
      //! [forceBoundary_Pressure]
      Cell offset(0, 0, 0);
      storage->transformBlockLocalToGlobalCell(offset, *block);

      for (auto cellIt = east.begin(); cellIt != east.end(); ++cellIt)
      {
         Cell globalCell = *cellIt + offset;
         const real_t y = real_c(globalCell[1]);

         real_t local_density =
            setup_.outflowPressure * (real_t(1.0) + real_t(0.01) * std::sin(real_t(2.0 * 3.1415926538) * y / H));

         handling->forceBoundary(PressureFlagUID, cellIt->x(), cellIt->y(), cellIt->z(),
                                 Pressure_T::LatticeDensity(local_density));
      }
      //! [forceBoundary_Pressure]
   }
   else if (setup_.outflowType == "Outlet")
   {
      //! [forceBoundary_Outlet]
      handling->forceBoundary(OutletFlagUID, east);
      //! [forceBoundary_Outlet]
   }
   else if (setup_.outflowType == "SimplePAB")
   {
      //! [forceBoundary_SimplePAB]
      handling->forceBoundary(SimplePABFlagUID, east);
      //! [forceBoundary_SimplePAB]
   }
   else
   {
      WALBERLA_ABORT("Please specify a valid outflow type.")
   }

   domainBB.yMin() -= ghost;
   domainBB.yMax() += ghost;

   // SOUTH - wall
   CellInterval south(domainBB.xMin(), domainBB.yMin(),
                      domainBB.zMin(), domainBB.xMax(),
                      domainBB.yMin(), domainBB.zMax());

   // NORTH - wall
   CellInterval north(domainBB.xMin(), domainBB.yMax(),
                      domainBB.zMin(), domainBB.xMax(),
                      domainBB.yMax(), domainBB.zMax());

   if (setup_.wallType == "NoSlip")
   {
      //! [forceBoundary_NoSlip]
      handling->forceBoundary(NoSlipFlagUID, south);
      //! [forceBoundary_NoSlip]
      handling->forceBoundary(NoSlipFlagUID, north);
   }
   else if (setup_.wallType == "FreeSlip")
   {
      //! [forceBoundary_FreeSlip]
      handling->forceBoundary(FreeSlipFlagUID, south);
      //! [forceBoundary_FreeSlip]
      handling->forceBoundary(FreeSlipFlagUID, north);
   }
   else
   {
      WALBERLA_ABORT("Please specify a valid wall type.")
   }

   //! [fillDomain]
   handling->fillWithDomain(domainBB);
   //! [fillDomain]

   return handling;
}

int main(int argc, char** argv)
{
   walberla::Environment walberlaEnv(argc, argv);

   auto blocks = blockforest::createUniformBlockGridFromConfig(walberlaEnv.config());

   // read parameters
   auto parameters = walberlaEnv.config()->getOneBlock("Parameters");

   const real_t omega = parameters.getParameter< real_t >("omega", real_c(1.4));
   const Vector3< real_t > initialVelocity =
      parameters.getParameter< Vector3< real_t > >("initialVelocity", Vector3< real_t >());
   const uint_t timesteps = parameters.getParameter< uint_t >("timesteps", uint_c(10));

   const double remainingTimeLoggerFrequency =
      parameters.getParameter< double >("remainingTimeLoggerFrequency", 3.0); // in seconds

   // create fields
   LatticeModel_T latticeModel = LatticeModel_T(lbm::collision_model::SRT(omega));
   BlockDataID pdfFieldID  = lbm::addPdfFieldToStorage(blocks, "pdf field", latticeModel, initialVelocity, real_t(1));
   BlockDataID flagFieldID = field::addFlagFieldToStorage< FlagField_T >(blocks, "flag field", FieldGhostLayers);

   // create and initialize boundary handling

   auto boundariesConfig = walberlaEnv.config()->getOneBlock("Boundaries");

   Setup setup;
   setup.wallType        = boundariesConfig.getParameter< std::string >("wallType", "NoSlip");
   setup.inflowType      = boundariesConfig.getParameter< std::string >("inflowType", "SimpleUBB");
   setup.outflowType     = boundariesConfig.getParameter< std::string >("outflowType", "SimplePressure");
   setup.inflowVelocity  = boundariesConfig.getParameter< Vector3< real_t > >("inflowVelocity", Vector3< real_t >());
   setup.outflowPressure = boundariesConfig.getParameter< real_t >("outflowPressure", real_t(1));

   setup.period = boundariesConfig.getParameter< real_t >("period", real_t(100));

   if (setup.inflowType == "ParserUBB") setup.parser = boundariesConfig.getBlock("Parser");

   setup.omega = omega;

   //! [timeTracker]
   std::shared_ptr< lbm::TimeTracker > timeTracker = std::make_shared< lbm::TimeTracker >();
   //! [timeTracker]

   //! [boundaryHandlingID]
   BlockDataID boundaryHandlingID = blocks->addStructuredBlockData< BoundaryHandling_T >(
      MyBoundaryHandling(flagFieldID, pdfFieldID, setup, timeTracker), "boundary handling");
   //! [boundaryHandlingID]

   // create time loop
   SweepTimeloop timeloop(blocks->getBlockStorage(), timesteps);

   // create communication for PdfField
   blockforest::communication::UniformBufferedScheme< CommunicationStencil_T > communication(blocks);
   communication.addPackInfo(make_shared< lbm::PdfFieldPackInfo< LatticeModel_T > >(pdfFieldID));

   // add LBM sweep and communication to time loop
   //! [boundarySweep]
   timeloop.add() << BeforeFunction(communication, "communication")
                  << Sweep(BoundaryHandling_T::getBlockSweep(boundaryHandlingID), "boundary handling");
   //! [boundarySweep]
   timeloop.add() << Sweep(
      makeSharedSweep(lbm::makeCellwiseSweep< LatticeModel_T, FlagField_T >(pdfFieldID, flagFieldID, FluidFlagUID)),
      "LB stream & collide");

   // increment time step counter
   //! [timeTracker_coupling]
   timeloop.addFuncAfterTimeStep(makeSharedFunctor(timeTracker), "time tracking");
   //! [timeTracker_coupling]

   // LBM stability check
   timeloop.addFuncAfterTimeStep(makeSharedFunctor(field::makeStabilityChecker< PdfField_T, FlagField_T >(
                                    walberlaEnv.config(), blocks, pdfFieldID, flagFieldID, FluidFlagUID)),
                                 "LBM stability check");

   // log remaining time
   timeloop.addFuncAfterTimeStep(timing::RemainingTimeLogger(timeloop.getNrOfTimeSteps(), remainingTimeLoggerFrequency),
                                 "remaining time logger");

   // add VTK output to time loop
   //   lbm::VTKOutput< LatticeModel_T, FlagField_T >::addToTimeloop(timeloop, blocks, walberlaEnv.config(), pdfFieldID,
   //                                                                flagFieldID, FluidFlagUID);

   auto vtkConfig = walberlaEnv.config()->getBlock("VTK");

   uint_t writeFrequency = vtkConfig.getBlock("fluid_field").getParameter< uint_t >("writeFrequency", uint_t(100));

   auto vtkOutput = vtk::createVTKOutput_BlockData(*blocks, "fluid_field", writeFrequency, FieldGhostLayers, false,
                                                   "vtk_out", "simulation_step", false, true, true, false, 0);

   auto velocityWriter  = std::make_shared< lbm::VelocityVTKWriter< LatticeModel_T, float > >(pdfFieldID, "velocity");
   auto flagFieldWriter = std::make_shared< field::VTKWriter< FlagField_T > >(flagFieldID, "flag field");

   vtkOutput->addCellDataWriter(velocityWriter);
   vtkOutput->addCellDataWriter(flagFieldWriter);

   timeloop.addFuncAfterTimeStep(vtk::writeFiles(vtkOutput), "VTKOutput");

   // create adaptors, so that the GUI also displays density and velocity
   // adaptors are like fields with the difference that they do not store values
   // but calculate the values based on other fields ( here the PdfField )
   field::addFieldAdaptor< lbm::Adaptor< LatticeModel_T >::Density >(blocks, pdfFieldID, "DensityAdaptor");
   field::addFieldAdaptor< lbm::Adaptor< LatticeModel_T >::VelocityVector >(blocks, pdfFieldID, "VelocityAdaptor");

   if (parameters.getParameter< bool >("useGui", false))
   {
      GUI gui(timeloop, blocks, argc, argv);
      lbm::connectToGui< LatticeModel_T >(gui);
      gui.run();
   }
   else
   {
      timeloop.run();
   }

   return EXIT_SUCCESS;
}
} // namespace walberla

int main(int argc, char** argv) { walberla::main(argc, argv); }
