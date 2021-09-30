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
//! \file GeneratedOutflowBC.cpp
//! \author Markus Holzer <markus.holzer@fau.de>
//! \brief Shear flow with dynamic UBB on the left (W) with linear flow profile from zero to u_max (flow in x-direction)
//!        along the y-direction. The upper wall (N) is a static UBB which applies u_max in x-direction.
//!        On the right (E) an outflow boundary is used and the lower wall (S) is a No-Slip wall.
//!        It is tested if the shear flow is correctly propagated through the entire domain after n timesteps.
//
//======================================================================================================================
#include "blockforest/Initialization.h"
#include "blockforest/communication/UniformBufferedScheme.h"

#include "core/Environment.h"
#include "core/timing/RemainingTimeLogger.h"

#include "field/AddToStorage.h"
#include "field/vtk/VTKWriter.h"

#include "geometry/InitBoundaryHandling.h"
#include "timeloop/SweepTimeloop.h"

// Generated Files
#include "GeneratedOutflowBC.h"

using namespace walberla;

using PackInfo_T  = lbm::GeneratedOutflowBC_PackInfo;
using flag_t      = walberla::uint8_t;
using FlagField_T = FlagField< flag_t >;

auto pdfFieldAdder = [](IBlock* const block, StructuredBlockStorage* const storage) {
   return new PdfField_T(storage->getNumberOfXCells(*block), storage->getNumberOfYCells(*block),
                         storage->getNumberOfZCells(*block), uint_t(1), field::fzyx,
                         make_shared< field::AllocateAligned< real_t, 64 > >());
};

////////////////////////////////////////////
// Linear Velocity Profile for left wall //
//////////////////////////////////////////

class ShearProfile
{
 public:

   ShearProfile( real_t inflow_velocity ) :
      inflow_velocity_( inflow_velocity ) {}

   Vector3< real_t > operator()( const Cell& pos, const shared_ptr< StructuredBlockForest >& SbF, IBlock& block ) const;

 private:

   const real_t inflow_velocity_;
}; // class ShearProfile

Vector3< real_t > ShearProfile::operator()( const Cell& pos, const shared_ptr< StructuredBlockForest >& SbF, IBlock& block ) const
{
   Cell globalCell;
   CellInterval domain = SbF->getDomainCellBB();
   real_t h_y          = domain.yMax() - domain.yMin();
   SbF->transformBlockLocalToGlobalCell(globalCell, block, pos);

   real_t u = inflow_velocity_ * (globalCell[1] / h_y);

   Vector3< real_t > result(u, 0.0, 0.0);
   return result;
}

//////////
// MAIN //
//////////

int main(int argc, char** argv)
{
   walberla::Environment walberlaEnv(argc, argv);

   auto blocks = blockforest::createUniformBlockGridFromConfig(walberlaEnv.config());

   // read parameters
   auto parameters = walberlaEnv.config()->getOneBlock("Parameters");

   const real_t omega     = parameters.getParameter< real_t >("omega", real_c(1.4));
   const real_t u_max     = parameters.getParameter< real_t >("u_max", real_t(0.05));
   const uint_t timesteps = parameters.getParameter< uint_t >("timesteps", uint_c(10));

   const double remainingTimeLoggerFrequency =
      parameters.getParameter< double >("remainingTimeLoggerFrequency", 3.0); // in seconds

   // create fields
   BlockDataID pdfFieldID     = blocks->addStructuredBlockData< PdfField_T >(pdfFieldAdder, "PDFs");
   BlockDataID velFieldID     = field::addToStorage< VelocityField_T >(blocks, "velocity", real_t(0), field::fzyx);
   BlockDataID densityFieldID = field::addToStorage< ScalarField_T >(blocks, "density", real_t(0), field::fzyx);

   BlockDataID flagFieldId = field::addFlagFieldToStorage< FlagField_T >(blocks, "flag field");

   pystencils::GeneratedOutflowBC_MacroSetter setterSweep(pdfFieldID, velFieldID);
   for (auto& block : *blocks)
      setterSweep(&block);

   // create and initialize boundary handling
   const FlagUID fluidFlagUID("Fluid");

   auto boundariesConfig = walberlaEnv.config()->getOneBlock("Boundaries");

   ShearProfile velocityCallback{u_max};
   std::function< Vector3< real_t >(const Cell&, const shared_ptr< StructuredBlockForest >&, IBlock&) >
      velocity_initialisation = velocityCallback;

   lbm::GeneratedOutflowBC_Dynamic_UBB ubb_dynamic(blocks, pdfFieldID, velocity_initialisation);
   lbm::GeneratedOutflowBC_Static_UBB ubb_static(blocks, pdfFieldID, u_max);
   lbm::GeneratedOutflowBC_NoSlip noSlip(blocks, pdfFieldID);
   lbm::GeneratedOutflowBC_Outflow outflow(blocks, pdfFieldID);

   geometry::initBoundaryHandling< FlagField_T >(*blocks, flagFieldId, boundariesConfig);
   geometry::setNonBoundaryCellsToDomain< FlagField_T >(*blocks, flagFieldId, fluidFlagUID);

   ubb_dynamic.fillFromFlagField< FlagField_T >(blocks, flagFieldId, FlagUID("UBB_Inflow"), fluidFlagUID);
   ubb_static.fillFromFlagField< FlagField_T >(blocks, flagFieldId, FlagUID("UBB_Wall"), fluidFlagUID);
   noSlip.fillFromFlagField< FlagField_T >(blocks, flagFieldId, FlagUID("NoSlip"), fluidFlagUID);
   outflow.fillFromFlagField< FlagField_T >(blocks, flagFieldId, FlagUID("Outflow"), fluidFlagUID);

   // create time loop
   SweepTimeloop timeloop(blocks->getBlockStorage(), timesteps);

   // create communication for PdfField
   blockforest::communication::UniformBufferedScheme< Stencil_T > communication(blocks);
   communication.addPackInfo(make_shared< PackInfo_T >(pdfFieldID));

   pystencils::GeneratedOutflowBC_Sweep UpdateSweep(densityFieldID, pdfFieldID, velFieldID, omega);

   // add LBM sweep and communication to time loop
   timeloop.add() << BeforeFunction(communication, "communication") << Sweep(noSlip, "noSlip boundary");
   timeloop.add() << Sweep(ubb_dynamic, "ubb inflow");
   timeloop.add() << Sweep(ubb_static, "ubb wall");
   timeloop.add() << Sweep(outflow, "outflow boundary");
   timeloop.add() << Sweep(UpdateSweep, "LB stream & collide");

   // log remaining time
   timeloop.addFuncAfterTimeStep(timing::RemainingTimeLogger(timeloop.getNrOfTimeSteps(), remainingTimeLoggerFrequency),
                                 "remaining time logger");

   // VTK Writer
   uint_t vtkWriteFrequency = parameters.getParameter< uint_t >("vtkWriteFrequency", 0);
   if (vtkWriteFrequency > 0)
   {
      auto vtkOutput = vtk::createVTKOutput_BlockData(*blocks, "GeneratedOutflowBC_VTK", vtkWriteFrequency, 0, false,
                                                      "vtk_out", "simulation_step", false, true, true, false, 0);

      auto velWriter     = make_shared< field::VTKWriter< VelocityField_T > >(velFieldID, "velocity");
      auto densityWriter = make_shared< field::VTKWriter< ScalarField_T > >(densityFieldID, "density");

      vtkOutput->addCellDataWriter(velWriter);
      vtkOutput->addCellDataWriter(densityWriter);

      timeloop.addFuncBeforeTimeStep(vtk::writeFiles(vtkOutput), "VTK Output");
   }
   timeloop.run();

   CellInterval domain = blocks->getDomainCellBB();
   real_t h_y          = domain.yMax() - domain.yMin();
   for (auto& block : *blocks)
   {
      auto velField = block.getData<VelocityField_T>(velFieldID);
      WALBERLA_FOR_ALL_CELLS_XYZ
      (
         velField,
   Cell globalCell;
         blocks->transformBlockLocalToGlobalCell(globalCell, block, Cell(x, y, z));
         WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(velField->get(x, y, z, 0), u_max * (globalCell[1] / h_y), 0.01)
      )
   }

   return EXIT_SUCCESS;
}
