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
//! \file GeneratedFreeSlip.cpp
//! \author Markus Holzer <markus.holzer@fau.de>
//! \brief Periodic flow in x direction with FreeSlip in Pipe setup. Checks if there are no gradients in the
//!        Velocity field in the end. Works for all streaming patterns (pull, push, aa, esotwist)
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
#include "GeneratedFreeSlip.h"

using namespace walberla;

using PackInfoEven_T  = lbm::GeneratedFreeSlip_PackInfoEven;
using PackInfoOdd_T  = lbm::GeneratedFreeSlip_PackInfoOdd;

using flag_t      = walberla::uint8_t;
using FlagField_T = FlagField< flag_t >;

auto pdfFieldAdder = [](IBlock* const block, StructuredBlockStorage* const storage) {
   return new PdfField_T(storage->getNumberOfXCells(*block), storage->getNumberOfYCells(*block),
                         storage->getNumberOfZCells(*block), uint_t(1), field::fzyx,
                         make_shared< field::AllocateAligned< real_t, 64 > >());
};

//////////
// MAIN //
//////////

class AlternatingBeforeFunction
{
 public:
   typedef std::function< void() > BeforeFunction;

   AlternatingBeforeFunction(BeforeFunction evenFunc, BeforeFunction oddFunc,
                             std::shared_ptr< lbm::TimestepTracker >& tracker)
      : tracker_(tracker), funcs_{ evenFunc, oddFunc } {};

   void operator()() { funcs_[tracker_->getCounter()](); }

 private:
   std::shared_ptr< lbm::TimestepTracker > tracker_;
   std::vector< BeforeFunction > funcs_;
};

int main(int argc, char** argv)
{
   walberla::Environment walberlaEnv(argc, argv);

   auto blocks = blockforest::createUniformBlockGridFromConfig(walberlaEnv.config());

   auto domainSetup                = walberlaEnv.config()->getOneBlock("DomainSetup");
   Vector3< uint_t > cellsPerBlock = domainSetup.getParameter< Vector3< uint_t > >("cellsPerBlock");

   // read parameters
   auto parameters = walberlaEnv.config()->getOneBlock("Parameters");

   const real_t omega     = parameters.getParameter< real_t >("omega", real_c(1.4));
   const uint_t timesteps = parameters.getParameter< uint_t >("timesteps", uint_c(10));

   const double remainingTimeLoggerFrequency =
      parameters.getParameter< double >("remainingTimeLoggerFrequency", 3.0); // in seconds

   // create fields
   BlockDataID pdfFieldID     = blocks->addStructuredBlockData< PdfField_T >(pdfFieldAdder, "PDFs");
   BlockDataID velFieldID     = field::addToStorage< VelocityField_T >(blocks, "velocity", real_t(0), field::fzyx);
   BlockDataID densityFieldID = field::addToStorage< ScalarField_T >(blocks, "density", real_t(1.0), field::fzyx);

   BlockDataID flagFieldId = field::addFlagFieldToStorage< FlagField_T >(blocks, "flag field");

   pystencils::GeneratedFreeSlip_MacroSetter setterSweep(densityFieldID, pdfFieldID, velFieldID);
   for (auto& block : *blocks)
      setterSweep(&block);

   // create and initialize boundary handling
   auto SpecialSetups = walberlaEnv.config()->getOneBlock("SpecialSetups");
   bool tubeSetup = SpecialSetups.getParameter< bool >("tubeSetup", false);
   bool slopedWall = SpecialSetups.getParameter< bool >("slopedWall", false);


   const FlagUID fluidFlagUID("Fluid");
   const FlagUID wallFlagUID("FreeSlip");

   auto boundariesConfig = walberlaEnv.config()->getOneBlock("Boundaries");

   lbm::GeneratedFreeSlip_NoSlip noSlip(blocks, pdfFieldID);
   lbm::GeneratedFreeSlip_FreeSlip freeSlip(blocks, pdfFieldID);

   geometry::initBoundaryHandling< FlagField_T >(*blocks, flagFieldId, boundariesConfig);

   if (tubeSetup || slopedWall)
   {
      for (auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt)
      {
         FlagField_T* flagField = blockIt->template getData< FlagField_T >(flagFieldId);
         auto wallFlag          = flagField->getOrRegisterFlag(wallFlagUID);

         auto My = blocks->getDomainCellBB().yMax() / 2.0;
         auto Mz = blocks->getDomainCellBB().zMax() / 2.0;

         WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ(flagField, {
            Cell globalCell = Cell(x, y, z);
            blocks->transformBlockLocalToGlobalCell(globalCell, *blockIt);

            if (tubeSetup)
            {
               real_t R = sqrt((globalCell[1] - My) * (globalCell[1] - My) + (globalCell[2] - Mz) * (globalCell[2] - Mz));
               if (R > real_c(cellsPerBlock[1]) * real_c(0.5)) addFlag(flagField->get(x, y, z), wallFlag);
            }
            if (slopedWall)
            {
               if (globalCell[2] - globalCell[1] > 0)
               {
                  addFlag(flagField->get(x, y, z), wallFlag);
               }
            }
         }) // WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ
      }
   }
   geometry::setNonBoundaryCellsToDomain< FlagField_T >(*blocks, flagFieldId, fluidFlagUID);

   noSlip.fillFromFlagField< FlagField_T >(blocks, flagFieldId, FlagUID("NoSlip"), fluidFlagUID);
   freeSlip.fillFromFlagField< FlagField_T >(blocks, flagFieldId, FlagUID("FreeSlip"), fluidFlagUID);

   // create time loop
   SweepTimeloop timeloop(blocks->getBlockStorage(), timesteps);

   // create communication for PdfField
   blockforest::communication::UniformBufferedScheme< Stencil_T > evenComm(blocks);
   evenComm.addPackInfo(make_shared< PackInfoEven_T >(pdfFieldID));

   blockforest::communication::UniformBufferedScheme< Stencil_T > oddComm(blocks);
   oddComm.addPackInfo(make_shared< PackInfoOdd_T >(pdfFieldID));

   // Timestep Tracking and Sweeps
   auto tracker = make_shared< lbm::TimestepTracker >(0);

   AlternatingBeforeFunction communication(evenComm, oddComm, tracker);

   lbm::GeneratedFreeSlip_Sweep UpdateSweep(densityFieldID, pdfFieldID, velFieldID, omega);

   // add LBM sweep and communication to time loop
   timeloop.add() << BeforeFunction(communication, "communication") << Sweep(noSlip.getSweep(tracker), "noSlip boundary");
   timeloop.add() << Sweep(freeSlip.getSweep(tracker), "freeSlip boundary");
   timeloop.add() << BeforeFunction(tracker->getAdvancementFunction(), "Timestep Advancement")
                  << Sweep(UpdateSweep.getSweep(tracker), "LB update rule");

   // log remaining time
   timeloop.addFuncAfterTimeStep(timing::RemainingTimeLogger(timeloop.getNrOfTimeSteps(), remainingTimeLoggerFrequency),
                                 "remaining time logger");

   // VTK Writer
   uint_t vtkWriteFrequency = parameters.getParameter< uint_t >("vtkWriteFrequency", 0);
   if (vtkWriteFrequency > 0)
   {
      auto vtkOutput = vtk::createVTKOutput_BlockData(*blocks, "GeneratedFreeSlip_VTK", vtkWriteFrequency, 1, false,
                                                      "vtk_out", "simulation_step", false, true, true, false, 0);

      auto velWriter     = make_shared< field::VTKWriter< VelocityField_T > >(velFieldID, "velocity");
      auto densityWriter = make_shared< field::VTKWriter< ScalarField_T > >(densityFieldID, "density");

      auto flagWriter = make_shared< field::VTKWriter< FlagField_T > >(flagFieldId, "flag");
      vtkOutput->addCellDataWriter(flagWriter);

      vtkOutput->addCellDataWriter(velWriter);
      vtkOutput->addCellDataWriter(densityWriter);

      timeloop.addFuncBeforeTimeStep(vtk::writeFiles(vtkOutput), "VTK Output");
   }
   timeloop.run();

   // simple check is only valid for tube
   if (tubeSetup)
   {
      for (auto& block : *blocks)
      {
         auto velField = block.getData< VelocityField_T >(velFieldID);
         for (cell_idx_t i = 0; i < cell_idx_c(cellsPerBlock[1] - 1); ++i)
         {
            real_t v1   = velField->get(cell_idx_c(cellsPerBlock[0] / 2), i, cell_idx_c(cellsPerBlock[2] / 2), 0);
            real_t v2   = velField->get(cell_idx_c(cellsPerBlock[0] / 2), i + 1, cell_idx_c(cellsPerBlock[2] / 2), 0);
            real_t grad = v2 - v1;
            // WALBERLA_LOG_DEVEL_VAR(grad)
            WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(grad, 0.0, 1e-16)
         }
      }
   }
   return EXIT_SUCCESS;
}
