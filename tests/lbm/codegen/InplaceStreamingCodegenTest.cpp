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
//! \file InplaceStreamingCodegenTest.cpp
//! \author Frederik Hennig <frederik.hennig@fau.de>
//
//======================================================================================================================

#include "InplaceStreamingCodegen.h"

#include "blockforest/all.h"

#include "core/Environment.h"

#include "field/all.h"

#include "geometry/InitBoundaryHandling.h"

#include "lbm/inplace_streaming/TimestepTracker.h"

#include "timeloop/SweepTimeloop.h"


using namespace walberla::lbmpy;

namespace walberla
{
typedef uint8_t boundaryFlag_t;
typedef FlagField< boundaryFlag_t > FlagField_T;

int main(int argc, char** argv)
{
   Environment walberlaEnv(argc, argv);
   auto blocks = blockforest::createUniformBlockGridFromConfig(walberlaEnv.config());

   auto params              = walberlaEnv.config()->getBlock("Parameters");
   uint_t timesteps         = params.getParameter< uint_t >("timesteps");
   uint_t vtkWriteFrequency = params.getParameter< uint_t >("vtkWriteFrequency");

   BlockDataID pullPdfFieldID    = field::addToStorage< PdfField_T >(blocks, "f_pull", real_c(0.0), field::fzyx);
   BlockDataID inplacePdfFieldID = field::addToStorage< PdfField_T >(blocks, "f_inplace", real_c(0.0), field::fzyx);

   BlockDataID pullVelocityFieldID = field::addToStorage< VelocityField_T >(blocks, "u_pull", real_c(0.0), field::fzyx);
   BlockDataID inplaceVelocityFieldID =
      field::addToStorage< VelocityField_T >(blocks, "u_inplace", real_c(0.0), field::fzyx);

   BlockDataID boundaryFlagFieldID = field::addFlagFieldToStorage< FlagField_T >(blocks, "boundaryFlagField", 1);

   // Set up boundary handling

   const FlagUID fluidFlagUID("Fluid");
   const FlagUID noslipFlagUID("NoSlip");
   const FlagUID ubbFlagUID("UBB");
   const FlagUID outflowFlagUID("Outflow");

   auto boundariesConfig = walberlaEnv.config()->getOneBlock("Boundaries");

   geometry::initBoundaryHandling< FlagField_T >(*blocks, boundaryFlagFieldID, boundariesConfig);
   geometry::setNonBoundaryCellsToDomain< FlagField_T >(*blocks, boundaryFlagFieldID, fluidFlagUID);

   PullNoSlip pullNoSlip(blocks, pullPdfFieldID);
   pullNoSlip.fillFromFlagField< FlagField_T >(blocks, boundaryFlagFieldID, noslipFlagUID, fluidFlagUID);

   PullUBB pullUBB(blocks, pullPdfFieldID);
   pullUBB.fillFromFlagField< FlagField_T >(blocks, boundaryFlagFieldID, ubbFlagUID, fluidFlagUID);

   PullOutflow pullOutflow(blocks, pullPdfFieldID);
   pullOutflow.fillFromFlagField< FlagField_T >(blocks, boundaryFlagFieldID, outflowFlagUID, fluidFlagUID);

   InPlaceNoSlip inplaceNoSlip(blocks, inplacePdfFieldID);
   inplaceNoSlip.fillFromFlagField< FlagField_T >(blocks, boundaryFlagFieldID, noslipFlagUID, fluidFlagUID);

   InPlaceUBB inplaceUBB(blocks, inplacePdfFieldID);
   inplaceUBB.fillFromFlagField< FlagField_T >(blocks, boundaryFlagFieldID, ubbFlagUID, fluidFlagUID);

   InPlaceOutflow inplaceOutflow(blocks, inplacePdfFieldID);
   inplaceOutflow.fillFromFlagField< FlagField_T >(blocks, boundaryFlagFieldID, outflowFlagUID, fluidFlagUID);

   // Init fields

   PullInit pullInit(pullPdfFieldID);
   InPlaceInit inplaceInit(inplacePdfFieldID);

   for (auto b = blocks->begin(); b != blocks->end(); b++)
   {
      pullInit(&(*b));
      inplaceInit(&(*b));
   }

   // Set up sweeps

   PullSweep pullLbmSweep(pullPdfFieldID, pullVelocityFieldID);
   InPlaceSweep inplaceLbmSweep(inplacePdfFieldID, inplaceVelocityFieldID);

   auto tracker = make_shared< lbm::TimestepTracker >(1);

   // Pull Timeloop

   SweepTimeloop pullTimeloop(blocks, timesteps);
   pullTimeloop.add() << Sweep(pullUBB);
   pullTimeloop.add() << Sweep(pullNoSlip);
   pullTimeloop.add() << Sweep(pullOutflow);
   pullTimeloop.add() << Sweep(pullLbmSweep);

   if (vtkWriteFrequency > 0)
   {
      pullTimeloop.addFuncAfterTimeStep(field::createVTKOutput< VelocityField_T, float >(
         pullVelocityFieldID, *blocks, "pull velocity", vtkWriteFrequency, uint_c(0), false, "pull_vtk"));
   }

   pullTimeloop.run();

   // In-Place Timeloop

   SweepTimeloop inplaceTimeloop(blocks, timesteps);
   inplaceTimeloop.add() << Sweep(inplaceUBB.getSweep(tracker));
   inplaceTimeloop.add() << Sweep(inplaceNoSlip.getSweep(tracker));
   inplaceTimeloop.add() << Sweep(inplaceOutflow.getSweep(tracker)) << AfterFunction(tracker->getAdvancementFunction());
   inplaceTimeloop.add() << Sweep(inplaceLbmSweep.getSweep(tracker));

   if (vtkWriteFrequency > 0)
   {
      inplaceTimeloop.addFuncAfterTimeStep(field::createVTKOutput< VelocityField_T, float >(
         inplaceVelocityFieldID, *blocks, "inplace velocity", vtkWriteFrequency, uint_c(0), false, "inplace_vtk"));
   }

   inplaceTimeloop.run();

   // Validate results

   for (auto b = blocks->begin(); b != blocks->end(); b++)
   {
      auto pullVelocityField    = b->getData< VelocityField_T >(pullVelocityFieldID);
      auto inplaceVelocityField = b->getData< VelocityField_T >(inplaceVelocityFieldID);

      CellInterval fieldSize = pullVelocityField->xyzSize();
      for (Cell c : fieldSize)
      {
         for (uint_t i = 0; i < pullVelocityField->fSize(); i++)
         {
            WALBERLA_CHECK_FLOAT_EQUAL(pullVelocityField->get(c, i), inplaceVelocityField->get(c, i))
         }
      }
   }

   return 0;
}
} // namespace walberla

int main(int argc, char** argv) { return walberla::main(argc, argv); }
