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
//! \file WallFunctionBounceTest.cpp
//! \author Helen Schottenhamml <helen.schottenhamml@ifpen.fr>
//! \brief Simple integration test for generated wall-function bounce (WFB). Unit tests that include physical
//!        feasibility checks are included in lbmpy; here, we only check for the successful integration in the
//!        waLBerla framework.
//
//======================================================================================================================
#include "core/math/AABB.h"

#include "core/Environment.h"
#include "core/timing/RemainingTimeLogger.h"

#include "field/AddToStorage.h"
#include "field/vtk/VTKWriter.h"

#include "geometry/InitBoundaryHandling.h"
#include "timeloop/SweepTimeloop.h"

// Generated Files
#include "WFBHeader.h"
#include "blockforest/Initialization.h"
#include "field/StabilityChecker.h"

using namespace walberla;

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

int main(int argc, char** argv)
{
   walberla::Environment walberlaEnv(argc, argv);

   const FlagUID fluidFlagUID("Fluid");
   const FlagUID wfbFlagUID("WFB");

   auto globalConfig = std::make_shared<config::Config>();
   auto & globalBlock = globalConfig->getWritableGlobalBlock();
   { // set up domain configuration
      auto & domainBlock = globalBlock.createBlock("DomainSetup");
      domainBlock.addParameter("blocks", "<1, 1, 1>");
      domainBlock.addParameter("cellsPerBlock", "<20, 10, 10>");
      domainBlock.addParameter("periodic", "<1, 0, 1>");
   }
   { // set up boundary configuration
      auto & boundariesBlock = globalBlock.createBlock("Boundaries");
      auto & wfbBlock = boundariesBlock.createBlock("Border");
      wfbBlock.addParameter("direction", "B");
      wfbBlock.addParameter("walldistance", "-1");
      wfbBlock.addParameter("flag", wfbFlagUID.toString());
   }
   { // set up stability checker config
      auto & checkerBlocker = globalBlock.createBlock("StabilityChecker");
      checkerBlocker.addParameter("checkFrequency", "1");
      checkerBlocker.addParameter("streamOutput", "1");
      checkerBlocker.addParameter("vtkOutput", "0");
   }

   auto blocks = blockforest::createUniformBlockGridFromConfig(globalConfig);

   // dummy values without physical importance
   const real_t omega{1.4};
   [[maybe_unused]] const real_t filterWidth{1e-3};
   const real_t targetUTau{1e-5};

   // create fields
   BlockDataID pdfFieldID     = blocks->addStructuredBlockData< PdfField_T >(pdfFieldAdder, "PDFs");
   BlockDataID velFieldID     = field::addToStorage< VelocityField_T >(blocks, "velocity", real_c(0.0), field::fzyx);

   BlockDataID flagFieldId = field::addFlagFieldToStorage< FlagField_T >(blocks, "flag field");

   auto boundariesConfig = globalConfig->getOneBlock("Boundaries");

#if defined FILTERED_VELOCITY
   lbm::WFB wfb(blocks, pdfFieldID, velFieldID, filterWidth, omega, targetUTau);
#else
   lbm::WFB wfb(blocks, pdfFieldID, velFieldID, omega, targetUTau);
#endif

   geometry::initBoundaryHandling< FlagField_T >(*blocks, flagFieldId, boundariesConfig);
   geometry::setNonBoundaryCellsToDomain< FlagField_T >(*blocks, flagFieldId, fluidFlagUID);

   wfb.fillFromFlagField< FlagField_T >(blocks, flagFieldId, wfbFlagUID, fluidFlagUID);

   // execute WFB once to ensure that it does something
   for (auto & block: *blocks)
      wfb(&block);

   // stability checker
   const auto stabilityChecker = field::makeStabilityChecker< PdfField_T >(blocks, pdfFieldID, 1, true, false, Set<SUID>::emptySet() );
   stabilityChecker->operator()();

   return EXIT_SUCCESS;
}
