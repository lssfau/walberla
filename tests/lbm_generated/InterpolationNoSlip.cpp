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
//! \file InterpolatedNoSlip.cpp
//! \author Markus Holzer <markus.holzer@fau.de>
//! \brief Couette flow driven by a UBB BC in Nord and wall boundary in the South. Remaining directions are periodic
//!        If Interpolation BC are used the distance of the southern wall can be controlled. The velocity in the
//!        first fluid cell is checked and compared with the velocity obtained from a default NoSlip BC.
//!        Depending on the set distance for the interpolation BCs the velocity is expected to be higher or lower
//
//======================================================================================================================
#include "blockforest/Initialization.h"
#include "blockforest/communication/UniformBufferedScheme.h"

#include "core/DataTypes.h"
#include "core/Environment.h"
#include "core/debug/TestSubsystem.h"
#include "core/logging/Initialization.h"
#include "core/math/Vector3.h"
#include "core/timing/RemainingTimeLogger.h"

#include "field/AddToStorage.h"
#include "field/FlagField.h"
#include "field/GhostLayerField.h"
#include "field/vtk/VTKWriter.h"

#include "geometry/InitBoundaryHandling.h"

#include "timeloop/SweepTimeloop.h"

#include "lbm_generated/communication/UniformGeneratedPdfPackInfo.h"
#include "lbm_generated/field/AddToStorage.h"
#include "lbm_generated/field/PdfField.h"

// include the generated header file. It includes all generated classes
#include "InterpolationNoSlipHeader.h"

using namespace walberla;
using namespace std::placeholders;

using StorageSpecification_T = lbm::InterpolationNoSlipStorageSpecification;
using Stencil_T              = StorageSpecification_T::Stencil;
using CommunicationStencil_T = StorageSpecification_T::CommunicationStencil;
using PdfField_T             = lbm_generated::PdfField< StorageSpecification_T >;
using PackInfo_T             = lbm_generated::UniformGeneratedPdfPackInfo< PdfField_T >;

using SweepCollection_T = lbm::InterpolationNoSlipSweepCollection;

using VectorField_T = GhostLayerField< real_t, StorageSpecification_T::Stencil::D >;
using ScalarField_T = GhostLayerField< real_t, 1 >;

using flag_t               = walberla::uint8_t;
using FlagField_T          = FlagField< flag_t >;
using BoundaryCollection_T = lbm::InterpolationNoSlipBoundaryCollection< FlagField_T >;

using blockforest::communication::UniformBufferedScheme;

class wallDistance
{
 public:
   wallDistance(const real_t q) : q_(q) {}

   real_t operator()(const Cell& fluidCell, const Cell& boundaryCell, const shared_ptr< StructuredBlockForest >& SbF,
                     IBlock& block) const;

 private:
   const real_t q_;
}; // class wallDistance

real_t wallDistance::operator()(const Cell& /*fluidCell*/, const Cell& /*boundaryCell*/,
                                const shared_ptr< StructuredBlockForest >& /*SbF*/, IBlock& /*block*/) const
{
   return q_;
}

int main(int argc, char** argv)
{
   debug::enterTestMode();
   walberla::Environment walberlaEnv(argc, argv);
   logging::configureLogging(walberlaEnv.config());

   auto blocks = blockforest::createUniformBlockGridFromConfig(walberlaEnv.config());

   auto domainSetup                = walberlaEnv.config()->getOneBlock("DomainSetup");
   Vector3< uint_t > cellsPerBlock = domainSetup.getParameter< Vector3< uint_t > >("cellsPerBlock");

   // read parameters
   auto parameters   = walberlaEnv.config()->getOneBlock("Parameters");
   const real_t omega        = parameters.getParameter< real_t >("omega", real_c(1.4));
   const real_t distanceWall = parameters.getParameter< real_t >("distanceWall", real_c(0.5));
   const uint_t timesteps    = parameters.getParameter< uint_t >("timesteps", uint_c(10)) + uint_c(1);

   WALBERLA_LOG_DEVEL_VAR(distanceWall)

   auto remainingTimeLoggerFrequency =
      parameters.getParameter< real_t >("remainingTimeLoggerFrequency", real_c(3.0)); // in seconds

   const StorageSpecification_T StorageSpec = StorageSpecification_T();
   BlockDataID pdfFieldId = lbm_generated::addPdfFieldToStorage(blocks, "pdf field", StorageSpec, uint_c(1));
   BlockDataID velFieldId = field::addToStorage< VectorField_T >(blocks, "Velocity", real_c(0.0), field::fzyx);

   BlockDataID flagFieldId = field::addFlagFieldToStorage< FlagField_T >(blocks, "flag field", uint_c(1));

   SweepCollection_T sweepCollection(blocks, pdfFieldId, velFieldId, omega);
   for (auto& block : *blocks)
   {
      sweepCollection.initialise(&block);
   }

   const FlagUID fluidFlagUID("Fluid");
   auto boundariesConfig = walberlaEnv.config()->getBlock("Boundaries");
   geometry::initBoundaryHandling< FlagField_T >(*blocks, flagFieldId, boundariesConfig);
   geometry::setNonBoundaryCellsToDomain< FlagField_T >(*blocks, flagFieldId, fluidFlagUID);

   const wallDistance wallDistanceCallback{ distanceWall };
   std::function< real_t(const Cell&, const Cell&, const shared_ptr< StructuredBlockForest >&, IBlock&) >
      wallDistanceFunctor = wallDistanceCallback;
   // For the BoundaryCollection a funcotr to calculate the wall distance for the Bouzidi NoSlip and for the QuadraticBB
   // have to be provided. In this test case we use the same function to calculate the wall distance
   BoundaryCollection_T boundaryCollection(blocks, flagFieldId, pdfFieldId, fluidFlagUID, omega, wallDistanceFunctor,
                                           wallDistanceFunctor);

   auto packInfo = std::make_shared< lbm_generated::UniformGeneratedPdfPackInfo< PdfField_T > >(pdfFieldId);
   UniformBufferedScheme< Stencil_T > communication(blocks);
   communication.addPackInfo(packInfo);

   SweepTimeloop timeloop(blocks->getBlockStorage(), timesteps);
   timeloop.add() << BeforeFunction(communication, "communication")
                  << Sweep(boundaryCollection.getSweep(BoundaryCollection_T::ALL), "Boundary Conditions");
   timeloop.add() << Sweep(sweepCollection.streamCollide(SweepCollection_T::ALL), "LBM StreamCollide");

   uint_t vtkWriteFrequency = parameters.getParameter< uint_t >("vtkWriteFrequency", 0);
   if (vtkWriteFrequency > 0)
   {
      auto vtkOutput = vtk::createVTKOutput_BlockData(*blocks, "InterpolationNoSlipVTK", vtkWriteFrequency, 0, false,
                                                      "vtk_out", "simulation_step", false, true, true, false, 0);

      auto velWriter = make_shared< field::VTKWriter< VectorField_T > >(velFieldId, "velocity");
      vtkOutput->addBeforeFunction([&]() {
         for (auto& block : *blocks)
         {
            sweepCollection.calculateMacroscopicParameters(&block);
         }
      });

      vtkOutput->addCellDataWriter(velWriter);
      timeloop.addFuncBeforeTimeStep(vtk::writeFiles(vtkOutput), "VTK Output");
   }

   if (remainingTimeLoggerFrequency > 0)
   {
      // log remaining time
      timeloop.addFuncAfterTimeStep(
         timing::RemainingTimeLogger(timeloop.getNrOfTimeSteps(), remainingTimeLoggerFrequency),
         "remaining time logger");
   }

   timeloop.run();

   // This is the velocity at the wall, when a NoSlip BC is used. This is similar to using an interpolation BC with a
   // wall distance of 0.5. This value can be obtained by either setting distanceWall to 0.5 in the Parameter file or
   // specifying the NoSlip BC at the southern boundary
   const real_t defaultNoSlipVelocity = real_c(0.0002);

   for (auto& block : *blocks)
   {
      sweepCollection.calculateMacroscopicParameters(&block);

      auto velField  = block.getData< VectorField_T >(velFieldId);
      auto velAtWall = velField->get(cell_idx_c(cellsPerBlock[0] / 2), 0, cell_idx_c(cellsPerBlock[2] / 2), 0);
      // WALBERLA_LOG_DEVEL_VAR(velAtWall)

      if (distanceWall > 0.49 && distanceWall < 0.51) { WALBERLA_CHECK_FLOAT_EQUAL(velAtWall, defaultNoSlipVelocity) }
      else if (distanceWall < 0.49) { WALBERLA_CHECK_GREATER(defaultNoSlipVelocity, velAtWall) }
      else if (distanceWall > 0.51) { WALBERLA_CHECK_LESS(defaultNoSlipVelocity, velAtWall) }
   }

   return EXIT_SUCCESS;
}
