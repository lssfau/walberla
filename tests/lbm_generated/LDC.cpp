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
//! \file LDC.cpp
//! \author Markus Holzer <markus.holzer@fau.de>
//
//======================================================================================================================

#include "blockforest/all.h"
#include "blockforest/communication/UniformBufferedScheme.h"

#include "core/all.h"

#include "domain_decomposition/all.h"

#include "field/all.h"
#include "geometry/all.h"
#include "timeloop/all.h"

#include "lbm_generated/communication/UniformGeneratedPdfPackInfo.h"
#include "lbm_generated/field/AddToStorage.h"
#include "lbm_generated/field/PdfField.h"

#include "lbm_generated/storage_specification/D3Q19StorageSpecification.h"
#include "lbm_generated/sweep_collection/D3Q19SRT.h"
#include "lbm_generated/boundary/D3Q19BoundaryCollection.h"


using namespace walberla;
using namespace std::placeholders;

using StorageSpecification_T = lbm::D3Q19StorageSpecification;
using Stencil_T              = StorageSpecification_T::Stencil;
using CommunicationStencil_T = StorageSpecification_T::CommunicationStencil;
using PdfField_T             = lbm_generated::PdfField< StorageSpecification_T >;

using SweepCollection_T = lbm::D3Q19SRT;

using VectorField_T = GhostLayerField< real_t, StorageSpecification_T::Stencil::D >;
using ScalarField_T = GhostLayerField< real_t, 1 >;

using flag_t      = walberla::uint8_t;
using FlagField_T = FlagField< flag_t >;
using BoundaryCollection_T = lbm::D3Q19BoundaryCollection< FlagField_T >;

using blockforest::communication::UniformBufferedScheme;

int main(int argc, char** argv)
{
   walberla::Environment walberlaEnv(argc, argv);
   mpi::MPIManager::instance()->useWorldComm();

   // read parameters
   auto parameters = walberlaEnv.config()->getOneBlock("Parameters");

   const real_t omega     = parameters.getParameter< real_t >("omega", real_c(1.4));
   const uint_t timesteps = parameters.getParameter< uint_t >("timesteps", uint_c(10)) + uint_c(1);

   const double remainingTimeLoggerFrequency =
      parameters.getParameter< double >("remainingTimeLoggerFrequency", real_c(3.0)); // in seconds

   auto blocks = blockforest::createUniformBlockGridFromConfig(walberlaEnv.config());

   StorageSpecification_T const StorageSpec = StorageSpecification_T();
   BlockDataID const pdfFieldId  = lbm_generated::addPdfFieldToStorage(blocks, "pdf field", StorageSpec, uint_c(1), field::fzyx);
   BlockDataID const velFieldId  = field::addToStorage< VectorField_T >(blocks, "Velocity", real_c(0.0), field::fzyx);
   BlockDataID const densityFieldId  = field::addToStorage< ScalarField_T >(blocks, "density", real_c(0.0), field::fzyx);
   BlockDataID const flagFieldId = field::addFlagFieldToStorage< FlagField_T >(blocks, "flag field", uint_c(1));

   const FlagUID fluidFlagUID("Fluid");

   auto boundariesConfig   = walberlaEnv.config()->getBlock("Boundaries");
   geometry::initBoundaryHandling< FlagField_T >(*blocks, flagFieldId, boundariesConfig);
   geometry::setNonBoundaryCellsToDomain< FlagField_T >(*blocks, flagFieldId, fluidFlagUID);

   BoundaryCollection_T boundaryCollection(blocks, flagFieldId, pdfFieldId, fluidFlagUID, real_c(1.0), real_c(0.05), real_c(0.0), real_c(0.0));
   SweepCollection_T sweepCollection(blocks, pdfFieldId, densityFieldId, velFieldId, omega);

   for (auto& block : *blocks)
   {
      sweepCollection.initialise(&block);
   }

   auto packInfo = std::make_shared<lbm_generated::UniformGeneratedPdfPackInfo< PdfField_T >>(pdfFieldId);
   UniformBufferedScheme< Stencil_T > communication(blocks);
   communication.addPackInfo(packInfo);

   SweepTimeloop timeLoop(blocks->getBlockStorage(), timesteps);

   timeLoop.add() << BeforeFunction(communication, "communication")
                  << Sweep(boundaryCollection.getSweep(BoundaryCollection_T::ALL), "Boundary Conditions");
   timeLoop.add() << Sweep(sweepCollection.streamCollide(SweepCollection_T::ALL), "LBM StreamCollide");
   //
   auto vtkWriteFrequency = parameters.getParameter< uint_t >("vtkWriteFrequency", 0);
   if (vtkWriteFrequency > 0)
   {
      auto vtkOutput = vtk::createVTKOutput_BlockData(*blocks, "ExampleVTK", vtkWriteFrequency, 0, false, "vtk_out",
                                                      "simulation_step", false, true, true, false, 0);

      auto velWriter = make_shared< field::VTKWriter< VectorField_T > >(velFieldId, "velocity");
      auto densWriter = make_shared< field::VTKWriter< ScalarField_T > >(densityFieldId, "density");
      vtkOutput->addBeforeFunction([&](){
      for (auto& block : *blocks)
      {
         sweepCollection.calculateMacroscopicParameters(&block);
      }
      });

      vtkOutput->addCellDataWriter(velWriter);
      vtkOutput->addCellDataWriter(densWriter);

      timeLoop.addFuncBeforeTimeStep(vtk::writeFiles(vtkOutput), "VTK Output");
   }

   // log remaining time
   timeLoop.addFuncAfterTimeStep(timing::RemainingTimeLogger(timeLoop.getNrOfTimeSteps(), remainingTimeLoggerFrequency),
                                 "remaining time logger");

   WALBERLA_LOG_INFO_ON_ROOT("Starting Simulation with " << timesteps << " timesteps")

   timeLoop.run();

   return EXIT_SUCCESS;
}
