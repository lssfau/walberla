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
//! \file ExcessMassDistributionFallbackTest.cpp
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//! \brief Test excess mass distribution in cases where the chosen model is not applicable.
//!
//! Test if fall back to other models works correctly with a simple two-dimensional 3x3 grid, where the center cell at
//! (1,1) is assumed to have converted from interface to liquid with excess mass 0.1.
//======================================================================================================================

#include "blockforest/Initialization.h"

#include "core/Environment.h"
#include "core/debug/TestSubsystem.h"

#include "field/FieldClone.h"

#include "lbm/blockforest/communication/SimpleCommunication.h"
#include "lbm/field/Adaptors.h"
#include "lbm/field/AddToStorage.h"
#include "lbm/free_surface/boundary/FreeSurfaceBoundaryHandling.h"
#include "lbm/free_surface/bubble_model/Geometry.h"
#include "lbm/free_surface/dynamics/ExcessMassDistributionModel.h"
#include "lbm/free_surface/dynamics/ExcessMassDistributionSweep.h"
#include "lbm/lattice_model/D2Q9.h"

#include "timeloop/SweepTimeloop.h"

namespace walberla
{
namespace free_surface
{
namespace ExcessMassDistributionFallbackTest
{
using LatticeModel_T = lbm::D2Q9< lbm::collision_model::SRT, true, lbm::force_model::None, 2 >;
using Stencil        = typename LatticeModel_T::Stencil;

using Communication_T = blockforest::SimpleCommunication< LatticeModel_T::CommunicationStencil >;

using PdfField_T    = lbm::PdfField< LatticeModel_T >;
using ScalarField_T = GhostLayerField< real_t, 1 >;
using VectorField_T = GhostLayerField< Vector3< real_t >, 1 >;

using flag_t                        = uint32_t;
using FlagField_T                   = FlagField< flag_t >;
using FreeSurfaceBoundaryHandling_T = FreeSurfaceBoundaryHandling< LatticeModel_T, FlagField_T, ScalarField_T >;

void runSimulation(const ExcessMassDistributionModel& excessMassDistributionModel,
                   const std::vector< Cell >& newInterfaceCells, const std::vector< Cell >& oldInterfaceCells,
                   const Vector3< real_t >& interfaceNormal)
{
   // define the domain size
   const Vector3< uint_t > numBlocks(uint_c(1));
   const Vector3< uint_t > domainSize(uint_c(3), uint_c(3), uint_c(1));
   const Vector3< uint_t > cellsPerBlock(domainSize[0] / numBlocks[0], domainSize[1] / numBlocks[1],
                                         domainSize[2] / numBlocks[2]);

   // create block grid
   const std::shared_ptr< StructuredBlockForest > blockForest =
      blockforest::createUniformBlockGrid(numBlocks[0], numBlocks[1], numBlocks[2],             // blocks
                                          cellsPerBlock[0], cellsPerBlock[1], cellsPerBlock[2], // cells
                                          real_c(1.0),                                          // dx
                                          true,                                                 // one block per process
                                          false, false, false);                                 // periodicity

   // create (dummy) lattice model
   LatticeModel_T latticeModel = LatticeModel_T(lbm::collision_model::SRT(real_c(1.8)));

   // add pdf field
   const BlockDataID pdfFieldID = lbm::addPdfFieldToStorage(blockForest, "PDF field", latticeModel,
                                                            Vector3< real_t >(real_c(0)), real_c(1), field::fzyx);

   // add normal field
   const BlockDataID normalFieldID = field::addToStorage< VectorField_T >(
      blockForest, "Normal field", Vector3< real_t >(real_c(0)), field::fzyx, uint_c(1));

   // add fill level field
   const BlockDataID fillFieldID =
      field::addToStorage< ScalarField_T >(blockForest, "Fill level field", real_c(1), field::fzyx, uint_c(1));

   // initialize fill levels and interface normal
   for (auto blockIt = blockForest->begin(); blockIt != blockForest->end(); ++blockIt)
   {
      ScalarField_T* const fillField   = blockIt->getData< ScalarField_T >(fillFieldID);
      VectorField_T* const normalField = blockIt->getData< VectorField_T >(normalFieldID);

      WALBERLA_FOR_ALL_CELLS(fillFieldIt, fillField, normalFieldIt, normalField, {
         // initialize interface cells
         for (const Cell& cell : newInterfaceCells)
         {
            if (fillFieldIt.cell() == cell) { *fillFieldIt = real_c(0.5); }
         }

         for (const Cell& cell : oldInterfaceCells)
         {
            if (fillFieldIt.cell() == cell) { *fillFieldIt = real_c(0.5); }
         }

         // this cell is assigned a fill level of 1.1 leading to an excess mass equivalent to a fill level of 0.1
         if (fillFieldIt.cell() == Cell(cell_idx_c(1), cell_idx_c(1), cell_idx_c(0))) { *fillFieldIt = real_c(1.1); }

         if (normalFieldIt.cell() == Cell(cell_idx_c(1), cell_idx_c(1), cell_idx_c(0)))
         {
            *normalFieldIt = interfaceNormal;
         }
      }) // WALBERLA_FOR_ALL_CELLS
   }

   // add boundary handling
   const std::shared_ptr< FreeSurfaceBoundaryHandling_T > freeSurfaceBoundaryHandling =
      std::make_shared< FreeSurfaceBoundaryHandling_T >(blockForest, pdfFieldID, fillFieldID);
   const BlockDataID flagFieldID                                      = freeSurfaceBoundaryHandling->getFlagFieldID();
   const typename FreeSurfaceBoundaryHandling_T::FlagInfo_T& flagInfo = freeSurfaceBoundaryHandling->getFlagInfo();
   freeSurfaceBoundaryHandling->initFlagsFromFillLevel();

   // initialize flags
   for (auto blockIt = blockForest->begin(); blockIt != blockForest->end(); ++blockIt)
   {
      FlagField_T* const flagField = blockIt->getData< FlagField_T >(flagFieldID);

      WALBERLA_FOR_ALL_CELLS(flagFieldIt, flagField, {
         // flag the cell with excess mass as newly converted to liquid
         if (flagFieldIt.cell() == Cell(cell_idx_c(1), cell_idx_c(1), cell_idx_c(0)))
         {
            field::addFlag(flagFieldIt, flagInfo.convertedFlag | flagInfo.convertToLiquidFlag);
         }

         // consider these cells to be newly-converted to interface
         for (const Cell& cell : newInterfaceCells)
         {
            if (flagFieldIt.cell() == cell) { field::addFlag(flagFieldIt, flagInfo.convertedFlag); }
         }
      }) // WALBERLA_FOR_ALL_CELLS
   }

   // initial communication
   Communication_T(blockForest, pdfFieldID, fillFieldID, flagFieldID, normalFieldID)();

   // create timeloop
   SweepTimeloop timeloop(blockForest, uint_c(1));

   if (excessMassDistributionModel.isEvenlyType())
   {
      const ExcessMassDistributionSweepInterfaceEvenly< LatticeModel_T, FlagField_T, ScalarField_T, VectorField_T >
         distributeMassSweep(excessMassDistributionModel, fillFieldID, flagFieldID, pdfFieldID, flagInfo);
      timeloop.add() << Sweep(distributeMassSweep, "Distribute excess mass");
   }
   else
   {
      if (excessMassDistributionModel.isWeightedType())
      {
         const ExcessMassDistributionSweepInterfaceWeighted< LatticeModel_T, FlagField_T, ScalarField_T, VectorField_T >
            distributeMassSweep(excessMassDistributionModel, fillFieldID, flagFieldID, pdfFieldID, flagInfo,
                                normalFieldID);
         timeloop.add() << Sweep(distributeMassSweep, "Distribute excess mass");
      }
   }

   timeloop.singleStep();

   // check if excess mass was distributed correctly; expected solutions were obtained manually
   for (auto blockIt = blockForest->begin(); blockIt != blockForest->end(); ++blockIt)
   {
      const ScalarField_T* const fillField = blockIt->getData< const ScalarField_T >(fillFieldID);

      WALBERLA_FOR_ALL_CELLS(fillFieldIt, fillField, {
         if (excessMassDistributionModel.getModelType() ==
             ExcessMassDistributionModel::ExcessMassModel::EvenlyOldInterface)
         {
            if (fillFieldIt.cell() == Cell(cell_idx_c(1), cell_idx_c(1), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(1), real_c(1e-4));
            }

            if (fillFieldIt.cell() == Cell(cell_idx_c(0), cell_idx_c(0), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(0.6), real_c(1e-4));
            }
         }

         if (excessMassDistributionModel.getModelType() ==
             ExcessMassDistributionModel::ExcessMassModel::EvenlyNewInterface)
         {
            if (fillFieldIt.cell() == Cell(cell_idx_c(1), cell_idx_c(1), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(1), real_c(1e-4));
            }

            if (fillFieldIt.cell() == Cell(cell_idx_c(0), cell_idx_c(0), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(0.6), real_c(1e-4));
            }
         }

         if (excessMassDistributionModel.getModelType() ==
                ExcessMassDistributionModel::ExcessMassModel::WeightedOldInterface &&
             !oldInterfaceCells.empty())
         {
            if (fillFieldIt.cell() == Cell(cell_idx_c(1), cell_idx_c(1), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(1), real_c(1e-4));
            }

            if (fillFieldIt.cell() == Cell(cell_idx_c(0), cell_idx_c(0), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(0.5), real_c(1e-4));
            }

            if (fillFieldIt.cell() == Cell(cell_idx_c(2), cell_idx_c(2), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(0.6), real_c(1e-4));
            }
         }

         if (excessMassDistributionModel.getModelType() ==
                ExcessMassDistributionModel::ExcessMassModel::WeightedNewInterface &&
             !newInterfaceCells.empty())
         {
            if (fillFieldIt.cell() == Cell(cell_idx_c(1), cell_idx_c(1), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(1), real_c(1e-4));
            }

            if (fillFieldIt.cell() == Cell(cell_idx_c(0), cell_idx_c(0), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(0.5), real_c(1e-4));
            }

            if (fillFieldIt.cell() == Cell(cell_idx_c(2), cell_idx_c(2), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(0.6), real_c(1e-4));
            }
         }

         if (excessMassDistributionModel.getModelType() ==
                ExcessMassDistributionModel::ExcessMassModel::WeightedOldInterface &&
             oldInterfaceCells.empty())
         {
            if (fillFieldIt.cell() == Cell(cell_idx_c(1), cell_idx_c(1), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(1), real_c(1e-4));
            }

            if (fillFieldIt.cell() == Cell(cell_idx_c(1), cell_idx_c(2), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(0.55), real_c(1e-4));
            }

            if (fillFieldIt.cell() == Cell(cell_idx_c(2), cell_idx_c(2), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(0.55), real_c(1e-4));
            }
         }

         if (excessMassDistributionModel.getModelType() ==
                ExcessMassDistributionModel::ExcessMassModel::WeightedNewInterface &&
             newInterfaceCells.empty())
         {
            if (fillFieldIt.cell() == Cell(cell_idx_c(1), cell_idx_c(1), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(1), real_c(1e-4));
            }

            if (fillFieldIt.cell() == Cell(cell_idx_c(1), cell_idx_c(2), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(0.55), real_c(1e-4));
            }

            if (fillFieldIt.cell() == Cell(cell_idx_c(2), cell_idx_c(2), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(0.55), real_c(1e-4));
            }
         }
      }) // WALBERLA_FOR_ALL_CELLS
   }

   MPIManager::instance()->resetMPI();
}

int main(int argc, char** argv)
{
   debug::enterTestMode();
   Environment env(argc, argv);

   Vector3< real_t > interfaceNormal(real_c(0));

   // (0,0) is the only interface cell (newly-converted) => EvenlyOldInterface must fall back to EvenlyNewInterface
   ExcessMassDistributionModel model = ExcessMassDistributionModel("EvenlyOldInterface");
   std::vector< Cell > newInterfaceCells{ Cell(cell_idx_c(0), cell_idx_c(0), cell_idx_c(0)) };
   std::vector< Cell > oldInterfaceCells{};
   WALBERLA_LOG_INFO_ON_ROOT("Testing model " << model.getFullModelSpecification());
   runSimulation(model, newInterfaceCells, oldInterfaceCells, interfaceNormal);

   // (0,0) is the only interface cell (not newly-converted) => EvenlyNewInterface must fall back to EvenlyOldInterface
   model             = ExcessMassDistributionModel("EvenlyNewInterface");
   newInterfaceCells = std::vector< Cell >{};
   oldInterfaceCells = std::vector< Cell >{ Cell(cell_idx_c(0), cell_idx_c(0), cell_idx_c(0)) };
   WALBERLA_LOG_INFO_ON_ROOT("Testing model " << model.getFullModelSpecification());
   runSimulation(model, newInterfaceCells, oldInterfaceCells, interfaceNormal);

   interfaceNormal = Vector3< real_t >(real_c(0.71), real_c(0.71), real_c(0));

   // (0,0) is old interface cell; (2,2) is newly-converted interface cell; interface normal points in direction (1,1)
   // => WeightedOldInterface must fall back to WeightedNewInterface
   model             = ExcessMassDistributionModel("WeightedOldInterface");
   newInterfaceCells = std::vector< Cell >{ Cell(cell_idx_c(2), cell_idx_c(2), cell_idx_c(0)) };
   oldInterfaceCells = std::vector< Cell >{ Cell(cell_idx_c(0), cell_idx_c(0), cell_idx_c(0)) };
   WALBERLA_LOG_INFO_ON_ROOT("Testing model " << model.getFullModelSpecification());
   runSimulation(model, newInterfaceCells, oldInterfaceCells, interfaceNormal);

   // (0,0) is newly-converted interface cell; (2,2) is old interface cell; interface normal points in direction (1,1)
   // => WeightedNewInterface must fall back to WeightedOldInterface
   model             = ExcessMassDistributionModel("WeightedNewInterface");
   newInterfaceCells = std::vector< Cell >{ Cell(cell_idx_c(0), cell_idx_c(0), cell_idx_c(0)) };
   oldInterfaceCells = std::vector< Cell >{ Cell(cell_idx_c(2), cell_idx_c(2), cell_idx_c(0)) };
   WALBERLA_LOG_INFO_ON_ROOT("Testing model " << model.getFullModelSpecification());
   runSimulation(model, newInterfaceCells, oldInterfaceCells, interfaceNormal);

   interfaceNormal = Vector3< real_t >(real_c(0), real_c(-1), real_c(0));

   // (1,2) and (2,2) are newly-converted interface cells; interface normal points in direction (0,-1)
   // => WeightedOldInterface must fall back to EvenlyAllInterface, as no interface cell is available in normal
   // direction
   model             = ExcessMassDistributionModel("WeightedOldInterface");
   newInterfaceCells = std::vector< Cell >{ Cell(cell_idx_c(1), cell_idx_c(2), cell_idx_c(0)),
                                            Cell(cell_idx_c(2), cell_idx_c(2), cell_idx_c(0)) };
   oldInterfaceCells = std::vector< Cell >{};
   WALBERLA_LOG_INFO_ON_ROOT("Testing model " << model.getFullModelSpecification());
   runSimulation(model, newInterfaceCells, oldInterfaceCells, interfaceNormal);

   // (1,2) and (2,2) are old interface cells; interface normal points in direction (0,-1)
   // => WeightedNewInterface must fall back to EvenlyAllInterface, as no interface cell is available in normal
   // direction
   model             = ExcessMassDistributionModel("WeightedNewInterface");
   newInterfaceCells = std::vector< Cell >{};
   oldInterfaceCells = std::vector< Cell >{ Cell(cell_idx_c(1), cell_idx_c(2), cell_idx_c(0)),
                                            Cell(cell_idx_c(2), cell_idx_c(2), cell_idx_c(0)) };
   WALBERLA_LOG_INFO_ON_ROOT("Testing model " << model.getFullModelSpecification());
   runSimulation(model, newInterfaceCells, oldInterfaceCells, interfaceNormal);

   return EXIT_SUCCESS;
}
} // namespace ExcessMassDistributionFallbackTest
} // namespace free_surface
} // namespace walberla

int main(int argc, char** argv) { return walberla::free_surface::ExcessMassDistributionFallbackTest::main(argc, argv); }