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
//! \file ExcessMassDistributionParallelTest.cpp
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//! \brief Test distribution of excess mass using two blocks (to verify also on a parallel environment).
//!
//! The tests and their expected solutions are also shown in the file
//! "/tests/lbm/free_surface/dynamics/ExcessMassDistributionParallelTest.odp".
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
namespace ExcessMassDistributionTest
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

void runSimulation(const ExcessMassDistributionModel& excessMassDistributionModel)
{
   // define the domain size
   const Vector3< uint_t > numBlocks(uint_c(2), uint_c(1), uint_c(1));
   const Vector3< uint_t > domainSize(uint_c(6), uint_c(3), uint_c(1));
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

   // add field for excess mass
   const BlockDataID excessMassFieldID =
      field::addToStorage< ScalarField_T >(blockForest, "Excess mass field", real_c(0), field::fzyx, uint_c(1));

   // add boundary handling
   const std::shared_ptr< FreeSurfaceBoundaryHandling_T > freeSurfaceBoundaryHandling =
      std::make_shared< FreeSurfaceBoundaryHandling_T >(blockForest, pdfFieldID, fillFieldID);
   const BlockDataID flagFieldID                                      = freeSurfaceBoundaryHandling->getFlagFieldID();
   const typename FreeSurfaceBoundaryHandling_T::FlagInfo_T& flagInfo = freeSurfaceBoundaryHandling->getFlagInfo();

   // initialize cells as in file "/tests/lbm/free_surface/dynamics/ExcessMassDistributionParallelTest.odp"
   for (auto blockIt = blockForest->begin(); blockIt != blockForest->end(); ++blockIt)
   {
      ScalarField_T* const fillField       = blockIt->getData< ScalarField_T >(fillFieldID);
      ScalarField_T* const excessMassField = blockIt->getData< ScalarField_T >(excessMassFieldID);
      FlagField_T* const flagField         = blockIt->getData< FlagField_T >(flagFieldID);
      PdfField_T* const pdfField           = blockIt->getData< PdfField_T >(pdfFieldID);
      VectorField_T* const normalField     = blockIt->getData< VectorField_T >(normalFieldID);

      WALBERLA_FOR_ALL_CELLS(
         fillFieldIt, fillField, excessMassFieldIt, excessMassField, flagFieldIt, flagField, pdfFieldIt, pdfField,
         normalFieldIt, normalField, {
            const Cell localCell = fillFieldIt.cell();

            // get global coordinate of this cell
            Cell globalCell;
            blockForest->transformBlockLocalToGlobalCell(globalCell, *blockIt, localCell);

            if (globalCell == Cell(cell_idx_c(0), cell_idx_c(0), cell_idx_c(0)))
            {
               *fillFieldIt       = real_c(1);
               *excessMassFieldIt = real_c(0.2);
               pdfField->setDensityAndVelocity(localCell, Vector3< real_t >(real_c(0)), real_c(1));
               field::addFlag(flagFieldIt, flagInfo.liquidFlag);
            }

            if (globalCell == Cell(cell_idx_c(1), cell_idx_c(0), cell_idx_c(0)))
            {
               *fillFieldIt = real_c(0.5);
               pdfField->setDensityAndVelocity(localCell, Vector3< real_t >(real_c(0)), real_c(1.5));
               field::addFlag(flagFieldIt, flagInfo.interfaceFlag | flagInfo.convertedFlag);
            }

            if (globalCell == Cell(cell_idx_c(2), cell_idx_c(0), cell_idx_c(0)))
            {
               *fillFieldIt = real_c(0.5);
               pdfField->setDensityAndVelocity(localCell, Vector3< real_t >(real_c(0)), real_c(1.4));
               field::addFlag(flagFieldIt, flagInfo.interfaceFlag);
            }

            if (globalCell == Cell(cell_idx_c(0), cell_idx_c(1), cell_idx_c(0)))
            {
               *fillFieldIt       = real_c(1);
               *excessMassFieldIt = real_c(0.1);
               pdfField->setDensityAndVelocity(localCell, Vector3< real_t >(real_c(0)), real_c(1));
               field::addFlag(flagFieldIt, flagInfo.liquidFlag);
            }

            if (globalCell == Cell(cell_idx_c(1), cell_idx_c(1), cell_idx_c(0)))
            {
               *fillFieldIt = real_c(1.1);
               pdfField->setDensityAndVelocity(localCell, Vector3< real_t >(real_c(0)), real_c(1));
               field::addFlag(flagFieldIt, flagInfo.liquidFlag | flagInfo.convertedFlag);
               *normalFieldIt = Vector3< real_t >(real_c(0.71), real_c(0.71), real_c(0));
            }

            if (globalCell == Cell(cell_idx_c(2), cell_idx_c(1), cell_idx_c(0)))
            {
               *fillFieldIt = real_c(0.5);
               pdfField->setDensityAndVelocity(localCell, Vector3< real_t >(real_c(0)), real_c(1.3));
               field::addFlag(flagFieldIt, flagInfo.interfaceFlag);
            }

            if (globalCell == Cell(cell_idx_c(0), cell_idx_c(2), cell_idx_c(0)))
            {
               *fillFieldIt = real_c(0);
               pdfField->setDensityAndVelocity(localCell, Vector3< real_t >(real_c(0)), real_c(1));
               field::addFlag(flagFieldIt, flagInfo.gasFlag);
            }

            if (globalCell == Cell(cell_idx_c(1), cell_idx_c(2), cell_idx_c(0)))
            {
               *fillFieldIt = real_c(0.5);
               pdfField->setDensityAndVelocity(localCell, Vector3< real_t >(real_c(0)), real_c(1.1));
               field::addFlag(flagFieldIt, flagInfo.interfaceFlag | flagInfo.convertedFlag);
            }

            if (globalCell == Cell(cell_idx_c(2), cell_idx_c(2), cell_idx_c(0)))
            {
               *fillFieldIt = real_c(0.5);
               pdfField->setDensityAndVelocity(localCell, Vector3< real_t >(real_c(0)), real_c(1.2));
               field::addFlag(flagFieldIt, flagInfo.interfaceFlag);
            }

            if (globalCell == Cell(cell_idx_c(3), cell_idx_c(0), cell_idx_c(0)))
            {
               *fillFieldIt = real_c(0.5);
               pdfField->setDensityAndVelocity(localCell, Vector3< real_t >(real_c(0)), real_c(0.5));
               field::addFlag(flagFieldIt, flagInfo.interfaceFlag | flagInfo.convertedFlag);
            }

            if (globalCell == Cell(cell_idx_c(4), cell_idx_c(0), cell_idx_c(0)))
            {
               *fillFieldIt = real_c(0.5);
               pdfField->setDensityAndVelocity(localCell, Vector3< real_t >(real_c(0)), real_c(0.6));
               field::addFlag(flagFieldIt, flagInfo.interfaceFlag | flagInfo.convertedFlag);
            }

            if (globalCell == Cell(cell_idx_c(5), cell_idx_c(0), cell_idx_c(0)))
            {
               *fillFieldIt       = real_c(1);
               *excessMassFieldIt = real_c(0.03);
               pdfField->setDensityAndVelocity(localCell, Vector3< real_t >(real_c(0)), real_c(1));
               field::addFlag(flagFieldIt, flagInfo.liquidFlag);
            }

            if (globalCell == Cell(cell_idx_c(3), cell_idx_c(1), cell_idx_c(0)))
            {
               *fillFieldIt = real_c(1.2);
               pdfField->setDensityAndVelocity(localCell, Vector3< real_t >(real_c(0)), real_c(0.95));
               field::addFlag(flagFieldIt, flagInfo.liquidFlag | flagInfo.convertedFlag);
               *normalFieldIt = Vector3< real_t >(real_c(0.71), real_c(0.71), real_c(0));
            }

            if (globalCell == Cell(cell_idx_c(4), cell_idx_c(1), cell_idx_c(0)))
            {
               *fillFieldIt = real_c(0.5);
               pdfField->setDensityAndVelocity(localCell, Vector3< real_t >(real_c(0)), real_c(0.7));
               field::addFlag(flagFieldIt, flagInfo.interfaceFlag | flagInfo.convertedFlag);
            }

            if (globalCell == Cell(cell_idx_c(5), cell_idx_c(1), cell_idx_c(0)))
            {
               *fillFieldIt       = real_c(1);
               *excessMassFieldIt = real_c(0.02);
               pdfField->setDensityAndVelocity(localCell, Vector3< real_t >(real_c(0)), real_c(1));
               field::addFlag(flagFieldIt, flagInfo.liquidFlag);
            }

            if (globalCell == Cell(cell_idx_c(3), cell_idx_c(2), cell_idx_c(0)))
            {
               *fillFieldIt = real_c(0.5);
               pdfField->setDensityAndVelocity(localCell, Vector3< real_t >(real_c(0)), real_c(0.9));
               field::addFlag(flagFieldIt, flagInfo.interfaceFlag);
            }

            if (globalCell == Cell(cell_idx_c(4), cell_idx_c(2), cell_idx_c(0)))
            {
               *fillFieldIt = real_c(0.5);
               pdfField->setDensityAndVelocity(localCell, Vector3< real_t >(real_c(0)), real_c(0.8));
               field::addFlag(flagFieldIt, flagInfo.interfaceFlag | flagInfo.convertedFlag);
            }

            if (globalCell == Cell(cell_idx_c(5), cell_idx_c(2), cell_idx_c(0)))
            {
               *fillFieldIt       = real_c(1);
               *excessMassFieldIt = real_c(0.01);
               pdfField->setDensityAndVelocity(localCell, Vector3< real_t >(real_c(0)), real_c(1));
               field::addFlag(flagFieldIt, flagInfo.liquidFlag);
            }
         }) // WALBERLA_FOR_ALL_CELLS
   }

   // initial communication
   Communication_T(blockForest, pdfFieldID, fillFieldID, flagFieldID, normalFieldID, excessMassFieldID)();

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
      else
      {
         if (excessMassDistributionModel.isEvenlyAllInterfaceFallbackLiquidType())
         {
            const ExcessMassDistributionSweepInterfaceAndLiquid< LatticeModel_T, FlagField_T, ScalarField_T,
                                                                 VectorField_T >
               distributeMassSweep(excessMassDistributionModel, fillFieldID, flagFieldID, pdfFieldID, flagInfo,
                                   excessMassFieldID);
            timeloop.add() << Sweep(distributeMassSweep, "Distribute excess mass");
         }
      }
   }

   timeloop.singleStep();

   // check if excess mass was distributed correctly; expected solutions were obtained manually, see file
   // "/tests/lbm/free_surface/dynamics/ExcessMassDistributionParallelTest.odp"
   for (auto blockIt = blockForest->begin(); blockIt != blockForest->end(); ++blockIt)
   {
      const ScalarField_T* const fillField       = blockIt->getData< const ScalarField_T >(fillFieldID);
      const ScalarField_T* const excessMassField = blockIt->getData< const ScalarField_T >(excessMassFieldID);

      WALBERLA_FOR_ALL_CELLS(fillFieldIt, fillField, excessMassFieldIt, excessMassField, {
         Cell globalCell;
         blockForest->transformBlockLocalToGlobalCell(globalCell, *blockIt, fillFieldIt.cell());

         if (excessMassDistributionModel.getModelType() ==
             ExcessMassDistributionModel::ExcessMassModel::EvenlyAllInterface)
         {
            // left block
            if (globalCell == Cell(cell_idx_c(0), cell_idx_c(0), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(1), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(1), cell_idx_c(0), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(0.513333333333333), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(2), cell_idx_c(0), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(0.53125), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(0), cell_idx_c(1), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(1), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(1), cell_idx_c(1), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(1), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(2), cell_idx_c(1), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(0.533653846153846), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(0), cell_idx_c(2), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(0), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(1), cell_idx_c(2), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(0.518181818181818), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(2), cell_idx_c(2), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(0.536458333333333), real_c(1e-4));
            }

            // right block
            if (globalCell == Cell(cell_idx_c(3), cell_idx_c(0), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(0.5475), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(4), cell_idx_c(0), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(0.539583333333333), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(5), cell_idx_c(0), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(1), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(3), cell_idx_c(1), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(1), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(4), cell_idx_c(1), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(0.533928571428571), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(5), cell_idx_c(1), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(1), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(3), cell_idx_c(2), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(0.526388888888889), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(4), cell_idx_c(2), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(0.5296875), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(5), cell_idx_c(2), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(1), real_c(1e-4));
            }
         }

         if (excessMassDistributionModel.getModelType() ==
             ExcessMassDistributionModel::ExcessMassModel::EvenlyOldInterface)
         {
            // left block
            if (globalCell == Cell(cell_idx_c(0), cell_idx_c(0), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(1), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(1), cell_idx_c(0), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(0.5), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(2), cell_idx_c(0), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(0.557738095238095), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(0), cell_idx_c(1), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(1), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(1), cell_idx_c(1), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(1), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(2), cell_idx_c(1), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(0.562179487179487), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(0), cell_idx_c(2), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(0), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(1), cell_idx_c(2), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(0.5), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(2), cell_idx_c(2), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(0.567361111111111), real_c(1e-4));
            }

            // right block
            if (globalCell == Cell(cell_idx_c(3), cell_idx_c(0), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(0.5), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(4), cell_idx_c(0), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(0.5), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(5), cell_idx_c(0), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(1), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(3), cell_idx_c(1), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(1), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(4), cell_idx_c(1), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(0.5), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(5), cell_idx_c(1), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(1), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(3), cell_idx_c(2), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(0.552777777777778), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(4), cell_idx_c(2), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(0.5), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(5), cell_idx_c(2), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(1), real_c(1e-4));
            }
         }

         if (excessMassDistributionModel.getModelType() ==
             ExcessMassDistributionModel::ExcessMassModel::EvenlyNewInterface)
         {
            // left block
            if (globalCell == Cell(cell_idx_c(0), cell_idx_c(0), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(1), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(1), cell_idx_c(0), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(0.533333333333333), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(2), cell_idx_c(0), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(0.5), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(0), cell_idx_c(1), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(1), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(1), cell_idx_c(1), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(1), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(2), cell_idx_c(1), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(0.5), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(0), cell_idx_c(2), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(0), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(1), cell_idx_c(2), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(0.545454545454546), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(2), cell_idx_c(2), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(0.5), real_c(1e-4));
            }

            // right block
            if (globalCell == Cell(cell_idx_c(3), cell_idx_c(0), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(0.595), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(4), cell_idx_c(0), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(0.579166666666667), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(5), cell_idx_c(0), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(1), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(3), cell_idx_c(1), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(1), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(4), cell_idx_c(1), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(0.567857142857143), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(5), cell_idx_c(1), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(1), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(3), cell_idx_c(2), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(0.5), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(4), cell_idx_c(2), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(0.559375), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(5), cell_idx_c(2), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(1), real_c(1e-4));
            }
         }

         if (excessMassDistributionModel.getModelType() ==
             ExcessMassDistributionModel::ExcessMassModel::WeightedAllInterface)
         {
            // left block
            if (globalCell == Cell(cell_idx_c(0), cell_idx_c(0), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(1), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(1), cell_idx_c(0), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(0.5), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(2), cell_idx_c(0), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(0.5), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(0), cell_idx_c(1), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(1), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(1), cell_idx_c(1), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(1), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(2), cell_idx_c(1), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(0.519230769230769), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(0), cell_idx_c(2), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(0), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(1), cell_idx_c(2), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(0.522727272727273), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(2), cell_idx_c(2), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(0.541666666666667), real_c(1e-4));
            }

            // right block
            if (globalCell == Cell(cell_idx_c(3), cell_idx_c(0), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(0.5), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(4), cell_idx_c(0), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(0.5), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(5), cell_idx_c(0), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(1), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(3), cell_idx_c(1), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(1), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(4), cell_idx_c(1), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(0.567857142857143), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(5), cell_idx_c(1), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(1), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(3), cell_idx_c(2), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(0.552777777777778), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(4), cell_idx_c(2), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(0.61875), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(5), cell_idx_c(2), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(1), real_c(1e-4));
            }
         }

         if (excessMassDistributionModel.getModelType() ==
             ExcessMassDistributionModel::ExcessMassModel::WeightedOldInterface)
         {
            // left block
            if (globalCell == Cell(cell_idx_c(0), cell_idx_c(0), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(1), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(1), cell_idx_c(0), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(0.5), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(2), cell_idx_c(0), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(0.5), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(0), cell_idx_c(1), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(1), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(1), cell_idx_c(1), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(1), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(2), cell_idx_c(1), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(0.525641025641026), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(0), cell_idx_c(2), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(0), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(1), cell_idx_c(2), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(0.5), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(2), cell_idx_c(2), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(0.555555555555556), real_c(1e-4));
            }

            // right block
            if (globalCell == Cell(cell_idx_c(3), cell_idx_c(0), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(0.5), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(4), cell_idx_c(0), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(0.5), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(5), cell_idx_c(0), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(1), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(3), cell_idx_c(1), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(1), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(4), cell_idx_c(1), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(0.5), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(5), cell_idx_c(1), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(1), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(3), cell_idx_c(2), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(0.711111111111111), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(4), cell_idx_c(2), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(0.5), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(5), cell_idx_c(2), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(1), real_c(1e-4));
            }
         }

         if (excessMassDistributionModel.getModelType() ==
             ExcessMassDistributionModel::ExcessMassModel::WeightedNewInterface)
         {
            // left block
            if (globalCell == Cell(cell_idx_c(0), cell_idx_c(0), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(1), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(1), cell_idx_c(0), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(0.5), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(2), cell_idx_c(0), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(0.5), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(0), cell_idx_c(1), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(1), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(1), cell_idx_c(1), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(1), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(2), cell_idx_c(1), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(0.5), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(0), cell_idx_c(2), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(0), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(1), cell_idx_c(2), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(0.590909090909091), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(2), cell_idx_c(2), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(0.5), real_c(1e-4));
            }

            // right block
            if (globalCell == Cell(cell_idx_c(3), cell_idx_c(0), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(0.5), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(4), cell_idx_c(0), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(0.5), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(5), cell_idx_c(0), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(1), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(3), cell_idx_c(1), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(1), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(4), cell_idx_c(1), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(0.59047619047619), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(5), cell_idx_c(1), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(1), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(3), cell_idx_c(2), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(0.5), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(4), cell_idx_c(2), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(0.658333333333333), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(5), cell_idx_c(2), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(1), real_c(1e-4));
            }
         }

         if (excessMassDistributionModel.getModelType() ==
             ExcessMassDistributionModel::ExcessMassModel::EvenlyAllInterfaceAndLiquid)
         {
            // left block
            if (globalCell == Cell(cell_idx_c(0), cell_idx_c(0), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(1), real_c(1e-4));
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*excessMassFieldIt, real_c(0.039285714285714), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(1), cell_idx_c(0), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(0.570634920634921), real_c(1e-4));
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*excessMassFieldIt, real_c(0), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(2), cell_idx_c(0), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(0.527168367346939), real_c(1e-4));
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*excessMassFieldIt, real_c(0), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(0), cell_idx_c(1), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(1), real_c(1e-4));
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*excessMassFieldIt, real_c(0.080952380952381), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(1), cell_idx_c(1), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(1), real_c(1e-4));
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*excessMassFieldIt, real_c(0.091666666666667), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(2), cell_idx_c(1), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(0.529258241758242), real_c(1e-4));
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*excessMassFieldIt, real_c(0), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(0), cell_idx_c(2), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(0), real_c(1e-4));
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*excessMassFieldIt, real_c(0), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(1), cell_idx_c(2), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(0.535714285714286), real_c(1e-4));
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*excessMassFieldIt, real_c(0), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(2), cell_idx_c(2), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(0.531696428571429), real_c(1e-4));
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*excessMassFieldIt, real_c(0), real_c(1e-4));
            }

            // right block
            if (globalCell == Cell(cell_idx_c(3), cell_idx_c(0), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(0.5475), real_c(1e-4));
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*excessMassFieldIt, real_c(0), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(4), cell_idx_c(0), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(0.562916666666667), real_c(1e-4));
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*excessMassFieldIt, real_c(0), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(5), cell_idx_c(0), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(1), real_c(1e-4));
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*excessMassFieldIt, real_c(0.004), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(3), cell_idx_c(1), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(1), real_c(1e-4));
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*excessMassFieldIt, real_c(0), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(4), cell_idx_c(1), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(0.558690476190476), real_c(1e-4));
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*excessMassFieldIt, real_c(0), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(5), cell_idx_c(1), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(1), real_c(1e-4));
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*excessMassFieldIt, real_c(0.013333333333333), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(3), cell_idx_c(2), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(0.526388888888889), real_c(1e-4));
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*excessMassFieldIt, real_c(0), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(4), cell_idx_c(2), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(0.538854166666667), real_c(1e-4));
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*excessMassFieldIt, real_c(0), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(5), cell_idx_c(2), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(1), real_c(1e-4));
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*excessMassFieldIt, real_c(0.004), real_c(1e-4));
            }
         }

         if (excessMassDistributionModel.getModelType() ==
             ExcessMassDistributionModel::ExcessMassModel::EvenlyAllInterfaceFallbackLiquid)
         {
            // left block
            if (globalCell == Cell(cell_idx_c(0), cell_idx_c(0), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(1), real_c(1e-4));
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*excessMassFieldIt, real_c(0), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(1), cell_idx_c(0), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(0.68), real_c(1e-4));
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*excessMassFieldIt, real_c(0), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(2), cell_idx_c(0), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(0.53125), real_c(1e-4));
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*excessMassFieldIt, real_c(0), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(0), cell_idx_c(1), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(1), real_c(1e-4));
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*excessMassFieldIt, real_c(0), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(1), cell_idx_c(1), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(1), real_c(1e-4));
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*excessMassFieldIt, real_c(0), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(2), cell_idx_c(1), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(0.533653846153846), real_c(1e-4));
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*excessMassFieldIt, real_c(0), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(0), cell_idx_c(2), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(0), real_c(1e-4));
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*excessMassFieldIt, real_c(0), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(1), cell_idx_c(2), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(0.563636363636364), real_c(1e-4));
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*excessMassFieldIt, real_c(0), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(2), cell_idx_c(2), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(0.536458333333333), real_c(1e-4));
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*excessMassFieldIt, real_c(0), real_c(1e-4));
            }

            // right block
            if (globalCell == Cell(cell_idx_c(3), cell_idx_c(0), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(0.5475), real_c(1e-4));
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*excessMassFieldIt, real_c(0), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(4), cell_idx_c(0), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(0.575694444444444), real_c(1e-4));
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*excessMassFieldIt, real_c(0), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(5), cell_idx_c(0), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(1), real_c(1e-4));
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*excessMassFieldIt, real_c(0), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(3), cell_idx_c(1), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(1), real_c(1e-4));
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*excessMassFieldIt, real_c(0), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(4), cell_idx_c(1), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(0.57202380952381), real_c(1e-4));
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*excessMassFieldIt, real_c(0), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(5), cell_idx_c(1), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(1), real_c(1e-4));
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*excessMassFieldIt, real_c(0), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(3), cell_idx_c(2), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(0.526388888888889), real_c(1e-4));
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*excessMassFieldIt, real_c(0), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(4), cell_idx_c(2), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(0.544270833333333), real_c(1e-4));
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*excessMassFieldIt, real_c(0), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(5), cell_idx_c(2), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(1), real_c(1e-4));
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*excessMassFieldIt, real_c(0), real_c(1e-4));
            }
         }

         if (excessMassDistributionModel.getModelType() ==
             ExcessMassDistributionModel::ExcessMassModel::EvenlyNewInterfaceFallbackLiquid)
         {
            // left block
            if (globalCell == Cell(cell_idx_c(0), cell_idx_c(0), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(1), real_c(1e-4));
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*excessMassFieldIt, real_c(0), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(1), cell_idx_c(0), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(0.7), real_c(1e-4));
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*excessMassFieldIt, real_c(0), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(2), cell_idx_c(0), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(0.5), real_c(1e-4));
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*excessMassFieldIt, real_c(0), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(0), cell_idx_c(1), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(1), real_c(1e-4));
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*excessMassFieldIt, real_c(0), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(1), cell_idx_c(1), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(1), real_c(1e-4));
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*excessMassFieldIt, real_c(0), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(2), cell_idx_c(1), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(0.5), real_c(1e-4));
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*excessMassFieldIt, real_c(0), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(0), cell_idx_c(2), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(0), real_c(1e-4));
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*excessMassFieldIt, real_c(0), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(1), cell_idx_c(2), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(0.590909090909091), real_c(1e-4));
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*excessMassFieldIt, real_c(0), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(2), cell_idx_c(2), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(0.5), real_c(1e-4));
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*excessMassFieldIt, real_c(0), real_c(1e-4));
            }

            // right block
            if (globalCell == Cell(cell_idx_c(3), cell_idx_c(0), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(0.595), real_c(1e-4));
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*excessMassFieldIt, real_c(0), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(4), cell_idx_c(0), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(0.615277777777778), real_c(1e-4));
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*excessMassFieldIt, real_c(0), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(5), cell_idx_c(0), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(1), real_c(1e-4));
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*excessMassFieldIt, real_c(0), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(3), cell_idx_c(1), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(1), real_c(1e-4));
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*excessMassFieldIt, real_c(0), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(4), cell_idx_c(1), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(0.605952380952381), real_c(1e-4));
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*excessMassFieldIt, real_c(0), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(5), cell_idx_c(1), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(1), real_c(1e-4));
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*excessMassFieldIt, real_c(0), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(3), cell_idx_c(2), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(0.5), real_c(1e-4));
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*excessMassFieldIt, real_c(0), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(4), cell_idx_c(2), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(0.573958333333333), real_c(1e-4));
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*excessMassFieldIt, real_c(0), real_c(1e-4));
            }

            if (globalCell == Cell(cell_idx_c(5), cell_idx_c(2), cell_idx_c(0)))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*fillFieldIt, real_c(1), real_c(1e-4));
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(*excessMassFieldIt, real_c(0), real_c(1e-4));
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

   ExcessMassDistributionModel model = ExcessMassDistributionModel("EvenlyAllInterface");
   WALBERLA_LOG_INFO_ON_ROOT("Testing model " << model.getFullModelSpecification());
   runSimulation(model);

   model = ExcessMassDistributionModel("EvenlyOldInterface");
   WALBERLA_LOG_INFO_ON_ROOT("Testing model " << model.getFullModelSpecification());
   runSimulation(model);

   model = ExcessMassDistributionModel("EvenlyNewInterface");
   WALBERLA_LOG_INFO_ON_ROOT("Testing model " << model.getFullModelSpecification());
   runSimulation(model);

   model = ExcessMassDistributionModel("WeightedAllInterface");
   WALBERLA_LOG_INFO_ON_ROOT("Testing model " << model.getFullModelSpecification());
   runSimulation(model);

   model = ExcessMassDistributionModel("WeightedOldInterface");
   WALBERLA_LOG_INFO_ON_ROOT("Testing model " << model.getFullModelSpecification());
   runSimulation(model);

   model = ExcessMassDistributionModel("WeightedNewInterface");
   WALBERLA_LOG_INFO_ON_ROOT("Testing model " << model.getFullModelSpecification());
   runSimulation(model);

   model = ExcessMassDistributionModel("EvenlyAllInterfaceAndLiquid");
   WALBERLA_LOG_INFO_ON_ROOT("Testing model " << model.getFullModelSpecification());
   runSimulation(model);

   model = ExcessMassDistributionModel("EvenlyAllInterfaceFallbackLiquid");
   WALBERLA_LOG_INFO_ON_ROOT("Testing model " << model.getFullModelSpecification());
   runSimulation(model);

   model = ExcessMassDistributionModel("EvenlyNewInterfaceFallbackLiquid");
   WALBERLA_LOG_INFO_ON_ROOT("Testing model " << model.getFullModelSpecification());
   runSimulation(model);

   return EXIT_SUCCESS;
}
} // namespace ExcessMassDistributionTest
} // namespace free_surface
} // namespace walberla

int main(int argc, char** argv) { return walberla::free_surface::ExcessMassDistributionTest::main(argc, argv); }