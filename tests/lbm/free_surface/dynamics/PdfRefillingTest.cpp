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
//! \file PdfRefillingTest.cpp
//! \ingroup lbm/free_surface/dynamics
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//! \author Michael Zikeli
//! \brief Test PDF refilling of gas cells that are converted to interface cells at the free surface boundary.
//!
//! For each test two helper functions are necessary, that are both passed in the main as callback functions
//!     - The initialization function which fills four sets (gasCells, interfaceCells, liquidCells, conversionCells).
//!     - The verification function, that checks whether the values are calculated correctly.
//! The domain is consists of 5x3x1 cells with periodic boundaries. Different scenarios are tested. The expected
//! solutions are obtained manually from the file "/tests/lbm/free_surface/dynamics/PdfRefillingTest.odp".
//======================================================================================================================

#include "blockforest/Initialization.h"

#include "core/Environment.h"
#include "core/debug/TestSubsystem.h"

#include "lbm/blockforest/communication/SimpleCommunication.h"
#include "lbm/field/Adaptors.h"
#include "lbm/field/AddToStorage.h"
#include "lbm/free_surface/FlagInfo.h"
#include "lbm/free_surface/boundary/FreeSurfaceBoundaryHandling.h"
#include "lbm/free_surface/dynamics/PdfRefillingModel.h"
#include "lbm/free_surface/dynamics/PdfRefillingSweep.h"
#include "lbm/free_surface/surface_geometry/SurfaceGeometryHandler.h"
#include "lbm/lattice_model/D2Q9.h"
#include "lbm/sweeps/CellwiseSweep.h"

#include "timeloop/SweepTimeloop.h"

#include <functional>
#include <iostream>

namespace walberla
{
namespace free_surface
{
namespace PdfRefillingTest
{
using ScalarField_T = GhostLayerField< real_t, 1 >;
using VectorField_T = GhostLayerField< Vector3< real_t >, 1 >;

using LatticeModel_T = lbm::D2Q9< lbm::collision_model::SRT, true, lbm::force_model::None, 2 >;
using Stencil_T      = typename LatticeModel_T::Stencil;

using Communication_T = blockforest::SimpleCommunication< LatticeModel_T::CommunicationStencil >;

using flag_t                        = uint32_t;
using FlagField_T                   = field::FlagField< flag_t >;
using FreeSurfaceBoundaryHandling_T = FreeSurfaceBoundaryHandling< LatticeModel_T, FlagField_T, ScalarField_T >;
using FlagInfo_T                    = FlagInfo< FlagField_T >;

using PdfField_T = lbm::PdfField< LatticeModel_T >;

using RefillingModel_T = PdfRefillingModel::RefillingModel;

using initCallback_T =
   std::function< void(std::set< Cell >&, std::set< Cell >&, std::set< Cell >&, std::set< Cell >&) >;
using verifyCallback_T = std::function< void(const PdfRefillingModel&, const Cell&, const PdfField_T* const) >;

void runSimulation(const PdfRefillingModel& pdfRefillingModel, const initCallback_T& fillInitializationSets,
                   const verifyCallback_T& verifyResults);

/***********************************************************************************************************************
 * Test: noRefillingCells
 *      | G | G | G | G | G |       | G | G | G | G | G |
 *      | G | G | G | G | G |  ==>  | G | G | I | G | G |
 *      | G | G | G | G | G |       | G | G | G | G | G |
 **********************************************************************************************************************/
void init_noRefillingCells(std::set< Cell >& gasCells, std::set< Cell >& interfaceCells, std::set< Cell >& liquidCells,
                           std::set< Cell >& conversionCells);
void verify_noRefillingCells(const PdfRefillingModel& pdfRefillingModel, const Cell& cell,
                             const PdfField_T* const pdfField);

/***********************************************************************************************************************
 * Test: fullRefillingCells
 *      | I | I | I | I | I |       | I | I | I | I | I |
 *      | I | G | I | I | I |  ==>  | I | I | I | I | I |
 *      | I | I | I | I | I |       | I | I | I | I | I |
 **********************************************************************************************************************/
void init_fullRefillingCells(std::set< Cell >& gasCells, std::set< Cell >& interfaceCells,
                             std::set< Cell >& liquidCells, std::set< Cell >& conversionCells);
void verify_fullRefillingCells(const PdfRefillingModel& pdfRefillingModel, const Cell& cell,
                               const PdfField_T* const pdfField);

/***********************************************************************************************************************
 * Test: someRefillingCells
 *      | G | I | I | I | I |       | I | I | I | I | I |
 *      | G | G | I | I | I |  ==>  | G | I | I | I | I |
 *      | G | G | G | I | I |       | G | G | I | I | I |
 **********************************************************************************************************************/
void init_someRefillingCells(std::set< Cell >& gasCells, std::set< Cell >& interfaceCells,
                             std::set< Cell >& liquidCells, std::set< Cell >& conversionCells);
void verify_someRefillingCells(PdfRefillingModel const& pdfRefillingModel, Cell const& cell,
                               PdfField_T const* const pdfField);

/***********************************************************************************************************************
 * Test: straightRefillingCells
 *      | G | G | I | I | I |       | G | G | I | I | I |
 *      | G | G | I | I | I |  ==>  | G | I | I | I | I |
 *      | G | G | I | I | I |       | G | G | I | I | I |
 **********************************************************************************************************************/
void init_straightRefillingCells(std::set< Cell >& gasCells, std::set< Cell >& interfaceCells,
                                 std::set< Cell >& liquidCells, std::set< Cell >& conversionCells);
void verify_straightRefillingCells(const PdfRefillingModel& pdfRefillingModel, const Cell& cell,
                                   const PdfField_T* const pdfField);

int main(int argc, char** argv)
{
   debug::enterTestMode();
   Environment env(argc, argv);

   for (auto modelType : PdfRefillingModel::getTypeIterator())
   {
      PdfRefillingModel model = PdfRefillingModel(modelType);

      WALBERLA_LOG_INFO_ON_ROOT("Testing model:" << model.getFullModelSpecification()
                                                 << " with setup \"noRefillingCells\":"
                                                 << "\n\t | G | G | G | G | G |       | G | G | G | G | G | "
                                                 << "\n\t | G | G | G | G | G |  ==>  | G | G | I | G | G | "
                                                 << "\n\t | G | G | G | G | G |       | G | G | G | G | G | \n");
      runSimulation(model, &init_noRefillingCells, &verify_noRefillingCells);

      WALBERLA_LOG_INFO_ON_ROOT("Testing model:" << model.getFullModelSpecification()
                                                 << " with setup \"fullRefillingCells\":"
                                                 << "\n\t | I | I | I | I | I |       | I | I | I | I | I | "
                                                 << "\n\t | I | G | I | I | I |  ==>  | I | I | I | I | I | "
                                                 << "\n\t | I | I | I | I | I |       | I | I | I | I | I | \n");
      runSimulation(model, init_fullRefillingCells, verify_fullRefillingCells);

      WALBERLA_LOG_INFO_ON_ROOT("Testing model:" << model.getFullModelSpecification()
                                                 << " with setup \"someRefillingCells\":"
                                                 << "\n\t | G | I | I | I | I |       | I | I | I | I | I | "
                                                 << "\n\t | G | G | I | I | I |  ==>  | G | I | I | I | I | "
                                                 << "\n\t | G | G | G | I | I |       | G | G | I | I | I | \n");
      runSimulation(model, &init_someRefillingCells, &verify_someRefillingCells);

      WALBERLA_LOG_INFO_ON_ROOT("Testing model:" << model.getFullModelSpecification()
                                                 << " with setup \"straightRefillingCells\":"
                                                 << "\n\t | G | G | I | I | I |       | G | G | I | I | I | "
                                                 << "\n\t | G | G | I | I | I |  ==>  | G | I | I | I | I | "
                                                 << "\n\t | G | G | I | I | I |       | G | G | I | I | I | \n");
      runSimulation(model, &init_straightRefillingCells, &verify_straightRefillingCells);
   }

   return EXIT_SUCCESS;
}

void init_noRefillingCells(std::set< Cell >& gasCells, std::set< Cell >& interfaceCells, std::set< Cell >& liquidCells,
                           std::set< Cell >& conversionCells)
{
   gasCells = std::set< Cell >(
      { Cell(cell_idx_c(0), cell_idx_c(0), cell_idx_c(0)), Cell(cell_idx_c(0), cell_idx_c(1), cell_idx_c(0)),
        Cell(cell_idx_c(0), cell_idx_c(2), cell_idx_c(0)), Cell(cell_idx_c(1), cell_idx_c(0), cell_idx_c(0)),
        Cell(cell_idx_c(1), cell_idx_c(1), cell_idx_c(0)), Cell(cell_idx_c(1), cell_idx_c(2), cell_idx_c(0)),
        Cell(cell_idx_c(2), cell_idx_c(0), cell_idx_c(0)), Cell(cell_idx_c(2), cell_idx_c(1), cell_idx_c(0)),
        Cell(cell_idx_c(2), cell_idx_c(2), cell_idx_c(0)), Cell(cell_idx_c(3), cell_idx_c(0), cell_idx_c(0)),
        Cell(cell_idx_c(3), cell_idx_c(1), cell_idx_c(0)), Cell(cell_idx_c(3), cell_idx_c(2), cell_idx_c(0)),
        Cell(cell_idx_c(4), cell_idx_c(0), cell_idx_c(0)), Cell(cell_idx_c(4), cell_idx_c(1), cell_idx_c(0)),
        Cell(cell_idx_c(4), cell_idx_c(2), cell_idx_c(0)) });

   interfaceCells = std::set< Cell >({});

   liquidCells = std::set< Cell >({});

   conversionCells = std::set< Cell >({ Cell(cell_idx_c(2), cell_idx_c(1), cell_idx_c(0)) });
}

void init_fullRefillingCells(std::set< Cell >& gasCells, std::set< Cell >& interfaceCells,
                             std::set< Cell >& liquidCells, std::set< Cell >& conversionCells)
{
   gasCells = std::set< Cell >({ Cell(cell_idx_c(1), cell_idx_c(1), cell_idx_c(0)) });

   interfaceCells = std::set< Cell >(
      { Cell(cell_idx_c(0), cell_idx_c(0), cell_idx_c(0)), Cell(cell_idx_c(0), cell_idx_c(1), cell_idx_c(0)),
        Cell(cell_idx_c(0), cell_idx_c(2), cell_idx_c(0)), Cell(cell_idx_c(1), cell_idx_c(0), cell_idx_c(0)),
        Cell(cell_idx_c(1), cell_idx_c(2), cell_idx_c(0)), Cell(cell_idx_c(2), cell_idx_c(0), cell_idx_c(0)),
        Cell(cell_idx_c(2), cell_idx_c(1), cell_idx_c(0)), Cell(cell_idx_c(2), cell_idx_c(2), cell_idx_c(0)),
        Cell(cell_idx_c(3), cell_idx_c(0), cell_idx_c(0)), Cell(cell_idx_c(3), cell_idx_c(1), cell_idx_c(0)),
        Cell(cell_idx_c(3), cell_idx_c(2), cell_idx_c(0)), Cell(cell_idx_c(4), cell_idx_c(0), cell_idx_c(0)),
        Cell(cell_idx_c(4), cell_idx_c(1), cell_idx_c(0)), Cell(cell_idx_c(4), cell_idx_c(2), cell_idx_c(0)) });

   liquidCells = std::set< Cell >({});

   conversionCells = std::set< Cell >({ Cell(cell_idx_c(1), cell_idx_c(1), cell_idx_c(0)) });
}

void init_someRefillingCells(std::set< Cell >& gasCells, std::set< Cell >& interfaceCells,
                             std::set< Cell >& liquidCells, std::set< Cell >& conversionCells)
{
   gasCells = std::set< Cell >(
      { Cell(cell_idx_c(0), cell_idx_c(0), cell_idx_c(0)), Cell(cell_idx_c(0), cell_idx_c(1), cell_idx_c(0)),
        Cell(cell_idx_c(1), cell_idx_c(0), cell_idx_c(0)), Cell(cell_idx_c(1), cell_idx_c(1), cell_idx_c(0)),
        Cell(cell_idx_c(0), cell_idx_c(2), cell_idx_c(0)), Cell(cell_idx_c(2), cell_idx_c(0), cell_idx_c(0)) });

   interfaceCells = std::set< Cell >(
      { Cell(cell_idx_c(2), cell_idx_c(2), cell_idx_c(0)), Cell(cell_idx_c(2), cell_idx_c(1), cell_idx_c(0)),
        Cell(cell_idx_c(1), cell_idx_c(2), cell_idx_c(0)), Cell(cell_idx_c(3), cell_idx_c(0), cell_idx_c(0)),
        Cell(cell_idx_c(3), cell_idx_c(1), cell_idx_c(0)), Cell(cell_idx_c(3), cell_idx_c(2), cell_idx_c(0)),
        Cell(cell_idx_c(4), cell_idx_c(0), cell_idx_c(0)), Cell(cell_idx_c(4), cell_idx_c(1), cell_idx_c(0)),
        Cell(cell_idx_c(4), cell_idx_c(2), cell_idx_c(0)) });

   liquidCells = std::set< Cell >({});

   conversionCells = std::set< Cell >({ Cell(cell_idx_c(0), cell_idx_c(2), cell_idx_c(0)),
                                        Cell(cell_idx_c(1), cell_idx_c(1), cell_idx_c(0)),
                                        Cell(cell_idx_c(2), cell_idx_c(0), cell_idx_c(0)) });
}

void init_straightRefillingCells(std::set< Cell >& gasCells, std::set< Cell >& interfaceCells,
                                 std::set< Cell >& liquidCells, std::set< Cell >& conversionCells)
{
   gasCells = std::set< Cell >(
      { Cell(cell_idx_c(0), cell_idx_c(0), cell_idx_c(0)), Cell(cell_idx_c(0), cell_idx_c(1), cell_idx_c(0)),
        Cell(cell_idx_c(0), cell_idx_c(2), cell_idx_c(0)), Cell(cell_idx_c(1), cell_idx_c(0), cell_idx_c(0)),
        Cell(cell_idx_c(1), cell_idx_c(1), cell_idx_c(0)), Cell(cell_idx_c(1), cell_idx_c(2), cell_idx_c(0)) });

   interfaceCells = std::set< Cell >(
      { Cell(cell_idx_c(2), cell_idx_c(0), cell_idx_c(0)), Cell(cell_idx_c(2), cell_idx_c(1), cell_idx_c(0)),
        Cell(cell_idx_c(2), cell_idx_c(2), cell_idx_c(0)), Cell(cell_idx_c(3), cell_idx_c(0), cell_idx_c(0)),
        Cell(cell_idx_c(3), cell_idx_c(1), cell_idx_c(0)), Cell(cell_idx_c(3), cell_idx_c(2), cell_idx_c(0)),
        Cell(cell_idx_c(4), cell_idx_c(0), cell_idx_c(0)), Cell(cell_idx_c(4), cell_idx_c(1), cell_idx_c(0)),
        Cell(cell_idx_c(4), cell_idx_c(2), cell_idx_c(0)) });

   liquidCells = std::set< Cell >({});

   conversionCells = std::set< Cell >({ Cell(cell_idx_c(1), cell_idx_c(1), cell_idx_c(0)) });
}

void verify_noRefillingCells(const PdfRefillingModel& pdfRefillingModel, const Cell& cell,
                             const PdfField_T* const pdfField)
{
   switch (pdfRefillingModel.getModelType())
   {
   case RefillingModel_T::EquilibriumRefilling:
   case RefillingModel_T::AverageRefilling:
   case RefillingModel_T::EquilibriumAndNonEquilibriumRefilling:
   case RefillingModel_T::ExtrapolationRefilling:
   case RefillingModel_T::GradsMomentsRefilling:
      WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, 0), real_c(0.444444444444444),
                                         real_c(1e-6)); // C, (0,0,0)
      WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, 1), real_c(0.111111111111111),
                                         real_c(1e-6)); // N, (0,1,0)
      WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, 2), real_c(0.111111111111111),
                                         real_c(1e-6)); // S, (0,-1,0)
      WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, 3), real_c(0.111111111111111),
                                         real_c(1e-6)); // W, (-1,0,0)
      WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, 4), real_c(0.111111111111111),
                                         real_c(1e-6)); // E, (1,0,0)
      WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, 5), real_c(0.027777777777778),
                                         real_c(1e-6)); // NW, (-1,1,0)
      WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, 6), real_c(0.027777777777778),
                                         real_c(1e-6)); // NE, (1,1,0)
      WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, 7), real_c(0.027777777777778),
                                         real_c(1e-6)); // SW, (-1,-1,0)
      WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, 8), real_c(0.027777777777778),
                                         real_c(1e-6)); // SE, (1,-1,0)
      break;
   default:
      WALBERLA_ABORT("The specified PDF refilling model " << pdfRefillingModel.getModelName()
                                                          << " is not available.\nSomething went wrong!");
   }
}

void verify_fullRefillingCells(const PdfRefillingModel& pdfRefillingModel, const Cell& cell,
                               const PdfField_T* const pdfField)
{
   switch (pdfRefillingModel.getModelType())
   {
   case RefillingModel_T::EquilibriumRefilling:
   case RefillingModel_T::EquilibriumAndNonEquilibriumRefilling:
   case RefillingModel_T::ExtrapolationRefilling:
      WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(0)), real_c(0.443017361111111),
                                         real_c(1e-6)); // C, (0,0,0)
      WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(1)), real_c(0.101584288194444),
                                         real_c(1e-6)); // N, (0,1,0)
      WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(2)), real_c(0.120750954861111),
                                         real_c(1e-6)); // S, (0,-1,0)
      WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(3)), real_c(0.099328038194445),
                                         real_c(1e-6)); // W, (-1,0,0)
      WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(4)), real_c(0.123494704861111),
                                         real_c(1e-6)); // E, (1,0,0)
      WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(5)), real_c(0.022800043402778),
                                         real_c(1e-6)); // NW, (-1,1,0)
      WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(6)), real_c(0.028320616319445),
                                         real_c(1e-6)); // NE, (1,1,0)
      WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(7)), real_c(0.027070616319445),
                                         real_c(1e-6)); // SW, (-1,-1,0)
      WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(8)), real_c(0.033633376736111),
                                         real_c(1e-6)); // SE, (1,-1,0)
      break;
   case RefillingModel_T::AverageRefilling:
      WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, 0), real_c(0.440902777777778),
                                         real_c(1e-6)); // C, (0,0,0)
      WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, 1), real_c(0.102829861111111),
                                         real_c(1e-6)); // N, (0,1,0)
      WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, 2), real_c(0.119496527777777),
                                         real_c(1e-6)); // S, (0,-1,0)
      WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, 3), real_c(0.097361111111111),
                                         real_c(1e-6)); // W, (-1,0,0)
      WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, 4), real_c(0.122777777777778),
                                         real_c(1e-6)); // E, (1,0,0)
      WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, 5), real_c(0.022803819444445),
                                         real_c(1e-6)); // NW, (-1,1,0)
      WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, 6), real_c(0.029470486111111),
                                         real_c(1e-6)); // NE, (1,1,0)
      WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, 7), real_c(0.030095486111111),
                                         real_c(1e-6)); // SW, (-1,-1,0)
      WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, 8), real_c(0.034262152777778),
                                         real_c(1e-6)); // SE, (1,-1,0)
      break;
   case RefillingModel_T::GradsMomentsRefilling:
      WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, 0), real_c(0.449190200617283),
                                         real_c(1e-6)); // C, (0,0,0)
      WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, 1), real_c(0.098497868441358),
                                         real_c(1e-6)); // N, (0,1,0)
      WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, 2), real_c(0.117664535108025),
                                         real_c(1e-6)); // S, (0,-1,0)
      WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, 3), real_c(0.100871248070988),
                                         real_c(1e-6)); // W, (-1,0,0)
      WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, 4), real_c(0.125037914737654),
                                         real_c(1e-6)); // E, (1,0,0)
      WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, 5), real_c(0.024343253279321),
                                         real_c(1e-6)); // NW, (-1,1,0)
      WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, 6), real_c(0.025234196566358),
                                         real_c(1e-6)); // NE, (1,1,0)
      WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, 7), real_c(0.023984196566358),
                                         real_c(1e-6)); // SW, (-1,-1,0)
      WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, 8), real_c(0.035176586612654),
                                         real_c(1e-6)); // SE, (1,-1,0)
      break;
   default:
      WALBERLA_ABORT("The specified PDF refilling model " << pdfRefillingModel.getModelName()
                                                          << " is not available.\nSomething went wrong!");
   }
}

void verify_someRefillingCells(const PdfRefillingModel& pdfRefillingModel, const Cell& cell,
                               const PdfField_T* const pdfField)
{
   switch (pdfRefillingModel.getModelType())
   {
   case RefillingModel_T::EquilibriumRefilling:
      if (cell.x() == cell_idx_c(1) && cell.y() == cell_idx_c(1) && cell.z() == cell_idx_c(0))
      {
         WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(0)), real_c(0.443777777777778),
                                            real_c(1e-6)); // C, (0,0,0)
         WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(1)), real_c(0.107661111111111),
                                            real_c(1e-6)); // N, (0,1,0)
         WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(2)), real_c(0.114327777777778),
                                            real_c(1e-6)); // S, (0,-1,0)
         WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(3)), real_c(0.101394444444444),
                                            real_c(1e-6)); // W, (-1,0,0)
         WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(4)), real_c(0.121394444444444),
                                            real_c(1e-6)); // E, (1,0,0)
         WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(5)), real_c(0.024602777777778),
                                            real_c(1e-6)); // NW, (-1,1,0)
         WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(6)), real_c(0.029452777777778),
                                            real_c(1e-6)); // NE, (1,1,0)
         WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(7)), real_c(0.026119444444445),
                                            real_c(1e-6)); // SW, (-1,-1,0)
         WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(8)), real_c(0.031269444444445),
                                            real_c(1e-6)); // SE, (1,-1,0)
      }
      else
      {
         if (cell.x() == cell_idx_c(0) && cell.y() == cell_idx_c(2) && cell.z() == cell_idx_c(0))
         {
            WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(0)), real_c(0.443506944444444),
                                               real_c(1e-6)); // C, (0,0,0)
            WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(1)), real_c(0.103629861111111),
                                               real_c(1e-6)); // N, (0,1,0)
            WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(2)), real_c(0.118629861111111),
                                               real_c(1e-6)); // S, (0,-1,0)
            WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(3)), real_c(0.101326736111111),
                                               real_c(1e-6)); // W, (-1,0,0)
            WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(4)), real_c(0.121326736111111),
                                               real_c(1e-6)); // E, (1,0,0)
            WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(5)), real_c(0.023688715277778),
                                               real_c(1e-6)); // NW, (-1,1,0)
            WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(6)), real_c(0.028351215277778),
                                               real_c(1e-6)); // NE, (1,1,0)
            WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(7)), real_c(0.027101215277778),
                                               real_c(1e-6)); // SW, (-1,-1,0)
            WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(8)), real_c(0.032438715277778),
                                               real_c(1e-6)); // SE, (1,-1,0)
         }
         else
         {
            if (cell.x() == cell_idx_c(2) && cell.y() == cell_idx_c(0) && cell.z() == cell_idx_c(0))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(0)), real_c(0.44262037037037),
                                                  real_c(1e-6)); // C,  (0,0,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(1)), real_c(0.104188425925926),
                                                  real_c(1e-6)); // N, (0,1,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(2)), real_c(0.117521759259259),
                                                  real_c(1e-6)); // S, (0,-1,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(3)), real_c(0.095712037037037),
                                                  real_c(1e-6)); // W, (-1,0,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(4)), real_c(0.127934259259259),
                                                  real_c(1e-6)); // E, (1,0,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(5)), real_c(0.022553009259259),
                                                  real_c(1e-6)); // NW, (-1,1,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(6)), real_c(0.030125231481482),
                                                  real_c(1e-6)); // NE, (1,1,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(7)), real_c(0.025403009259259),
                                                  real_c(1e-6)); // SW, (-1,-1,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(8)), real_c(0.033941898148148),
                                                  real_c(1e-6)); // SE, (1,-1,0)
            }
            else { WALBERLA_ABORT("The cells that need refilling are not set correctly."); }
         }
      }
      break;
   case RefillingModel_T::AverageRefilling:
      if (cell.x() == cell_idx_c(1) && cell.y() == cell_idx_c(1) && cell.z() == cell_idx_c(0))
      {
         WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, 0), real_c(0.441666666666667),
                                            real_c(1e-6)); // C, (0,0,0)
         WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, 1), real_c(0.110416666666666),
                                            real_c(1e-6)); // N, (0,1,0)
         WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, 2), real_c(0.110416666666666),
                                            real_c(1e-6)); // S, (0,-1,0)
         WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, 3), real_c(0.095833333333333),
                                            real_c(1e-6)); // W, (-1,0,0)
         WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, 4), real_c(0.119166666666667),
                                            real_c(1e-6)); // E, (1,0,0)
         WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, 5), real_c(0.023958333333333),
                                            real_c(1e-6)); // NW, (-1,1,0)
         WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, 6), real_c(0.032291666666667),
                                            real_c(1e-6)); // NE, (1,1,0)
         WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, 7), real_c(0.033958333333333),
                                            real_c(1e-6)); // SW, (-1,-1,0)
         WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, 8), real_c(0.032291666666667),
                                            real_c(1e-6)); // SE, (1,-1,0)
      }
      else if (cell.x() == cell_idx_c(0) && cell.y() == cell_idx_c(2) && cell.z() == cell_idx_c(0))
      {
         WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, 0), real_c(0.441111111111111),
                                            real_c(1e-6)); // C,  (0,0,0)
         WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, 1), real_c(0.099340277777778),
                                            real_c(1e-6)); // N, (0,1,0)
         WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, 2), real_c(0.116840277777778),
                                            real_c(1e-6)); // S, (0,-1,0)
         WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, 3), real_c(0.098715277777778),
                                            real_c(1e-6)); // W, (-1,0,0)
         WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, 4), real_c(0.108715277777778),
                                            real_c(1e-6)); // E, (1,0,0)
         WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, 5), real_c(0.022256944444444),
                                            real_c(1e-6)); // NW, (-1,1,0)
         WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, 6), real_c(0.042881944444445),
                                            real_c(1e-6)); // NE, (1,1,0)
         WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, 7), real_c(0.035381944444445),
                                            real_c(1e-6)); // SW, (-1,-1,0)
         WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, 8), real_c(0.034756944444445),
                                            real_c(1e-6)); // SE, (1,-1,0)
      }
      else if (cell.x() == cell_idx_c(2) && cell.y() == cell_idx_c(0) && cell.z() == cell_idx_c(0))
      {
         WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, 0), real_c(0.44),
                                            real_c(1e-6)); // C, (0,0,0)
         WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, 1), real_c(0.102708333333333),
                                            real_c(1e-6)); // N, (0,1,0)
         WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, 2), real_c(0.109375),
                                            real_c(1e-6)); // S, (0,-1,0)
         WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, 3), real_c(0.092847222222222),
                                            real_c(1e-6)); // W, (-1,0,0)
         WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, 4), real_c(0.126736111111111),
                                            real_c(1e-6)); // E, (1,0,0)
         WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, 5), real_c(0.021805555555556),
                                            real_c(1e-6)); // NW, (-1,1,0)
         WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, 6), real_c(0.035694444444445),
                                            real_c(1e-6)); // NE, (1,1,0)
         WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, 7), real_c(0.035138888888889),
                                            real_c(1e-6)); // SW, (-1,-1,0)
         WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, 8), real_c(0.035694444444445),
                                            real_c(1e-6)); // SE, (1,-1,0)
      }
      else { WALBERLA_ABORT("The cells that need refilling are not set correctly!"); }
      break;
   case RefillingModel_T::EquilibriumAndNonEquilibriumRefilling:
      if (cell.x() == cell_idx_c(1) && cell.y() == cell_idx_c(1) && cell.z() == cell_idx_c(0))
      {
         WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(0)), real_c(0.452777777777777),
                                            real_c(1e-6)); // C, (0,0,0)
         WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(1)), real_c(0.144811111111111),
                                            real_c(1e-6)); // N, (0,1,0)
         WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(2)), real_c(0.101477777777778),
                                            real_c(1e-6)); // S, (0,-1,0)
         WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(3)), real_c(0.051994444444444),
                                            real_c(1e-6)); // W, (-1,0,0)
         WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(4)), real_c(0.111994444444444),
                                            real_c(1e-6)); // E, (1,0,0)
         WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(5)), real_c(0.021427777777778),
                                            real_c(1e-6)); // NW, (-1,1,0)
         WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(6)), real_c(0.020377777777778),
                                            real_c(1e-6)); // NE, (1,1,0)
         WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(7)), real_c(0.062044444444445),
                                            real_c(1e-6)); // SW, (-1,-1,0)
         WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(8)), real_c(0.033094444444445),
                                            real_c(1e-6)); // SE, (1,-1,0)
      }
      else
      {
         if (cell.x() == cell_idx_c(0) && cell.y() == cell_idx_c(2) && cell.z() == cell_idx_c(0))
         {
            WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(0)), real_c(0.443506944444444),
                                               real_c(1e-6)); // C, (0,0,0)
            WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(1)), real_c(0.103629861111111),
                                               real_c(1e-6)); // N, (0,1,0)
            WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(2)), real_c(0.118629861111111),
                                               real_c(1e-6)); // S, (0,-1,0)
            WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(3)), real_c(0.101326736111111),
                                               real_c(1e-6)); // W, (-1,0,0)
            WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(4)), real_c(0.121326736111111),
                                               real_c(1e-6)); // E, (1,0,0)
            WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(5)), real_c(0.023688715277778),
                                               real_c(1e-6)); // NW, (-1,1,0)
            WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(6)), real_c(0.028351215277778),
                                               real_c(1e-6)); // NE, (1,1,0)
            WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(7)), real_c(0.027101215277778),
                                               real_c(1e-6)); // SW, (-1,-1,0)
            WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(8)), real_c(0.032438715277778),
                                               real_c(1e-6)); // SE, (1,-1,0)
         }
         else
         {
            if (cell.x() == cell_idx_c(2) && cell.y() == cell_idx_c(0) && cell.z() == cell_idx_c(0))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(0)), real_c(0.43522037037037),
                                                  real_c(1e-6)); // C, (0,0,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(1)), real_c(0.112788425925926),
                                                  real_c(1e-6)); // N, (0,1,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(2)), real_c(0.046121759259259),
                                                  real_c(1e-6)); // S, (0,-1,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(3)), real_c(0.078962037037037),
                                                  real_c(1e-6)); // W, (-1,0,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(4)), real_c(0.181184259259259),
                                                  real_c(1e-6)); // E, (1,0,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(5)), real_c(0.021353009259259),
                                                  real_c(1e-6)); // NW, (-1,1,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(6)), real_c(0.008175231481481),
                                                  real_c(1e-6)); // NE, (1,1,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(7)), real_c(0.078453009259259),
                                                  real_c(1e-6)); // SW, (-1,-1,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(8)), real_c(0.037741898148148),
                                                  real_c(1e-6)); // SE, (1,-1,0)
            }
            else { WALBERLA_ABORT("The cells that need refilling are not set correctly!"); }
         }
      }
      break;
   case RefillingModel_T::ExtrapolationRefilling:
      if (cell.x() == cell_idx_c(1) && cell.y() == cell_idx_c(1) && cell.z() == cell_idx_c(0))
      {
         WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, 0), real_c(0.452777777777777),
                                            real_c(1e-6)); // C, (0,0,0)
         WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, 1), real_c(0.141527777777778),
                                            real_c(1e-6)); // N, (0,1,0)
         WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, 2), real_c(0.104861111111111),
                                            real_c(1e-6)); // S, (0,-1,0)
         WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, 3), real_c(0.128611111111111),
                                            real_c(1e-6)); // W, (-1,0,0)
         WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, 4), real_c(0.035277777777778),
                                            real_c(1e-6)); // E, (1,0,0)
         WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, 5), real_c(0.037986111111111),
                                            real_c(1e-6)); // NW, (-1,1,0)
         WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, 6), real_c(0.002152777777778),
                                            real_c(1e-6)); // NE, (1,1,0)
         WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, 7), real_c(0.083819444444445),
                                            real_c(1e-6)); // SW, (-1,-1,0)
         WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, 8), real_c(0.012986111111111),
                                            real_c(1e-6)); // SE, (1,-1,0)
      }
      else
      {
         if (cell.x() == cell_idx_c(0) && cell.y() == cell_idx_c(2) && cell.z() == cell_idx_c(0))
         {
            WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(0)), real_c(0.443506944444444),
                                               real_c(1e-6)); // C, (0,0,0)
            WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(1)), real_c(0.103629861111111),
                                               real_c(1e-6)); // N, (0,1,0)
            WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(2)), real_c(0.118629861111111),
                                               real_c(1e-6)); // S, (0,-1,0)
            WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(3)), real_c(0.101326736111111),
                                               real_c(1e-6)); // W, (-1,0,0)
            WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(4)), real_c(0.121326736111111),
                                               real_c(1e-6)); // E, (1,0,0)
            WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(5)), real_c(0.023688715277778),
                                               real_c(1e-6)); // NW, (-1,1,0)
            WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(6)), real_c(0.028351215277778),
                                               real_c(1e-6)); // NE, (1,1,0)
            WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(7)), real_c(0.027101215277778),
                                               real_c(1e-6)); // SW, (-1,-1,0)
            WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(8)), real_c(0.032438715277778),
                                               real_c(1e-6)); // SE, (1,-1,0)
         }
         else
         {
            if (cell.x() == cell_idx_c(2) && cell.y() == cell_idx_c(0) && cell.z() == cell_idx_c(0))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(0)), real_c(0.429444444444444),
                                                  real_c(1e-6)); // C, (0,0,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(1)), real_c(0.076527777777778),
                                                  real_c(1e-6)); // N, (0,1,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(2)), real_c(0.083194444444444),
                                                  real_c(1e-6)); // S, (0,-1,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(3)), real_c(0.066111111111111),
                                                  real_c(1e-6)); // W, (-1,0,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(4)), real_c(0.196111111111111),
                                                  real_c(1e-6)); // E, (1,0,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(5)), real_c(0.011319444444444),
                                                  real_c(1e-6)); // NW, (-1,1,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(6)), real_c(0.001319444444444),
                                                  real_c(1e-6)); // NE, (1,1,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(7)), real_c(0.082986111111111),
                                                  real_c(1e-6)); // SW, (-1,-1,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(8)), real_c(0.052986111111111),
                                                  real_c(1e-6)); // SE, (1,-1,0)
            }
            else { WALBERLA_ABORT("The cells that need refilling are not set correctly!"); }
         }
      }
      break;
   case RefillingModel_T::GradsMomentsRefilling:
      if (cell.x() == cell_idx_c(1) && cell.y() == cell_idx_c(1) && cell.z() == cell_idx_c(0))
      {
         WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, 0), real_c(0.46353086419753),
                                            real_c(1e-6)); // C, (0,0,0)
         WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, 1), real_c(0.110747530864197),
                                            real_c(1e-6)); // N, (0,1,0)
         WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, 2), real_c(0.117414197530864),
                                            real_c(1e-6)); // S, (0,-1,0)
         WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, 3), real_c(0.09336975308642),
                                            real_c(1e-6)); // W, (-1,0,0)
         WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, 4), real_c(0.11336975308642),
                                            real_c(1e-6)); // E, (1,0,0)
         WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, 5), real_c(0.023522530864198),
                                            real_c(1e-6)); // NW, (-1,1,0)
         WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, 6), real_c(0.02559475308642),
                                            real_c(1e-6)); // NE, (1,1,0)
         WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, 7), real_c(0.022261419753087),
                                            real_c(1e-6)); // SW, (-1,-1,0)
         WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, 8), real_c(0.030189197530864),
                                            real_c(1e-6)); // SE, (1,-1,0)
      }
      else
      {
         if (cell.x() == cell_idx_c(0) && cell.y() == cell_idx_c(2) && cell.z() == cell_idx_c(0))
         {
            WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(0)), real_c(0.443506944444444),
                                               real_c(1e-6)); // C, (0,0,0)
            WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(1)), real_c(0.103629861111111),
                                               real_c(1e-6)); // N, (0,1,0)
            WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(2)), real_c(0.118629861111111),
                                               real_c(1e-6)); // S, (0,-1,0)
            WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(3)), real_c(0.101326736111111),
                                               real_c(1e-6)); // W, (-1,0,0)
            WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(4)), real_c(0.121326736111111),
                                               real_c(1e-6)); // E, (1,0,0)
            WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(5)), real_c(0.022994270833333),
                                               real_c(1e-6)); // NW, (-1,1,0)
            WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(6)), real_c(0.029045659722222),
                                               real_c(1e-6)); // NE, (1,1,0)
            WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(7)), real_c(0.027795659722222),
                                               real_c(1e-6)); // SW, (-1,-1,0)
            WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(8)), real_c(0.031744270833333),
                                               real_c(1e-6)); // SE, (1,-1,0)
         }
         else
         {
            if (cell.x() == cell_idx_c(2) && cell.y() == cell_idx_c(0) && cell.z() == cell_idx_c(0))
            {
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(0)), real_c(0.454143004115226),
                                                  real_c(1e-6)); // C, (0,0,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(1)), real_c(0.104291306584362),
                                                  real_c(1e-6)); // N, (0,1,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(2)), real_c(0.117624639917695),
                                                  real_c(1e-6)); // S, (0,-1,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(3)), real_c(0.092728497942387),
                                                  real_c(1e-6)); // W, (-1,0,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(4)), real_c(0.124950720164609),
                                                  real_c(1e-6)); // E, (1,0,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(5)), real_c(0.02389045781893),
                                                  real_c(1e-6)); // NW, (-1,1,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(6)), real_c(0.025907124485597),
                                                  real_c(1e-6)); // NE, (1,1,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(7)), real_c(0.021184902263375),
                                                  real_c(1e-6)); // SW, (-1,-1,0)
               WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(8)), real_c(0.035279346707819),
                                                  real_c(1e-6)); // SE, (1,-1,0)
            }
            else { WALBERLA_ABORT("The cells that need refilling are not set correctly!"); }
         }
      }
      break;
   default:
      WALBERLA_ABORT("The specified PDF refilling model " << pdfRefillingModel.getModelName()
                                                          << " is not available.\nSomething went wrong!");
   }
}

void verify_straightRefillingCells(const PdfRefillingModel& pdfRefillingModel, const Cell& cell,
                                   const PdfField_T* const pdfField)
{
   switch (pdfRefillingModel.getModelType())
   {
   case RefillingModel_T::EquilibriumRefilling:
      if (cell.x() == cell_idx_c(1) && cell.y() == cell_idx_c(1) && cell.z() == cell_idx_c(0))
      {
         WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(0)), real_c(0.444259259259259),
                                            real_c(1e-6)); // C, (0,0,0)
         WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(1)), real_c(0.107781481481481),
                                            real_c(1e-6)); // N, (0,1,0)
         WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(2)), real_c(0.114448148148148),
                                            real_c(1e-6)); // S, (0,-1,0)
         WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(3)), real_c(0.106709259259259),
                                            real_c(1e-6)); // W, (-1,0,0)
         WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(4)), real_c(0.115598148148148),
                                            real_c(1e-6)); // E, (1,0,0)
         WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(5)), real_c(0.025889814814815),
                                            real_c(1e-6)); // NW, (-1,1,0)
         WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(6)), real_c(0.02804537037037),
                                            real_c(1e-6)); // NE, (1,1,0)
         WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(7)), real_c(0.027489814814815),
                                            real_c(1e-6)); // SW, (-1,-1,0)
         WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(8)), real_c(0.029778703703704),
                                            real_c(1e-6)); // SE, (1,-1,0)
      }
      else { WALBERLA_ABORT("The cells that need refilling are not set correctly!"); }
      break;
   case RefillingModel_T::AverageRefilling:
      if (cell.x() == cell_idx_c(1) && cell.y() == cell_idx_c(1) && cell.z() == cell_idx_c(0))
      {
         WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, 0), real_c(0.442222222222222),
                                            real_c(1e-6)); // C, (0,0,0)
         WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, 1), real_c(0.110555555555555),
                                            real_c(1e-6)); // N, (0,1,0)
         WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, 2), real_c(0.110555555555555),
                                            real_c(1e-6)); // S, (0,-1,0)
         WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, 3), real_c(0.101111111111111),
                                            real_c(1e-6)); // W, (-1,0,0)
         WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, 4), real_c(0.113333333333333),
                                            real_c(1e-6)); // E, (1,0,0)
         WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, 5), real_c(0.025277777777778),
                                            real_c(1e-6)); // NW, (-1,1,0)
         WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, 6), real_c(0.030833333333333),
                                            real_c(1e-6)); // NE, (1,1,0)
         WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, 7), real_c(0.035277777777778),
                                            real_c(1e-6)); // SW, (-1,-1,0)
         WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, 8), real_c(0.030833333333333),
                                            real_c(1e-6)); // SE, (1,-1,0)
      }
      else { WALBERLA_ABORT("The cells that need refilling are not set correctly!"); }
      break;
   case RefillingModel_T::EquilibriumAndNonEquilibriumRefilling:
      if (cell.x() == cell_idx_c(1) && cell.y() == cell_idx_c(1) && cell.z() == cell_idx_c(0))
      {
         WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(0)), real_c(0.449659259259259),
                                            real_c(1e-6)); // C, (0,0,0)
         WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(1)), real_c(0.10168148148148),
                                            real_c(1e-6)); // N, (0,1,0)
         WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(2)), real_c(0.188348148148148),
                                            real_c(1e-6)); // S, (0,-1,0)
         WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(3)), real_c(0.121459259259259),
                                            real_c(1e-6)); // W, (-1,0,0)
         WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(4)), real_c(0.060348148148148),
                                            real_c(1e-6)); // E, (1,0,0)
         WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(5)), real_c(0.027489814814815),
                                            real_c(1e-6)); // NW, (-1,1,0)
         WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(6)), real_c(0.050095370370371),
                                            real_c(1e-6)); // NE, (1,1,0)
         WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(7)), real_c(-0.025460185185185),
                                            real_c(1e-6)); // SW, (-1,-1,0)
         WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(8)), real_c(0.026378703703704),
                                            real_c(1e-6)); // SE, (1,-1,0)
      }
      else { WALBERLA_ABORT("The cells that need refilling are not set correctly!"); }
      break;
   case RefillingModel_T::ExtrapolationRefilling:
      if (cell.x() == cell_idx_c(1) && cell.y() == cell_idx_c(1) && cell.z() == cell_idx_c(0))
      {
         WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(0)), real_c(0.441111111111112),
                                            real_c(1e-6)); // C, (0,0,0)
         WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(1)), real_c(0.128194444444443),
                                            real_c(1e-6)); // N, (0,1,0)
         WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(2)), real_c(0.154861111111111),
                                            real_c(1e-6)); // S, (0,-1,0)
         WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(3)), real_c(0.094861111111111),
                                            real_c(1e-6)); // W, (-1,0,0)
         WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(4)), real_c(0.098194444444445),
                                            real_c(1e-6)); // E, (1,0,0)
         WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(5)), real_c(0.025694444444445),
                                            real_c(1e-6)); // NW, (-1,1,0)
         WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(6)), real_c(0.069027777777778),
                                            real_c(1e-6)); // NE, (1,1,0)
         WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(7)), real_c(-0.037638888888889),
                                            real_c(1e-6)); // SW, (-1,-1,0)
         WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(8)), real_c(0.025694444444444),
                                            real_c(1e-6)); // SE, (1,-1,0)
      }
      else { WALBERLA_ABORT("The cells that need refilling are not set correctly!"); }
      break;
   case RefillingModel_T::GradsMomentsRefilling:
      if (cell.x() == cell_idx_c(1) && cell.y() == cell_idx_c(1) && cell.z() == cell_idx_c(0))
      {
         WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(0)), real_c(0.465658436213991),
                                            real_c(1e-6)); // C, (0,0,0)
         WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(1)), real_c(0.113131275720165),
                                            real_c(1e-6)); // N, (0,1,0)
         WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(2)), real_c(0.119797942386831),
                                            real_c(1e-6)); // S, (0,-1,0)
         WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(3)), real_c(0.096009670781893),
                                            real_c(1e-6)); // W, (-1,0,0)
         WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(4)), real_c(0.104898559670782),
                                            real_c(1e-6)); // E, (1,0,0)
         WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(5)), real_c(0.023677880658436),
                                            real_c(1e-6)); // NW, (-1,1,0)
         WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(6)), real_c(0.024907510288066),
                                            real_c(1e-6)); // NE, (1,1,0)
         WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(7)), real_c(0.02435195473251),
                                            real_c(1e-6)); // SW, (-1,-1,0)
         WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(pdfField->get(cell, cell_idx_c(8)), real_c(0.027566769547325),
                                            real_c(1e-6)); // SE, (1,-1,0)
      }
      else { WALBERLA_ABORT("The cells that need refilling are not set correctly!"); }
      break;
   default:
      WALBERLA_ABORT("The specified PDF refilling model " << pdfRefillingModel.getModelName()
                                                          << " is not available.\nSomething went wrong!");
   }
}

void runSimulation(const PdfRefillingModel& pdfRefillingModel, const initCallback_T& fillInitializationSets,
                   const verifyCallback_T& verifyResults)
{
   // define the domain size (5x3x1)
   const Vector3< uint_t > numBlocks(uint_c(1));
   const Vector3< uint_t > domainSize(uint_c(5), uint_c(3), uint_c(1));
   const Vector3< uint_t > cellsPerBlock(domainSize[0] / numBlocks[0], domainSize[1] / numBlocks[1],
                                         domainSize[2] / numBlocks[2]);

   // create block grid
   const std::shared_ptr< StructuredBlockForest > blockForest =
      blockforest::createUniformBlockGrid(numBlocks[0], numBlocks[1], numBlocks[2],             // blocks
                                          cellsPerBlock[0], cellsPerBlock[1], cellsPerBlock[2], // cells
                                          real_c(1.0),                                          // dx
                                          true,                                                 // one block per process
                                          true, true, false);                                   // periodicity

   // relaxation rate (used by GradsMomentsRefilling)
   real_t omega = real_c(1.8);

   // create lattice model
   LatticeModel_T latticeModel = LatticeModel_T(lbm::collision_model::SRT(omega));

   // add PDF field (and a copy that stores the original PDF field)
   const BlockDataID pdfFieldID =
      lbm::addPdfFieldToStorage(blockForest, "PDF source field", latticeModel, uint_c(1), field::fzyx);

   // field to store PDFs before refilling, i.e., to test if non-refilling cells are modified
   const BlockDataID pdfOrgFieldID =
      lbm::addPdfFieldToStorage(blockForest, "PDF destination field", latticeModel, uint_c(1), field::fzyx);

   // add fill level field (MUST be initialized with 1, i.e., fluid everywhere for this test; otherwise the fluid
   // flag is not detected below by initFlagsFromFillLevel())
   BlockDataID fillFieldID =
      field::addToStorage< ScalarField_T >(blockForest, "Fill level field", real_c(1), field::fzyx, uint_c(1));

   // define the initial properties of the cells ( gas, liquid, interface )
   std::set< Cell > gasCells;
   std::set< Cell > interfaceCells;
   std::set< Cell > liquidCells;
   std::set< Cell > conversionCells;
   fillInitializationSets(gasCells, interfaceCells, liquidCells, conversionCells);

   real_t initDensity = real_c(1.0);

   for (auto blockIt = blockForest->begin(); blockIt != blockForest->end(); ++blockIt)
   {
      ScalarField_T* const fillField = blockIt->getData< ScalarField_T >(fillFieldID);
      PdfField_T* const pdfField     = blockIt->getData< PdfField_T >(pdfFieldID);
      PdfField_T* pdfOrgField        = blockIt->getData< PdfField_T >(pdfOrgFieldID);

      std::vector< std::pair< Cell, Vector3< real_t > > > velocities = {
         { Cell{ cell_idx_c(0), cell_idx_c(0), cell_idx_c(0) },
           Vector3< real_t >{ real_c(0.05), real_c(-0.10), real_c(0.00) } },
         { Cell{ cell_idx_c(1), cell_idx_c(0), cell_idx_c(0) },
           Vector3< real_t >{ real_c(0.00), real_c(-0.05), real_c(0.00) } },
         { Cell{ cell_idx_c(2), cell_idx_c(0), cell_idx_c(0) },
           Vector3< real_t >{ real_c(0.00), real_c(0.00), real_c(0.00) } },
         { Cell{ cell_idx_c(3), cell_idx_c(0), cell_idx_c(0) },
           Vector3< real_t >{ real_c(0.05), real_c(-0.10), real_c(0.00) } },
         { Cell{ cell_idx_c(4), cell_idx_c(0), cell_idx_c(0) },
           Vector3< real_t >{ real_c(0.00), real_c(-0.05), real_c(0.00) } },
         { Cell{ cell_idx_c(0), cell_idx_c(1), cell_idx_c(0) },
           Vector3< real_t >{ real_c(0.10), real_c(-0.05), real_c(0.00) } },
         { Cell{ cell_idx_c(1), cell_idx_c(1), cell_idx_c(0) },
           Vector3< real_t >{ real_c(0.05), real_c(-0.10), real_c(0.00) } },
         { Cell{ cell_idx_c(2), cell_idx_c(1), cell_idx_c(0) },
           Vector3< real_t >{ real_c(0.10), real_c(0.00), real_c(0.00) } },
         { Cell{ cell_idx_c(3), cell_idx_c(1), cell_idx_c(0) },
           Vector3< real_t >{ real_c(0.10), real_c(-0.05), real_c(0.00) } },
         { Cell{ cell_idx_c(4), cell_idx_c(1), cell_idx_c(0) },
           Vector3< real_t >{ real_c(0.05), real_c(-0.10), real_c(0.00) } },
         { Cell{ cell_idx_c(0), cell_idx_c(2), cell_idx_c(0) },
           Vector3< real_t >{ real_c(0.05), real_c(0.00), real_c(0.00) } },
         { Cell{ cell_idx_c(1), cell_idx_c(2), cell_idx_c(0) },
           Vector3< real_t >{ real_c(0.05), real_c(0.00), real_c(0.00) } },
         { Cell{ cell_idx_c(2), cell_idx_c(2), cell_idx_c(0) },
           Vector3< real_t >{ real_c(0.00), real_c(0.00), real_c(0.00) } },
         { Cell{ cell_idx_c(3), cell_idx_c(2), cell_idx_c(0) },
           Vector3< real_t >{ real_c(0.05), real_c(0.00), real_c(0.00) } },
         { Cell{ cell_idx_c(4), cell_idx_c(2), cell_idx_c(0) },
           Vector3< real_t >{ real_c(0.05), real_c(0.00), real_c(0.00) } }
      };

      for (auto velocity = velocities.begin(); velocity != velocities.end(); velocity++)
      {
         pdfField->setDensityAndVelocity(velocity->first, velocity->second, initDensity);
         pdfOrgField->setDensityAndVelocity(velocity->first, velocity->second, initDensity);
      }

      // modify the PDFs to achieve inequality for the equalAndNonEqualRefilling
      for (auto d = Stencil_T::begin(); d != Stencil_T::end(); ++d)
      {
         using D = stencil::Direction;
         D dir   = d.direction();
         Cell cell1(cell_idx_c(2), cell_idx_c(2), cell_idx_c(0));
         Cell cell2(cell_idx_c(3), cell_idx_c(0), cell_idx_c(0));
         Cell cell3(cell_idx_c(4), cell_idx_c(1), cell_idx_c(0));
         Cell cell4(cell_idx_c(3), cell_idx_c(1), cell_idx_c(0));
         Cell cell5(cell_idx_c(4), cell_idx_c(2), cell_idx_c(0));
         real_t PDFOffset(real_c(0.03));

         if (dir == D::SW)
         {
            pdfField->get(cell1, dir) += PDFOffset;
            pdfOrgField->get(cell1, dir) += PDFOffset;
            pdfField->get(cell3, dir) += PDFOffset;
            pdfOrgField->get(cell3, dir) += PDFOffset;
            pdfField->get(cell4, dir) += PDFOffset;
            pdfOrgField->get(cell4, dir) += PDFOffset;
         }
         if (dir == D::E)
         {
            pdfField->get(cell1, dir) -= PDFOffset;
            pdfOrgField->get(cell1, dir) -= PDFOffset;
            pdfField->get(cell3, dir) -= PDFOffset;
            pdfOrgField->get(cell3, dir) -= PDFOffset;
            pdfField->get(cell5, dir) -= PDFOffset;
            pdfOrgField->get(cell5, dir) -= PDFOffset;
         }
         if (dir == D::S)
         {
            pdfField->get(cell2, dir) -= PDFOffset;
            pdfOrgField->get(cell2, dir) -= PDFOffset;
            pdfField->get(cell3, dir) -= PDFOffset;
            pdfOrgField->get(cell3, dir) -= PDFOffset;
            pdfField->get(cell4, dir) -= PDFOffset;
            pdfOrgField->get(cell4, dir) -= PDFOffset;
         }
         if (dir == D::NE)
         {
            pdfField->get(cell2, dir) += PDFOffset;
            pdfOrgField->get(cell2, dir) += PDFOffset;
            pdfField->get(cell3, dir) += PDFOffset;
            pdfOrgField->get(cell3, dir) += PDFOffset;
            pdfField->get(cell5, dir) += PDFOffset;
            pdfOrgField->get(cell5, dir) += PDFOffset;
         }
      }

      pdfOrgField = pdfField->clone();

      // initialize fill level
      WALBERLA_FOR_ALL_CELLS(fillFieldIt, fillField, {
         if (gasCells.find(fillFieldIt.cell()) != gasCells.end()) { *fillFieldIt = real_c(0.0); }
         else
         {
            if (interfaceCells.find(fillFieldIt.cell()) != interfaceCells.end()) { *fillFieldIt = real_c(0.5); }
            else
            {
               if (liquidCells.find(fillFieldIt.cell()) != liquidCells.end()) { *fillFieldIt = real_c(1.0); }
            }
         }
      }) // WALBERLA_FOR_ALL_CELLS
   }

   // add boundary handling
   const std::shared_ptr< FreeSurfaceBoundaryHandling_T > freeSurfaceBoundaryHandling =
      std::make_shared< FreeSurfaceBoundaryHandling_T >(blockForest, pdfFieldID, fillFieldID);
   const BlockDataID flagFieldID          = freeSurfaceBoundaryHandling->getFlagFieldID();
   const FlagInfo< FlagField_T > flagInfo = freeSurfaceBoundaryHandling->getFlagInfo();
   freeSurfaceBoundaryHandling->initFlagsFromFillLevel();

   // create timeloop
   uint_t timesteps = uint_c(1);
   SweepTimeloop timeloop(blockForest, timesteps);

   Communication_T(blockForest, pdfFieldID, pdfOrgFieldID, flagFieldID)();

   using geometryHandler   = SurfaceGeometryHandler< LatticeModel_T, FlagField_T, ScalarField_T, VectorField_T >;
   auto geometryHandlerPtr = std::make_shared< geometryHandler >(blockForest, freeSurfaceBoundaryHandling, fillFieldID,
                                                                 "FiniteDifferenceMethod", false, false, real_c(0.0));

   for (auto blockIt = blockForest->begin(); blockIt != blockForest->end(); ++blockIt)
   {
      ScalarField_T* const fillField = blockIt->getData< ScalarField_T >(fillFieldID);
      FlagField_T* const flagField   = blockIt->getData< FlagField_T >(flagFieldID);

      // convert a cell from gas to interface
      WALBERLA_FOR_ALL_CELLS(fillFieldIt, fillField, flagFieldIt, flagField, {
         if (conversionCells.find(flagFieldIt.cell()) != conversionCells.end())
         {
            // fill one cell entirely and create an interface cell from it; the fill level must be 1.1 to obtain the
            // normal as used in the manual computations of the expected solutions
            (*fillFieldIt) = real_c(1.1);

            // mark cell to be reinitialized
            field::removeFlag(flagFieldIt, flagInfo.gasFlag);
            field::addFlag(flagFieldIt, flagInfo.interfaceFlag);
            field::addFlag(flagFieldIt, flagInfo.convertedFlag);
            field::addFlag(flagFieldIt, flagInfo.convertFromGasToInterfaceFlag);
         }
      }) // WALBERLA_FOR_ALL_CELLS

      // perform refilling
      switch (pdfRefillingModel.getModelType())
      { // the scope for each "case" is required since variables are defined within "case"
      case PdfRefillingModel::RefillingModel::EquilibriumRefilling: {
         EquilibriumRefillingSweep< LatticeModel_T, FlagField_T > equilibriumRefillingSweep(pdfFieldID, flagFieldID,
                                                                                            flagInfo, true);
         equilibriumRefillingSweep(blockIt.get());
         break;
      }

      case PdfRefillingModel::RefillingModel::AverageRefilling: {
         AverageRefillingSweep< LatticeModel_T, FlagField_T > averageRefillingSweep(pdfFieldID, flagFieldID, flagInfo,
                                                                                    true);
         averageRefillingSweep(blockIt.get());
         break;
      }

      case PdfRefillingModel::RefillingModel::EquilibriumAndNonEquilibriumRefilling: {
         EquilibriumAndNonEquilibriumRefillingSweep< LatticeModel_T, FlagField_T, ScalarField_T, VectorField_T >
            equilibriumAndNonEquilibriumRefillingSweep(pdfFieldID, flagFieldID, fillFieldID, flagInfo, uint_c(3), true);
         equilibriumAndNonEquilibriumRefillingSweep(blockIt.get());
         break;
      }

      case PdfRefillingModel::RefillingModel::ExtrapolationRefilling: {
         ExtrapolationRefillingSweep< LatticeModel_T, FlagField_T, ScalarField_T, VectorField_T >
            extrapolationRefillingSweep(pdfFieldID, flagFieldID, fillFieldID, flagInfo, uint_c(3), true);
         extrapolationRefillingSweep(blockIt.get());
         break;
      }

      case PdfRefillingModel::RefillingModel::GradsMomentsRefilling: {
         GradsMomentsRefillingSweep< LatticeModel_T, FlagField_T > gradsMomentsRefillingSweep(pdfFieldID, flagFieldID,
                                                                                              flagInfo, omega, true);
         gradsMomentsRefillingSweep(blockIt.get());
         break;
      }
      default:
         WALBERLA_ABORT("The specified pdf refilling model is not available.");
      }
   }

   for (auto blockIt = blockForest->begin(); blockIt != blockForest->end(); ++blockIt)
   {
      const PdfField_T* const pdfField     = blockIt->getData< const PdfField_T >(pdfFieldID);
      const PdfField_T* const pdfOrgField  = blockIt->getData< const PdfField_T >(pdfOrgFieldID);
      const FlagField_T* const flagField   = blockIt->getData< const FlagField_T >(flagFieldID);
      const ScalarField_T* const fillField = blockIt->getData< const ScalarField_T >(fillFieldID);

      WALBERLA_FOR_ALL_CELLS(pdfFieldIt, pdfField, pdfFieldOrgIt, pdfOrgField, flagFieldIt, flagField, fillFieldIt,
                             fillField, {
                                if (conversionCells.find(flagFieldIt.cell()) == conversionCells.end())
                                {
                                   // check cells that were not converted if their PDFs changed
                                   WALBERLA_CHECK_FLOAT_EQUAL(pdfFieldIt[0], pdfFieldOrgIt[0]); // C, (0,0,0)
                                   WALBERLA_CHECK_FLOAT_EQUAL(pdfFieldIt[1], pdfFieldOrgIt[1]); // N, (0,1,0)
                                   WALBERLA_CHECK_FLOAT_EQUAL(pdfFieldIt[2], pdfFieldOrgIt[2]); // S, (0,-1,0)
                                   WALBERLA_CHECK_FLOAT_EQUAL(pdfFieldIt[3], pdfFieldOrgIt[3]); // W, (-1,0,0)
                                   WALBERLA_CHECK_FLOAT_EQUAL(pdfFieldIt[4], pdfFieldOrgIt[4]); // E, (1,0,0)
                                   WALBERLA_CHECK_FLOAT_EQUAL(pdfFieldIt[5], pdfFieldOrgIt[5]); // NW, (-1,1,0)
                                   WALBERLA_CHECK_FLOAT_EQUAL(pdfFieldIt[6], pdfFieldOrgIt[6]); // NE, (1,1,0)
                                   WALBERLA_CHECK_FLOAT_EQUAL(pdfFieldIt[7], pdfFieldOrgIt[7]); // SW, (-1,-1,0)
                                   WALBERLA_CHECK_FLOAT_EQUAL(pdfFieldIt[8], pdfFieldOrgIt[8]); // SE, (1,-1,0)
                                }
                                else
                                {
                                   // check cells converted from gas to interface
                                   verifyResults(pdfRefillingModel, pdfFieldIt.cell(), pdfField);
                                }
                             }) // WALBERLA_FOR_ALL_CELLS
   }

   MPIManager::instance()->resetMPI();
}

} // namespace PdfRefillingTest
} // namespace free_surface
} // namespace walberla

int main(int argc, char** argv) { return walberla::free_surface::PdfRefillingTest::main(argc, argv); }