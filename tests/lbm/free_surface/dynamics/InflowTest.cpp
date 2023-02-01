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
//! \file InflowTest.cpp
//! \ingroup lbm/free_surface/dynamics
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//! \brief Test inflow boundary condition.
//!
//! Set inflow boundaries and initialize gas everywhere. After performing one time step, it is evaluated whether gas
//! cells have been converted to interface and initialized correctly according to the neighboring inflow cells.
//======================================================================================================================

#include "blockforest/Initialization.h"
#include "blockforest/communication/UniformBufferedScheme.h"

#include "core/Environment.h"
#include "core/debug/TestSubsystem.h"

#include "field/AddToStorage.h"
#include "field/adaptors/AdaptorCreators.h"
#include "field/communication/PackInfo.h"

#include "lbm/field/Adaptors.h"
#include "lbm/field/AddToStorage.h"
#include "lbm/free_surface//FlagInfo.h"
#include "lbm/free_surface/boundary/FreeSurfaceBoundaryHandling.h"
#include "lbm/free_surface/dynamics/StreamReconstructAdvectSweep.h"
#include "lbm/free_surface/dynamics/SurfaceDynamicsHandler.h"
#include "lbm/free_surface/surface_geometry/SurfaceGeometryHandler.h"
#include "lbm/lattice_model/CollisionModel.h"
#include "lbm/lattice_model/D3Q19.h"
#include "lbm/lattice_model/ForceModel.h"

#include "timeloop/SweepTimeloop.h"

#include <algorithm>

namespace walberla
{
namespace free_surface
{
namespace InflowTest
{
// define types
using Flag_T         = uint32_t;
using ScalarField_T  = GhostLayerField< real_t, 1 >;
using VectorField_T  = GhostLayerField< Vector3< real_t >, 1 >;
using FlagField_T    = FlagField< Flag_T >;
using LatticeModel_T = lbm::D3Q19< lbm::collision_model::SRT, true, lbm::force_model::None, 2 >;

void testInflow()
{
   // define types
   using Stencil_T                     = typename LatticeModel_T::Stencil;
   using PdfField_T                    = lbm::PdfField< LatticeModel_T >;
   using FreeSurfaceBoundaryHandling_T = FreeSurfaceBoundaryHandling< LatticeModel_T, FlagField_T, ScalarField_T >;

   // define the domain size
   const Vector3< uint_t > numBlocks(uint_c(1));
   const Vector3< uint_t > domainSize(uint_c(10), uint_c(3), uint_c(1));
   const Vector3< uint_t > cellsPerBlock(domainSize[0] / numBlocks[0], domainSize[1] / numBlocks[1],
                                         domainSize[2] / numBlocks[2]);

   // create block grid
   const std::shared_ptr< StructuredBlockForest > blockForest =
      blockforest::createUniformBlockGrid(numBlocks[0], numBlocks[1], numBlocks[2],             // blocks
                                          cellsPerBlock[0], cellsPerBlock[1], cellsPerBlock[2], // cells
                                          real_c(1.0),                                          // dx
                                          true,                                                 // one block per process
                                          false, false, true);                                  // periodicity

   real_t relaxRate = real_c(0.51);

   // create lattice model with omega=0.51
   LatticeModel_T latticeModel(relaxRate);

   // add fields
   BlockDataID pdfFieldID = lbm::addPdfFieldToStorage(blockForest, "PDF field", latticeModel, field::fzyx);
   BlockDataID fillFieldID =
      field::addToStorage< ScalarField_T >(blockForest, "Fill levels", real_c(0.0), field::fzyx, uint_c(2));
   BlockDataID forceDensityFieldID = field::addToStorage< VectorField_T >(
      blockForest, "Force density field", Vector3< real_t >(real_c(0)), field::fzyx, uint_c(1));
   BlockDataID densityAdaptor = field::addFieldAdaptor< typename lbm::Adaptor< LatticeModel_T >::Density >(
      blockForest, pdfFieldID, "DensityAdaptor");
   BlockDataID velocityAdaptor = field::addFieldAdaptor< typename lbm::Adaptor< LatticeModel_T >::VelocityVector >(
      blockForest, pdfFieldID, "VelocityAdaptor");

   // add boundary handling
   const std::shared_ptr< FreeSurfaceBoundaryHandling_T > freeSurfaceBoundaryHandling =
      std::make_shared< FreeSurfaceBoundaryHandling_T >(blockForest, pdfFieldID, fillFieldID);
   const BlockDataID flagFieldID                                      = freeSurfaceBoundaryHandling->getFlagFieldID();
   const typename FreeSurfaceBoundaryHandling_T::FlagInfo_T& flagInfo = freeSurfaceBoundaryHandling->getFlagInfo();

   // set inflow boundary conditions in some cells of the western domain border
   freeSurfaceBoundaryHandling->setInflowInCell(Cell(cell_idx_c(0), cell_idx_c(-1), cell_idx_c(0)),
                                                Vector3< real_t >(real_c(0), real_c(0.01), real_c(0)));
   freeSurfaceBoundaryHandling->setInflowInCell(Cell(cell_idx_c(-1), cell_idx_c(0), cell_idx_c(0)),
                                                Vector3< real_t >(real_c(0.01), real_c(0), real_c(0)));
   freeSurfaceBoundaryHandling->setInflowInCell(Cell(cell_idx_c(-1), cell_idx_c(1), cell_idx_c(0)),
                                                Vector3< real_t >(real_c(0.02), real_c(0), real_c(0)));
   freeSurfaceBoundaryHandling->setInflowInCell(Cell(cell_idx_c(-1), cell_idx_c(2), cell_idx_c(0)),
                                                Vector3< real_t >(real_c(0.03), real_c(0.01), real_c(0)));
   freeSurfaceBoundaryHandling->setInflowInCell(Cell(cell_idx_c(-1), cell_idx_c(3), cell_idx_c(0)),
                                                Vector3< real_t >(real_c(0.04), real_c(0), real_c(0)));
   freeSurfaceBoundaryHandling->setInflowInCell(Cell(cell_idx_c(-1), cell_idx_c(4), cell_idx_c(0)),
                                                Vector3< real_t >(real_c(0.05), real_c(0.02), real_c(0)));

   // these inflow cells should not have any influence as their velocity direction points away from the neighboring gas
   // cells
   freeSurfaceBoundaryHandling->setInflowInCell(Cell(cell_idx_c(-1), cell_idx_c(5), cell_idx_c(0)),
                                                Vector3< real_t >(real_c(-0.06), real_c(0), real_c(0)));
   freeSurfaceBoundaryHandling->setInflowInCell(Cell(cell_idx_c(-1), cell_idx_c(6), cell_idx_c(0)),
                                                Vector3< real_t >(real_c(-0.07), real_c(0), real_c(0)));

   freeSurfaceBoundaryHandling->initFlagsFromFillLevel();

   // add bubble model
   auto bubbleModel = std::make_shared< bubble_model::BubbleModel< Stencil_T > >(blockForest, true);
   bubbleModel->initFromFillLevelField(fillFieldID);
   bubbleModel->setAtmosphere(Cell(cell_idx_c(1), cell_idx_c(1), cell_idx_c(0)), real_c(1));

   // create timeloop
   SweepTimeloop timeloop(blockForest, uint_c(1));

   // add communication
   blockforest::communication::UniformBufferedScheme< typename LatticeModel_T::Stencil > comm(blockForest);
   comm.addPackInfo(std::make_shared< field::communication::PackInfo< PdfField_T > >(pdfFieldID));
   comm.addPackInfo(std::make_shared< field::communication::PackInfo< ScalarField_T > >(fillFieldID));
   comm.addPackInfo(
      std::make_shared< field::communication::PackInfo< FlagField_T > >(freeSurfaceBoundaryHandling->getFlagFieldID()));

   // communicate
   comm();

   // add surface geometry handler
   SurfaceGeometryHandler< LatticeModel_T, FlagField_T, ScalarField_T, VectorField_T > geometryHandler(
      blockForest, freeSurfaceBoundaryHandling, fillFieldID, "FiniteDifferenceMethod", false, false, real_c(0));

   ConstBlockDataID curvatureFieldID = geometryHandler.getConstCurvatureFieldID();
   ConstBlockDataID normalFieldID    = geometryHandler.getConstNormalFieldID();

   geometryHandler.addSweeps(timeloop);

   // add boundary handling for standard boundaries and free surface boundaries
   SurfaceDynamicsHandler< LatticeModel_T, FlagField_T, ScalarField_T, VectorField_T > dynamicsHandler(
      blockForest, pdfFieldID, flagFieldID, fillFieldID, forceDensityFieldID, normalFieldID, curvatureFieldID,
      freeSurfaceBoundaryHandling, bubbleModel, "NormalBasedKeepCenter", "EquilibriumRefilling", "EvenlyNewInterface",
      relaxRate, Vector3< real_t >(real_c(0)), real_c(0), false, real_c(1e-3), real_c(1e-1));

   dynamicsHandler.addSweeps(timeloop);

   timeloop.singleStep();

   // evaluate if inflow boundary has generated the correct interface cells
   for (auto blockIt = blockForest->begin(); blockIt != blockForest->end(); ++blockIt)
   {
      const ScalarField_T* const fillField = blockIt->getData< const ScalarField_T >(fillFieldID);
      const typename lbm::Adaptor< LatticeModel_T >::Density* const densityField =
         blockIt->getData< const typename lbm::Adaptor< LatticeModel_T >::Density >(densityAdaptor);
      const typename lbm::Adaptor< LatticeModel_T >::VelocityVector* const velocityField =
         blockIt->getData< const typename lbm::Adaptor< LatticeModel_T >::VelocityVector >(velocityAdaptor);
      const FlagField_T* const flagField =
         blockIt->getData< const FlagField_T >(freeSurfaceBoundaryHandling->getFlagFieldID());

      WALBERLA_FOR_ALL_CELLS(fillFieldIt, fillField, velocityFieldIt, velocityField, densityFieldIt, densityField,
                             flagFieldIt, flagField, {
                                if (flagFieldIt.cell() == Cell(cell_idx_c(0), cell_idx_c(0), cell_idx_c(0)))
                                {
                                   WALBERLA_CHECK(flagInfo.isInterface(flagFieldIt));

                                   // velocity must be the average from inflow cells (-1,0,0) and (1,-1,0)
                                   WALBERLA_CHECK_FLOAT_EQUAL((*velocityFieldIt)[0], real_c(0.005), real_c(1e-15));
                                   WALBERLA_CHECK_FLOAT_EQUAL((*velocityFieldIt)[1], real_c(0.005), real_c(1e-15));
                                   WALBERLA_CHECK_FLOAT_EQUAL((*velocityFieldIt)[2], real_c(0), real_c(1e-15));
                                   WALBERLA_CHECK_FLOAT_EQUAL(*densityFieldIt, real_c(1), real_c(1e-15));
                                   continue;
                                }

                                if (flagFieldIt.cell() == Cell(cell_idx_c(0), cell_idx_c(1), cell_idx_c(0)))
                                {
                                   WALBERLA_CHECK(flagInfo.isInterface(flagFieldIt));

                                   // velocity must be identical to the inflow cell (-1,1,0); no other inflow cell
                                   // should influence this cell
                                   WALBERLA_CHECK_FLOAT_EQUAL((*velocityFieldIt)[0], real_c(0.02), real_c(1e-15));
                                   WALBERLA_CHECK_FLOAT_EQUAL((*velocityFieldIt)[1], real_c(0), real_c(1e-15));
                                   WALBERLA_CHECK_FLOAT_EQUAL((*velocityFieldIt)[2], real_c(0), real_c(1e-15));
                                   WALBERLA_CHECK_FLOAT_EQUAL(*densityFieldIt, real_c(1), real_c(1e-15));
                                   continue;
                                }

                                if (flagFieldIt.cell() == Cell(cell_idx_c(0), cell_idx_c(2), cell_idx_c(0)))
                                {
                                   WALBERLA_CHECK(flagInfo.isInterface(flagFieldIt));

                                   // velocity must be identical to the inflow cell (-1,2,0); no other inflow cell
                                   // should influence this cell
                                   WALBERLA_CHECK_FLOAT_EQUAL((*velocityFieldIt)[0], real_c(0.03), real_c(1e-15));
                                   WALBERLA_CHECK_FLOAT_EQUAL((*velocityFieldIt)[1], real_c(0.01), real_c(1e-15));
                                   WALBERLA_CHECK_FLOAT_EQUAL((*velocityFieldIt)[2], real_c(0), real_c(1e-15));
                                   WALBERLA_CHECK_FLOAT_EQUAL(*densityFieldIt, real_c(1), real_c(1e-15));
                                   continue;
                                }

                                if (flagFieldIt.cell() == Cell(cell_idx_c(0), cell_idx_c(3), cell_idx_c(0)))
                                {
                                   WALBERLA_CHECK(flagInfo.isInterface(flagFieldIt));

                                   // velocity must be the average from inflow cells (-1,2,0) and (-1,3,0)
                                   WALBERLA_CHECK_FLOAT_EQUAL((*velocityFieldIt)[0], real_c(0.035), real_c(1e-15));
                                   WALBERLA_CHECK_FLOAT_EQUAL((*velocityFieldIt)[1], real_c(0.005), real_c(1e-15));
                                   WALBERLA_CHECK_FLOAT_EQUAL((*velocityFieldIt)[2], real_c(0), real_c(1e-15));
                                   WALBERLA_CHECK_FLOAT_EQUAL(*densityFieldIt, real_c(1), real_c(1e-15));
                                   continue;
                                }

                                if (flagFieldIt.cell() == Cell(cell_idx_c(0), cell_idx_c(4), cell_idx_c(0)))
                                {
                                   WALBERLA_CHECK(flagInfo.isInterface(flagFieldIt));

                                   // velocity must be identical to inflow cell (-1,4,0); no other inflow cell
                                   // should influence this cell
                                   WALBERLA_CHECK_FLOAT_EQUAL((*velocityFieldIt)[0], real_c(0.05), real_c(1e-15));
                                   WALBERLA_CHECK_FLOAT_EQUAL((*velocityFieldIt)[1], real_c(0.02), real_c(1e-15));
                                   WALBERLA_CHECK_FLOAT_EQUAL((*velocityFieldIt)[2], real_c(0), real_c(1e-15));
                                   WALBERLA_CHECK_FLOAT_EQUAL(*densityFieldIt, real_c(1), real_c(1e-15));
                                   WALBERLA_LOG_DEVEL_VAR(*velocityFieldIt);
                                   WALBERLA_LOG_DEVEL_VAR(flagInfo.isInterface(flagFieldIt));
                                   continue;
                                }

                                if (flagFieldIt.cell() == Cell(cell_idx_c(0), cell_idx_c(5), cell_idx_c(0)))
                                {
                                   // cell must be converted due to velocity vector pointing from inflow
                                   // cell(-1,4,0) to this cell
                                   WALBERLA_CHECK(flagInfo.isInterface(flagFieldIt));

                                   // velocity must be identical to inflow cell (-1,4,0); no other inflow cell
                                   // should influence this cell
                                   WALBERLA_CHECK_FLOAT_EQUAL((*velocityFieldIt)[0], real_c(0.05), real_c(1e-15));
                                   WALBERLA_CHECK_FLOAT_EQUAL((*velocityFieldIt)[1], real_c(0.02), real_c(1e-15));
                                   WALBERLA_CHECK_FLOAT_EQUAL((*velocityFieldIt)[2], real_c(0), real_c(1e-15));
                                   WALBERLA_CHECK_FLOAT_EQUAL(*densityFieldIt, real_c(1), real_c(1e-15));
                                   continue;
                                }

                                // cell(0,6,0)
                                WALBERLA_CHECK(flagInfo.isGas(flagFieldIt));
                             }) // WALBERLA_FOR_ALL_CELLS
   }

   MPIManager::instance()->resetMPI();
}

int main(int argc, char** argv)
{
   debug::enterTestMode();
   Environment env(argc, argv);

   testInflow();

   return EXIT_SUCCESS;
}

} // namespace InflowTest
} // namespace free_surface
} // namespace walberla

int main(int argc, char** argv) { return walberla::free_surface::InflowTest::main(argc, argv); }