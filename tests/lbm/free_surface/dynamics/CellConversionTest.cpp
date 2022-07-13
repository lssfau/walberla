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
//! \file CellConversionTest.cpp
//! \ingroup lbm/free_surface/dynamics
//! \author Martin Bauer
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//! \brief Test cell conversions on a flat surface by checking mass conservation and density balance.
//!
//! Initialize a pool of liquid in half of the domain and a box-shaped atmosphere (density 1.1) bubble in the remaining
//! cells. All sweeps from SurfaceDynamicsHandler are performed. Due to the higher density in the gas, the interface
//! cells that separate liquid and gas are emptied and their fill level must become negative. These cells are converted
//! to liquid, while former fluid cells must become interface cells. After this process, mass must be conserved and the
//! fluid density must balance the atmosphere density.
//
//======================================================================================================================

#include "blockforest/Initialization.h"
#include "blockforest/communication/UniformBufferedScheme.h"

#include "core/Environment.h"
#include "core/debug/TestSubsystem.h"

#include "field/AddToStorage.h"
#include "field/communication/PackInfo.h"

#include "lbm/field/Adaptors.h"
#include "lbm/field/AddToStorage.h"
#include "lbm/free_surface/FlagInfo.h"
#include "lbm/free_surface/boundary/FreeSurfaceBoundaryHandling.h"
#include "lbm/free_surface/dynamics/SurfaceDynamicsHandler.h"
#include "lbm/lattice_model/CollisionModel.h"
#include "lbm/lattice_model/D3Q19.h"
#include "lbm/lattice_model/ForceModel.h"

#include "timeloop/SweepTimeloop.h"

#include <limits>

namespace walberla
{
namespace free_surface
{
namespace CellConversionTest
{
// define types
using flag_t = uint32_t;

// define fields
using ScalarField_T = GhostLayerField< real_t, 1 >;
using VectorField_T = GhostLayerField< Vector3< real_t >, 1 >;
using FlagField_T   = FlagField< flag_t >;

// compute the total mass of the fluid (in all interface and liquid cells)
template< typename LatticeModel_T, typename FlagField_T >
real_t computeTotalMass(
   const std::weak_ptr< StructuredBlockForest >& blockForest, ConstBlockDataID pdfField, ConstBlockDataID fillField,
   ConstBlockDataID flagField,
   const typename FreeSurfaceBoundaryHandling< LatticeModel_T, FlagField_T, ScalarField_T >::FlagInfo_T& flags);

// compute the minimum and maximum density in all interface and liquid cells
template< typename LatticeModel_T, typename FlagField_T >
void computeMinMaxDensity(
   const std::weak_ptr< StructuredBlockForest >& blockForest, ConstBlockDataID pdfField, ConstBlockDataID flagField,
   const typename FreeSurfaceBoundaryHandling< LatticeModel_T, FlagField_T, ScalarField_T >::FlagInfo_T& flags,
   real_t& minDensity, real_t& maxDensity);

template< typename LatticeModel_T >
void testCellConversion()
{
   // define types
   using Stencil_T  = typename LatticeModel_T::Stencil;
   using PdfField_T = lbm::PdfField< LatticeModel_T >;

   using FreeSurfaceBoundaryHandling_T = FreeSurfaceBoundaryHandling< LatticeModel_T, FlagField_T, ScalarField_T >;

   // define the domain size
   const Vector3< uint_t > numBlocks(uint_c(1));
   const Vector3< uint_t > domainSize(uint_c(3), uint_c(10), uint_c(2));
   const Vector3< uint_t > cellsPerBlock(domainSize[0] / numBlocks[0], domainSize[1] / numBlocks[1],
                                         domainSize[2] / numBlocks[2]);

   // create block grid
   const std::shared_ptr< StructuredBlockForest > blockForest =
      blockforest::createUniformBlockGrid(numBlocks[0], numBlocks[1], numBlocks[2],             // blocks
                                          cellsPerBlock[0], cellsPerBlock[1], cellsPerBlock[2], // cells
                                          real_c(1.0),                                          // dx
                                          true,                                                 // one block per process
                                          true, false, true);                                   // periodicity

   real_t relaxRate = real_c(0.51);

   // create lattice model with omega=0.51
   LatticeModel_T latticeModel(relaxRate);

   // add fields
   BlockDataID pdfFieldID = lbm::addPdfFieldToStorage(blockForest, "PDF field", latticeModel, field::fzyx);
   BlockDataID fillFieldID =
      field::addToStorage< ScalarField_T >(blockForest, "Fill levels", real_c(1.0), field::fzyx, uint_c(1));
   BlockDataID forceFieldID = field::addToStorage< VectorField_T >(
      blockForest, "Force field", Vector3< real_t >(real_c(0)), field::fzyx, uint_c(1));
   BlockDataID curvatureFieldID =
      field::addToStorage< ScalarField_T >(blockForest, "Curvature", real_c(0.0), field::fzyx, uint_c(1));
   BlockDataID normalFieldID = field::addToStorage< VectorField_T >(
      blockForest, "Normals", Vector3< real_t >(real_c(0), real_c(1), real_c(0)), field::fzyx, uint_c(1));

   // add boundary handling
   const std::shared_ptr< FreeSurfaceBoundaryHandling_T > freeSurfaceBoundaryHandling =
      std::make_shared< FreeSurfaceBoundaryHandling_T >(blockForest, pdfFieldID, fillFieldID);
   BlockDataID flagFieldID                                            = freeSurfaceBoundaryHandling->getFlagFieldID();
   const typename FreeSurfaceBoundaryHandling_T::FlagInfo_T& flagInfo = freeSurfaceBoundaryHandling->getFlagInfo();

   // add box-shaped gas bubble (occupies about half of the domain in y-direction)
   AABB box    = blockForest->getDomain();
   auto newMin = box.min() + Vector3< real_t >(real_c(0), real_c(0.5) * box.ySize() + real_c(1 - 0.02), real_c(0));
   box.initMinMaxCorner(newMin, box.max());
   freeSurfaceBoundaryHandling->addFreeSurfaceObject(box);

   // set no slip boundary conditions at the southern and northern domain borders
   freeSurfaceBoundaryHandling->setNoSlipAtBorder(stencil::S);
   freeSurfaceBoundaryHandling->setNoSlipAtBorder(stencil::N);
   freeSurfaceBoundaryHandling->initFlagsFromFillLevel();

   // add bubble model
   auto bubbleModel = std::make_shared< bubble_model::BubbleModel< Stencil_T > >(blockForest, false);
   bubbleModel->initFromFillLevelField(fillFieldID);

   // get some cell located inside the bubble (used as representative cell to set bubble to atmosphere bubble type)
   Cell cellInBubble;
   cellInBubble.x() = cell_idx_c(box.xMin() + real_c(0.5) * box.xSize());
   cellInBubble.y() = cell_idx_c(box.yMin() + real_c(0.5) * box.ySize());
   cellInBubble.z() = cell_idx_c(box.zMin() + real_c(0.5) * box.zSize());

   // set bubble to atmosphere with constant density higher than the initial fluid density
   const real_t atmDensity = real_c(1.1);
   bubbleModel->setAtmosphere(cellInBubble, atmDensity);

   // create timeloop; 400 time steps are required (for very low omega of 0.51) to ensure that fluid density has
   // stabilized
   SweepTimeloop timeloop(blockForest, uint_c(400));

   // add communication
   blockforest::communication::UniformBufferedScheme< Stencil_T > comm(blockForest);
   comm.addPackInfo(std::make_shared< field::communication::PackInfo< PdfField_T > >(pdfFieldID));
   comm.addPackInfo(std::make_shared< field::communication::PackInfo< ScalarField_T > >(fillFieldID));
   comm.addPackInfo(std::make_shared< field::communication::PackInfo< FlagField_T > >(flagFieldID));

   // communicate
   comm();

   // add various sweeps for surface dynamics
   SurfaceDynamicsHandler< LatticeModel_T, FlagField_T, ScalarField_T, VectorField_T > dynamicsHandler(
      blockForest, pdfFieldID, flagFieldID, fillFieldID, forceFieldID, normalFieldID, curvatureFieldID,
      freeSurfaceBoundaryHandling, bubbleModel, "NormalBasedKeepCenter", "EquilibriumRefilling", "EvenlyNewInterface",
      relaxRate, Vector3< real_t >(real_c(0)), real_c(0), false, false, real_c(1e-3), real_c(1e-1));
   dynamicsHandler.addSweeps(timeloop);

   real_t initialMass =
      computeTotalMass< LatticeModel_T, FlagField_T >(blockForest, pdfFieldID, fillFieldID, flagFieldID, flagInfo);

   timeloop.run();

   // in the absence of forces, fluid density must balance atmosphere density
   real_t rhoMin = real_c(0);
   real_t rhoMax = real_c(0);
   computeMinMaxDensity< LatticeModel_T, FlagField_T >(blockForest, pdfFieldID, flagFieldID, flagInfo, rhoMin, rhoMax);
   WALBERLA_CHECK_FLOAT_EQUAL(atmDensity, rhoMin);
   WALBERLA_CHECK_FLOAT_EQUAL(atmDensity, rhoMax);

   // mass must be conserved
   real_t finalMass =
      computeTotalMass< LatticeModel_T, FlagField_T >(blockForest, pdfFieldID, fillFieldID, flagFieldID, flagInfo);
   WALBERLA_CHECK_FLOAT_EQUAL(initialMass, finalMass);

   MPIManager::instance()->resetMPI();
}

// compute the total mass of the fluid (in all interface and liquid cells)
template< typename LatticeModel_T, typename FlagField_T >
real_t computeTotalMass(
   const std::weak_ptr< StructuredBlockForest >& blockForestPtr, ConstBlockDataID pdfFieldID,
   ConstBlockDataID fillFieldID, ConstBlockDataID flagFieldID,
   const typename FreeSurfaceBoundaryHandling< LatticeModel_T, FlagField_T, ScalarField_T >::FlagInfo_T& flagInfo)
{
   const auto blockForest = blockForestPtr.lock();
   WALBERLA_CHECK_NOT_NULLPTR(blockForest);

   real_t mass = real_c(0);

   for (auto blockIt = blockForest->begin(); blockIt != blockForest->end(); ++blockIt)
   {
      const lbm::PdfField< LatticeModel_T >* const pdfField =
         blockIt->getData< const lbm::PdfField< LatticeModel_T > >(pdfFieldID);
      const ScalarField_T* const fillField = blockIt->getData< const ScalarField_T >(fillFieldID);
      const FlagField_T* const flagField   = blockIt->getData< const FlagField_T >(flagFieldID);

      // iterate over all interface and liquid cells and compute the total mass of fluid in the domain
      WALBERLA_FOR_ALL_CELLS_OMP(fillFieldIt, fillField, pdfFieldIt, pdfField, flagFieldIt, flagField,
                                 omp parallel for schedule(static) reduction(+:mass),
                                 {
         const real_t rho = lbm::getDensity< LatticeModel_T >(pdfField->latticeModel(), pdfFieldIt);
         if (flagInfo.isInterface(flagFieldIt)) { mass += rho * (*fillFieldIt); }
         else
         {
            if (flagInfo.isLiquid(flagFieldIt)) { mass += rho; }
         }
                                 }); // WALBERLA_FOR_ALL_CELLS_OMP
   }

   WALBERLA_LOG_RESULT("Total current mass " << mass << ".");

   return mass;
}

// compute the minimum and maximum density in all interface and liquid cells
template< typename LatticeModel_T, typename FlagField_T >
void computeMinMaxDensity(
   const std::weak_ptr< StructuredBlockForest >& blockForestPtr, ConstBlockDataID pdfFieldID,
   ConstBlockDataID flagFieldID,
   const typename FreeSurfaceBoundaryHandling< LatticeModel_T, FlagField_T, ScalarField_T >::FlagInfo_T& flags,
   real_t& minDensity, real_t& maxDensity)
{
   const auto blockForest = blockForestPtr.lock();
   WALBERLA_CHECK_NOT_NULLPTR(blockForest);

   minDensity = std::numeric_limits< real_t >::max();
   maxDensity = std::numeric_limits< real_t >::min();

   for (auto blockIt = blockForest->begin(); blockIt != blockForest->end(); ++blockIt)
   {
      const lbm::PdfField< LatticeModel_T >* const pdfField =
         blockIt->getData< const lbm::PdfField< LatticeModel_T > >(pdfFieldID);
      const FlagField_T* const flagField = blockIt->getData< const FlagField_T >(flagFieldID);

      // iterate over all interface and liquid cells and find the minimum and maximum density; explicitly avoid OpenMP
      // (problematic to reduce max and min)
      WALBERLA_FOR_ALL_CELLS_OMP(flagFieldIt, flagField, pdfFieldIt, pdfField, omp critical, {
         if (flags.isInterface(flagFieldIt) || flags.isLiquid(flagFieldIt))
         {
            const real_t rho = lbm::getDensity< LatticeModel_T >(pdfField->latticeModel(), pdfFieldIt);
            minDensity       = std::min(rho, minDensity);
            maxDensity       = std::max(rho, maxDensity);
         }
      }) // WALBERLA_FOR_ALL_CELLS
   }
}

int main(int argc, char** argv)
{
   debug::enterTestMode();
   Environment env(argc, argv);

   WALBERLA_LOG_INFO("Testing with D2Q9 stencil.")
   testCellConversion< walberla::lbm::D2Q9< walberla::lbm::collision_model::SRT, true > >();

   WALBERLA_LOG_INFO("Testing with D3Q19 stencil.")
   testCellConversion< walberla::lbm::D3Q19< walberla::lbm::collision_model::SRT, true > >();

   WALBERLA_LOG_INFO("Testing with D3Q27 stencil.")
   testCellConversion< walberla::lbm::D3Q27< walberla::lbm::collision_model::SRT, true > >();

   return EXIT_SUCCESS;
}
} // namespace CellConversionTest
} // namespace free_surface
} // namespace walberla

int main(int argc, char** argv) { return walberla::free_surface::CellConversionTest::main(argc, argv); }