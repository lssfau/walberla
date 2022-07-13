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
//! \file AdvectionTest.cpp
//! \ingroup lbm/free_surface/dynamics
//! \author Martin Bauer
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//! \brief Test mass advection on a flat surface.
//!
//! Initialize a pool of liquid in half of the domain and a box-shaped atmosphere (density 1.01) bubble in the remaining
//! cells. Only the StreamReconstructAdvectSweep is performed. Due to the higher density in the gas, the interface cells
//! that separate liquid and gas are emptied and their fill level must become negative. These cells must be marked
//! for conversion and the fluid density must balance the atmosphere density.
//
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
#include "lbm/free_surface/dynamics/PdfReconstructionModel.h"
#include "lbm/free_surface/dynamics/StreamReconstructAdvectSweep.h"
#include "lbm/lattice_model/CollisionModel.h"
#include "lbm/lattice_model/D3Q19.h"
#include "lbm/lattice_model/ForceModel.h"
#include "lbm/sweeps/CellwiseSweep.h"
#include "lbm/sweeps/SweepWrappers.h"

#include "timeloop/SweepTimeloop.h"

namespace walberla
{
namespace free_surface
{
namespace AdvectionTest
{
// define types
using Flag_T        = uint32_t;
using ScalarField_T = GhostLayerField< real_t, 1 >;
using VectorField_T = GhostLayerField< Vector3< real_t >, 1 >;
using FlagField_T   = FlagField< Flag_T >;

template< typename LatticeModel_T >
void testAdvection()
{
   // define types
   using Stencil_T                     = typename LatticeModel_T::Stencil;
   using PdfField_T                    = lbm::PdfField< LatticeModel_T >;
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

   // create lattice model with omega=0.51
   LatticeModel_T latticeModel(real_c(0.51));

   // add fields
   BlockDataID pdfFieldID = lbm::addPdfFieldToStorage(blockForest, "PDF field", latticeModel, field::fzyx);
   BlockDataID fillFieldID =
      field::addToStorage< ScalarField_T >(blockForest, "Fill levels", real_c(1.0), field::fzyx, uint_c(1));
   BlockDataID curvatureFieldID =
      field::addToStorage< ScalarField_T >(blockForest, "Curvature", real_c(0.0), field::fzyx, uint_c(1));
   BlockDataID normalFieldID = field::addToStorage< VectorField_T >(
      blockForest, "Normals", Vector3< real_t >(real_c(0), real_c(1), real_c(0)), field::fzyx, uint_c(1));

   BlockDataID densityAdaptor = field::addFieldAdaptor< typename lbm::Adaptor< LatticeModel_T >::Density >(
      blockForest, pdfFieldID, "DensityAdaptor");

   // add boundary handling
   const std::shared_ptr< FreeSurfaceBoundaryHandling_T > freeSurfaceBoundaryHandling =
      std::make_shared< FreeSurfaceBoundaryHandling_T >(blockForest, pdfFieldID, fillFieldID);
   const typename FreeSurfaceBoundaryHandling_T::FlagInfo_T& flagInfo = freeSurfaceBoundaryHandling->getFlagInfo();

   // add box-shaped gas bubble (occupies about half of the domain in y-direction)
   AABB box    = blockForest->getDomain();
   auto newMin = box.min() + Vector3< real_t >(real_c(0), real_c(0.5) * box.ySize() + real_c(0.01), real_c(0));
   box.initMinMaxCorner(newMin, box.max());
   freeSurfaceBoundaryHandling->addFreeSurfaceObject(box);

   // set no slip boundary conditions at the southern and northern domain borders
   freeSurfaceBoundaryHandling->setNoSlipAtBorder(stencil::S);
   freeSurfaceBoundaryHandling->setNoSlipAtBorder(stencil::N);
   freeSurfaceBoundaryHandling->initFlagsFromFillLevel();

   // add bubble model
   bubble_model::BubbleModel< Stencil_T > bubbleModel(blockForest, false);
   bubbleModel.initFromFillLevelField(fillFieldID);

   // get some cell located inside the bubble (used as representative cell to set bubble to atmosphere bubble type)
   Cell cellInBubble;
   cellInBubble.x() = cell_idx_c(box.xMin() + real_c(0.5) * box.xSize());
   cellInBubble.y() = cell_idx_c(box.yMin() + real_c(0.5) * box.ySize());
   cellInBubble.z() = cell_idx_c(box.zMin() + real_c(0.5) * box.zSize());

   // set bubble to atmosphere with constant density higher than the initial fluid density
   const real_t atmDensity = real_c(1.01);
   bubbleModel.setAtmosphere(cellInBubble, atmDensity);

   // create timeloop; 300 time steps are required (for very low omega of 0.51) to ensure that fluid density has
   // stabilized
   SweepTimeloop timeloop(blockForest, uint_c(300));

   // add communication
   blockforest::communication::UniformBufferedScheme< typename LatticeModel_T::Stencil > comm(blockForest);
   comm.addPackInfo(std::make_shared< field::communication::PackInfo< PdfField_T > >(pdfFieldID));
   comm.addPackInfo(std::make_shared< field::communication::PackInfo< ScalarField_T > >(fillFieldID));
   comm.addPackInfo(
      std::make_shared< field::communication::PackInfo< FlagField_T > >(freeSurfaceBoundaryHandling->getFlagFieldID()));

   // communicate
   comm();

   const PdfReconstructionModel pdfRecModel = PdfReconstructionModel("NormalBasedKeepCenter");

   // add free surface boundary sweep for
   // - reconstruction of PDFs in interface cells
   // - advection of mass
   // - marking interface cells for conversion
   // - update bubble volumes
   StreamReconstructAdvectSweep< LatticeModel_T, typename FreeSurfaceBoundaryHandling_T::BoundaryHandling_T,
                                 FlagField_T, typename FreeSurfaceBoundaryHandling_T::FlagInfo_T, ScalarField_T,
                                 VectorField_T, false >
      streamReconstructAdvectSweep(real_c(0), freeSurfaceBoundaryHandling->getHandlingID(), fillFieldID,
                                   freeSurfaceBoundaryHandling->getFlagFieldID(), pdfFieldID, normalFieldID,
                                   curvatureFieldID, flagInfo, &bubbleModel, pdfRecModel, false, real_c(1e-3),
                                   real_c(1e-1));
   timeloop.add() << Sweep(streamReconstructAdvectSweep);

   // add boundary handling sweep
   timeloop.add() << BeforeFunction(comm) << Sweep(freeSurfaceBoundaryHandling->getBoundarySweep());

   // add LBM collision sweep
   auto lbmSweep = lbm::makeCellwiseSweep< LatticeModel_T, FlagField_T >(
      pdfFieldID, freeSurfaceBoundaryHandling->getFlagFieldID(), flagIDs::liquidInterfaceFlagIDs);
   timeloop.add() << Sweep(lbm::makeCollideSweep(lbmSweep));

   timeloop.run();

   // evaluate
   for (auto blockIt = blockForest->begin(); blockIt != blockForest->end(); ++blockIt)
   {
      const ScalarField_T* const fillField = blockIt->getData< const ScalarField_T >(fillFieldID);
      const typename lbm::Adaptor< LatticeModel_T >::Density* const densityField =
         blockIt->getData< const typename lbm::Adaptor< LatticeModel_T >::Density >(densityAdaptor);
      const FlagField_T* const flagField =
         blockIt->getData< const FlagField_T >(freeSurfaceBoundaryHandling->getFlagFieldID());

      WALBERLA_FOR_ALL_CELLS(fillFieldIt, fillField, densityFieldIt, densityField, flagFieldIt, flagField, {
         if (flagInfo.isInterface(flagFieldIt))
         {
            // fill level in interface cells must be negative
            WALBERLA_CHECK_LESS(*fillFieldIt, real_c(0));

            // due to negative fill level, these cells must be marked for conversion to gas
            WALBERLA_CHECK(isFlagSet(flagFieldIt, flagInfo.convertToGasFlag));
         }

         // in the absence of forces, fluid density must balance atmosphere density
         if (flagInfo.isLiquid(flagFieldIt)) { WALBERLA_CHECK_FLOAT_EQUAL(*densityFieldIt, atmDensity); }
      }) // WALBERLA_FOR_ALL_CELLS
   }

   MPIManager::instance()->resetMPI();
}

int main(int argc, char** argv)
{
   debug::enterTestMode();
   Environment env(argc, argv);

   WALBERLA_LOG_INFO("Testing with D2Q9 stencil.");
   testAdvection< lbm::D2Q9< lbm::collision_model::SRT, true, lbm::force_model::None, 2 > >();

   WALBERLA_LOG_INFO("Testing with D3Q19 stencil.");
   testAdvection< lbm::D3Q19< lbm::collision_model::SRT, true, lbm::force_model::None, 2 > >();

   WALBERLA_LOG_INFO("Testing with D3Q27 stencil.");
   testAdvection< lbm::D3Q27< lbm::collision_model::SRT, true, lbm::force_model::None, 2 > >();

   return EXIT_SUCCESS;
}

} // namespace AdvectionTest
} // namespace free_surface
} // namespace walberla

int main(int argc, char** argv) { return walberla::free_surface::AdvectionTest::main(argc, argv); }