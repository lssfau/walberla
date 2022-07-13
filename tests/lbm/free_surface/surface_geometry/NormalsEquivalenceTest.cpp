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
//! \file NormalsEquivalenceTest.cpp
//! \ingroup lbm/free_surface/surface_geometry
//! \author Matthias Markl <matthias.markl@fau.de>
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//! \brief Test if the NormalSweep (unrolled version) is equal to a loop-based computation of the interface normal.
//!
//! Test (with respect to a maximum error of 1e-13) if the explicit normal computation (implemented in NormalSweep)
//! gives the same result as the loop-based computation of the interface normal. The former method is faster and does
//! not loop over all D3Q27 directions but calculates the interface normal explicitly. See function computeNormal() in
//! src/lbm/free_surface/surface_geometry/NormalSweep.impl.h
//
//======================================================================================================================

#include "blockforest/Initialization.h"

#include "core/Environment.h"
#include "core/debug/TestSubsystem.h"

#include "field/AddToStorage.h"

#include "lbm/blockforest/communication/SimpleCommunication.h"
#include "lbm/free_surface/FlagInfo.h"
#include "lbm/free_surface/surface_geometry/NormalSweep.h"

#include "stencil/D3Q27.h"

#include "timeloop/SweepTimeloop.h"

#include <cstdlib>

namespace walberla
{
namespace free_surface
{
namespace NormalsEquivalenceTest
{
// define types
using Stencil_T = stencil::D3Q27;
using flag_t    = uint32_t;

// define fields
using ScalarField_T = GhostLayerField< real_t, 1 >;
using VectorField_T = GhostLayerField< Vector3< real_t >, 1 >;
using FlagField_T   = FlagField< flag_t >;

int main(int argc, char** argv)
{
   debug::enterTestMode();
   Environment env(argc, argv);

   // define the domain size
   const Vector3< uint_t > numBlocks(uint_c(1));
   const Vector3< uint_t > domainSize(uint_c(32), uint_c(32), uint_c(32));
   const Vector3< uint_t > cellsPerBlock(domainSize[0] / numBlocks[0], domainSize[1] / numBlocks[1],
                                         domainSize[2] / numBlocks[2]);

   // create block grid
   const std::shared_ptr< StructuredBlockForest > blockForest =
      blockforest::createUniformBlockGrid(numBlocks[0], numBlocks[1], numBlocks[2],             // blocks
                                          cellsPerBlock[0], cellsPerBlock[1], cellsPerBlock[2], // cells
                                          real_c(1.0),                                          // dx
                                          true,                                                 // one block per process
                                          true, true, true);                                    // periodicity

   // add fields
   BlockDataID flagFieldID = field::addFlagFieldToStorage< FlagField_T >(blockForest, "Flags", uint_c(2));
   BlockDataID fillFieldID =
      field::addToStorage< ScalarField_T >(blockForest, "Fill levels", real_c(0), field::fzyx, uint_c(1));
   BlockDataID normalFieldLoopID = field::addToStorage< VectorField_T >(
      blockForest, "Normals loop", Vector3< real_t >(real_c(0)), field::fzyx, uint_c(1));
   BlockDataID normalFieldExplID = field::addToStorage< VectorField_T >(
      blockForest, "Normals explicit", Vector3< real_t >(real_c(0)), field::fzyx, uint_c(1));

   // add boundary handling
   const auto flagInfo = FlagInfo< FlagField_T >(Set< FlagUID >(), Set< FlagUID >());
   for (auto blockIt = blockForest->begin(); blockIt != blockForest->end(); ++blockIt)
   {
      flagInfo.registerFlags(blockIt->getData< FlagField_T >(flagFieldID));
   }
   WALBERLA_ASSERT(flagInfo.isConsistentAcrossBlocksAndProcesses(blockForest, flagFieldID));

   // add communication
   blockforest::SimpleCommunication< stencil::D3Q27 > commFill(blockForest, fillFieldID);
   blockforest::SimpleCommunication< stencil::D3Q27 > commFlag(blockForest, flagFieldID);

   // initialize flag field (only interface cells) and fill levels (random values)
   std::srand(static_cast< unsigned int >(std::time(nullptr)));
   for (auto blockIt = blockForest->begin(); blockIt != blockForest->end(); ++blockIt)
   {
      ScalarField_T* const fillField = blockIt->getData< ScalarField_T >(fillFieldID);
      FlagField_T* const flagField   = blockIt->getData< FlagField_T >(flagFieldID);

      // set whole flag field to interface
      flagField->set(flagInfo.interfaceFlag);

      // set random fill levels in fill field
      WALBERLA_FOR_ALL_CELLS(fillFieldIt, fillField, {
         *fillFieldIt = real_c(std::rand()) / real_c(std::numeric_limits< int >::max());
      }) // WALBERLA_FOR_ALL_CELLS
   }

   // communicate fill level field, and flag field
   commFill();
   commFlag();

   // loop-based normal computation, i.e., old and slower version of NormalSweep operator()
   auto normalsFunc = [&](IBlock* block) {
      // get fields
      VectorField_T* const normalField     = block->getData< VectorField_T >(normalFieldLoopID);
      const ScalarField_T* const fillField = block->getData< const ScalarField_T >(fillFieldID);
      const FlagField_T* const flagField   = block->getData< const FlagField_T >(flagFieldID);

      WALBERLA_FOR_ALL_CELLS(normalFieldIt, normalField, fillFieldIt, fillField, flagFieldIt, flagField, {
         bool computeNormalInCell = flagInfo.isInterface(flagFieldIt);

         Vector3< real_t >& normal = *normalFieldIt;

         if (computeNormalInCell)
         {
            if (!isFlagInNeighborhood< stencil::D3Q27 >(flagFieldIt, flagInfo.obstacleFlagMask))
            {
               normal.set(real_c(0), real_c(0), real_c(0));

               // loop over all directions to compute interface normal; compare this loop with the explicit (non-loop)
               // computation in function computeNormal() in src/lbm/free_surface/surface_geometry/NormalSweep.impl.h
               for (auto d = stencil::D3Q27::beginNoCenter(); d != stencil::D3Q27::end(); ++d)
               {
                  const real_t weightedFill =
                     real_c(stencil::gaussianMultipliers[d.toIdx()]) * fillFieldIt.neighbor(*d);
                  normal[0] += real_c(d.cx()) * weightedFill;
                  normal[1] += real_c(d.cy()) * weightedFill;
                  normal[2] += real_c(d.cz()) * weightedFill;
               }
            }

            // normalize and negate normal (to make it point from gas to liquid)
            normal = real_c(-1) * normal.getNormalizedOrZero();
         }
         else { normal.set(real_c(0), real_c(0), real_c(0)); }
      }) // WALBERLA_FOR_ALL_CELLS
   };

   real_t err_max  = real_c(0);
   real_t err_mean = real_c(0);

#ifdef WALBERLA_DOUBLE_ACCURACY
   real_t tolerance = real_c(1e-13);
#else
   real_t tolerance = real_c(1e-4);
#endif

   auto evaluateError = [&](IBlock* block) {
      const VectorField_T* const loopField = block->getData< const VectorField_T >(normalFieldLoopID);
      const VectorField_T* const explField = block->getData< const VectorField_T >(normalFieldExplID);

      err_mean = real_c(0);
      err_max  = real_c(0);

      WALBERLA_FOR_ALL_CELLS(loopFieldIt, loopField, explFieldIt, explField, {
         const Vector3< real_t > diff = *loopFieldIt - *explFieldIt;
         const real_t length          = diff.length();

         if (length > tolerance)
         {
            WALBERLA_ABORT("Unequal normals at " << loopFieldIt.cell() << " diff=" << diff << ", |diff|=" << length);
         }
         err_mean += length;
         err_max = std::max(err_max, length);
      }) // WALBERLA_FOR_ALL_CELLS
   };

   // create timeloop
   SweepTimeloop timeloop(blockForest, uint_c(1));

   // add regular sweep for computing interface normals
   NormalSweep< Stencil_T, FlagField_T, ScalarField_T, VectorField_T > normalsSweep(
      normalFieldExplID, fillFieldID, flagFieldID, flagIDs::interfaceFlagID, flagIDs::liquidInterfaceGasFlagIDs,
      flagIDs::gasFlagID, false, false, false, false);
   timeloop.add() << Sweep(normalsSweep, "Normal sweep (explicit)");

   // add sweeps
   timeloop.add() << Sweep(normalsFunc, "Normal loop-based");
   timeloop.add() << Sweep(evaluateError, "Error computation");

   WcTimingPool timeloopTiming;
   timeloop.run(timeloopTiming);
   // timeloopTiming.logResultOnRoot();

   WALBERLA_LOG_RESULT("Max Error: " << err_max);
   WALBERLA_LOG_RESULT("Mean Error: " << err_mean / real_c(domainSize[0] * domainSize[1] * domainSize[2]));

   return EXIT_SUCCESS;
}
} // namespace NormalsEquivalenceTest
} // namespace free_surface
} // namespace walberla

int main(int argc, char** argv) { return walberla::free_surface::NormalsEquivalenceTest::main(argc, argv); }
