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
//! \file SplitDetectionTest.cpp
//! \ingroup lbm/free_surface/bubble_model
//! \author Martin Bauer
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//! \brief Test bubble split detection with 2D and 3D bubbles.
//
//======================================================================================================================

#include "core/Environment.h"
#include "core/debug/TestSubsystem.h"

#include "lbm/free_surface/bubble_model/BubbleModel.h"

#include "stencil/D2Q9.h"
#include "stencil/D3Q19.h"
#include "stencil/D3Q27.h"

namespace walberla
{
namespace free_surface
{
namespace SplitDetectionTest
{
using namespace bubble_model;

// class derived from BubbleModel to access its protected members for testing purposes
template< typename Stencil_T >
class BubbleModelTest : public BubbleModel< Stencil_T >
{
 public:
   static bool checkForSplit(BubbleField_T* bf, const Cell& cell, BubbleID prevID)
   {
      return BubbleModel< Stencil_T >::checkForSplit(bf, cell, prevID);
   }
}; // class BubbleModelTest

void initSlice(BubbleField_T* bf, cell_idx_t z, const BubbleID array3x3[])
{
   for (cell_idx_t y = cell_idx_c(0); y < cell_idx_c(3); ++y)
   {
      for (cell_idx_t x = cell_idx_c(0); x < cell_idx_c(3); ++x)
      {
         bf->get(x, y, z) = array3x3[3 * (2 - y) + x];
      }
   }
}

void test2D_notConnected()
{
   // create a 3x3x3 bubble field (without ghost layer) and initialize with invalid IDs
   BubbleField_T bubbleField(uint_c(3), uint_c(3), uint_c(3), uint_c(0), INVALID_BUBBLE_ID);

   // initialize a 2D slice of the bubble field (disconnected bubble)
   const BubbleID N      = INVALID_BUBBLE_ID;
   const BubbleID init[] = { 2, 2, 2, N, N, N, 2, 2, 2 };

   initSlice(&bubbleField, cell_idx_c(1), init);

   WALBERLA_CHECK(BubbleModelTest< stencil::D2Q9 >::checkForSplit(
                     &bubbleField, Cell(cell_idx_c(1), cell_idx_c(1), cell_idx_c(1)), 2) == true);
}

void test2D_connected()
{
   // create a 3x3x3 bubble field (without ghost layer) and initialize with invalid IDs
   BubbleField_T bubbleField(uint_c(3), uint_c(3), uint_c(3), uint_c(0), INVALID_BUBBLE_ID);

   // initialize a 2D slice of the bubble field (connected bubble)
   const BubbleID N      = INVALID_BUBBLE_ID;
   const BubbleID init[] = { 2, 2, 2, 2, N, N, 2, 2, 2 };

   initSlice(&bubbleField, cell_idx_c(1), init);

   WALBERLA_CHECK(BubbleModelTest< stencil::D2Q9 >::checkForSplit(
                     &bubbleField, Cell(cell_idx_c(1), cell_idx_c(1), cell_idx_c(1)), 2) == false);
}

void test3D_connected()
{
   // create a 3x3x3 bubble field (without ghost layer) and initialize with invalid IDs
   BubbleField_T bubbleField(uint_c(3), uint_c(3), uint_c(3), uint_c(0), INVALID_BUBBLE_ID);

   // initialize whole bubble field (connected bubble)
   const BubbleID N    = INVALID_BUBBLE_ID;
   const BubbleID z0[] = { 2, 2, 2, N, N, N, 2, 2, 2 };
   const BubbleID z1[] = { N, 2, N, N, N, N, N, N, N };
   const BubbleID z2[] = { N, 2, N, N, 2, N, N, 2, N };

   initSlice(&bubbleField, cell_idx_c(0), z0);
   initSlice(&bubbleField, cell_idx_c(0), z1);
   initSlice(&bubbleField, cell_idx_c(0), z2);

   WALBERLA_CHECK(BubbleModelTest< stencil::D3Q19 >::checkForSplit(
                     &bubbleField, Cell(cell_idx_c(1), cell_idx_c(1), cell_idx_c(1)), 2) == false);
   WALBERLA_CHECK(BubbleModelTest< stencil::D3Q27 >::checkForSplit(
                     &bubbleField, Cell(cell_idx_c(1), cell_idx_c(1), cell_idx_c(1)), 2) == false);
}

void test3D_notConnected()
{
   // create a 3x3x3 bubble field (without ghost layer) and initialize with invalid IDs
   BubbleField_T bubbleField(uint_c(3), uint_c(3), uint_c(3), uint_c(0), INVALID_BUBBLE_ID);

   // initialize whole bubble field (disconnected bubble)
   const BubbleID N    = INVALID_BUBBLE_ID;
   const BubbleID z0[] = { 2, 2, 2, N, N, N, 2, 2, 2 };
   const BubbleID z1[] = { N, N, N, N, N, N, N, N, N };
   const BubbleID z2[] = { N, 2, N, N, 2, N, N, 2, N };

   initSlice(&bubbleField, 0, z0);
   initSlice(&bubbleField, 1, z1);
   initSlice(&bubbleField, 2, z2);

   WALBERLA_CHECK(BubbleModelTest< stencil::D3Q19 >::checkForSplit(
                     &bubbleField, Cell(cell_idx_c(1), cell_idx_c(1), cell_idx_c(1)), 2) == true);
   WALBERLA_CHECK(BubbleModelTest< stencil::D3Q27 >::checkForSplit(
                     &bubbleField, Cell(cell_idx_c(1), cell_idx_c(1), cell_idx_c(1)), 2) == true);
}

int main(int argc, char** argv)
{
   debug::enterTestMode();
   Environment env(argc, argv);

   test2D_notConnected();
   test2D_connected();

   test3D_connected();
   test3D_notConnected();

   return EXIT_SUCCESS;
}
} // namespace SplitDetectionTest
} // namespace free_surface
} // namespace walberla

int main(int argc, char** argv) { return walberla::free_surface::SplitDetectionTest::main(argc, argv); }