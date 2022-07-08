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
//! \file RegionalFloodFillTest.cpp
//! \ingroup lbm/free_surface/bubble_model
//! \author Martin Bauer
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//! \brief Test flood fill algorithm with D3Q7 and D3Q19 stencil.
//
//======================================================================================================================

#include "lbm/free_surface/bubble_model/RegionalFloodFill.h"

#include "core/Environment.h"
#include "core/debug/TestSubsystem.h"

#include "field/Printers.h"

#include "stencil/D3Q19.h"
#include "stencil/D3Q7.h"

namespace walberla
{
namespace free_surface
{
namespace RegionalFloodFillTest
{
int main(int argc, char** argv)
{
   debug::enterTestMode();
   Environment walberlaEnv(argc, argv);

   // create 3x3 field with 1 ghost layer (initialized with 0)
   GhostLayerField< int, 1 > field(3, 3, 1, 1, 0);

   // initial values of the field; the flood fill starting cell is marked with /**/;
   // attention: the y coordinate is upside down in this array
   const int initValues[5][5] = {
      { 1, 1, 1, 1, 1 }, { 0, 0 /**/, 0, 1, 0 }, { 1, 0, 0, 1, 0 }, { 1, 0, 0, 1, 0 }, { 1, 1, 1, 0, 0 },
   };

   // test scenario: detect the connection from (1,0) and (0,2) with starting cell (0,0)
   for (cell_idx_t y = cell_idx_c(-1); y < cell_idx_c(4); ++y)
   {
      for (cell_idx_t x = cell_idx_c(-1); x < cell_idx_c(4); ++x)
      {
         field(x, y, cell_idx_c(0)) = initValues[y + 1][x + 1];
      }
   }

   // print the initialized field (for debugging purposes)
   //   std::cout << "Initialized field:" << std::endl;
   //   field::printSlice(std::cout, field, 2, 0);

   // connection should not be found since search neighborhood (2) is too small
   bubble_model::RegionalFloodFill< int, stencil::D3Q19 > neigh2_D3Q19(
      &field, Cell(cell_idx_c(0), cell_idx_c(0), cell_idx_c(0)), stencil::S, 1, cell_idx_c(2));
   WALBERLA_CHECK_EQUAL(neigh2_D3Q19.connected(stencil::NW), false);

   // connection should be found since search neighborhood (3) is large enough
   bubble_model::RegionalFloodFill< int, stencil::D3Q19 > neigh3_D3Q19(
      &field, Cell(cell_idx_c(0), cell_idx_c(0), cell_idx_c(0)), stencil::S, 1, cell_idx_c(3));
   WALBERLA_CHECK_EQUAL(neigh3_D3Q19.connected(stencil::NW), true);

   // connection should be found since search neighborhood (3) is large enough
   bubble_model::RegionalFloodFill< int, stencil::D3Q7 > neigh3_D3Q7(
      &field, Cell(cell_idx_c(0), cell_idx_c(0), cell_idx_c(0)), stencil::S, 1, cell_idx_c(3));
   WALBERLA_CHECK_EQUAL(neigh3_D3Q7.connected(stencil::NW), false);

   return EXIT_SUCCESS;
}
} // namespace RegionalFloodFillTest
} // namespace free_surface
} // namespace walberla

int main(int argc, char** argv) { return walberla::free_surface::RegionalFloodFillTest::main(argc, argv); }