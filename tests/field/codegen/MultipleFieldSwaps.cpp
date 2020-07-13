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
//! \file MultipleFieldSwaps.cpp
//! \author Markus Holzer <markus.holzer@fau.de>
//
//======================================================================================================================

#include "MultipleFieldSwaps.h"
#include "blockforest/Initialization.h"
#include "blockforest/communication/UniformDirectScheme.h"

#include "core/Environment.h"
#include "core/debug/TestSubsystem.h"

#include "field/AddToStorage.h"

#include "stencil/D2Q9.h"


using namespace walberla;

typedef GhostLayerField<double,1> ScalarField;
void testMultipleFieldSwaps()
{
   uint_t xSize = 5;
   uint_t ySize = 5;
   // Create blocks
   shared_ptr< StructuredBlockForest > blocks = blockforest::createUniformBlockGrid (
      uint_t(1) , uint_t(1),  uint_t(1),  // number of blocks in x,y,z direction
      xSize, ySize, uint_t(1),            // how many cells per block (x,y,z)
      real_t(1),                          // dx: length of one cell in physical coordinates
      false,                              // one block per process - "false" means all blocks to one process
      true, true, true );                 // full periodicity


   BlockDataID fieldID_1 = field::addToStorage<ScalarField>(blocks, "Field_1", real_t(1.0));
   BlockDataID fieldID_2 = field::addToStorage<ScalarField>(blocks, "Field_2", real_t(1.0));
   BlockDataID fieldID_3 = field::addToStorage<ScalarField>(blocks, "Field_3", real_t(1.0));

   pystencils::MultipleFieldSwaps kernel(fieldID_1, fieldID_2, fieldID_3);

   for (auto &block: *blocks)
      kernel(&block);

   for (auto& block : *blocks)
   {
      auto field_1 = block.getData< ScalarField >(fieldID_1);
      auto field_2 = block.getData< ScalarField >(fieldID_2);
      auto field_3 = block.getData< ScalarField >(fieldID_3);
      // clang-format off
      WALBERLA_FOR_ALL_CELLS_XYZ(field_1, Cell globalCell;
         blocks->transformBlockLocalToGlobalCell(globalCell, block, Cell(x, y, z));
         WALBERLA_CHECK_FLOAT_EQUAL(field_1->get(x, y, z), real_t(2.0))
         WALBERLA_CHECK_FLOAT_EQUAL(field_2->get(x, y, z), real_t(2.0))
         WALBERLA_CHECK_FLOAT_EQUAL(field_3->get(x, y, z), real_t(2.0))
      )
      // clang-format on
   }
}



int main( int argc, char ** argv )
{
   mpi::Environment env( argc, argv );
   debug::enterTestMode();

   testMultipleFieldSwaps();

   return EXIT_SUCCESS;
}