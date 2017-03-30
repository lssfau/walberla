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
//! \file 01_BlocksAndFields.cpp
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#include "blockforest/Initialization.h"
#include "core/Environment.h"
#include "field/Field.h"
#include "gui/Gui.h"
#include "timeloop/SweepTimeloop.h"


using namespace walberla;


Field<real_t,1> * createFields( IBlock * const block, StructuredBlockStorage * const storage )
{
   return new Field<real_t,1> ( storage->getNumberOfXCells( *block ), // number of cells in x direction for this block
                                storage->getNumberOfYCells( *block ), // number of cells in y direction for this block
                                storage->getNumberOfZCells( *block ), // number of cells in z direction for this block
                                real_c(0) );                          // initial value
}


int main( int argc, char ** argv )
{
   walberla::Environment env( argc, argv );

   // create blocks
   shared_ptr< StructuredBlockForest > blocks = blockforest::createUniformBlockGrid(
             uint_c( 3), uint_c(2), uint_c( 4), // number of blocks in x,y,z direction
             uint_c(10), uint_c(8), uint_c(12), // how many cells per block (x,y,z)
             real_c(0.5),                       // dx: length of one cell in physical coordinates
             false,                             // one block per process? - "false" means all blocks to one process
             false, false, false );             // no periodicity

   // add a field to all blocks
   blocks->addStructuredBlockData< Field<real_t,1> >( &createFields, "My Field" );

   SweepTimeloop timeloop( blocks, uint_c(1) );

   GUI gui( timeloop, blocks, argc, argv );
   gui.run();

   return EXIT_SUCCESS;
}
