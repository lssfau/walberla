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
//! \file 02_Sweeps.cpp
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#include "blockforest/Initialization.h"
#include "core/Environment.h"
#include "field/Field.h"
#include "gui/Gui.h"
#include "timeloop/SweepTimeloop.h"

#include <functional>


namespace walberla {

// some arbitrary value for our bogus algorithm
const int ARBITRARY_VALUE = 424242;



Field<real_t,1> * createFields( IBlock * const block, StructuredBlockStorage * const storage )
{
   return new Field<real_t,1> ( storage->getNumberOfXCells( *block ), // number of cells in x direction for this block
                                storage->getNumberOfYCells( *block ), // number of cells in y direction for this block
                                storage->getNumberOfZCells( *block ), // number of cells in z direction for this block
                                real_c(0) );                          // initial value
}


// function sweep
void simpleSweep( IBlock * block, BlockDataID fieldBlockDataID )
{
   // retrieve the field from the block
   auto field = block->getData< Field<real_t,1> >( fieldBlockDataID );

   // some bogus "algorithm"
   for( auto iter = field->begin(); iter != field->end(); ++iter )
   {
      if( *iter > real_c(ARBITRARY_VALUE) )
         *iter /= real_c(2);
      else
         *iter *= real_c(2);
   }
}


// the same "algorithm", but now as Class Sweep
class SimpleSweep
{
public:

   SimpleSweep( BlockDataID fieldID )
      : fieldID_( fieldID )
   {}

   void operator()( IBlock * block )
   {
      // the fieldID is now a member! (was set by the constructor)
      auto field = block->getData< Field<real_t,1> >( fieldID_ );

      // some bogus "algorithm"
      for( auto iter = field->begin(); iter != field->end(); ++iter )
      {
         if( *iter > real_c(ARBITRARY_VALUE) )
            *iter /= real_c(2);
         else
            *iter *= real_c(2);
      }
   }

private:

   BlockDataID fieldID_;
};



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

   // add a field to all blocks - and store the returned block data ID which is needed to access the field
   BlockDataID fieldID = blocks->addStructuredBlockData< Field<real_t,1> >( &createFields, "My Field" );

   // iterate all blocks and initialize with random values
   for( auto blockIterator = blocks->begin(); blockIterator != blocks->end(); ++blockIterator )
   {
      IBlock & currentBlock = *blockIterator;

      // get the field stored on the current block
      Field<real_t,1> * field = currentBlock.getData< Field<real_t,1> >( fieldID );

      // iterate the field and set random values
      for( auto iter = field->begin(); iter != field->end(); ++iter )
         *iter = real_c( rand() % ARBITRARY_VALUE );
   }

   // create time loop
   const uint_t numberOfTimesteps = uint_c(10); // number of time steps for non-gui runs
   SweepTimeloop timeloop( blocks, numberOfTimesteps );

   // registering the function sweep
   auto pointerToTwoArgFunction = & simpleSweep;
   auto pointerToOneArgFunction = std::bind( pointerToTwoArgFunction, std::placeholders::_1, fieldID );
   timeloop.add() << Sweep( pointerToOneArgFunction, "BogusAlgorithm" );

   // registering the class sweep
   timeloop.add() << Sweep( SimpleSweep(fieldID), "BogusAlgorithmButNowAsFunctor" );

   // two sweeps were registered, so both are executed in each time step

   GUI gui( timeloop, blocks, argc, argv );
   gui.run();

   return EXIT_SUCCESS;
}
} // namespace walberla

int main( int argc, char* argv[] )
{
  return walberla::main( argc, argv );
}