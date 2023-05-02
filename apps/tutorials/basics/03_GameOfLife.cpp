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
//! \file 03_GameOfLife.cpp
//! \author Martin Bauer <martin.bauer@fau.de>
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#include "blockforest/Initialization.h"
#include "blockforest/communication/UniformBufferedScheme.h"

#include "core/Environment.h"

#include "field/AddToStorage.h"
#include "field/communication/PackInfo.h"

#include "geometry/initializer/ScalarFieldFromGrayScaleImage.h"
#include "geometry/structured/GrayScaleImage.h"

#include "gui/Gui.h"

#include "stencil/D2Q9.h"

#include "timeloop/SweepTimeloop.h"


namespace walberla {


class GameOfLifeSweep
{
public:

   GameOfLifeSweep( BlockDataID fieldID )
      : fieldID_( fieldID )
   {}

   void operator()( IBlock * block );

private:

   BlockDataID fieldID_;
};

void GameOfLifeSweep::operator()( IBlock * block )
{
   ///////////////////////////////////
   // version using field iterators //
   ///////////////////////////////////

   /*
   auto field = block->getData< GhostLayerField<real_t,1> >( fieldID_ );

   shared_ptr< GhostLayerField<real_t,1> > copy( field->clone() );

   auto fieldIter = field->begin();
   auto  copyIter =  copy->begin();

   while( fieldIter != field->end() )
   {
      // count number of living neighbors
      int liveNeighbors = 0;
      for( auto dir = stencil::D2Q9::beginNoCenter(); dir != stencil::D2Q9::end(); ++dir )
         if( copyIter.neighbor( *dir ) > real_c(0.5) )
            ++liveNeighbors;

      // cell dies because of under- or over-population
      if( liveNeighbors < 2 || liveNeighbors > 3 )
         *fieldIter = real_c(0);
      if( liveNeighbors == 3 ) // cell comes alive
         *fieldIter = real_c(1);

      ++fieldIter;
      ++copyIter;
   }
   */

   /////////////////////////////////////////
   // version using field iterators macro //
   /////////////////////////////////////////

   /*
   auto field = block->getData< GhostLayerField<real_t,1> >( fieldID_ );

   shared_ptr< GhostLayerField<real_t,1> > copy( field->clone() );

   WALBERLA_FOR_ALL_CELLS( fieldIter, field, copyIter, copy,

      // count number of living neighbors
      int liveNeighbors = 0;
      for( auto dir = stencil::D2Q9::beginNoCenter(); dir != stencil::D2Q9::end(); ++dir )
         if( copyIter.neighbor( *dir ) > real_c(0.5) )
            ++liveNeighbors;

      // cell dies because of under- or over-population
      if( liveNeighbors < 2 || liveNeighbors > 3 )
         *fieldIter = real_c(0);
      if( liveNeighbors == 3 ) // cell comes alive
         *fieldIter = real_c(1);

   ) // WALBERLA_FOR_ALL_CELLS
   */

   ////////////////////////////
   // version using xyz loop //
   ////////////////////////////

   /*
   auto field = block->getData< GhostLayerField<real_t,1> >( fieldID_ );

   shared_ptr< GhostLayerField<real_t,1> > copy( field->clone() );

   auto xyz = field->xyzSize();

   for( cell_idx_t z = xyz.zMin(); z <= xyz.zMax(); ++z ) {
      for( cell_idx_t y = xyz.yMin(); y <= xyz.yMax(); ++y ) {
         for( cell_idx_t x = xyz.xMin(); x <= xyz.xMax(); ++x )
         {
            // count number of living neighbors
            int liveNeighbors = 0;
            for( auto dir = stencil::D2Q9::beginNoCenter(); dir != stencil::D2Q9::end(); ++dir )
               if( copy->getNeighbor( x, y, z, *dir ) > real_c(0.5) )
                  ++liveNeighbors;

            // cell dies because of under- or over-population
            if( liveNeighbors < 2 || liveNeighbors > 3 )
               field->get(x,y,z) = real_c(0);
            if( liveNeighbors == 3 ) // cell comes alive
               field->get(x,y,z) = real_c(1);
         }
      }
   }
   */

   //////////////////////////////////
   // version using xyz loop macro //
   //////////////////////////////////

   auto field = block->getData< GhostLayerField<real_t,1> >( fieldID_ );

   shared_ptr< GhostLayerField<real_t,1> > copy( field->clone() );

   WALBERLA_FOR_ALL_CELLS_XYZ( field,

      // count number of living neighbors
      int liveNeighbors = 0;
      for( auto dir = stencil::D2Q9::beginNoCenter(); dir != stencil::D2Q9::end(); ++dir )
         if( copy->getNeighbor( x, y, z, *dir ) > real_c(0.5) )
            ++liveNeighbors;

      // cell dies because of under- or over-population
      if( liveNeighbors < 2 || liveNeighbors > 3 )
         field->get(x,y,z) = real_c(0);
      if( liveNeighbors == 3 ) // cell comes alive
         field->get(x,y,z) = real_c(1);

   ) // WALBERLA_FOR_ALL_CELLS_XYZ
}



int main( int argc, char ** argv )
{
   walberla::Environment env( argc, argv );

   geometry::GrayScaleImage image( "GosperGliderGun.png" );

   // create blocks
   shared_ptr< StructuredBlockForest > blocks = blockforest::createUniformBlockGrid(
            uint_c(1) ,              uint_c(2),                           uint_c(1), // number of blocks in x,y,z direction
            image.size( uint_c(0) ), image.size( uint_c(1) ) / uint_c(2), uint_c(1), // how many cells per block (x,y,z)
            real_c(1),                                                               // dx: length of one cell in physical coordinates
            false,                                                                   // one block per process? - "false" means all blocks to one process
            false, false, false );                                                   // no periodicity

   using ScalarField = GhostLayerField<real_t, 1>;
   BlockDataID fieldID = field::addToStorage<ScalarField>( blocks,      // block storage
                                                           "My Field",  // name
                                                           real_c(0),   // initial value
                                                           field::fzyx, // layout (not relevant for scalar fields)
                                                           uint_c(1)    // number of ghost layers
                                                           );

   // create the scheme ...
   blockforest::communication::UniformBufferedScheme<stencil::D2Q9> myCommScheme( blocks );
   // ... and add a PackInfo that packs/unpacks our field
   myCommScheme.addPackInfo( make_shared< field::communication::PackInfo<ScalarField> >( fieldID  ) );

   // initializing the field from an image
   using geometry::initializer::ScalarFieldFromGrayScaleImage;
   ScalarFieldFromGrayScaleImage fieldInitializer( *blocks, fieldID );
   fieldInitializer.init( image, uint_c(2), false );

   // create time loop
   const uint_t numberOfTimesteps = uint_c(10); // number of timesteps for non-gui runs
   SweepTimeloop timeloop( blocks, numberOfTimesteps );

   // registering the sweep
   timeloop.add() << BeforeFunction( myCommScheme, "Communication" )
                  << Sweep( GameOfLifeSweep(fieldID), "GameOfLifeSweep" );

   GUI gui ( timeloop, blocks, argc, argv );
   gui.run();

   return EXIT_SUCCESS;
}
} // namespace walberla

int main( int argc, char* argv[] )
{
  return walberla::main( argc, argv );
}