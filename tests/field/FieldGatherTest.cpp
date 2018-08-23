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
//! \file CollectTest.h
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================


#include "core/debug/TestSubsystem.h"
#include "core/Environment.h"

#include "blockforest/Initialization.h"

#include "field/Gather.h"
#include "field/AddToStorage.h"


namespace walberla {


int main( int argc, char ** argv )
{
   debug::enterTestMode();
   walberla::Environment walberlaEnv( argc, argv );

   bool oneBlockPerProcess = MPIManager::instance()->numProcesses() != 1;
   auto  blocks = blockforest::createUniformBlockGrid( 1, 3, 1,       // blocks in x,y,z
                                                       10u, 20u, 1u,  // nr of cells per block
                                                       1.0,           // dx
                                                       oneBlockPerProcess
                                                       );
   typedef GhostLayerField<cell_idx_t,3> MyField;
   BlockDataID fieldID = field::addToStorage<MyField>( blocks, "Field" );


   for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt )
   {
      MyField * field = blockIt->getData<MyField>( fieldID );

      for( auto cellIt = field->beginXYZ(); cellIt != field->end(); ++cellIt )
      {
         Cell globalCell = cellIt.cell();
         blocks->transformBlockLocalToGlobalCell( globalCell, *blockIt );
         cellIt.getF( 0 ) = globalCell[0];
         cellIt.getF( 1 ) = globalCell[1];
         cellIt.getF( 2 ) = globalCell[2];
      }
   }


   MyField gatheredField(0,0,0,0);
   CellInterval boundingBox = blocks->getDomainCellBB();
   boundingBox.min()[ 1 ] = 10;
   boundingBox.max()[ 1 ] = 29;

   auto targetRank = MPIManager::instance()->numProcesses() -1;
   field::gather<MyField>( gatheredField, blocks, fieldID, boundingBox, targetRank );

   WALBERLA_EXCLUSIVE_WORLD_SECTION( targetRank )
   {
      for( auto cellIt = gatheredField.beginXYZ(); cellIt != gatheredField.end(); ++cellIt )
      {
         WALBERLA_CHECK_EQUAL( cellIt.getF(0) - boundingBox.min()[0], cellIt.cell()[0] );
         WALBERLA_CHECK_EQUAL( cellIt.getF(1) - boundingBox.min()[1], cellIt.cell()[1] );
         WALBERLA_CHECK_EQUAL( cellIt.getF(2) - boundingBox.min()[2], cellIt.cell()[2] );
      }
   }


   return 0;
}
}

int main( int argc, char ** argv )
{
   return walberla::main(argc,argv);
}
