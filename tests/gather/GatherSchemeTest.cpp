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
//! \file GatherSchemeTest.cpp
//! \ingroup gather
//! \author Martin Bauer <martin.bauer@fau.de>
//! \brief Tests MPIGatherScheme, FileGatherScheme and CellGatherPackInfo
//
//======================================================================================================================

#include "gather/CellGatherPackInfo.h"
#include "gather/FileGatherScheme.h"
#include "gather/GnuPlotGraphWriter.h"
#include "gather/MPIGatherScheme.h"

#include "blockforest/Initialization.h"

#include "core/Environment.h"
#include "core/cell/CellInterval.h"
#include "core/debug/TestSubsystem.h"

#include "field/GhostLayerField.h"

#include <iostream>


namespace walberla {

typedef GhostLayerField<cell_idx_t, 1> ScalarField;



/**
 * Test Description:
 *
 * Domain Setup:
 *
 * -------------------------------------
 * |        |        |        |        |
 * |        |     lllllllll   |        |
 * |        |        |        |        |
 * -------------------------------------
 * |        |        |        |        |
 * |  lllllllllll    |        |        |
 * |        |        |        |        |
 * -------------------------------------
 *     P0        P1       P2       P3
 *
 * y
 * |
 * -- x
 *
 * Data along two lines is gathered, each line is a separate PackInfo.
 *    - 4 processes, 2 blocks per process
 *    - there are processes that send data for one, two and no PackInfo
 *    - one scalar field, the value is the same as global x coordinate
 *      for easier test validation
 *
 */


ScalarField * createField( IBlock* const block, StructuredBlockStorage* const storage )
{
   auto f = new  ScalarField (storage->getNumberOfXCells( *block ),
                              storage->getNumberOfYCells( *block ),
                              storage->getNumberOfZCells( *block ),
                              1);

   for ( auto i = f->beginWithGhostLayer(); i != f->end(); ++i )
   {
      Cell globalCell;
      storage->transformBlockLocalToGlobalCell( globalCell, *block, i.cell() );
      *i = globalCell.x();
   }
   return f;
}


class CheckingDataProcessor : public gather::DataProcessor
{
public:
   CheckingDataProcessor( const CellInterval & ci )
        : cellInterval_ ( ci )
    {}

   void process(const std::vector<std::vector<real_t> > & data) override
   {
      uint_t counter =0;
      for( auto it = cellInterval_.begin(); it != cellInterval_.end(); ++it )
      {
         WALBERLA_CHECK_FLOAT_EQUAL( real_c( it->x() ), data[counter][1] );
         counter++;
      }
   }

private:
   CellInterval cellInterval_;
};



int main(int argc, char ** argv )
{
   walberla::Environment env ( argc, argv );
   debug::enterTestMode();

   try {

   auto mpiManager = MPIManager::instance();



   using blockforest::createUniformBlockGrid;


   WALBERLA_CHECK_EQUAL( mpiManager->numProcesses(), 4 );
   auto  blocks = createUniformBlockGrid( 4, 2, 1,       // blocks in x,y,z
                                          10u, 10u, 10u, // nr of cells per block
                                          1.0,           // dx
                                          uint_t(4), uint_t(1), uint_t(1)    // nr of processes in x,y,z
                                          );
   /*
   auto blocks = createUniformBlockGrid( 4,2,1,
                                         10u,10u,10u,
                                         1.0,
                                         false,              // one block per process
                                         false,false,false,  // periodicity
                                         false );
   */

   BlockDataID fieldID = blocks->addStructuredBlockData<ScalarField>( &createField, "MyField" );


   using namespace gather;


   typedef CellGatherPackInfo<ScalarField, CellInterval > CellGatherPI;

   // PackInfo across process P0 and P1
   CellInterval lowerLeftLine (   5, 5, 5,
                                 14, 5, 5 );
   auto dp1 = make_shared<CheckingDataProcessor>( lowerLeftLine );
   auto lowerLeftPI = make_shared<CellGatherPI> ( blocks, fieldID, lowerLeftLine, dp1 );


   // PackInfo across process P1 and P2
   CellInterval upperRightLine ( 15,15, 5,
                                 24,15, 5 );
   auto dp2 = make_shared<CheckingDataProcessor>( upperRightLine );
   auto upperRightPI = make_shared<CellGatherPI> ( blocks, fieldID, upperRightLine, dp2 );



   // Test MPI Gather
   MPIGatherScheme mpiGatherScheme( blocks->getBlockStorage() );
   mpiGatherScheme.addPackInfo( lowerLeftPI  );
   mpiGatherScheme.addPackInfo( upperRightPI );

   for(int i=0; i< 10; ++i )
      mpiGatherScheme();


   // Test File Gather
   FileGatherScheme fileGatherScheme( blocks->getBlockStorage() );
   fileGatherScheme.addPackInfo( lowerLeftPI  );
   fileGatherScheme.addPackInfo( upperRightPI );

   for(int i=0; i< 10; ++i )
      fileGatherScheme();

   }
   catch( std::exception & e )
   {
      std::cout<<  "Caught Exception: " <<std::endl;
      std::cout<<  e.what() <<std::endl;
      std::cout<<  "Description end.. rethrowing..." <<std::endl;
      throw e;
   }
   return 0;
}
} // namespace walberla

int main( int argc, char* argv[] )
{
  return walberla::main( argc, argv );
}