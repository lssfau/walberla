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
//! \file FieldPackInfoTest.cpp
//! \ingroup field
//! \author Martin Bauer <martin.bauer@fau.de>
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//! \brief Tests if a Field is correctly packed into buffers
//
//======================================================================================================================

#include "field/AddToStorage.h"
#include "field/GhostLayerField.h"
#include "field/communication/PackInfo.h"
#include "field/communication/UniformPullReductionPackInfo.h"

#include "blockforest/Initialization.h"

#include "core/debug/TestSubsystem.h"
#include "core/mpi/MPIManager.h"

#include <cstring>


namespace walberla {


void testScalarField( IBlock * block, BlockDataID fieldId )
{
   GhostLayerField<int,1> & field = *(block->getData<GhostLayerField<int,1> > (fieldId));
   field.setWithGhostLayer( 0 );

   WALBERLA_CHECK_EQUAL(field.xSize(), 2);
   WALBERLA_CHECK_EQUAL(field.ySize(), 2);
   WALBERLA_CHECK_EQUAL(field.zSize(), 2);

   // initialize the bottom boundary
   field(0,0,0) = 1;
   field(0,1,0) = 2;
   field(1,0,0) = 3;
   field(1,1,0) = 4;

   // -------------- Local Communication Test ----------------------

   // communicate periodic from bottom to top
   field::communication::PackInfo< GhostLayerField<int,1> > pi (fieldId);
   pi.communicateLocal( block, block, stencil::B );

   WALBERLA_CHECK_EQUAL ( field(0,0,+2), 1 );
   WALBERLA_CHECK_EQUAL ( field(0,1,+2), 2 );
   WALBERLA_CHECK_EQUAL ( field(1,0,+2), 3 );
   WALBERLA_CHECK_EQUAL ( field(1,1,+2), 4 );

   // -------------- Buffer Communication Test ---------------------

   // Reset
   field(0,0,2) = 0;
   field(0,1,2) = 0;
   field(1,0,2) = 0;
   field(1,1,2) = 0;

   mpi::GenericSendBuffer<> sendBuf;
   pi.packData( block, stencil::B, sendBuf );

   // Manually copy over the send to the receive buffer
   mpi::GenericRecvBuffer<> recvBuf;
   recvBuf.resize( sendBuf.size() );
   memcpy( recvBuf.ptr(), sendBuf.ptr(), sendBuf.size()* sizeof(mpi::GenericSendBuffer<>::ElementType) );

   pi.unpackData( block, stencil::T, recvBuf );

   WALBERLA_CHECK_EQUAL ( field(0,0,+2), 1 );
   WALBERLA_CHECK_EQUAL ( field(0,1,+2), 2 );
   WALBERLA_CHECK_EQUAL ( field(1,0,+2), 3 );
   WALBERLA_CHECK_EQUAL ( field(1,1,+2), 4 );
}

void testScalarFieldPullReduction( IBlock * block, BlockDataID fieldId )
{
   GhostLayerField<int,1> & field = *(block->getData<GhostLayerField<int,1> > (fieldId));
   field.setWithGhostLayer( 0 );

   WALBERLA_CHECK_EQUAL(field.xSize(), 2);
   WALBERLA_CHECK_EQUAL(field.ySize(), 2);
   WALBERLA_CHECK_EQUAL(field.zSize(), 2);

   // initialize the bottom ghost layer cells
   field(0,0,-1) = 1;
   field(0,1,-1) = 2;
   field(1,0,-1) = 3;
   field(1,1,-1) = 4;

   // initialize the top interior cells
   field(0,0,1) = 1;
   field(0,1,1) = 1;
   field(1,0,1) = 1;
   field(1,1,1) = 1;

   // communicate periodic from bottom to top with uniform pull scheme
   field::communication::UniformPullReductionPackInfo<std::plus, GhostLayerField<int,1> > pi1 (fieldId);
   pi1.communicateLocal( block, block, stencil::B );

   // check values in top ghost layer
   WALBERLA_CHECK_EQUAL ( field(0,0,2), 0 );
   WALBERLA_CHECK_EQUAL ( field(0,1,2), 0 );
   WALBERLA_CHECK_EQUAL ( field(1,0,2), 0 );
   WALBERLA_CHECK_EQUAL ( field(1,1,2), 0 );

   // check values in top interior cells
   WALBERLA_CHECK_EQUAL ( field(0,0,1), 2 );
   WALBERLA_CHECK_EQUAL ( field(0,1,1), 3 );
   WALBERLA_CHECK_EQUAL ( field(1,0,1), 4 );
   WALBERLA_CHECK_EQUAL ( field(1,1,1), 5 );

   // communicate periodic from top to bottom with standard form to sync ghost layers
   field::communication::PackInfo< GhostLayerField<int,1> > pi2 (fieldId);
   pi2.communicateLocal( block, block, stencil::T );

   // check values in bottom ghost layer
   WALBERLA_CHECK_EQUAL ( field(0,0,-1), 2 );
   WALBERLA_CHECK_EQUAL ( field(0,1,-1), 3 );
   WALBERLA_CHECK_EQUAL ( field(1,0,-1), 4 );
   WALBERLA_CHECK_EQUAL ( field(1,1,-1), 5 );

   // check values in top interior cells
   WALBERLA_CHECK_EQUAL ( field(0,0,1), 2 );
   WALBERLA_CHECK_EQUAL ( field(0,1,1), 3 );
   WALBERLA_CHECK_EQUAL ( field(1,0,1), 4 );
   WALBERLA_CHECK_EQUAL ( field(1,1,1), 5 );

}

int main(int argc, char **argv)
{
   using blockforest::createUniformBlockGrid;

   debug::enterTestMode();
   MPIManager::instance()->initializeMPI(&argc,&argv);


   // Create a BlockForest with 2x2x2 cells per block
   uint_t processes = uint_c( MPIManager::instance()->numProcesses() );
   auto blocks = createUniformBlockGrid(processes,1 ,1, //blocks
                                        2,2,2,          //cells
                                        1,              //dx
                                        false,          //one block per process
                                        true,true,true);//periodicity

   // Create a Field with the same number of cells as the block
   BlockDataID scalarFieldId = field::addToStorage<GhostLayerField<int,1> > ( blocks, "ScalarField" );

   for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt ) // block loop
      testScalarField( &(*blockIt), scalarFieldId );

   // Create a BlockForest with 8x8x8 cells per block
   blocks = createUniformBlockGrid(processes,1 ,1, //blocks
                                   8,8,8,          //cells
                                   1,              //dx
                                   false,          //one block per process
                                   true,true,true);//periodicity

   // Create a Field with one quarter as many cells per dimension, i.e. a field with the same size as the one above
   auto getSize = []( const shared_ptr< StructuredBlockStorage > &, IBlock * const ) {
      return Vector3<uint_t>(2,2,2);
   };
   scalarFieldId = field::addToStorage<GhostLayerField<int,1> > ( blocks, "ScalarField", getSize );

   for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt ) // block loop
      testScalarField( &(*blockIt), scalarFieldId );

   for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt ) // block loop
      testScalarFieldPullReduction( &(*blockIt), scalarFieldId );

   return 0;
}
} // namespace walberla

int main( int argc, char* argv[] )
{
  return walberla::main( argc, argv );
}