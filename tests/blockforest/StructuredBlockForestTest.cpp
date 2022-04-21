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
//! \file StructuredBlockForestTest.cpp
//! \ingroup blockforest
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#include "blockforest/SetupBlockForest.h"
#include "blockforest/StructuredBlockForest.h"

#include "core/debug/TestSubsystem.h"
#include "core/mpi/Environment.h"

#include <random>

namespace walberla {
namespace blockforest {



static void rootBlockExclusionFunction( std::vector<uint8_t>& excludeBlock, const SetupBlockForest::RootBlockAABB& aabb ) {

   for( uint_t i = 0; i != excludeBlock.size(); ++i )
   {
      AABB bb = aabb(i);
      if( bb.contains( real_c(25), real_c(25), real_c(25) ) )
         excludeBlock[i] = 1;
   }
}



static void refinementSelectionFunction( SetupBlockForest& forest ) {

   SetupBlock* block = forest.getRootBlock( 7 );

   if( !block->hasChildren() )
      block->setMarker( true );
}



static void workloadMemorySUIDAssignmentFunction( SetupBlockForest& forest ) {

   std::vector< SetupBlock* > blocks;
   forest.getBlocks( blocks );

   for( uint_t i = 0; i != blocks.size(); ++i ) {
      blocks[i]->setMemory( 1.0 );
      blocks[i]->setWorkload( 1.0 );
   }
}



static uint_t* blockdata23( const IBlock* const /*block*/ ) {

   return new uint_t(23);
}

static int* blockdata42( const IBlock* const /*block*/ ) {

   return new int(42);
}



static uint_t* blockdataZ( const IBlock* const block, const StructuredBlockStorage* const storage ) {

   return new uint_t( storage->getNumberOfZCells(*block) );
}

static uint_t* blockdataY( const IBlock* const block, const StructuredBlockStorage* const storage ) {

   return new uint_t( storage->getNumberOfYCells(*block) );
}



static void test() {

   // SetupBlockForest (determine domain decomposition & perform process distribution)

   SetupBlockForest sforest;

   sforest.addRootBlockExclusionFunction( rootBlockExclusionFunction );
   sforest.addRefinementSelectionFunction( refinementSelectionFunction );
   sforest.addWorkloadMemorySUIDAssignmentFunction( workloadMemorySUIDAssignmentFunction );

   AABB domain( 0, 0, 0, 100, 100, 100 );
   sforest.init( domain, 2, 2, 2, true, false, false );

// sforest.writeVTKOutput( "SetupBlockForest_0" );

   sforest.assignAllBlocksToRootProcess();

// sforest.writeVTKOutput( "SetupBlockForest_1" );

   // setup StructuredBlockForest

   StructuredBlockForest forest( make_shared< BlockForest >( 0, sforest, true ), 50, 60, 70 );
   forest.createCellBoundingBoxes();

   BlockDataID bdid23 = forest.addBlockData( "23" ) << BlockDataCreator< uint_t >( blockdata23 );
   BlockDataID bdid42 = forest.addBlockData< int >( blockdata42, "42" );

   BlockDataID bdidZ = forest.addStructuredBlockData( "Z" ) << StructuredBlockDataCreator< uint_t >( blockdataZ );
   BlockDataID bdidY = forest.addStructuredBlockData< uint_t >( blockdataY, "Y" );

   // block data

   for( auto block = forest.begin(); block != forest.end(); ++block ) {
      WALBERLA_CHECK_EQUAL( *(block->getData< uint_t >( bdid23 )), 23 )
      WALBERLA_CHECK_EQUAL( *(block->getData< int    >( bdid42 )), 42 )
      WALBERLA_CHECK_EQUAL( *(block->getData< uint_t >( bdidZ  )), 70 )
      WALBERLA_CHECK_EQUAL( *(block->getData< uint_t >( bdidY  )), 60 )
   }

   // getDomainCellBB

   CellInterval bb = forest.getDomainCellBB();
   WALBERLA_CHECK_EQUAL( bb.xMin(),   0 )
   WALBERLA_CHECK_EQUAL( bb.yMin(),   0 )
   WALBERLA_CHECK_EQUAL( bb.zMin(),   0 )
   WALBERLA_CHECK_EQUAL( bb.xMax(),  99 )
   WALBERLA_CHECK_EQUAL( bb.yMax(), 119 )
   WALBERLA_CHECK_EQUAL( bb.zMax(), 139 )

   bb = forest.getDomainCellBB(1);
   WALBERLA_CHECK_EQUAL( bb.xMin(),   0 )
   WALBERLA_CHECK_EQUAL( bb.yMin(),   0 )
   WALBERLA_CHECK_EQUAL( bb.zMin(),   0 )
   WALBERLA_CHECK_EQUAL( bb.xMax(), 199 )
   WALBERLA_CHECK_EQUAL( bb.yMax(), 239 )
   WALBERLA_CHECK_EQUAL( bb.zMax(), 279 )

   // getNumberOf*Cells

   WALBERLA_CHECK_EQUAL( forest.getNumberOfXCells(), 100 )
   WALBERLA_CHECK_EQUAL( forest.getNumberOfYCells(), 120 )
   WALBERLA_CHECK_EQUAL( forest.getNumberOfZCells(), 140 )

   WALBERLA_CHECK_EQUAL( forest.getNumberOfXCells(1), 200 )
   WALBERLA_CHECK_EQUAL( forest.getNumberOfYCells(1), 240 )
   WALBERLA_CHECK_EQUAL( forest.getNumberOfZCells(1), 280 )

   WALBERLA_CHECK_EQUAL( forest.getNumberOfCells(0), 100 )
   WALBERLA_CHECK_EQUAL( forest.getNumberOfCells(1), 120 )
   WALBERLA_CHECK_EQUAL( forest.getNumberOfCells(2), 140 )

   WALBERLA_CHECK_EQUAL( forest.getNumberOfCells(0,1), 200 )
   WALBERLA_CHECK_EQUAL( forest.getNumberOfCells(1,1), 240 )
   WALBERLA_CHECK_EQUAL( forest.getNumberOfCells(2,1), 280 )

   // dx/dy/dz

   WALBERLA_CHECK_FLOAT_EQUAL( forest.dx(), 100.0 / 100.0 )
   WALBERLA_CHECK_FLOAT_EQUAL( forest.dy(), 100.0 / 120.0 )
   WALBERLA_CHECK_FLOAT_EQUAL( forest.dz(), 100.0 / 140.0 )

   WALBERLA_CHECK_FLOAT_EQUAL( forest.dx(1), 100.0 / 200.0 )
   WALBERLA_CHECK_FLOAT_EQUAL( forest.dy(1), 100.0 / 240.0 )
   WALBERLA_CHECK_FLOAT_EQUAL( forest.dz(1), 100.0 / 280.0 )

   // mapToPeriodicDomain

   Cell cell( cell_idx_c(42), cell_idx_c(234), cell_idx_c(567) );
   forest.mapToPeriodicDomain( cell );

   WALBERLA_CHECK_EQUAL( cell.x(),  42 )
   WALBERLA_CHECK_EQUAL( cell.y(), 234 )
   WALBERLA_CHECK_EQUAL( cell.z(), 567 )

   cell.x() = cell_idx_c(-123);
   forest.mapToPeriodicDomain( cell );

   WALBERLA_CHECK_EQUAL( cell.x(),  77 )
   WALBERLA_CHECK_EQUAL( cell.y(), 234 )
   WALBERLA_CHECK_EQUAL( cell.z(), 567 )

   cell.x() = cell_idx_c(-223);
   forest.mapToPeriodicDomain( cell, 1 );

   WALBERLA_CHECK_EQUAL( cell.x(), 177 )
   WALBERLA_CHECK_EQUAL( cell.y(), 234 )
   WALBERLA_CHECK_EQUAL( cell.z(), 567 )

   cell.x() = cell_idx_c(123);
   forest.mapToPeriodicDomain( cell );

   WALBERLA_CHECK_EQUAL( cell.x(),  23 )
   WALBERLA_CHECK_EQUAL( cell.y(), 234 )
   WALBERLA_CHECK_EQUAL( cell.z(), 567 )

   cell.x() = cell_idx_c(123);
   forest.mapToPeriodicDomain( cell, 1 );

   WALBERLA_CHECK_EQUAL( cell.x(), 123 )
   WALBERLA_CHECK_EQUAL( cell.y(), 234 )
   WALBERLA_CHECK_EQUAL( cell.z(), 567 )

   cell.x() = cell_idx_c(223);
   forest.mapToPeriodicDomain( cell, 1 );

   WALBERLA_CHECK_EQUAL( cell.x(),  23 )
   WALBERLA_CHECK_EQUAL( cell.y(), 234 )
   WALBERLA_CHECK_EQUAL( cell.z(), 567 )

   // getCell

   forest.getCell( cell, real_c(23.5), real_c(49.9), real_c(50.1) );

   WALBERLA_CHECK_EQUAL( cell.x(), 23 )
   WALBERLA_CHECK_EQUAL( cell.y(), 59 )
   WALBERLA_CHECK_EQUAL( cell.z(), 70 )

   forest.getCell( cell, real_c(-23.5), real_c(49.9), real_c(50.1) );

   WALBERLA_CHECK_EQUAL( cell.x(), -24 )
   WALBERLA_CHECK_EQUAL( cell.y(),  59 )
   WALBERLA_CHECK_EQUAL( cell.z(),  70 )

   forest.getCell( cell, real_c(123.5), real_c(49.9), real_c(50.1) );

   WALBERLA_CHECK_EQUAL( cell.x(), 123 )
   WALBERLA_CHECK_EQUAL( cell.y(),  59 )
   WALBERLA_CHECK_EQUAL( cell.z(),  70 )

   forest.getCell( cell, real_c(23.4), real_c(49.9), real_c(50.1), 1 );

   WALBERLA_CHECK_EQUAL( cell.x(),  46 )
   WALBERLA_CHECK_EQUAL( cell.y(), 119 )
   WALBERLA_CHECK_EQUAL( cell.z(), 140 )

   forest.getCell( cell, real_c(-23.4), real_c(49.9), real_c(50.1), 1 );

   WALBERLA_CHECK_EQUAL( cell.x(), -47 )
   WALBERLA_CHECK_EQUAL( cell.y(), 119 )
   WALBERLA_CHECK_EQUAL( cell.z(), 140 )

   forest.getCell( cell, real_c(123.4), real_c(49.9), real_c(50.1), 1 );

   WALBERLA_CHECK_EQUAL( cell.x(), 246 )
   WALBERLA_CHECK_EQUAL( cell.y(), 119 )
   WALBERLA_CHECK_EQUAL( cell.z(), 140 )

   forest.getCell( cell, real_c(-1.1), real_c(0), real_c(0) );
   WALBERLA_CHECK_EQUAL( cell.x(), -2 )

   forest.getCell( cell, real_c(-1.0), real_c(0), real_c(0) );
   WALBERLA_CHECK_EQUAL( cell.x(), -1 )

   forest.getCell( cell, real_c(-0.5), real_c(0), real_c(0) );
   WALBERLA_CHECK_EQUAL( cell.x(), -1 )

   forest.getCell( cell, real_c(-0.0), real_c(0), real_c(0) );
   WALBERLA_CHECK_EQUAL( cell.x(),  0 )

   forest.getCell( cell,  real_c(0.0), real_c(0), real_c(0) );
   WALBERLA_CHECK_EQUAL( cell.x(),  0 )

   forest.getCell( cell,  real_c(0.5), real_c(0), real_c(0) );
   WALBERLA_CHECK_EQUAL( cell.x(),  0 )

   forest.getCell( cell,  real_c(1.0), real_c(0), real_c(0) );
   WALBERLA_CHECK_EQUAL( cell.x(),  1 )

   forest.getCell( cell,  real_c(1.1), real_c(0), real_c(0) );
   WALBERLA_CHECK_EQUAL( cell.x(),  1 )

   // getCellCenter

   real_t x, y, z;
   forest.getCellCenter( x, y, z, Cell(23,59,70) );

   WALBERLA_CHECK_FLOAT_EQUAL( x, 23.5 )
   WALBERLA_CHECK_FLOAT_EQUAL( y, ( 59 + 0.5 ) * 100.0 / 120.0 )
   WALBERLA_CHECK_FLOAT_EQUAL( z, ( 70 + 0.5 ) * 100.0 / 140.0 )

   forest.getCellCenter( x, y, z, Cell(-24,49,72) );

   WALBERLA_CHECK_FLOAT_EQUAL( x, -23.5 )
   WALBERLA_CHECK_FLOAT_EQUAL( y, ( 49 + 0.5 ) * 100.0 / 120.0 )
   WALBERLA_CHECK_FLOAT_EQUAL( z, ( 72 + 0.5 ) * 100.0 / 140.0 )

   forest.getCellCenter( x, y, z, Cell(123,47,62) );

   WALBERLA_CHECK_FLOAT_EQUAL( x, 123.5 )
   WALBERLA_CHECK_FLOAT_EQUAL( y, ( 47 + 0.5 ) * 100.0 / 120.0 )
   WALBERLA_CHECK_FLOAT_EQUAL( z, ( 62 + 0.5 ) * 100.0 / 140.0 )

   forest.getCellCenter( x, y, z, Cell(46,119,140), 1 );

   WALBERLA_CHECK_FLOAT_EQUAL( x, 23.25 )
   WALBERLA_CHECK_FLOAT_EQUAL( y, ( 119 + 0.5 ) * 100.0 / 240.0 )
   WALBERLA_CHECK_FLOAT_EQUAL( z, ( 140 + 0.5 ) * 100.0 / 280.0 )

   forest.getCellCenter( x, y, z, Cell(-47,109,142), 1 );

   WALBERLA_CHECK_FLOAT_EQUAL( x, -23.25 )
   WALBERLA_CHECK_FLOAT_EQUAL( y, ( 109 + 0.5 ) * 100.0 / 240.0 )
   WALBERLA_CHECK_FLOAT_EQUAL( z, ( 142 + 0.5 ) * 100.0 / 280.0 )

   forest.getCellCenter( x, y, z, Cell(246,19,1400), 1 );

   WALBERLA_CHECK_FLOAT_EQUAL( x, 123.25 )
   WALBERLA_CHECK_FLOAT_EQUAL( y, ( 19   + 0.5 ) * 100.0 / 240.0 )
   WALBERLA_CHECK_FLOAT_EQUAL( z, ( 1400 + 0.5 ) * 100.0 / 280.0 )

   // getCellBBFromAABB & getCellBBFromCellAlignedAABB

   forest.getCellBBFromAABB( bb, AABB(real_c(2.0),real_c(0.0),real_c(0.0),real_c(23.0),real_c(25.0),real_c(50.0)) );

   WALBERLA_CHECK_EQUAL( bb, CellInterval(Cell(2,0,0),Cell(22,29,69)) )

   forest.getCellBBFromCellAlignedAABB( bb, AABB(real_c(2.0),real_c(0.0),real_c(0.0),real_c(23.0),real_c(25.0),real_c(50.0)) );

   WALBERLA_CHECK_EQUAL( bb, CellInterval(Cell(2,0,0),Cell(22,29,69)) )

   WALBERLA_CHECK( forest.isCellAlignedAABB( AABB(real_c(2.0),real_c(0.0),real_c(0.0),real_c(23.0),real_c(25.0),real_c(50.0)) ) )

   forest.getCellBBFromAABB( bb, AABB(real_c(2.3),real_c(0.12),real_c(0.42),real_c(23.1),real_c(25.1),real_c(50.1)) );

   WALBERLA_CHECK_EQUAL( bb, CellInterval(Cell(2,0,0),Cell(23,30,70)) )

   WALBERLA_CHECK( !forest.isCellAlignedAABB( AABB(real_c(2.3),real_c(0.12),real_c(0.42),real_c(23.1),real_c(25.1),real_c(50.1)) ) )

   forest.getCellBBFromAABB( bb, AABB(real_c(2.5),real_c(0.0),real_c(0.0),real_c(23.5),real_c(25.0),real_c(50.0)), 1 );

   WALBERLA_CHECK_EQUAL( bb, CellInterval(Cell(5,0,0),Cell(46,59,139)) )

   forest.getCellBBFromCellAlignedAABB( bb, AABB(real_c(2.5),real_c(0.0),real_c(0.0),real_c(23.5),real_c(25.0),real_c(50.0)), 1 );

   WALBERLA_CHECK_EQUAL( bb, CellInterval(Cell(5,0,0),Cell(46,59,139)) )

   WALBERLA_CHECK( forest.isCellAlignedAABB( AABB(real_c(2.5),real_c(0.0),real_c(0.0),real_c(23.5),real_c(25.0),real_c(50.0)), 1 ) )
   WALBERLA_CHECK( !forest.isCellAlignedAABB( AABB(real_c(2.5),real_c(0.0),real_c(0.0),real_c(23.5),real_c(25.0),real_c(50.0)) ) )

   forest.getCellBBFromAABB( bb, AABB(real_c(2.7),real_c(0.12),real_c(0.22),real_c(23.6),real_c(25.1),real_c(50.1)), 1 );

   WALBERLA_CHECK_EQUAL( bb, CellInterval(Cell(5,0,0),Cell(47,60,140)) )

   WALBERLA_CHECK( !forest.isCellAlignedAABB( AABB(real_c(2.7),real_c(0.12),real_c(0.22),real_c(23.6),real_c(25.1),real_c(50.1)), 1 ) )

   // getAABBFromCellBB

   AABB aabb;

   forest.getAABBFromCellBB( aabb, CellInterval(Cell(3,0,0),Cell(46,29,139)) );

   WALBERLA_CHECK_EQUAL( aabb, AABB(real_c(3),real_c(0),real_c(0),real_c(47),real_c(25),real_c(100)) )

   forest.getAABBFromCellBB( aabb, CellInterval(Cell(3,0,0),Cell(46,59,279)), 1 );

   WALBERLA_CHECK_EQUAL( aabb, AABB(real_c(1.5),real_c(0),real_c(0),real_c(23.5),real_c(25),real_c(100)) )

   // getBlock & getBlockCellBB

   const IBlock* block = forest.getBlock( Cell(-1,23,42) );
   WALBERLA_CHECK_NULLPTR( block )

   block = forest.getBlock( Cell(1,23,42) );
   WALBERLA_CHECK_NULLPTR( block )

   block = forest.getBlock( Cell(50,23,42) );
   WALBERLA_CHECK_NOT_NULLPTR( block )
   bb = forest.getBlockCellBB( *block );
   WALBERLA_CHECK_EQUAL( bb, CellInterval(Cell(50,0,0),Cell(99,59,69)) )

   block = forest.getBlock( Cell(49,63,42) );
   WALBERLA_CHECK_NOT_NULLPTR( block )
   bb = forest.getBlockCellBB( *block );
   WALBERLA_CHECK_EQUAL( bb, CellInterval(Cell(0,60,0),Cell(49,119,69)) )

   block = forest.getBlock( Cell(52,63,42) );
   WALBERLA_CHECK_NOT_NULLPTR( block )
   bb = forest.getBlockCellBB( *block );
   WALBERLA_CHECK_EQUAL( bb, CellInterval(Cell(50,60,0),Cell(99,119,69)) )

   block = forest.getBlock( Cell(1,23,72) );
   WALBERLA_CHECK_NOT_NULLPTR( block )
   bb = forest.getBlockCellBB( *block );
   WALBERLA_CHECK_EQUAL( bb, CellInterval(Cell(0,0,70),Cell(49,59,139)) )

   block = forest.getBlock( Cell(51,63,72) );
   WALBERLA_CHECK_NULLPTR( block )

   ConstBlockDataID bid = forest.getBlockCellBBId();

   block = forest.getBlock( Cell(23,42,240), 1 );
   WALBERLA_CHECK_NULLPTR( block )

   forest.getCellCenter( x, y, z, Cell(23,42,240), 1 );
   block = forest.getBlock( x, y, z );
   WALBERLA_CHECK_NOT_NULLPTR( block )

   block = forest.getBlock( Cell(100,179,200), 1 );
   WALBERLA_CHECK_NOT_NULLPTR( block )
   bb = *( block->getData< CellInterval >( bid ) );
   WALBERLA_CHECK_EQUAL( bb, CellInterval(Cell(100,120,140),Cell(149,179,209)) )

   block = forest.getBlock( Cell(155,230,279), 1 );
   WALBERLA_CHECK_NOT_NULLPTR( block )
   bb = *( block->getData< CellInterval >( bid ) );
   WALBERLA_CHECK_EQUAL( bb, CellInterval(Cell(150,180,210),Cell(199,239,279)) )

   std::mt19937 eng( 23 );
   for( int i = 0; i < 23; ++i )
   {
      const Vector3<real_t> globalPoint = forest.getDomain().getScaled( real_t(1.25) ).randomPoint( eng );

      for( auto blockIt = forest.begin(); blockIt != forest.end(); ++blockIt )
      {
          Vector3<real_t> localPoint0, localPoint1;

          forest.transformGlobalToBlockLocal( localPoint0, *blockIt, globalPoint );

          localPoint1 = globalPoint;
          forest.transformGlobalToBlockLocal( localPoint1, *blockIt );

          WALBERLA_CHECK_EQUAL( localPoint0, localPoint1 )

          Vector3<real_t> globalPoint0, globalPoint1;

          forest.transformBlockLocalToGlobal( globalPoint0, *blockIt, localPoint0 );

          globalPoint1 = localPoint0;
          forest.transformBlockLocalToGlobal( globalPoint1, *blockIt );

          WALBERLA_CHECK_FLOAT_EQUAL( globalPoint[0], globalPoint0[0] )
          WALBERLA_CHECK_FLOAT_EQUAL( globalPoint[1], globalPoint0[1] )
          WALBERLA_CHECK_FLOAT_EQUAL( globalPoint[2], globalPoint0[2] )
          WALBERLA_CHECK_FLOAT_EQUAL( globalPoint[0], globalPoint1[0] )
          WALBERLA_CHECK_FLOAT_EQUAL( globalPoint[1], globalPoint1[1] )
          WALBERLA_CHECK_FLOAT_EQUAL( globalPoint[2], globalPoint1[2] )
      }
   }
}



} // namespace blockforest
} // namespace walberla



int main( int argc, char** argv )
{
   walberla::debug::enterTestMode();

   walberla::mpi::Environment env( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

   walberla::blockforest::test();

   return EXIT_SUCCESS;
}
