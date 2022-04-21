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
//! \file BlockForestTest.cpp
//! \ingroup blockforest
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#include "blockforest/BlockForest.h"
#include "blockforest/SetupBlockForest.h"
#include "core/debug/TestSubsystem.h"
#include "core/math/Random.h"
#include "core/mpi/Environment.h"
#include "core/mpi/MPIManager.h"

#include <vector>


namespace walberla {
namespace blockforest {



class RootBlockExclusion { // excludes root block (0,0,0)

public:

   RootBlockExclusion( const AABB& domain, const uint_t xSize, const uint_t ySize, const uint_t zSize ) :
      x_( domain.xMin() + ( domain.xSize() / real_c( xSize ) ) / real_c(2) ),
      y_( domain.yMin() + ( domain.ySize() / real_c( ySize ) ) / real_c(2) ),
      z_( domain.zMin() + ( domain.zSize() / real_c( zSize ) ) / real_c(2) ) {}

   void operator()( std::vector<uint8_t>& excludeBlock, const SetupBlockForest::RootBlockAABB& aabb )
   {
      for( uint_t i = 0; i != excludeBlock.size(); ++i )
      {
         AABB bb = aabb(i);
         if( bb.contains( x_, y_, z_ ) )
            excludeBlock[i] = 1;
      }
   }

private:

   const real_t x_;
   const real_t y_;
   const real_t z_;
};



static void refinementSelectionFunctionRandom( SetupBlockForest& forest ) {

   std::vector< SetupBlock* > blocks;
   forest.getBlocks( blocks, forest.getDepth() );

   const uint_t max = blocks.size() >> 3;

   for( uint_t i = 0; i != max; ++i ) {
      SetupBlock* const block = blocks[ math::intRandom( uint_t(0), uint_c( blocks.size()-1 ) ) ];
      if( block->getLevel() < 4 ) block->setMarker( true );
   }
}



static void workloadMemorySUIDAssignmentFunction( SetupBlockForest& forest ) {

   std::vector< SetupBlock* > blocks;
   forest.getBlocks( blocks );

   std::vector< SUID > suids;
   for( uint_t i = 0; i != forest.getNumberOfLevels(); ++i ) {
      std::ostringstream oss;
      oss << "Level_" << i;
      suids.emplace_back( oss.str(), false );
   }

   for( uint_t i = 0; i != blocks.size(); ++i ) {
      SetupBlock* block = blocks[i];
      block->setMemory( numeric_cast< memory_t >( 1.0 ) );
      block->setWorkload( numeric_cast< workload_t >( std::pow( 2.0, numeric_cast<double>(block->getLevel()) ) ) );
      if( block->getLevel() > 0 )
         block->addState( suids[ block->getLevel() ] );
   }
}



static memory_t communicationCalculationFunction( const SetupBlock* const a, const SetupBlock* const b ) {

   uint_t faces[] = { 4, 10, 12, 13, 15, 21 };

   for( uint_t i = 0; i != 6; ++i ) {
      for( uint_t j = 0; j != a->getNeighborhoodSectionSize(faces[i]); ++j )
         if( a->getNeighbor(faces[i],j) == b )
            return 1000.0;
   }

   return 1.0;
}



static uint_t* blockdata23( const IBlock* const /*block*/ ) {

   return new uint_t(23);
}



static int* blockdata42( const IBlock* const /*block*/ ) {

   return new int(42);
}



static uint_t* blockdata5( const IBlock* const /*block*/ ) {

   return new uint_t(5);
}



class Base
{
public:
   virtual ~Base() = default;
   bool operator==( const Base& /*rhs*/ ) const { return true; }
           uint_t override() const { return 1; }
   virtual uint_t func()     const { return 2; }
};

class Derived : public Base
{
public:
   uint_t override() const { return 10; }
   uint_t func()     const override { return 20; }
};



static Base* blockdataBase( const IBlock* const /*block*/ )
{
   return new Base();
}

static Derived* blockdataDerived( const IBlock* const /*block*/ )
{
   return new Derived();
}



class SecondBase
{
public:
   virtual ~SecondBase() = default;
   bool operator==( const SecondBase& /*rhs*/ ) const { return true; }
   uint_t override() const { return 100; }
};

class Multi : public Base, public SecondBase
{
public:
   bool operator==( const Multi& /*rhs*/ ) const { return true; }
   uint_t override() const { return 1000; }
   uint_t func()     const override { return 2000; }
};



static Multi* blockdataMulti( const IBlock* const /*block*/ )
{
   return new Multi();
}



static void test() {

   SetupBlockForest sforest;

   sforest.addRefinementSelectionFunction( refinementSelectionFunctionRandom );
   sforest.addWorkloadMemorySUIDAssignmentFunction( workloadMemorySUIDAssignmentFunction );

   real_t xmin = math::realRandom( real_c(-100), real_c(100) );
   real_t xmax = math::realRandom( xmin + real_c(10), real_c(120) );
   real_t ymin = math::realRandom( real_c(-100), real_c(100) );
   real_t ymax = math::realRandom( ymin + real_c(10), real_c(120) );
   real_t zmin = math::realRandom( real_c(-100), real_c(100) );
   real_t zmax = math::realRandom( zmin + real_c(10), real_c(120) );

   AABB domain( xmin, ymin, zmin, xmax, ymax, zmax );
#ifdef NDEBUG
   const uint_t xSize = math::intRandom( uint_t(4), uint_t(8) );
   const uint_t ySize = math::intRandom( uint_t(4), uint_t(8) );
   const uint_t zSize = math::intRandom( uint_t(4), uint_t(8) );
#else
   const uint_t xSize = uint_c(3);
   const uint_t ySize = uint_c(3);
   const uint_t zSize = uint_c(3);
#endif

   sforest.addRootBlockExclusionFunction( RootBlockExclusion( domain, xSize, ySize, zSize ) );

   sforest.init( domain, xSize, ySize, zSize, math::boolRandom(), math::boolRandom(), math::boolRandom() );

   //sforest.writeVTKOutput( "SetupBlockForest_0" );

   std::vector< SetupBlock* > blocks;
   sforest.getBlocks( blocks );

   const uint_t numberOfBlocks = blocks.size();

   memory_t maxMemory = 0.0;
   for( uint_t i = 0; i != numberOfBlocks; ++i )
      maxMemory += blocks[i]->getMemory();

   GlobalLoadBalancing::MetisConfiguration< SetupBlock > metisConfig( true, false, communicationCalculationFunction );

   sforest.calculateProcessDistribution_Default( uint_c(4), maxMemory, "hilbert", 10, false, metisConfig );

   //sforest.writeVTKOutput( "SetupBlockForest_1" );

   WALBERLA_ROOT_SECTION()
   {
      sforest.saveToFile( "blockforest.sav" );

      // endian independent serialization

#ifdef NDEBUG
      const uint_t size = 100000;
#else
      const uint_t size = 1000;
#endif
      std::vector< real_t  > send( size,  real_c(0) );
      std::vector< uint8_t > buffer( size * ( sizeof(real_t) + 3 ) );

      realToByteArray( real_c(0), buffer, 0 );
      for( uint_t i = 1; i != size; ++i )
      {
         send[i] = math::realRandom( real_c(-1e23), real_c(1e23) );
         realToByteArray( send[i], buffer, i * ( sizeof(real_t) + 3 ) );
      }

      for( uint_t i = 0; i != size; ++i )
      {
         const real_t recv = byteArrayToReal< real_t >( buffer, i * ( sizeof(real_t) + 3 ) );
         WALBERLA_CHECK( realIsIdentical( recv, send[i] ) )
      }
   }

   WALBERLA_MPI_WORLD_BARRIER()

   BlockForest bforest( uint_c( MPIManager::instance()->rank() ), sforest, true );

   BlockForest fforest_1( uint_c( MPIManager::instance()->rank() ), "blockforest.sav", false, true ); // all processes read the same file
   WALBERLA_CHECK_EQUAL( bforest, fforest_1 )

   BlockForest fforest_2( uint_c( MPIManager::instance()->rank() ), "blockforest.sav", true, true ); // only the root process reads the file
   WALBERLA_CHECK_EQUAL( bforest, fforest_2 )

   SUID level1( "Level_1" );
   SUID level2( "Level_2" );
   SUID level3( "Level_3" );
   SUID level4( "Level_4" );

   BlockDataID data1 = bforest.addBlockData( "blockdata23_42" ) << BlockDataCreator< uint_t >( blockdata23, "blockdata23", level3 )
                                                                << BlockDataCreator<    int >( blockdata42, "blockdata42", level4 );

   BlockDataID data2 = bforest.addBlockData< uint_t >( blockdata5, "blockdata5", Set<SUID>::emptySet(), level2 );

   ConstBlockDataID data3 = bforest.addBlockData< Base >( blockdataBase, "blockdataBase" );
   ConstBlockDataID data4 = bforest.addBlockData< Derived >( blockdataDerived, "blockdataDerived" );

   BlockDataID data5 = bforest.addBlockData< Multi >( blockdataMulti, "blockdataMulti" );

   for( BlockStorage::const_iterator it = bforest.begin(); it != bforest.end(); ++it ) {
      const Block* block = static_cast< const Block* >( it.get() );

      if( block->getLevel() == 0 ) {

         WALBERLA_CHECK( block->getState().empty() )

         WALBERLA_CHECK_NULLPTR( block->getData< uint_t >( data1 ) )
         WALBERLA_CHECK_NOT_NULLPTR( block->getData< uint_t >( data2 ) )
         WALBERLA_CHECK_EQUAL( *(block->getData< uint_t >( data2 )), uint_c(5) )

         WALBERLA_CHECK( !block->isBlockDataAllocated( data1 ) )
         WALBERLA_CHECK( block->isBlockDataAllocated( data2 ) )
         WALBERLA_CHECK( block->isBlockDataAllocated( data3 ) )
         WALBERLA_CHECK( block->isBlockDataAllocated( data4 ) )
         WALBERLA_CHECK( block->isBlockDataAllocated( data5 ) )

         WALBERLA_CHECK( block->isDataOfType<uint_t>( data2 ) )
         WALBERLA_CHECK( !block->isDataOfType<int>( data2 ) )
         WALBERLA_CHECK( block->isDataClassOrSubclassOf<uint_t>( data2 ) )
         WALBERLA_CHECK( !block->isDataClassOrSubclassOf<BlockDataID>( data2 ) )
      }
      else if( block->getLevel() == 1 ) {

         WALBERLA_CHECK_EQUAL( block->getState().size(), uint_c(1) )

         WALBERLA_CHECK_NULLPTR( block->getData< uint_t >( data1 ) )
         WALBERLA_CHECK_NOT_NULLPTR( block->getData< uint_t >( data2 ) )
         WALBERLA_CHECK_EQUAL( *(block->getData< uint_t >( data2 )), uint_c(5) )

         WALBERLA_CHECK( !block->isBlockDataAllocated( data1 ) )
         WALBERLA_CHECK( block->isBlockDataAllocated( data2 ) )
         WALBERLA_CHECK( block->isBlockDataAllocated( data3 ) )
         WALBERLA_CHECK( block->isBlockDataAllocated( data4 ) )
         WALBERLA_CHECK( block->isBlockDataAllocated( data5 ) )

         WALBERLA_CHECK( block->isDataOfType<uint_t>( data2 ) )
         WALBERLA_CHECK( !block->isDataOfType<SUID>( data2 ) )
         WALBERLA_CHECK( block->isDataClassOrSubclassOf<uint_t>( data2 ) )
         WALBERLA_CHECK( !block->isDataClassOrSubclassOf<double>( data2 ) )
      }
      else if( block->getLevel() == 2 ) {

         WALBERLA_CHECK_EQUAL( block->getState().size(), uint_c(1) )

         WALBERLA_CHECK_NULLPTR( block->getData< uint_t >( data1 ) )
         WALBERLA_CHECK_NULLPTR( block->getData< uint_t >( data2 ) )

         WALBERLA_CHECK( !block->isBlockDataAllocated( data1 ) )
         WALBERLA_CHECK( !block->isBlockDataAllocated( data2 ) )
         WALBERLA_CHECK( block->isBlockDataAllocated( data3 ) )
         WALBERLA_CHECK( block->isBlockDataAllocated( data4 ) )
         WALBERLA_CHECK( block->isBlockDataAllocated( data5 ) )
      }
      else if( block->getLevel() == 3 ) {

         WALBERLA_CHECK_EQUAL( block->getState().size(), uint_c(1) )

         WALBERLA_CHECK_NOT_NULLPTR( block->getData< uint_t >( data1 ) )
         WALBERLA_CHECK_EQUAL( *(block->getData< uint_t >( data1 )), uint_c(23) )
         WALBERLA_CHECK_NOT_NULLPTR( block->getData< uint_t >( data2 ) )
         WALBERLA_CHECK_EQUAL( *(block->getData< uint_t >( data2 )), uint_c( 5) )

         WALBERLA_CHECK( block->isBlockDataAllocated( data1 ) )
         WALBERLA_CHECK( block->isBlockDataAllocated( data2 ) )
         WALBERLA_CHECK( block->isBlockDataAllocated( data3 ) )
         WALBERLA_CHECK( block->isBlockDataAllocated( data4 ) )
         WALBERLA_CHECK( block->isBlockDataAllocated( data5 ) )

         WALBERLA_CHECK( block->isDataOfSameType( data1, data2 ) )
      }
      else if( block->getLevel() == 4 ) {

         WALBERLA_CHECK_EQUAL( block->getState().size(), uint_c(1) )

         WALBERLA_CHECK_NOT_NULLPTR( block->getData< int >( data1 ) )
         WALBERLA_CHECK_EQUAL( *(block->getData< int >( data1 )), uint_c(42) )
         WALBERLA_CHECK_NOT_NULLPTR( block->getData< uint_t >( data2 ) )
         WALBERLA_CHECK_EQUAL( *(block->getData< uint_t >( data2 )), uint_c( 5) )

         WALBERLA_CHECK( block->isBlockDataAllocated( data1 ) )
         WALBERLA_CHECK( block->isBlockDataAllocated( data2 ) )
         WALBERLA_CHECK( block->isBlockDataAllocated( data3 ) )
         WALBERLA_CHECK( block->isBlockDataAllocated( data4 ) )
         WALBERLA_CHECK( block->isBlockDataAllocated( data5 ) )

         WALBERLA_CHECK( !block->isDataOfSameType( data1, data2 ) )

         const Base* base = block->getData< Base >( data3 );

         WALBERLA_CHECK_EQUAL( base->override(), 1 )
         WALBERLA_CHECK_EQUAL( base->func(), 2 )

         WALBERLA_CHECK( block->isDataOfType<Base>( data3 ) )
         WALBERLA_CHECK( !block->isDataOfType<Derived>( data3 ) )

         WALBERLA_CHECK( block->isDataClassOrSubclassOf<Base>( data3 ) )
         WALBERLA_CHECK( !block->isDataClassOrSubclassOf<Derived>( data3 ) )
         WALBERLA_CHECK( !block->isDataClassOrSubclassOf<int>( data3 ) )

         WALBERLA_CHECK( !block->isDataSubclassOf<Base>( data3 ) )
         WALBERLA_CHECK( !block->isDataSubclassOf<Derived>( data3 ) )
         WALBERLA_CHECK( !block->isDataSubclassOf<int>( data3 ) )

                        base    = block->getData< Base >( data4 );
         const Derived* derived = block->getData< Derived >( data4 );

         WALBERLA_CHECK_EQUAL( base->override(), 1 )
         WALBERLA_CHECK_EQUAL( base->func(), 20 )

         WALBERLA_CHECK_EQUAL( derived->override(), 10 )
         WALBERLA_CHECK_EQUAL( derived->func(), 20 )

         WALBERLA_CHECK( !block->isDataOfType<Base>( data4 ) )
         WALBERLA_CHECK( block->isDataOfType<Derived>( data4 ) )

         WALBERLA_CHECK( block->isDataClassOrSubclassOf<Base>( data4 ) )
         WALBERLA_CHECK( block->isDataClassOrSubclassOf<Derived>( data4 ) )
         WALBERLA_CHECK( !block->isDataClassOrSubclassOf<int>( data4 ) )

         WALBERLA_CHECK( block->isDataSubclassOf<Base>( data4 ) )
         WALBERLA_CHECK( !block->isDataSubclassOf<Derived>( data4 ) )
         WALBERLA_CHECK( !block->isDataSubclassOf<int>( data4 ) )

                           base  = block->getData< Base >( data5 );
         const SecondBase* sBase = block->getData< SecondBase >( data5 );
         const Multi*      multi = block->getData< Multi >( data5 );

         WALBERLA_CHECK( reinterpret_cast< const void* >( base ) != reinterpret_cast< const void* >( sBase ) )
         WALBERLA_CHECK( ( reinterpret_cast< const void* >(  base ) != reinterpret_cast< const void* >( multi ) ) ||
                ( reinterpret_cast< const void* >( sBase ) != reinterpret_cast< const void* >( multi ) ) )

         WALBERLA_CHECK_EQUAL( base->override(), 1 )
         WALBERLA_CHECK_EQUAL( base->func(), 2000 )

         WALBERLA_CHECK_EQUAL( sBase->override(), 100 )

         WALBERLA_CHECK_EQUAL( multi->override(), 1000 )
         WALBERLA_CHECK_EQUAL( multi->func(), 2000 )

         WALBERLA_CHECK( !block->isDataOfType<Base>( data5 ) )
         WALBERLA_CHECK( !block->isDataOfType<SecondBase>( data5 ) )
         WALBERLA_CHECK( block->isDataOfType<Multi>( data5 ) )

         WALBERLA_CHECK( block->isDataClassOrSubclassOf<Base>( data5 ) )
         WALBERLA_CHECK( block->isDataClassOrSubclassOf<SecondBase>( data5 ) )
         WALBERLA_CHECK( block->isDataClassOrSubclassOf<Multi>( data5 ) )
         WALBERLA_CHECK( !block->isDataClassOrSubclassOf<Derived>( data5 ) )

         WALBERLA_CHECK( block->isDataSubclassOf<Base>( data5 ) )
         WALBERLA_CHECK( block->isDataSubclassOf<SecondBase>( data5 ) )
         WALBERLA_CHECK( !block->isDataSubclassOf<Multi>( data5 ) )
         WALBERLA_CHECK( !block->isDataSubclassOf<Derived>( data5 ) )
      }
   }

   std::vector< shared_ptr< IBlockID > > blockIds;
   bforest.getAllBlocks( blockIds );

   std::set< BlockID > setupBlockSet;
   std::set< BlockID > forestBlockSet;

   for( auto block = sforest.begin(); block != sforest.end(); ++block ) {
      WALBERLA_CHECK_EQUAL( setupBlockSet.find( block->getId() ), setupBlockSet.end() )
      setupBlockSet.insert( block->getId() );
   }

   for( auto blockId = blockIds.begin(); blockId != blockIds.end(); ++blockId ) {
      const BlockID id = *dynamic_cast< BlockID* >( blockId->get() );
      WALBERLA_CHECK_EQUAL( forestBlockSet.find( id ), forestBlockSet.end() )
      forestBlockSet.insert( id );
   }

   WALBERLA_CHECK_EQUAL( setupBlockSet, forestBlockSet )
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
