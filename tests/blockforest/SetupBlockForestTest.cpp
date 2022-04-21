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
//! \file SetupBlockForestTest.cpp
//! \ingroup blockforest
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#include "blockforest/SetupBlockForest.h"
#include "core/debug/Debug.h"
#include "core/debug/TestSubsystem.h"
#include "core/math/Random.h"
#include "core/mpi/Environment.h"

#include <set>
#include <vector>


namespace walberla {
namespace blockforest {



static void refinementSelectionFunctionAll( SetupBlockForest& forest ) {

   std::vector< SetupBlock* > blocks;
   forest.getBlocks( blocks );

   for( uint_t i = 0; i != blocks.size(); ++i )
      if( blocks[i]->getLevel() < 2 ) blocks[i]->setMarker( true );
}



static void refinementSelectionFunctionRandom( SetupBlockForest& forest ) {

   // 20x20x20 = 8000 = 1 1111 0100 0000 = 13

#ifndef WALBERLA_BLOCKFOREST_PRIMITIVE_BLOCKID
   const uint_t maxLevel = 25;
#else
   WALBERLA_ASSERT_GREATER_EQUAL( math::UINT_BITS, 32 )
   const uint_t maxLevel = ( math::UINT_BITS - 14 ) / 3;
#endif

   std::vector< SetupBlock* > blocks;
   forest.getBlocks( blocks, forest.getDepth() );

   const uint_t max = blocks.size() >> 3;

   for( uint_t i = 0; i != max; ++i ) {
      SetupBlock* const block = blocks[ math::intRandom( uint_t(0), uint_c( blocks.size()-1 ) ) ];
      if( block->getLevel() < maxLevel ) block->setMarker( true );
   }
}



static void checkCollectorConsistency( SetupBlockForest& forest ) {

   std::vector< SetupBlock* > blocks;
   forest.getBlocks( blocks );

   std::vector< SetupBlock* > itBlocks;
   for( SetupBlockForest::iterator it = forest.begin(); it != forest.end(); ++it )
      itBlocks.push_back( it.get() );

   WALBERLA_CHECK_EQUAL( blocks.size(), itBlocks.size() )
   for( uint_t i = 0; i != blocks.size(); ++i )
      WALBERLA_CHECK( blocks[i] == itBlocks[i] )

   std::set< SetupBlock* > baseSet;
   for( uint_t i = 0; i != blocks.size(); ++i )
      baseSet.insert( blocks[i] );
   WALBERLA_CHECK_EQUAL( baseSet.size(), blocks.size() )

   std::vector< SetupBlock* > hilbertBlocks;
   forest.getHilbertOrder( hilbertBlocks );
   WALBERLA_CHECK_EQUAL( hilbertBlocks.size(), blocks.size() )

   std::set< SetupBlock* > hilbertSet;
   for( uint_t i = 0; i != hilbertBlocks.size(); ++i )
      hilbertSet.insert( hilbertBlocks[i] );
   WALBERLA_CHECK_EQUAL( hilbertSet.size(), hilbertBlocks.size() )

   std::set< SetupBlock* >::iterator baseIterator    = baseSet.begin();
   std::set< SetupBlock* >::iterator hilbertIterator = hilbertSet.begin();

   while( baseIterator != baseSet.end() ) {
      WALBERLA_CHECK( *baseIterator == *hilbertIterator )
      ++baseIterator;
      ++hilbertIterator;
   }

   std::vector< SetupBlock* > aabbBlocks;
   forest.getBlocksOverlappedByAABB( aabbBlocks, forest.getDomain() );
   WALBERLA_CHECK_EQUAL( aabbBlocks.size(), blocks.size() )

   std::set< SetupBlock* > aabbSet;
   for( uint_t i = 0; i != aabbBlocks.size(); ++i )
      aabbSet.insert( aabbBlocks[i] );
   WALBERLA_CHECK_EQUAL( aabbSet.size(), aabbBlocks.size() )

                                     baseIterator = baseSet.begin();
   std::set< SetupBlock* >::iterator aabbIterator = aabbSet.begin();

   while( baseIterator != baseSet.end() ) {
      WALBERLA_CHECK( *baseIterator == *aabbIterator )
      ++baseIterator;
      ++aabbIterator;
   }
}



static void checkNeighborhoodConsistency( const SetupBlockForest& forest ) {

   std::vector< const SetupBlock* > blocks;
   forest.getBlocks( blocks );

   const int blockssize = int_c( blocks.size() );
#ifdef _OPENMP
   #pragma omp parallel for schedule(static)
#endif
   for( int i = 0; i < blockssize; ++i ) {

      const SetupBlock* const block = blocks[uint_c(i)];

      std::vector< real_t > neighborhoodSectionBlockCenters;

      for( uint_t n = 0; n != 26; ++n ) {

         std::vector< bool > hit( block->getNeighborhoodSectionSize(n), false );

         constructNeighborhoodSectionBlockCenters( n, block->getAABB(), neighborhoodSectionBlockCenters );

         WALBERLA_CHECK_EQUAL( neighborhoodSectionBlockCenters.size() % 3, uint_c(0) )

         for( uint_t p = 0; p != neighborhoodSectionBlockCenters.size(); p += 3 ) {

            real_t x = neighborhoodSectionBlockCenters[p];
            real_t y = neighborhoodSectionBlockCenters[p+1];
            real_t z = neighborhoodSectionBlockCenters[p+2];

            // treat periodicity
            if( x <  forest.getDomain().xMin() && forest.isXPeriodic() ) x = forest.getDomain().xMax() - forest.getDomain().xMin() + x;
            if( x >= forest.getDomain().xMax() && forest.isXPeriodic() ) x = forest.getDomain().xMin() - forest.getDomain().xMax() + x;
            if( y <  forest.getDomain().yMin() && forest.isYPeriodic() ) y = forest.getDomain().yMax() - forest.getDomain().yMin() + y;
            if( y >= forest.getDomain().yMax() && forest.isYPeriodic() ) y = forest.getDomain().yMin() - forest.getDomain().yMax() + y;
            if( z <  forest.getDomain().zMin() && forest.isZPeriodic() ) z = forest.getDomain().zMax() - forest.getDomain().zMin() + z;
            if( z >= forest.getDomain().zMax() && forest.isZPeriodic() ) z = forest.getDomain().zMin() - forest.getDomain().zMax() + z;

            bool noHit = true;
            for( uint_t c = 0; c != block->getNeighborhoodSectionSize(n) && noHit; ++c ) {
               if( block->getNeighbor(n,c)->getAABB().contains(x,y,z) ) {
                  hit[c] = true;
                  noHit = false;
               }
            }

            // either one neighbor must be hit OR the block is located at the border of the (non-periodic) simulation domain
            if( noHit )
               WALBERLA_CHECK( forest.getBlock(x,y,z) == nullptr )
         }

         // every neighbor must be hit by at least one point
         for( uint_t c = 0; c != block->getNeighborhoodSectionSize(n); ++c )
            WALBERLA_CHECK( hit[c] )

         neighborhoodSectionBlockCenters.clear();
      }
   }
}



static void test() {

   for( uint_t i = 0; i < 5; ++i ) {

      SetupBlockForest forest;

      forest.addRefinementSelectionFunction( refinementSelectionFunctionAll );

      real_t xmin = math::realRandom( real_c(-100), real_c(100) );
      real_t xmax = math::realRandom( xmin + real_c(10), real_c(120) );
      real_t ymin = math::realRandom( real_c(-100), real_c(100) );
      real_t ymax = math::realRandom( ymin + real_c(10), real_c(120) );
      real_t zmin = math::realRandom( real_c(-100), real_c(100) );
      real_t zmax = math::realRandom( zmin + real_c(10), real_c(120) );

      AABB domain( xmin, ymin, zmin, xmax, ymax, zmax );
      forest.init( domain, math::intRandom( uint_t(5), uint_t(20) ), math::intRandom( uint_t(5), uint_t(20) ), math::intRandom( uint_t(5), uint_t(20) ),
                           math::boolRandom(), math::boolRandom(), math::boolRandom() );

      checkNeighborhoodConsistency( forest );
      checkCollectorConsistency( forest );
   }

   for( uint_t i = 0; i < 5; ++i ) {

      SetupBlockForest forest;

      forest.addRefinementSelectionFunction( refinementSelectionFunctionRandom );

      real_t xmin = math::realRandom( real_c(-100), real_c(100) );
      real_t xmax = math::realRandom( xmin + real_c(10), real_c(120) );
      real_t ymin = math::realRandom( real_c(-100), real_c(100) );
      real_t ymax = math::realRandom( ymin + real_c(10), real_c(120) );
      real_t zmin = math::realRandom( real_c(-100), real_c(100) );
      real_t zmax = math::realRandom( zmin + real_c(10), real_c(120) );

      AABB domain( xmin, ymin, zmin, xmax, ymax, zmax );
      forest.init( domain, math::intRandom( uint_t(5), uint_t(20) ), math::intRandom( uint_t(5), uint_t(20) ), math::intRandom( uint_t(5), uint_t(20) ),
                           math::boolRandom(), math::boolRandom(), math::boolRandom() );

      checkNeighborhoodConsistency( forest );
      checkCollectorConsistency( forest );
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
