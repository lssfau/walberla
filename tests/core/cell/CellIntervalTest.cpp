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
//! \file CellIntervalTest.cpp
//! \ingroup core
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//! \brief Unit test for class walberla::CellInterval
//
//======================================================================================================================

#include "core/cell/CellInterval.h"
#include "core/debug/TestSubsystem.h"

#include <random>
#include <iostream>
#include <iterator>


using namespace walberla;

typedef std::mersenne_twister_engine< walberla::uint32_t, 32, 351, 175, 19, 0xccab8ee7, 11, 0xffffffff, 7, 0x31b6ab00, 15, 0xffe50000, 17, 0xa37d3c92 > mt11213b;

CellInterval makeRandomInterval(uint_t maxSize)
{
   static mt11213b rng;
   std::uniform_int_distribution<cell_idx_t> dist( std::numeric_limits<cell_idx_t>::min(), std::numeric_limits<cell_idx_t>::max() - cell_idx_c( maxSize ) );
   std::uniform_int_distribution<uint_t> dist2( uint_t(0), maxSize );

   cell_idx_t xMin = dist(rng);
   cell_idx_t yMin = dist(rng);
   cell_idx_t zMin = dist(rng);

   cell_idx_t xMax = xMin + cell_idx_c( dist2(rng) );
   cell_idx_t yMax = yMin + cell_idx_c( dist2(rng) );
   cell_idx_t zMax = zMin + cell_idx_c( dist2(rng) );

   return CellInterval(xMin, yMin, zMin, xMax, yMax, zMax);
}

CellInterval makeRandomEmptyInterval(uint_t maxSize)
{
   static mt11213b rng;
   std::uniform_int_distribution<cell_idx_t> dist( std::numeric_limits<cell_idx_t>::min() + cell_idx_c( maxSize ), std::numeric_limits<cell_idx_t>::max() - cell_idx_c( maxSize ) );
   std::uniform_int_distribution<uint_t> dist2( uint_t(1), maxSize );
   std::uniform_int_distribution<uint_t> dist3( uint_t(0), uint_t(1) );

   cell_idx_t xMin = dist(rng);
   cell_idx_t yMin = dist(rng);
   cell_idx_t zMin = dist(rng);


   bool xNegative = false;
   bool yNegative = false;
   bool zNegative = false;

   while( !( xNegative || yNegative || zNegative) )
   {
      // determine which component(s) of the CI should be negative (one must be to make the CI empty)
      xNegative = dist3( rng ) == 1;
      yNegative = dist3( rng ) == 1;
      zNegative = dist3( rng ) == 1;
   }

   cell_idx_t xMax = xMin + ( xNegative ? - cell_idx_c( dist2(rng) ) : cell_idx_c( dist2(rng) ) );
   cell_idx_t yMax = yMin + ( yNegative ? - cell_idx_c( dist2(rng) ) : cell_idx_c( dist2(rng) ) );
   cell_idx_t zMax = zMin + ( zNegative ? - cell_idx_c( dist2(rng) ) : cell_idx_c( dist2(rng) ) );

   return CellInterval(xMin, yMin, zMin, xMax, yMax, zMax);
}

Cell makeRandomCell()
{
   static mt11213b rng;
   std::uniform_int_distribution<cell_idx_t> dist( std::numeric_limits<cell_idx_t>::min(), std::numeric_limits<cell_idx_t>::max() );
   return Cell( dist(rng), dist(rng), dist(rng) );
}

inline void testCI( const CellInterval & ci )
{
   WALBERLA_CHECK( !ci.empty() );

   Cell minCorner = ci.min();
   Cell maxCorner = ci.max();

   WALBERLA_CHECK_EQUAL( minCorner.x(), ci.xMin() );
   WALBERLA_CHECK_EQUAL( minCorner.y(), ci.yMin() );
   WALBERLA_CHECK_EQUAL( minCorner.z(), ci.zMin() );
   WALBERLA_CHECK_EQUAL( maxCorner.x(), ci.xMax() );
   WALBERLA_CHECK_EQUAL( maxCorner.y(), ci.yMax() );
   WALBERLA_CHECK_EQUAL( maxCorner.z(), ci.zMax() );

   WALBERLA_CHECK_EQUAL( ci.positiveIndicesOnly(), minCorner.positiveIndicesOnly() && maxCorner.positiveIndicesOnly() );

   WALBERLA_CHECK_EQUAL( ci, ci );
   WALBERLA_CHECK( !( ci != ci ) );
   WALBERLA_CHECK( ci.overlaps(ci) );

   WALBERLA_CHECK_EQUAL( ci.xSize(), ci.xMax() - ci.xMin() + 1 );
   WALBERLA_CHECK_EQUAL( ci.ySize(), ci.yMax() - ci.yMin() + 1 );
   WALBERLA_CHECK_EQUAL( ci.zSize(), ci.zMax() - ci.zMin() + 1 );

   WALBERLA_CHECK_EQUAL( ci.numCells(), ci.xSize() * ci.ySize() * ci.zSize() );

   CellInterval ci_copied( ci );
   CellInterval ci_assigned;
   ci_assigned = ci;

   WALBERLA_CHECK_EQUAL( ci, ci_copied );
   WALBERLA_CHECK_EQUAL( ci, ci_assigned );

   std::stringstream ss;
   const int forty_two = 42;
   const int twenty_three = 23;
   ss << forty_two << ci << twenty_three;
   CellInterval rci;
   int i0, i1;
   ss >> i0 >> rci >> i1;
   WALBERLA_CHECK_EQUAL( ci, rci );
   WALBERLA_CHECK_EQUAL( i0, forty_two );
   WALBERLA_CHECK_EQUAL( i1, twenty_three );
   WALBERLA_CHECK( !ss.bad() );
   WALBERLA_CHECK( ss.eof() );
}

void testIterators( const CellInterval & ci )
{
   uint_t ctr = 0;
   auto it = ci.begin();
   for( cell_idx_t z = ci.zMin(); z <= ci.zMax(); ++z)
      for( cell_idx_t y = ci.yMin(); y <= ci.yMax(); ++y)
         for( cell_idx_t x = ci.xMin(); x <= ci.xMax(); ++x, ++it)
         {
            WALBERLA_CHECK( ci.contains( *it ) );

            WALBERLA_CHECK_EQUAL( x, it->x() );
            WALBERLA_CHECK_EQUAL( y, it->y() );
            WALBERLA_CHECK_EQUAL( z, it->z() );

            ++ctr;
         }

   WALBERLA_CHECK_EQUAL( it, ci.end() );
   WALBERLA_CHECK_EQUAL( ctr, ci.numCells() );
   WALBERLA_CHECK_EQUAL( std::distance( ci.begin(), ci.end() ), ci.numCells() );

   ctr = 0;
   auto it2 = ci.end();
   --it2;
   for( cell_idx_t z = ci.zMax(); z >= ci.zMin(); --z)
   {
      for (cell_idx_t y = ci.yMax(); y >= ci.yMin(); --y)
      {
         for( cell_idx_t x = ci.xMax(); x >= ci.xMin(); --x, --it2)
         {
            WALBERLA_CHECK( ci.contains( *it2 ) );

            WALBERLA_CHECK_EQUAL( x, it2->x() );
            WALBERLA_CHECK_EQUAL( y, it2->y() );
            WALBERLA_CHECK_EQUAL( z, it2->z() );

            ++ctr;
         }
      }
   }
   WALBERLA_CHECK_EQUAL(ctr, ci.numCells());
   WALBERLA_CHECK_EQUAL(it2, --ci.begin());
}

void testEmptyCI( const CellInterval & ci )
{
   WALBERLA_CHECK( ci.empty() );

   for(int i = 0; i < 100; ++i)
      WALBERLA_CHECK( !ci.contains(makeRandomCell()) );

   for(int i = 0; i < 100; ++i)
      WALBERLA_CHECK( !ci.overlaps(makeRandomInterval(1000000)) );

   std::stringstream ss;
   const int forty_two = 42;
   const int twenty_three = 23;
   ss << forty_two << ci << twenty_three;
   CellInterval rci;
   int i0, i1;
   ss >> i0 >> rci >> i1;
   WALBERLA_CHECK_EQUAL( ci, rci );
   WALBERLA_CHECK_EQUAL( i0, forty_two );
   WALBERLA_CHECK_EQUAL( i1, twenty_three );
   WALBERLA_CHECK( !ss.bad() );
   WALBERLA_CHECK( ss.eof() );

   size_t numIterations = 0;
   for( auto it = ci.begin(); it != ci.end() && numIterations < 1; ++it )
      ++numIterations;

   WALBERLA_CHECK_EQUAL( numIterations, 0 );
}

int main( int /*argc*/, char** /*argv*/ ) {

   debug::enterTestMode();

   WALBERLA_CHECK( CellInterval().empty() );

   for(int i = 0; i < 1000; ++i)
   {
      CellInterval emptyCI( makeRandomEmptyInterval(10000000) );
      testEmptyCI(emptyCI);
   }

   for(int i = 0; i < 1000; ++i)
   {
      CellInterval ci( makeRandomInterval(10000000) );
      testCI(ci);
   }

   for(int i = 0; i < 1000; ++i)
   {
      CellInterval emptyCI( makeRandomEmptyInterval(100) );
      testEmptyCI(emptyCI);
   }

   for(int i = 0; i < 1000; ++i)
   {
      CellInterval ci( makeRandomInterval(100) );
      testCI(ci);
   }

   for(int i = 0; i < 10; ++i)
   {
      CellInterval ci( makeRandomInterval(100) );
      testIterators(ci);
   }

   return 0;
}
