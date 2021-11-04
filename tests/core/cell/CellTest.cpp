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
//! \file CellTest.cpp
//! \ingroup core
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//! \brief Unit test for class walberla::Cell
//
//======================================================================================================================

#include "core/cell/Cell.h"
#include "core/math/Vector2.h"
#include "core/math/Vector3.h"
#include "core/debug/TestSubsystem.h"

#include <random>
#include <unordered_set>


using namespace walberla;

void testCell( cell_idx_t x, cell_idx_t y, cell_idx_t z )
{
   Cell c(x,y,z);

   WALBERLA_CHECK_EQUAL( c.x(), x );
   WALBERLA_CHECK_EQUAL( c.y(), y );
   WALBERLA_CHECK_EQUAL( c.z(), z );

   WALBERLA_CHECK_EQUAL( c[0], x );
   WALBERLA_CHECK_EQUAL( c[1], y );
   WALBERLA_CHECK_EQUAL( c[2], z );

   WALBERLA_CHECK_EQUAL( c.positiveIndicesOnly(), x >= cell_idx_t(0) && y >= cell_idx_t(0) && z >= cell_idx_t(0) );

   Cell sum = c + c;
   WALBERLA_CHECK_EQUAL( sum.x(), x + x );
   WALBERLA_CHECK_EQUAL( sum.y(), y + y );
   WALBERLA_CHECK_EQUAL( sum.z(), z + z );
   Cell tmp(c);
   tmp += c;
   WALBERLA_CHECK_EQUAL( tmp, sum );

   Cell difference = c - c;
   WALBERLA_CHECK_EQUAL( difference.x(), cell_idx_t(0) );
   WALBERLA_CHECK_EQUAL( difference.y(), cell_idx_t(0) );
   WALBERLA_CHECK_EQUAL( difference.z(), cell_idx_t(0) );
   tmp = c;
   tmp -= c;
   WALBERLA_CHECK_EQUAL( tmp, difference);

   WALBERLA_CHECK_EQUAL( c, +c );
   WALBERLA_CHECK_EQUAL( Cell(0,0,0) - c, -c );

   WALBERLA_CHECK( c == c );
   WALBERLA_CHECK( !(c != c) );
   WALBERLA_CHECK( !(c < c) );

   Cell cellCopied( c );
   Cell cellAssigned;
   cellAssigned = c;

   WALBERLA_CHECK_EQUAL( c, cellCopied );
   WALBERLA_CHECK_EQUAL( c, cellAssigned );

   std::stringstream ss;
   const int forty_two = 42;
   const int twenty_three = 23;
   ss << forty_two << c << twenty_three;
   Cell rc;
   int i0, i1;
   ss >> i0 >> rc >> i1;
   WALBERLA_CHECK_EQUAL( c, rc );
   WALBERLA_CHECK_EQUAL( i0, forty_two );
   WALBERLA_CHECK_EQUAL( i1, twenty_three );
   WALBERLA_CHECK( !ss.bad() );
   WALBERLA_CHECK( ss.eof() );
}

void testBinaryOperators( const Cell & c0, const Cell & c1 )
{
   Cell sum = c0 + c1;
   WALBERLA_CHECK_EQUAL( sum.x(), c0.x() + c1.x() );
   WALBERLA_CHECK_EQUAL( sum.y(), c0.y() + c1.y() );
   WALBERLA_CHECK_EQUAL( sum.z(), c0.z() + c1.z() );
   Cell tmp = c0;
   tmp += c1;
   WALBERLA_CHECK_EQUAL( tmp, sum );

   Cell difference = c0 - c1;
   WALBERLA_CHECK_EQUAL( difference.x(), c0.x() - c1.x() );
   WALBERLA_CHECK_EQUAL( difference.y(), c0.y() - c1.y() );
   WALBERLA_CHECK_EQUAL( difference.z(), c0.z() - c1.z() );
   tmp = c0;
   tmp -= c1;
   WALBERLA_CHECK_EQUAL( tmp, difference);

   WALBERLA_CHECK_EQUAL( c0 == c1, c0.x() == c1.x() && c0.y() == c1.y() && c0.z() == c1.z() );
   WALBERLA_CHECK_EQUAL( c0 != c1, c0.x() != c1.x() || c0.y() != c1.y() || c0.z() != c1.z() );

   bool isLess = ( c0.z() < c1.z() || ( c0.z() == c1.z() && c0.y() < c1.y() ) || ( c0.z() == c1.z() && c0.y() == c1.y() && c0.x() < c1.x() ) );


   WALBERLA_CHECK_EQUAL( c0 < c1, isLess );

   std::stringstream ss;
   const int forty_two = 42;
   const int twenty_three = 23;
   const int thirteen = 13;
   ss << forty_two << c0 << twenty_three << c1 << thirteen;
   Cell rc0, rc1;
   int i0, i1, i2;
   ss >> i0 >> rc0 >> i1 >> rc1 >> i2;
   WALBERLA_CHECK_EQUAL( c0, rc0 );
   WALBERLA_CHECK_EQUAL( c1, rc1 );
   WALBERLA_CHECK_EQUAL( i0, forty_two );
   WALBERLA_CHECK_EQUAL( i1, twenty_three );
   WALBERLA_CHECK_EQUAL( i2, thirteen );
   WALBERLA_CHECK( !ss.bad() );
   WALBERLA_CHECK( ss.eof() );
}

void testHashAlgorithm()
{
   auto const hasher = std::hash< walberla::Cell >();

   // check hash concatenates individual elements
   std::size_t const prefix = hasher(walberla::Cell{15, 6, 0});
   cell_idx_t const max_z = (sizeof(std::size_t) >= 8)? 2<<21 : 2<<10;
   std::size_t mismatches = 0;
   for( cell_idx_t z = 0; z < max_z; z += 5 )
   {
      auto const cell = walberla::Cell{15, 6, z};
      auto const expected_hash = prefix + static_cast<std::size_t>(z);
      if( hasher(cell) != expected_hash )
      {
         mismatches++;
      }
   }
   WALBERLA_CHECK_EQUAL( mismatches, 0 );

   // check hash collisions (use a small block size to limit memory footprint)
   cell_idx_t const block_size = 128;
   cell_idx_t const ghost_layer = 8;
   std::unordered_set<std::size_t> keys{};
   std::size_t collisions = 0;
   for( auto x = -ghost_layer; x < block_size + ghost_layer; ++x )
   {
      for( auto y = -ghost_layer; y < block_size + ghost_layer; ++y )
      {
         for( auto z = -ghost_layer; z < block_size + ghost_layer; ++z )
         {
            auto const cell = walberla::Cell{x, y, z};
            auto const hash = hasher(cell);
            if (keys.count(hash))
            {
               collisions++;
            }
            else
            {
               keys.emplace(hash);
            }
         }
      }
   }
   WALBERLA_CHECK_EQUAL( collisions, 0 );

   // check hash matches with Vector2 and Vector3
   auto const hasher2 = std::hash< walberla::Vector2<int> >();
   auto const hasher3 = std::hash< walberla::Vector3<int> >();
   
   auto const cell = walberla::Cell{15, 6, 42};
   auto const vec3 = walberla::Vector3<int>{15, 6, 42};
   WALBERLA_CHECK_EQUAL( hasher(cell), hasher3(vec3) );
   auto const vec2 = walberla::Vector2<int>{15, 6};
   auto const vec3b = walberla::Vector3<int>{15, 6, 0};
   WALBERLA_CHECK_EQUAL( hasher2(vec2), hasher3(vec3b) );
}


int main( int /*argc*/, char** /*argv*/ ) {

   debug::enterTestMode();

   typedef std::mersenne_twister_engine< walberla::uint32_t, 32, 351, 175, 19, 0xccab8ee7, 11, 0xffffffff, 7, 0x31b6ab00, 15, 0xffe50000, 17, 0xa37d3c92 > mt11213b;
   mt11213b rng;
   std::uniform_int_distribution<cell_idx_t> dist( std::numeric_limits<cell_idx_t>::min(), std::numeric_limits<cell_idx_t>::max() );

   for(int i = 0; i < 100000; ++i)
   {
      cell_idx_t x0 = dist(rng);
      cell_idx_t y0 = dist(rng);
      cell_idx_t z0 = dist(rng);
      cell_idx_t x1 = dist(rng);
      cell_idx_t y1 = dist(rng);
      cell_idx_t z1 = dist(rng);

      testCell(x0, y0, z0);
      testCell(x1, y1, z1);

      testBinaryOperators( Cell(x0, y0, z0), Cell(x1, y1, z1) );
   }

   testHashAlgorithm();

   return 0;
}
