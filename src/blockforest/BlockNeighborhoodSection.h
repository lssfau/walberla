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
//! \file BlockNeighborhoodSection.h
//! \ingroup blockforest
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/DataTypes.h"
#include "core/debug/Debug.h"
#include "stencil/Directions.h"

#include <array>

namespace walberla {
namespace blockforest {



inline uint_t getBlockNeighborhoodSectionIndex( const int x, const int y, const int z )
{
   WALBERLA_ASSERT_LESS_EQUAL( -1, x ) WALBERLA_ASSERT_LESS_EQUAL( x, 1 )
   WALBERLA_ASSERT_LESS_EQUAL( -1, y ) WALBERLA_ASSERT_LESS_EQUAL( y, 1 )
   WALBERLA_ASSERT_LESS_EQUAL( -1, z ) WALBERLA_ASSERT_LESS_EQUAL( z, 1 )

   WALBERLA_ASSERT( !( x == 0 && y == 0 && z == 0 ) )

   const uint_t index = uint_c( (z+1) * 9 + (y+1) * 3 + x+1 );

   if( index > uint_t(13) )
      return index - uint_t(1);
   return index;
}



inline uint_t getBlockNeighborhoodSectionIndex( const uint_t x, const uint_t y, const uint_t z )
{
   WALBERLA_ASSERT_LESS( x, 3 )
   WALBERLA_ASSERT_LESS( y, 3 )
   WALBERLA_ASSERT_LESS( z, 3 )

   WALBERLA_ASSERT( !( x == 1 && y == 1 && z == 1 ) )

   const uint_t index = z * uint_t(9) + y * uint_t(3) + x;

   if( index > uint_t(13) )
      return index - uint_t(1);
   return index;
}



inline uint_t getBlockNeighborhoodSectionIndex( const stencil::Direction & d )
{
   WALBERLA_ASSERT_UNEQUAL( d, stencil::C )
   WALBERLA_ASSERT( !( stencil::cx[d] == 0 && stencil::cy[d] == 0 && stencil::cz[d] == 0 ) )

   return getBlockNeighborhoodSectionIndex( stencil::cx[d], stencil::cy[d], stencil::cz[d] );
}



inline uint_t getBlockMaxNeighborhoodSectionSize( const uint_t sectionIndex )
{
   WALBERLA_ASSERT_LESS( sectionIndex, uint_t(26) )

   // faces
   if( sectionIndex == uint_t( 4) || sectionIndex == uint_t(10) || sectionIndex == uint_t(12) ||
       sectionIndex == uint_t(13) || sectionIndex == uint_t(15) || sectionIndex == uint_t(21) )
      return uint_t(4);
   // corners
   if( sectionIndex == uint_t( 0) || sectionIndex == uint_t( 2) || sectionIndex == uint_t( 6) || sectionIndex == uint_t( 8) ||
       sectionIndex == uint_t(17) || sectionIndex == uint_t(19) || sectionIndex == uint_t(23) || sectionIndex == uint_t(25) )
      return uint_t(1);
   // edges
   return uint_t(2);
}



inline const std::array<uint_t, 6> & getFaceNeighborhoodSectionIndices()
{
   static std::array<uint_t, 6> faces{ { uint_t(4), uint_t(10), uint_t(12), uint_t(13), uint_t(15), uint_t(21) } };
   return faces;
}


inline const std::array<uint_t,12> & getEdgeNeighborhoodSectionIndices()
{
   static std::array<uint_t,12> edges{{ uint_t( 1), uint_t( 3), uint_t( 5), uint_t( 7), uint_t( 9), uint_t(11),
                                        uint_t(14), uint_t(16), uint_t(18), uint_t(20), uint_t(22), uint_t(24) }};
   return edges;
}


inline const std::array<uint_t, 8> & getCornerNeighborhoodSectionIndices()
{
   static std::array<uint_t, 8> corners {{ uint_t( 0), uint_t( 2), uint_t( 6), uint_t( 8),
                                           uint_t(17), uint_t(19), uint_t(23), uint_t(25) }};
   return corners;
}


} // namespace blockforest
} // namespace walberla
