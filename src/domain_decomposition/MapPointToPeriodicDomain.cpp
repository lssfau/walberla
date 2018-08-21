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
//! \file MapPointToPeriodicDomain.cpp
//! \ingroup domain_decomposition
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#include "MapPointToPeriodicDomain.h"
#include "core/debug/Debug.h"



namespace walberla {
namespace domain_decomposition {



void mapPointToPeriodicDomain( const std::array< bool, 3 > & periodic, const AABB & domain, real_t & x, real_t & y, real_t & z )
{
   if( periodic[0] )
   {
      const real_t domainSize = domain.xMax() - domain.xMin();
      if( x < domain.xMin() )
      {
         real_t shift = ( domain.xMin() - x ) - ( std::floor( ( domain.xMin() - x ) / domainSize ) * domainSize );
         shift = std::max( shift, real_t(0) );
         x = domain.xMax() - shift;
         if( isIdentical( x, domain.xMax() ) || x < domain.xMin() )
            x = std::nextafter(domain.xMax(), domain.xMin());
      }
      else if( x >= domain.xMax() )
      {
         real_t shift = ( x - domain.xMax() ) - ( std::floor( ( x - domain.xMax() ) / domainSize ) * domainSize );
         shift = std::max( shift, real_t(0) );
         x = domain.xMin() + shift;
         if( x >= domain.xMax() )
            x = domain.xMin();
      }
   }

   if( periodic[1] )
   {
      const real_t domainSize = domain.yMax() - domain.yMin();
      if( y < domain.yMin() )
      {
         real_t shift = ( domain.yMin() - y ) - ( std::floor( ( domain.yMin() - y ) / domainSize ) * domainSize );
         shift = std::max( shift, real_t(0) );
         y = domain.yMax() - shift;
         if( isIdentical( y, domain.yMax() ) || y < domain.yMin() )
            y = std::nextafter(domain.yMax(), domain.yMin());
      }
      else if( y >= domain.yMax() )
      {
         real_t shift = ( y - domain.yMax() ) - ( std::floor( ( y - domain.yMax() ) / domainSize ) * domainSize );
         shift = std::max( shift, real_t(0) );
         y = domain.yMin() + shift;
         if( y >= domain.yMax() )
            y = domain.yMin();
      }
   }

   if( periodic[2] )
   {
      const real_t domainSize = domain.zMax() - domain.zMin();
      if( z < domain.zMin() )
      {
         real_t shift = ( domain.zMin() - z ) - ( std::floor( ( domain.zMin() - z ) / domainSize ) * domainSize );
         shift = std::max( shift, real_t(0) );
         z = domain.zMax() - shift;
         if( isIdentical( z, domain.zMax() ) || z < domain.zMin() )
            z = std::nextafter(domain.zMax(), domain.zMin());
      }
      else if( z >= domain.zMax() )
      {
         real_t shift = ( z - domain.zMax() ) - ( std::floor( ( z - domain.zMax() ) / domainSize ) * domainSize );
         shift = std::max( shift, real_t(0) );
         z = domain.zMin() + shift;
         if( z >= domain.zMax() )
            z = domain.zMin();
      }
   }

   WALBERLA_ASSERT( !periodic[0] || ( domain.xMin() <= x && x < domain.xMax() ) );
   WALBERLA_ASSERT( !periodic[1] || ( domain.yMin() <= y && y < domain.yMax() ) );
   WALBERLA_ASSERT( !periodic[2] || ( domain.zMin() <= z && z < domain.zMax() ) );
}



} // namespace domain_decomposition
} // namespace walberla
