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
//! \file SetupBlock.cpp
//! \ingroup blockforest
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#include "SetupBlock.h"

#include <set>


namespace walberla {
namespace blockforest {



void SetupBlock::assembleNeighborhood() {

   neighborhood_.clear();

   std::set< SetupBlock* > neighborhood;

   for( uint_t n = 0; n != 26; ++n )
      for( uint_t i = 0; i != neighborhoodSection_[n].size(); ++i )
         if( neighborhood.insert( neighborhoodSection_[n][i] ).second )
            neighborhood_.push_back( neighborhoodSection_[n][i] );
}



void SetupBlock::split() {

   WALBERLA_ASSERT( children_.empty() );

   const real_t xMid = ( aabb_.xMin() + aabb_.xMax() ) / real_c(2);
   const real_t yMid = ( aabb_.yMin() + aabb_.yMax() ) / real_c(2);
   const real_t zMid = ( aabb_.zMin() + aabb_.zMax() ) / real_c(2);   
   
   for( uint_t c = 0; c != 8; ++c )
   {
      BlockID childId( Id_, c );

      children_.push_back( new SetupBlock( this, childId,
         ( ( c & 1 ) ? xMid : aabb_.xMin() ), // xmin (incl.)
         ( ( c & 2 ) ? yMid : aabb_.yMin() ), // ymin (incl.)
         ( ( c & 4 ) ? zMid : aabb_.zMin() ), // zmin (incl.)
         ( ( c & 1 ) ? aabb_.xMax() : xMid ), // xmax (excl.)
         ( ( c & 2 ) ? aabb_.yMax() : yMid ), // ymax (excl.)
         ( ( c & 4 ) ? aabb_.zMax() : zMid ), // zmax (excl.)
         level_ + 1 ) );
   }
}



} // namespace blockforest
} // namespace walberla
