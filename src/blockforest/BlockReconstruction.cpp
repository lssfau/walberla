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
//! \file BlockReconstruction.cpp
//! \ingroup blockforest
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#include "BlockReconstruction.h"


namespace walberla {
namespace blockforest {



uint_t BlockReconstruction::reconstructAABB( AABB& aabb, const BlockID& blockId, const AABB& domain, const uint_t xSize, const uint_t ySize,
                                             const uint_t zSize, const uint_t treeIdDigits ) {

   WALBERLA_ASSERT_GREATER_EQUAL( blockId.getUsedBits(), treeIdDigits );
   WALBERLA_ASSERT_EQUAL( ( blockId.getUsedBits() - treeIdDigits ) % 3, 0 );

   BlockID id( blockId );

   const uint_t levels = ( id.getUsedBits() - treeIdDigits ) / 3;

   std::vector< uint_t > branchId( levels );

   for( uint_t i = levels; i-- != 0; ) {
      branchId[i] = id.getBranchId();
      id.removeBranchId();
   }

   uint_t index = id.getTreeIndex();

   const uint_t z = index / ( xSize * ySize );
            index = index % ( xSize * ySize );
   const uint_t y = index / xSize;
   const uint_t x = index % xSize;

   const real_t xw = ( domain.xMax() - domain.xMin() ) / static_cast< real_t >( xSize );
   const real_t yw = ( domain.yMax() - domain.yMin() ) / static_cast< real_t >( ySize );
   const real_t zw = ( domain.zMax() - domain.zMin() ) / static_cast< real_t >( zSize );

   real_t xMin =                                    domain.xMin() + static_cast< real_t >(  x  ) * xw;
   real_t xMax = ( x+1 == xSize ) ? domain.xMax() : domain.xMin() + static_cast< real_t >( x+1 ) * xw;

   real_t yMin =                                    domain.yMin() + static_cast< real_t >(  y  ) * yw;
   real_t yMax = ( y+1 == ySize ) ? domain.yMax() : domain.yMin() + static_cast< real_t >( y+1 ) * yw;

   real_t zMin =                                    domain.zMin() + static_cast< real_t >(  z  ) * zw;
   real_t zMax = ( z+1 == zSize ) ? domain.zMax() : domain.zMin() + static_cast< real_t >( z+1 ) * zw;

   for( uint_t i = 0; i != levels; ++i ) {

      const real_t xMid = ( xMin + xMax ) / real_c(2);
      const real_t yMid = ( yMin + yMax ) / real_c(2);
      const real_t zMid = ( zMin + zMax ) / real_c(2);

      xMin = ( branchId[i] & 1 ) ? xMid : xMin;
      xMax = ( branchId[i] & 1 ) ? xMax : xMid;

      yMin = ( branchId[i] & 2 ) ? yMid : yMin;
      yMax = ( branchId[i] & 2 ) ? yMax : yMid;

      zMin = ( branchId[i] & 4 ) ? zMid : zMin;
      zMax = ( branchId[i] & 4 ) ? zMax : zMid;
   }

   aabb.initMinMaxCorner( xMin, yMin, zMin, xMax, yMax, zMax );

   return levels;
}



} // namespace blockforest
} // namespace walberla
