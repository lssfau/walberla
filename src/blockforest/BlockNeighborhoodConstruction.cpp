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
//! \file BlockNeighborhoodConstruction.cpp
//! \ingroup blockforest
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#include "BlockNeighborhoodConstruction.h"

#include "core/debug/Debug.h"


namespace walberla {
namespace blockforest {



inline static void pushBackNeighborhoodSectionBlockCenters( const AABB& blockAABB, const real_t x, const real_t y, const real_t z,
                                                            std::vector< real_t >& centers ) {

   centers.push_back( blockAABB.xMin() + x * ( blockAABB.xMax() - blockAABB.xMin() ) ); // x
   centers.push_back( blockAABB.yMin() + y * ( blockAABB.yMax() - blockAABB.yMin() ) ); // y
   centers.push_back( blockAABB.zMin() + z * ( blockAABB.zMax() - blockAABB.zMin() ) ); // z
}



void constructNeighborhoodSectionBlockCenters( uint_t sectionIndex, const AABB& blockAABB, std::vector< real_t >& centers ) {

   switch( sectionIndex ) {

   // bottom slice

   case  0: // -1 -1 -1

      pushBackNeighborhoodSectionBlockCenters( blockAABB, -0.25, -0.25, -0.25, centers );
      break;

   case  1: //  0 -1 -1

      pushBackNeighborhoodSectionBlockCenters( blockAABB, +0.25, -0.25, -0.25, centers );
      pushBackNeighborhoodSectionBlockCenters( blockAABB, +0.75, -0.25, -0.25, centers );
      break;

   case  2: //  1 -1 -1

      pushBackNeighborhoodSectionBlockCenters( blockAABB, +1.25, -0.25, -0.25, centers );
      break;

   case  3: // -1  0 -1

      pushBackNeighborhoodSectionBlockCenters( blockAABB, -0.25, +0.25, -0.25, centers );
      pushBackNeighborhoodSectionBlockCenters( blockAABB, -0.25, +0.75, -0.25, centers );
      break;

   case  4: //  0  0 -1

      pushBackNeighborhoodSectionBlockCenters( blockAABB, +0.25, +0.25, -0.25, centers );
      pushBackNeighborhoodSectionBlockCenters( blockAABB, +0.75, +0.25, -0.25, centers );
      pushBackNeighborhoodSectionBlockCenters( blockAABB, +0.25, +0.75, -0.25, centers );
      pushBackNeighborhoodSectionBlockCenters( blockAABB, +0.75, +0.75, -0.25, centers );
      break;

   case  5: //  1  0 -1

      pushBackNeighborhoodSectionBlockCenters( blockAABB, +1.25, +0.25, -0.25, centers );
      pushBackNeighborhoodSectionBlockCenters( blockAABB, +1.25, +0.75, -0.25, centers );
      break;

   case  6: // -1  1 -1

      pushBackNeighborhoodSectionBlockCenters( blockAABB, -0.25, +1.25, -0.25, centers );
      break;

   case  7: //  0  1 -1

      pushBackNeighborhoodSectionBlockCenters( blockAABB, +0.25, +1.25, -0.25, centers );
      pushBackNeighborhoodSectionBlockCenters( blockAABB, +0.75, +1.25, -0.25, centers );
      break;

   case  8: //  1  1 -1

      pushBackNeighborhoodSectionBlockCenters( blockAABB, +1.25, +1.25, -0.25, centers );
      break;

   // middle slice

   case  9: // -1 -1  0

      pushBackNeighborhoodSectionBlockCenters( blockAABB, -0.25, -0.25, +0.25, centers );
      pushBackNeighborhoodSectionBlockCenters( blockAABB, -0.25, -0.25, +0.75, centers );
      break;

   case 10: //  0 -1  0

      pushBackNeighborhoodSectionBlockCenters( blockAABB, +0.25, -0.25, +0.25, centers );
      pushBackNeighborhoodSectionBlockCenters( blockAABB, +0.75, -0.25, +0.25, centers );
      pushBackNeighborhoodSectionBlockCenters( blockAABB, +0.25, -0.25, +0.75, centers );
      pushBackNeighborhoodSectionBlockCenters( blockAABB, +0.75, -0.25, +0.75, centers );
      break;

   case 11: //  1 -1  0

      pushBackNeighborhoodSectionBlockCenters( blockAABB, +1.25, -0.25, +0.25, centers );
      pushBackNeighborhoodSectionBlockCenters( blockAABB, +1.25, -0.25, +0.75, centers );
      break;

   case 12: // -1  0  0

      pushBackNeighborhoodSectionBlockCenters( blockAABB, -0.25, +0.25, +0.25, centers );
      pushBackNeighborhoodSectionBlockCenters( blockAABB, -0.25, +0.75, +0.25, centers );
      pushBackNeighborhoodSectionBlockCenters( blockAABB, -0.25, +0.25, +0.75, centers );
      pushBackNeighborhoodSectionBlockCenters( blockAABB, -0.25, +0.75, +0.75, centers );
      break;

   // 0 0 0

   case 13: //  1  0  0

      pushBackNeighborhoodSectionBlockCenters( blockAABB, +1.25, +0.25, +0.25, centers );
      pushBackNeighborhoodSectionBlockCenters( blockAABB, +1.25, +0.75, +0.25, centers );
      pushBackNeighborhoodSectionBlockCenters( blockAABB, +1.25, +0.25, +0.75, centers );
      pushBackNeighborhoodSectionBlockCenters( blockAABB, +1.25, +0.75, +0.75, centers );
      break;

   case 14: // -1  1  0

      pushBackNeighborhoodSectionBlockCenters( blockAABB, -0.25, +1.25, +0.25, centers );
      pushBackNeighborhoodSectionBlockCenters( blockAABB, -0.25, +1.25, +0.75, centers );
      break;

   case 15: //  0  1  0

      pushBackNeighborhoodSectionBlockCenters( blockAABB, +0.25, +1.25, +0.25, centers );
      pushBackNeighborhoodSectionBlockCenters( blockAABB, +0.75, +1.25, +0.25, centers );
      pushBackNeighborhoodSectionBlockCenters( blockAABB, +0.25, +1.25, +0.75, centers );
      pushBackNeighborhoodSectionBlockCenters( blockAABB, +0.75, +1.25, +0.75, centers );
      break;

   case 16: //  1  1  0

      pushBackNeighborhoodSectionBlockCenters( blockAABB, +1.25, +1.25, +0.25, centers );
      pushBackNeighborhoodSectionBlockCenters( blockAABB, +1.25, +1.25, +0.75, centers );
      break;

   // top slice

   case 17: // -1 -1  1

      pushBackNeighborhoodSectionBlockCenters( blockAABB, -0.25, -0.25, +1.25, centers );
      break;

   case 18: //  0 -1  1

      pushBackNeighborhoodSectionBlockCenters( blockAABB, +0.25, -0.25, +1.25, centers );
      pushBackNeighborhoodSectionBlockCenters( blockAABB, +0.75, -0.25, +1.25, centers );
      break;

   case 19: //  1 -1  1

      pushBackNeighborhoodSectionBlockCenters( blockAABB, +1.25, -0.25, +1.25, centers );
      break;

   case 20: // -1  0  1

      pushBackNeighborhoodSectionBlockCenters( blockAABB, -0.25, +0.25, +1.25, centers );
      pushBackNeighborhoodSectionBlockCenters( blockAABB, -0.25, +0.75, +1.25, centers );
      break;

   case 21: //  0  0  1

      pushBackNeighborhoodSectionBlockCenters( blockAABB, +0.25, +0.25, +1.25, centers );
      pushBackNeighborhoodSectionBlockCenters( blockAABB, +0.75, +0.25, +1.25, centers );
      pushBackNeighborhoodSectionBlockCenters( blockAABB, +0.25, +0.75, +1.25, centers );
      pushBackNeighborhoodSectionBlockCenters( blockAABB, +0.75, +0.75, +1.25, centers );
      break;

   case 22: //  1  0  1

      pushBackNeighborhoodSectionBlockCenters( blockAABB, +1.25, +0.25, +1.25, centers );
      pushBackNeighborhoodSectionBlockCenters( blockAABB, +1.25, +0.75, +1.25, centers );
      break;

   case 23: // -1  1  1

      pushBackNeighborhoodSectionBlockCenters( blockAABB, -0.25, +1.25, +1.25, centers );
      break;

   case 24: //  0  1  1

      pushBackNeighborhoodSectionBlockCenters( blockAABB, +0.25, +1.25, +1.25, centers );
      pushBackNeighborhoodSectionBlockCenters( blockAABB, +0.75, +1.25, +1.25, centers );
      break;

   case 25: //  1  1  1

      pushBackNeighborhoodSectionBlockCenters( blockAABB, +1.25, +1.25, +1.25, centers );
      break;

   default:
      WALBERLA_ASSERT( false )
   }
}



} // namespace blockforest
} // namespace walberla
