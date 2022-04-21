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
//! \file CellSet.cpp
//! \ingroup core
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#include "CellSet.h"

#include "core/debug/Debug.h"


namespace walberla {
namespace cell {


/// Complexity is O(N), where N == this->size()
CellInterval CellSet::boundingBox() const {

   WALBERLA_ASSERT( !empty() )

   Set<Cell>::const_iterator beginIt = Set<Cell>::begin();
   Set<Cell>::const_iterator endIt   = Set<Cell>::end();

   CellInterval interval( beginIt->x(), beginIt->y(), beginIt->z(), beginIt->x(), beginIt->y(), beginIt->z() );

   for( Set<Cell>::const_iterator cellIt = ++beginIt; cellIt != endIt; ++cellIt ) {

      if( cellIt->x() < interval.xMin() ) interval.xMin() = cellIt->x();
      if( cellIt->y() < interval.yMin() ) interval.yMin() = cellIt->y();
      if( cellIt->z() < interval.zMin() ) interval.zMin() = cellIt->z();

      if( cellIt->x() > interval.xMax() ) interval.xMax() = cellIt->x();
      if( cellIt->y() > interval.yMax() ) interval.yMax() = cellIt->y();
      if( cellIt->z() > interval.zMax() ) interval.zMax() = cellIt->z();
   }

   return interval;
}


} // namespace cell
} // namespace walberla
