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
//! \file CellVector.cpp
//! \ingroup core
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//
//======================================================================================================================

#include "CellVector.h"

#include <algorithm>
#include <ostream>
#include <stdexcept>


namespace walberla {
namespace cell {



//**********************************************************************************************************************
/**
 * \brief   Removes duplicate entries.
 *
 * During the operation the CellVector is sorted. Iterators get invalidated.
 *
 * \return  The number of duplicate entries removed.
 */
//**********************************************************************************************************************
CellVector::difference_type CellVector::removeDuplicates()
{
   std::sort(cells_.begin(), cells_.end());

   iterator newEnd = std::unique(cells_.begin(), cells_.end());

   const difference_type numDuplicates = std::distance(newEnd, cells_.end());

   cells_.erase(newEnd, cells_.end());

   return numDuplicates;
}



//**********************************************************************************************************************
/**
 * \brief   Calculates the axis-aligned bounding box of the cell vector.
 *
 * Complexity is O(N), where N == this->size()
 *
 * \exception  std::domain_error Thrown when the CellVector is empty.
 *
 * \return  CellInterval that stores the bounding box.
 */
//**********************************************************************************************************************
CellInterval CellVector::boundingBox() const
{
   if(empty())
      throw std::domain_error("Calculation of bounding box from empty CellVector is not possible!");

   CellInterval interval(cells_.front().x(), cells_.front().y(), cells_.front().z(), cells_.front().x(), cells_.front().y(), cells_.front().z());

   for(CellVector::const_iterator cellIt = cells_.begin() + 1; cellIt != cells_.end(); ++cellIt)
   {
      if(cellIt->x() < interval.xMin())
         interval.xMin() = cellIt->x();
      if(cellIt->y() < interval.yMin())
         interval.yMin() = cellIt->y();
      if(cellIt->z() < interval.zMin())
         interval.zMin() = cellIt->z();

      if(cellIt->x() > interval.xMax())
         interval.xMax() = cellIt->x();
      if(cellIt->y() > interval.yMax())
         interval.yMax() = cellIt->y();
      if(cellIt->z() > interval.zMax())
         interval.zMax() = cellIt->z();
   }
   return interval;
}



//**********************************************************************************************************************
/** \brief   Output stream operator for CellIntervals
 *
 * The Cell vector is serialized in the Form "[2]{(1 1 1) (1 1 2)}" for a CellVector with 2 entries
 *
 * \param[in,out] os    output stream
 * \param[in]  cells    CellVector to serialize
 *
 * \return  the modified output stream.
 */
//**********************************************************************************************************************
std::ostream & operator<<(std::ostream & os, const CellVector & cells)
{
   os << "[" << cells.size() << "]{";

   if(!cells.empty())
   {
      CellVector::const_iterator cellIt;
      for( cellIt = cells.begin(); cellIt != cells.end() - 1; ++cellIt)
         os << *cellIt << " ";
      os << *cellIt;
   }

   os << "}";

   return os;
}



} // namespace cell
} // namespace walberla
