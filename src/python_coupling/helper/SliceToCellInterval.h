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
//! \file SliceToCellInterval.h
//! \ingroup python_coupling
//! \author Martin Bauer <martin.bauer@fau.de>
//! \author Markus Holzer <markus.holzer@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/cell/CellInterval.h"

#include "domain_decomposition/StructuredBlockStorage.h"

#include "python_coupling/PythonWrapper.h"
namespace py = pybind11;
namespace walberla
{
namespace python_coupling
{
namespace internal
{
inline cell_idx_t normalizeIdx(py::object pyIndex, uint_t coordinateSize)
{
   cell_idx_t index;

   try{
       index = pyIndex.cast<cell_idx_t>();
   }
   catch (std::exception &){
       try {
           auto test = pyIndex.cast<real_t>();
           index = cell_idx_c( test * real_t( coordinateSize) );
       }
       catch (std::exception &) {
           throw py::cast_error("Incompatible index data type");
       }
   }

   if (index < 0)
      return cell_idx_c(coordinateSize) + 1 + index;
   else
      return index;
}

} // namespace internal

//*******************************************************************************************************************
/*! Creates a CellInterval as subset from the complete domain-cell-bounding-box based on a Python slice
 *
 *     Example: Python Slice: [ :, 3, -1 ]  and a domain size of ( 3,4,5 )
 *                 - x coordinate is the complete valid x-range indicated by the semicolon: i.e. [0,3)
 *                 - y coordinate is just a normal index i.e. the range from [3,4)
 *                 - z coordiante is the first valid coordinate from the back [4,5)
 *
 *     Python slices are tuples with slice classes as entry. Each slice has start, stop and step.
 *     Steps are not supported since they can not be encoded in a CellInterval
 */
//*******************************************************************************************************************
inline CellInterval globalPythonSliceToCellInterval(const shared_ptr< StructuredBlockStorage >& blocks,
                                                    py::tuple indexTuple)
{
   using internal::normalizeIdx;

   CellInterval bounds = blocks->getDomainCellBB();

   if (len(indexTuple) != 3)
   {
      throw py::index_error("Slice needs three components");
   }

   CellInterval interval;
   for (uint_t i = 0; i < 3; ++i)
   {
      if (!py::isinstance< py::slice >(indexTuple[i]))
      {
         cell_idx_t idx    = normalizeIdx(indexTuple[i], uint_c(bounds.max()[i]));
         interval.min()[i] = idx;
         interval.max()[i] = idx;
      }
      else if (py::isinstance< py::slice >(indexTuple[i]))
      {
         py::slice s = py::cast< py::slice >(indexTuple[i]);
         // Min
         if ( py::isinstance< py::none >(s.attr("start")) )
            interval.min()[i] = bounds.min()[i];
         else
            interval.min()[i] = normalizeIdx( s.attr("start"), uint_c( bounds.min()[i] ) );

         // Max
         if ( py::isinstance< py::none >(s.attr("stop")) )
            interval.max()[i] = bounds.max()[i];
         else
            interval.max()[i] = normalizeIdx( s.attr("stop"), uint_c( bounds.max()[i] ) );
      }
   }
   return interval;
}

} // namespace python_coupling
} // namespace walberla
