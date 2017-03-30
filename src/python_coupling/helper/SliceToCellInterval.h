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
//
//======================================================================================================================

#pragma once

#include "python_coupling/PythonWrapper.h"
#include "core/cell/CellInterval.h"
#include "domain_decomposition/StructuredBlockStorage.h"
#include <boost/python/slice.hpp>

namespace walberla {
namespace python_coupling {


   namespace internal
   {
      inline cell_idx_t normalizeIdx( boost::python::object pyIndex, uint_t coordinateSize )
      {
         using namespace boost::python;

         cell_idx_t index;
         if ( extract<cell_idx_t>( pyIndex ).check() )
            index = extract<cell_idx_t> ( pyIndex );
         else if ( extract<double>( pyIndex ).check() )
            index = cell_idx_c( double ( extract<double>( pyIndex ) ) * double( coordinateSize) );
         else {
            PyErr_SetString( PyExc_IndexError, "Incompatible index data type" );
            throw error_already_set();
         }

         if ( index < 0 )
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
   inline CellInterval globalPythonSliceToCellInterval( const shared_ptr<StructuredBlockStorage> & blocks,
                                                        boost::python::tuple indexTuple )
   {
      using namespace boost::python;
      using internal::normalizeIdx;

      CellInterval bounds = blocks->getDomainCellBB();

      if ( len(indexTuple) != 3 )
      {
         PyErr_SetString( PyExc_IndexError, "Slice needs three components" );
         throw error_already_set();
      }

      CellInterval interval;
      for( uint_t i=0; i<3; ++i )
      {
         if( ! extract< slice >(indexTuple[i]).check() )
         {
            cell_idx_t idx = normalizeIdx( indexTuple[i], uint_c( bounds.max()[i] ) );
            interval.min()[i] = idx;
            interval.max()[i] = idx;
         }
         else if ( extract< slice >(indexTuple[i]).check() )
         {
             slice s = extract< slice >(indexTuple[i]);

             // Min
             if ( s.start() == object() )
                interval.min()[i] = bounds.min()[i];
             else
                interval.min()[i] = normalizeIdx( s.start(), uint_c( bounds.max()[i] ) );

             // Max
             if ( s.stop() == object() )
                interval.max()[i] = bounds.max()[i];
             else
                interval.max()[i] = normalizeIdx( s.stop(), uint_c( bounds.max()[i] ) );

             if ( s.step() != object() ) {
                PyErr_SetString( PyExc_IndexError, "Steps in slice not supported." );
                throw error_already_set();
             }
         }
      }
      return interval;
   }



   //*******************************************************************************************************************
   /*! Creates a CellInterval based on a Python Slice as subset of a field
   *
   *   Similar to globalPythonSliceToCellInterval() with the following additional features:
   *     - slice may have a forth component: [ :, 3, -1, 'g' ] with the only valid entry 'g' for ghost layers
   *     - if this ghost layer marker is present, coordinate 0 addresses the outermost ghost layer, otherwise the
   *       first inner cell is addressed
   */
   //*******************************************************************************************************************
   template<typename Field_T>
   CellInterval localPythonSliceToCellInterval( const Field_T & field,
                                                boost::python::tuple indexTuple )
   {
      using namespace boost::python;
      using internal::normalizeIdx;

      bool withGhostLayer=false;

      if ( len(indexTuple) != 3 )
      {
         if ( len(indexTuple) == 4 )
         {
            std::string marker =  extract<std::string>( indexTuple[3]);
            if ( marker == std::string("g") )
               withGhostLayer = true;
            else
            {
               PyErr_SetString( PyExc_IndexError, "Unknown marker in slice" );
               throw error_already_set();
            }
         }
         else
         {
            PyErr_SetString( PyExc_IndexError, "Slice needs three components ( + optional ghost layer marker )" );
            throw error_already_set();
         }
      }

      cell_idx_t gl = cell_idx_c( field.nrOfGhostLayers() );

      CellInterval bounds;;
      if ( withGhostLayer )
      {
         bounds =  field.xyzSizeWithGhostLayer();
         bounds.shift( gl,gl,gl );
      }
      else
         bounds = field.xyzSize();


      CellInterval interval;

      for( uint_t i=0; i<3; ++i )
      {
         if( !  extract< slice >(indexTuple[i]).check()  )
         {
            interval.min()[i] = normalizeIdx( indexTuple[i], uint_c(bounds.max()[i]) );
            interval.max()[i] = normalizeIdx( indexTuple[i], uint_c(bounds.max()[i]) );
         }
         else
         {
             slice s = extract< slice >(indexTuple[i]);

             // Min
             if ( s.start() == object() )
                interval.min()[i] = bounds.min()[i];
             else
                interval.min()[i] = normalizeIdx( s.start(), uint_c(bounds.max()[i]) );

             // Max
             if ( s.stop() == object() )
                interval.max()[i] = bounds.max()[i];
             else
                interval.max()[i] = normalizeIdx(  s.stop(), uint_c(bounds.max()[i]) );

             if ( s.step() != object() ) {
                PyErr_SetString( PyExc_IndexError, "Steps in slice not supported." );
                throw error_already_set();
             }
         }
      }

      if ( withGhostLayer )
         interval.shift( -gl,-gl,-gl );

      // Range check
      if ( ! field.xyzAllocSize().contains( interval ) ) {
         PyErr_SetString( PyExc_IndexError, "Index out of bounds." );
         throw error_already_set();
      }

      return interval;
   }




} // namespace python_coupling
} // namespace walberla


