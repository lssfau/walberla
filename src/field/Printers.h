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
//! \file Printers.h
//! \ingroup field
//! \author Martin Bauer <martin.bauer@fau.de>
//! \brief Console output for Fields
//
//======================================================================================================================

#pragma once

#include "Field.h"
#include "FlagField.h"
#include "GhostLayerField.h"
#include "core/DataTypes.h"

#include <iomanip>
#include <iostream>


namespace walberla {
namespace field {

   /****************************************************************************************************************//**
    *  Prints a slice along a coordinate axis
    *
    * \param os          output stream
    * \param field       the field to print
    * \param sliceCoord  dimension to slice through, 0=x, 1=y, 2=z
    * \param sliceValue  fixed value of the sliceCoordinate
    * \param f           fixed f value
   ********************************************************************************************************************/
   template<typename T, uint_t fs>
   std::ostream &  printSlice( std::ostream & os,
                               const Field<T,fs> & field,
                               int sliceCoord,
                               cell_idx_t sliceValue,
                               cell_idx_t f=0 );


   /****************************************************************************************************************//**
   * Overload of printSlice for GhostLayerFields
   ********************************************************************************************************************/
   template<typename T, uint_t fs>
   std::ostream &  printSlice( std::ostream & os,
                               const GhostLayerField<T,fs> & field,
                               int sliceCoord,
                               cell_idx_t sliceValue,
                               cell_idx_t f=0 );



   /****************************************************************************************************************//**
   * Overload of printSlice for FlagFields
   ********************************************************************************************************************/
   template<typename T>
   std::ostream &  printSlice( std::ostream & os,
                               const FlagField<T> & field,
                               int sliceCoord,
                               cell_idx_t sliceValue,
                               cell_idx_t f=0 );


} // namespace field
} // namespace walberla


#include "Printers.impl.h"

