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
//! \file CellArray.h
//! \ingroup core
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "Cell.h"
#include "CellSet.h"
#include "CellVector.h"
#include "core/Array.h"


namespace walberla {
namespace cell {


/// An array of cells

class CellArray : public Array<Cell> {

public:

   CellArray() : Array<Cell>() {}

   CellArray( const CellVector& cells )
   {
      this->array_ = cells.empty() ? nullptr : new Cell[ cells.size() ];
      this->size_  = cells.size();

      for( uint_t i = 0; i != this->size_; ++i )
         this->array_[i] = cells[i];
   }

   CellArray( const CellSet& cells )
   {
      this->array_ = cells.empty() ? nullptr : new Cell[ cells.size() ];
      this->size_  = cells.size();

      uint_t i = 0;
      for( auto cell = cells.begin(); cell != cells.end(); ++cell, ++i )
         this->array_[i] = *cell;
   }

}; // class CellArray



} // namespace cell

using cell::CellArray;

} // namespace walberla
