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
//! \file ChainedFilter.h
//! \ingroup vtk
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/DataTypes.h"
#include "core/cell/CellSet.h"
#include "domain_decomposition/StructuredBlockStorage.h"

#include <functional>


namespace walberla {
namespace vtk {



class ChainedFilter {

public:

   using CellFilter = std::function<void (CellSet &, const IBlock &, const StructuredBlockStorage &, const uint_t)>;

   void addFilter( const CellFilter& filter ) { filters_.push_back( filter ); }

   void operator()( CellSet& filteredCells, const IBlock& block, const StructuredBlockStorage& storage, const uint_t ghostLayers = uint_t(0) ) const
   {
      if( !filters_.empty() )
         filters_[0]( filteredCells, block, storage, ghostLayers );
      for( uint_t i = 1; i < filters_.size(); ++i )
      {
         CellSet cells;
         filters_[i]( cells, block, storage, ghostLayers );
         filteredCells &= cells;
      }
      // result (in terms of set theory): intersection of all filters
   }

private:

   std::vector< CellFilter > filters_;

}; // class ChainedFilter



} // namespace vtk
} // namespace walberla


