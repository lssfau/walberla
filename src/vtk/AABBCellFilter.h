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
//! \file AABBCellFilter.h
//! \ingroup vtk
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/DataTypes.h"
#include "core/cell/CellSet.h"
#include "domain_decomposition/StructuredBlockStorage.h"


namespace walberla {
namespace vtk {



//**********************************************************************************************************************
/*!
*   \brief Implementation of a cell filter based on an axis-aligned bounding box
*
*   All cells that are included in or intersect with the given axis-aligned bounding box are selected by this filter.
*/
//**********************************************************************************************************************

class AABBCellFilter {

public:

   AABBCellFilter( const AABB & aabb ) : aabb_( aabb ) {}

   //*******************************************************************************************************************
   /*!
   *   All cells (in block local cell coordinates!) that are included in or intersect with the axis-aligned bounding
   *   box are inserted into the set "filteredCells", thereby marking these cells. For more details on cell filters
   *   and including cells in and excluding cells from the VTK output, see the documentation of class VTKOutput.
   */
   //*******************************************************************************************************************
   void operator()( CellSet& filteredCells, const IBlock& block, const StructuredBlockStorage& storage, const uint_t ghostLayers = uint_t(0) ) const
   {
      CellInterval cellBB;
      storage.getCellBBFromAABB( cellBB, aabb_, storage.getLevel(block) );
      storage.transformGlobalToBlockLocalCellInterval( cellBB, block );

      const cell_idx_t start = cell_idx_c(-1) * cell_idx_c(ghostLayers);

      for( cell_idx_t z = start; z != cell_idx_c( storage.getNumberOfZCells( block ) ) + cell_idx_c(ghostLayers); ++z )
         for( cell_idx_t y = start; y != cell_idx_c( storage.getNumberOfYCells( block ) ) + cell_idx_c(ghostLayers); ++y )
            for( cell_idx_t x = start; x != cell_idx_c( storage.getNumberOfXCells( block ) ) + cell_idx_c(ghostLayers); ++x )
               if( cellBB.contains(x,y,z) ) filteredCells.insert(x,y,z);
   }

private:

   const AABB aabb_;

}; // class AABBCellFilter



} // namespace vtk
} // namespace walberla

