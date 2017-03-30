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
//! \file CellBBCellFilter.h
//! \ingroup vtk
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/DataTypes.h"
#include "core/cell/CellInterval.h"
#include "core/cell/CellSet.h"
#include "domain_decomposition/StructuredBlockStorage.h"


namespace walberla {
namespace vtk {



//**********************************************************************************************************************
/*!
*   \brief Implementation of a cell filter based on a cell bounding box
*
*   All cells that are contained within the given cell bounding box are selected by this filter. The cell bounding box
*   is always given in global cell coordinates and only acts on blocks on level "level_" (default=0).
*/
//**********************************************************************************************************************

class CellBBCellFilter {

public:

   CellBBCellFilter( const CellInterval& cellBB, const uint_t level = 0 ) : cellBB_( cellBB ), level_( level ) {}

   //*******************************************************************************************************************
   /*!
   *   All cells (in block local cell coordinates!) that are contained within the given cell bounding box (this cell
   *   bounding box is given in global cell coordinates!) are inserted into the set "filteredCells", thereby marking
   *   these cells. For more details on cell filters and including cells in and excluding cells from the
   *   VTK output, see the documentation of class VTKOutput.
   *   Attention: If the specified level "level_" does not match with the level the block resides on, no cells are
   *              selected.
   */
   //*******************************************************************************************************************
   void operator()( CellSet& filteredCells, const IBlock& block, const StructuredBlockStorage& storage, const uint_t ghostLayers ) const
   {
      if( storage.getLevel(block) != level_ )
         return;

      CellInterval cellBB;
      storage.transformGlobalToBlockLocalCellInterval( cellBB, block, cellBB_ );

      const cell_idx_t start = cell_idx_c(-1) * cell_idx_c(ghostLayers);

      for( cell_idx_t z = start; z != cell_idx_c( storage.getNumberOfZCells( block ) ) + cell_idx_c(ghostLayers); ++z )
         for( cell_idx_t y = start; y != cell_idx_c( storage.getNumberOfYCells( block ) ) + cell_idx_c(ghostLayers); ++y )
            for( cell_idx_t x = start; x != cell_idx_c( storage.getNumberOfXCells( block ) ) + cell_idx_c(ghostLayers); ++x )
               if( cellBB.contains(x,y,z) ) filteredCells.insert(x,y,z);
   }

private:

   const CellInterval cellBB_;
   const uint_t level_;

}; // class CellBBCellFilter



} // namespace vtk
} // namespace walberla

