//==============================================================================================================================================================
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
//! \file GhostRegions.h
//! \ingroup field
//! \author Paulo Carvalho <prcjunior@inf.ufpr.br>
//! \author Martin Bauer <martin.bauer@fau.de>
//
//==============================================================================================================================================================

#pragma once

#include "core/cell/CellInterval.h"

namespace walberla {
namespace field {


//*******************************************************************************************************************
/*!\brief Constructs CellInterval containing the ghost region in the specified direction
 *
 * \param f   field
 * \param d   direction of the ghost layer  For W,E,N,S,T,B   a slice is returned
 *                                          for NW, NE, ..    an edge is returned
 *                                          for TBE, TNW, ... a corner ( single cell ) is returned
 * \param thickness how many ghost layer to return, if thickness < nrOfGhostLayers() the ghostlayers
 *                  nearest to the domain are returned
 * \param fullSlice  if true also the ghost cells in directions orthogonal to d are contained in the
 *                   returned cell interval. Example (d=W ): if fullSlice then also the ghost cells in N-S and T-B
 *                   are included.
 * \return CellInterval describing the ghost layer
 *******************************************************************************************************************/
template< typename GhostLayerField_T>
CellInterval getGhostRegion( const GhostLayerField_T & f, stencil::Direction d,
                             cell_idx_t thickness, bool fullSlice )
{
   const cell_idx_t sizeArr [] = { cell_idx_c( f.xSize() ),
                                   cell_idx_c( f.ySize() ),
                                   cell_idx_c( f.zSize() )};

   WALBERLA_ASSERT_GREATER( thickness, 0 );
   WALBERLA_ASSERT_LESS_EQUAL( uint_c(thickness), f.nrOfGhostLayers() );
   const cell_idx_t ghosts = cell_idx_c ( thickness );

   cell_idx_t fullSliceInc = fullSlice ? cell_idx_c( f.nrOfGhostLayers() ) : 0;

   CellInterval ci;
   for( uint_t dim = 0; dim< 3; ++dim )
      switch ( stencil::c[dim][d] )
      {
         case -1: ci.min()[dim] =     -ghosts;     ci.max()[dim] =         0                   - 1; break;
         case  0: ci.min()[dim] = -fullSliceInc;   ci.max()[dim] =  sizeArr[dim]+fullSliceInc  - 1; break;
         case  1: ci.min()[dim] =   sizeArr[dim];  ci.max()[dim] =  sizeArr[dim]+ghosts        - 1; break;
      }
   return ci;
}



//*******************************************************************************************************************
/*!\brief Constructs CellInterval containing the last slice before the ghost layer begins
 *        for a given direction.
 *
 * \param f   field
 * \param d   direction of the border .     For W,E,N,S,T,B   a slice is returned
 *                                          for NW, NE, ..    an edge is returned
 *                                          for TBE, TNW, ... a corner ( single cell ) is returned
 * \param thickness  how many slices to return
 * \param fullSlice  if true also the ghost cells in directions orthogonal to d are contained in the \
 *                   returned cell interval. Example (d=W ): if fullSlice then also the ghost layer in N-S and T-B
 *                   are included, otherwise only inner cells are returned
 * \return CellInterval describing the slice before ghost layer
 *******************************************************************************************************************/
template< typename GhostLayerField_T>
CellInterval getSliceBeforeGhostLayer(const GhostLayerField_T & f, stencil::Direction d,
                                      cell_idx_t thickness, bool fullSlice )
{
   WALBERLA_ASSERT_GREATER( thickness, 0 );

   const cell_idx_t sizeArr [] = { cell_idx_c( f.xSize() ),
                                   cell_idx_c( f.ySize() ),
                                   cell_idx_c( f.zSize() )};


   CellInterval ci;
   cell_idx_t fullSliceInc = fullSlice ? cell_idx_c( f.nrOfGhostLayers()) : 0;
   for( uint_t dim = 0; dim< 3; ++dim )
      switch ( stencil::c[dim][d] )
      {
         case -1: ci.min()[dim] =                      0;  ci.max()[dim] =     thickness              - 1; break;
         case  0: ci.min()[dim] =          -fullSliceInc;  ci.max()[dim] =  sizeArr[dim] +fullSliceInc- 1; break;
         case  1: ci.min()[dim] = sizeArr[dim]-thickness;  ci.max()[dim] =  sizeArr[dim]              - 1; break;
      }
   return ci;
}


} // namespace field
} // namespace walberla
