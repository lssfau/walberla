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
//! \file Neumann.h
//! \ingroup pde
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/cell/CellInterval.h"
#include "domain_decomposition/StructuredBlockStorage.h"
#include "stencil/D3Q6.h"



namespace walberla {
namespace pde {



//**********************************************************************************************************************
/*!
*   \brief Class for setting Neumann domain boundaries for 1st and 2nd derivatives
*
*   ATTENTION: only works for cell-based setups where domain boundaries are located exactly in-between two cells !
*/
//**********************************************************************************************************************

template< typename PdeField >
class NeumannDomainBoundary
{
public:

   NeumannDomainBoundary( StructuredBlockStorage & blocks, const BlockDataID & fieldId ) :
      blocks_( blocks ), fieldId_( fieldId )
   {
      for( uint_t i = 0; i != stencil::D3Q6::Size; ++i )
      {
         includeBoundary_[i] = true;
         order_[i] = uint_t(1);
         value_[i] = real_t(0);
      }
      dx_[ stencil::D3Q6::idx[ stencil::W ] ] = blocks.dx();
      dx_[ stencil::D3Q6::idx[ stencil::E ] ] = blocks.dx();
      dx_[ stencil::D3Q6::idx[ stencil::S ] ] = blocks.dy();
      dx_[ stencil::D3Q6::idx[ stencil::N ] ] = blocks.dy();
      dx_[ stencil::D3Q6::idx[ stencil::B ] ] = blocks.dz();
      dx_[ stencil::D3Q6::idx[ stencil::T ] ] = blocks.dz();
   }

   void includeBoundary( const stencil::Direction & direction ) { includeBoundary_[stencil::D3Q6::idx[direction]] = true; }
   void excludeBoundary( const stencil::Direction & direction ) { includeBoundary_[stencil::D3Q6::idx[direction]] = false; }

   void setOrder( const uint_t order ) { for( uint_t i = 0; i != stencil::D3Q6::Size; ++i ) order_[i] = order; }
   void setOrder( const stencil::Direction & direction, const uint_t order ) { order_[stencil::D3Q6::idx[direction]] = order; }

   void setValue( const real_t value ) { for( uint_t i = 0; i != stencil::D3Q6::Size; ++i ) value_[i] = value; }
   void setValue( const stencil::Direction & direction, const real_t value ) { value_[stencil::D3Q6::idx[direction]] = value; }

   void setDx( const real_t dx ) { for( uint_t i = 0; i != stencil::D3Q6::Size; ++i ) dx_[i] = dx; }
   void setDx( const stencil::Direction & direction, const real_t dx ) { dx_[stencil::D3Q6::idx[direction]] = dx; }

   void operator()();

protected:

   void apply( PdeField * p, const CellInterval & interval, const cell_idx_t cx, const cell_idx_t cy, const cell_idx_t cz,
               const uint_t order, const real_t value, const real_t dx ) const;



   StructuredBlockStorage & blocks_;
   BlockDataID fieldId_;

   bool includeBoundary_[ stencil::D3Q6::Size ];
   uint_t order_[ stencil::D3Q6::Size ];

   real_t value_[ stencil::D3Q6::Size ];
   real_t dx_[ stencil::D3Q6::Size ];

}; // class NeumannDomainBoundary



template< typename PdeField >
void NeumannDomainBoundary< PdeField >::operator()()
{
   for( auto block = blocks_.begin(); block != blocks_.end(); ++block )
   {
      PdeField * p = block->template getData< PdeField >( fieldId_ );

      if( includeBoundary_[ stencil::D3Q6::idx[ stencil::W ] ] && blocks_.atDomainXMinBorder( *block ) )
      {
         apply( p, CellInterval( cell_idx_t(-1), cell_idx_t(0), cell_idx_t(0),
                                 cell_idx_t(-1), cell_idx_c(p->ySize()) - cell_idx_t(1), cell_idx_c(p->zSize()) - cell_idx_t(1) ),
                cell_idx_t(1), cell_idx_t(0), cell_idx_t(0),
                order_[ stencil::D3Q6::idx[ stencil::W ] ], value_[ stencil::D3Q6::idx[ stencil::W ] ], dx_[ stencil::D3Q6::idx[ stencil::W ] ] );
      }
      if( includeBoundary_[ stencil::D3Q6::idx[ stencil::E ] ] && blocks_.atDomainXMaxBorder( *block ) )
      {
         apply( p, CellInterval( cell_idx_c(p->xSize()), cell_idx_t(0), cell_idx_t(0),
                                 cell_idx_c(p->xSize()), cell_idx_c(p->ySize()) - cell_idx_t(1), cell_idx_c(p->zSize()) - cell_idx_t(1) ),
                cell_idx_t(-1), cell_idx_t(0), cell_idx_t(0),
                order_[ stencil::D3Q6::idx[ stencil::E ] ], value_[ stencil::D3Q6::idx[ stencil::E ] ], dx_[ stencil::D3Q6::idx[ stencil::E ] ] );
      }

      if( includeBoundary_[ stencil::D3Q6::idx[ stencil::S ] ] && blocks_.atDomainYMinBorder( *block ) )
      {
         apply( p, CellInterval( cell_idx_t(0), cell_idx_t(-1), cell_idx_t(0),
                                 cell_idx_c(p->xSize()) - cell_idx_t(1), cell_idx_t(-1), cell_idx_c(p->zSize()) - cell_idx_t(1) ),
                cell_idx_t(0), cell_idx_t(1), cell_idx_t(0),
                order_[ stencil::D3Q6::idx[ stencil::S ] ], value_[ stencil::D3Q6::idx[ stencil::S ] ], dx_[ stencil::D3Q6::idx[ stencil::S ] ] );
      }
      if( includeBoundary_[ stencil::D3Q6::idx[ stencil::N ] ] && blocks_.atDomainYMaxBorder( *block ) )
      {
         apply( p, CellInterval( cell_idx_t(0), cell_idx_c(p->ySize()), cell_idx_t(0),
                                 cell_idx_c(p->xSize()) - cell_idx_t(1), cell_idx_c(p->ySize()), cell_idx_c(p->zSize()) - cell_idx_t(1) ),
                cell_idx_t(0), cell_idx_t(-1), cell_idx_t(0),
                order_[ stencil::D3Q6::idx[ stencil::N ] ], value_[ stencil::D3Q6::idx[ stencil::N ] ], dx_[ stencil::D3Q6::idx[ stencil::N ] ] );
      }

      if( includeBoundary_[ stencil::D3Q6::idx[ stencil::B ] ] && blocks_.atDomainZMinBorder( *block ) )
      {
         apply( p, CellInterval( cell_idx_t(0), cell_idx_t(0), cell_idx_t(-1),
                                 cell_idx_c(p->xSize()) - cell_idx_t(1), cell_idx_c(p->ySize()) - cell_idx_t(1), cell_idx_t(-1) ),
                cell_idx_t(0), cell_idx_t(0), cell_idx_t(1),
                order_[ stencil::D3Q6::idx[ stencil::B ] ], value_[ stencil::D3Q6::idx[ stencil::B ] ], dx_[ stencil::D3Q6::idx[ stencil::B ] ] );
      }
      if( includeBoundary_[ stencil::D3Q6::idx[ stencil::T ] ] && blocks_.atDomainZMaxBorder( *block ) )
      {
         apply( p, CellInterval( cell_idx_t(0), cell_idx_t(0), cell_idx_c(p->zSize()),
                                 cell_idx_c(p->xSize()) - cell_idx_t(1), cell_idx_c(p->ySize()) - cell_idx_t(1), cell_idx_c(p->zSize()) ),
                cell_idx_t(0), cell_idx_t(0), cell_idx_t(-1),
                order_[ stencil::D3Q6::idx[ stencil::T ] ], value_[ stencil::D3Q6::idx[ stencil::T ] ], dx_[ stencil::D3Q6::idx[ stencil::T ] ] );
      }
   }
}



template< typename PdeField >
void NeumannDomainBoundary< PdeField >::apply( PdeField * p, const CellInterval & interval,
                                               const cell_idx_t cx, const cell_idx_t cy, const cell_idx_t cz,
                                               const uint_t order, const real_t value, const real_t dx ) const
{
   if( order == uint_t(1) )
   {
      if( isIdentical( value, real_t(0) ) )
      {
         WALBERLA_FOR_ALL_CELLS_IN_INTERVAL_XYZ( interval,
            p->get(x,y,z) = p->get( x + cx, y + cy, z + cz );  // (dp / dx) == 0 _on_ the boundary
         )
      }
      else
      {
         const real_t vdx = value * dx;
         WALBERLA_FOR_ALL_CELLS_IN_INTERVAL_XYZ( interval,
            p->get(x,y,z) = vdx + p->get( x + cx, y + cy, z + cz );  // (dp / dx) == value _on_ the boundary
         )
      }
   }
   else
   {
      WALBERLA_ASSERT_EQUAL( order, uint_t(2) );

      if( isIdentical( value, real_t(0) ) )
      {      
         WALBERLA_FOR_ALL_CELLS_IN_INTERVAL_XYZ( interval,

            const real_t pBoundary = p->get( x + cx, y + cy, z + cz );
            const real_t pInner    = p->get( x + cell_idx_t(2) * cx, y + cell_idx_t(2) * cy, z + cell_idx_t(2) * cz );

            const real_t boundaryValue = pBoundary + real_c(0.5) * ( pBoundary - pInner ); // extrapolation of value _on_ the boundary

            p->get(x,y,z) = real_t(2) * boundaryValue - pBoundary;  // (d^2 p / dx^2) == 0 _on_ the boundary
         )
      }
      else
      {
         const real_t vdx = value * real_c(0.25) * dx;
         WALBERLA_FOR_ALL_CELLS_IN_INTERVAL_XYZ( interval,

            const real_t pBoundary = p->get( x + cx, y + cy, z + cz );
            const real_t pInner    = p->get( x + cell_idx_t(2) * cx, y + cell_idx_t(2) * cy, z + cell_idx_t(2) * cz );

            const real_t boundaryValue = pBoundary + real_c(0.5) * ( pBoundary - pInner ); // extrapolation of value _on_ the boundary

            p->get(x,y,z) = vdx + real_t(2) * boundaryValue - pBoundary;  // (d^2 p / dx^2) == value _on_ the boundary
         )
      }
   }
}



} // namespace pde
} // namespace walberla
