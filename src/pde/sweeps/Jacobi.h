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
//! \file Jacobi.h
//! \ingroup pde
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "StencilFieldSweepBase.h"
#include "stencil/Directions.h"



namespace walberla {
namespace pde {



template< typename Stencil_T >
class Jacobi : public StencilFieldSweepBase< Stencil_T >
{
public:

   typedef typename StencilFieldSweepBase< Stencil_T >::Field_T         Field_T;
   typedef typename StencilFieldSweepBase< Stencil_T >::StencilField_T  StencilField_T;

   // block has NO dst u field
   Jacobi( const BlockDataID & uFieldId, const BlockDataID & fFieldId, const BlockDataID & stencilFieldId ) :
      StencilFieldSweepBase< Stencil_T >( uFieldId, fFieldId, stencilFieldId ) {}

   // every block has a dedicated dst u field
   Jacobi( const BlockDataID & src, const BlockDataID & dst, const BlockDataID & fFieldId, const BlockDataID & stencilFieldId ) :
      StencilFieldSweepBase< Stencil_T >( src, dst, fFieldId, stencilFieldId ) {}

   void operator()( IBlock * const block );
};



template< typename Stencil_T >
void Jacobi< Stencil_T >::operator()( IBlock * const block )
{
   Field_T * sf( nullptr );
   Field_T * df( nullptr );
   Field_T * ff( nullptr );
   StencilField_T * stencil( nullptr );
   this->getFields( block, sf, df, ff, stencil );

   WALBERLA_ASSERT_GREATER_EQUAL( sf->nrOfGhostLayers(), 1 );
   
   WALBERLA_FOR_ALL_CELLS_XYZ( sf,

      df->get(x,y,z) = ff->get(x,y,z);

      for( auto dir = Stencil_T::beginNoCenter(); dir != Stencil_T::end(); ++dir )
         df->get(x,y,z) -= stencil->get( x, y, z, dir.toIdx() ) * sf->getNeighbor( x, y, z, *dir );

      df->get(x,y,z) /= stencil->get( x, y, z, Stencil_T::idx[stencil::C] );

   ) // WALBERLA_FOR_ALL_CELLS_XYZ

   sf->swapDataPointers( df );
}



} // namespace pde
} // namespace walberla
