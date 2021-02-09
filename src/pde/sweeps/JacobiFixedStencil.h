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
//! \file JacobiFixedStencil.h
//! \ingroup pde
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "StencilSweepBase.h"
#include "field/iterators/IteratorMacros.h"
#include "stencil/Directions.h"



namespace walberla {
namespace pde {



template< typename Stencil_T >
class JacobiFixedStencil : public StencilSweepBase< Stencil_T >
{
public:

   typedef typename StencilSweepBase< Stencil_T >::Field_T Field_T;

   // block has NO dst u field
   JacobiFixedStencil( const BlockDataID & uFieldId, const BlockDataID & fFieldId, const std::vector< real_t > & weights ) :
      StencilSweepBase< Stencil_T >( uFieldId, fFieldId, weights ) {}

   // every block has a dedicated dst u field
   JacobiFixedStencil( const BlockDataID & src, const BlockDataID & dst, const BlockDataID & fFieldId, const std::vector< real_t > & weights ) :
      StencilSweepBase< Stencil_T >( src, dst, fFieldId, weights ) {}

   void operator()( IBlock * const block );
};



template< typename Stencil_T >
void JacobiFixedStencil< Stencil_T >::operator()( IBlock * const block )
{
   Field_T * sf( nullptr );
   Field_T * df( nullptr );
   Field_T * ff( nullptr );
   this->getFields( block, sf, df, ff );

   WALBERLA_ASSERT_GREATER_EQUAL( sf->nrOfGhostLayers(), 1 );

   // stencil weights
   real_t weights[ Stencil_T::Size ];
   for( auto dir = Stencil_T::beginNoCenter(); dir != Stencil_T::end(); ++dir )
      weights[ dir.toIdx() ] = this->w( dir.toIdx() );
   weights[ Stencil_T::idx[ stencil::C ] ] = real_t(1) / this->w( Stencil_T::idx[ stencil::C ] ); // center already inverted here!
 
   WALBERLA_FOR_ALL_CELLS_XYZ( sf,
   
      df->get(x,y,z) = ff->get(x,y,z);

      for( auto dir = Stencil_T::beginNoCenter(); dir != Stencil_T::end(); ++dir )
         df->get(x,y,z) -= weights[ dir.toIdx() ] * sf->getNeighbor(x,y,z,*dir);

      df->get(x,y,z) *= weights[ Stencil_T::idx[ stencil::C ] ];   
   
   ) // WALBERLA_FOR_ALL_CELLS_XYZ

   sf->swapDataPointers( df );
}



} // namespace pde
} // namespace walberla
