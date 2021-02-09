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
//! \file SOR.h
//! \ingroup pde
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "StencilFieldSweepBase.h"
#include "core/math/Uint.h"
#include "stencil/Directions.h"

#include <map>
#include <functional>


namespace walberla {
namespace pde {



template< typename Stencil_T >
class SOR : public StencilFieldSweepBase< Stencil_T >
{
public:

   typedef typename StencilFieldSweepBase< Stencil_T >::Field_T         Field_T;
   typedef typename StencilFieldSweepBase< Stencil_T >::StencilField_T  StencilField_T;

   SOR( const shared_ptr< domain_decomposition::StructuredBlockStorage > & blocks,
                    const BlockDataID & uFieldId, const BlockDataID & fFieldId, const BlockDataID & stencilFieldId, const real_t omega ) :
      StencilFieldSweepBase< Stencil_T >( uFieldId, fFieldId, stencilFieldId ), blocks_( blocks ), omega_( omega ) {}

   void operator()( IBlock * const block ) const { WALBERLA_ABORT( "You are not allowed to use class 'SOR' as a standard sweep!\n"
                                                                   "Use the member functions 'getRedSweep' and 'getBlackSweep' instead." ); }

   void update( IBlock * const block, const bool rb );

   std::function< void ( IBlock * const ) > getRedSweep()
   {
      return std::bind( &SOR::update, this, std::placeholders::_1, true );
   }

   std::function< void ( IBlock * const ) > getBlackSweep()
   {
      return std::bind( &SOR::update, this, std::placeholders::_1, false );
   }

private:

   shared_ptr< domain_decomposition::StructuredBlockStorage > blocks_;
   real_t omega_;
};



template< typename Stencil_T >
void SOR< Stencil_T >::update( IBlock * const block, const bool rb )
{
#ifndef NDEBUG
   for( auto dir = Stencil_T::beginNoCenter(); dir != Stencil_T::end(); ++dir )
      WALBERLA_ASSERT( realIsIdentical( dir.length(), real_t(1) ) );
#endif

   Field_T * uf( nullptr );
   Field_T * ff( nullptr );
   StencilField_T * stencil( nullptr );
   this->getFields( block, uf, ff, stencil );

   WALBERLA_ASSERT_GREATER_EQUAL( uf->nrOfGhostLayers(), 1 );

   const cell_idx_t zero = cell_idx_t(0);
   const cell_idx_t one  = cell_idx_t(1);

   const real_t omegaInv = real_t(1) - omega_;

   Cell cell(zero,zero,zero);
   blocks_->transformBlockLocalToGlobalCell( cell, *block );

   WALBERLA_FOR_ALL_CELLS_YZ( uf,

      Cell c( cell );
      c.y() += y;
      c.z() += z;

      const cell_idx_t xBegin = ( (((c.x() & one) + (c.y() & one) + (c.z() & one)) & one) == zero ) ? (rb ? zero : one) : (rb ? one : zero);

      const cell_idx_t xSize = cell_idx_c( uf->xSize() );
      for( cell_idx_t x = xBegin; x < xSize; x += cell_idx_t(2) )
      {
         real_t value = ff->get(x,y,z);

         for( auto dir = Stencil_T::beginNoCenter(); dir != Stencil_T::end(); ++dir )
            value -= stencil->get( x, y, z, dir.toIdx() ) * uf->getNeighbor(x,y,z,*dir);

         value /= stencil->get( x, y, z, Stencil_T::idx[stencil::C] );

         uf->get(x,y,z) = omegaInv * uf->get(x,y,z) + omega_ * value;
      }

   ) // WALBERLA_FOR_ALL_CELLS_YZ
}



} // namespace pde
} // namespace walberla
