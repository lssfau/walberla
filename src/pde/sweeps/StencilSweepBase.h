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
//! \file StencilSweepBase.h
//! \ingroup pde
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "SweepBase.h"

#include <vector>



namespace walberla {
namespace pde {



template< typename Stencil_T >
class StencilSweepBase : public SweepBase
{
public:

   typedef SweepBase::Field_T Field_T;

   // block has NO dst u field
   StencilSweepBase( const BlockDataID & uFieldId, const BlockDataID & fFieldId, const std::vector< real_t > & weights ) :
      SweepBase( uFieldId, fFieldId )
   {
      WALBERLA_ASSERT_EQUAL( weights.size(), Stencil_T::Size );
      for( uint_t i = uint_t(0); i < Stencil_T::Size; ++i )
         w_[i] = weights[i];
   }

   // every block has a dedicated dst u field
   StencilSweepBase( const BlockDataID & src, const BlockDataID & dst, const BlockDataID & fFieldId, const std::vector< real_t > & weights ) :
      SweepBase( src, dst, fFieldId )
   {
      WALBERLA_ASSERT_EQUAL( weights.size(), Stencil_T::Size );
      for( uint_t i = uint_t(0); i < Stencil_T::Size; ++i )
         w_[i] = weights[i];
   }

protected:

   inline real_t w( const uint_t i ) const { WALBERLA_ASSERT_LESS( i, Stencil_T::Size ); return w_[i]; }



   real_t w_[ Stencil_T::Size ];
};



} // namespace pde
} // namespace walberla
