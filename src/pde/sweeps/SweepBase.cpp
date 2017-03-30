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
//! \file SweepBase.cpp
//! \ingroup pde
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#include "field/iterators/IteratorMacros.h"

#include "SweepBase.h"



namespace walberla {
namespace pde {



SweepBase::Field_T * SweepBase::getDstField( IBlock * const block, Field_T * const src )
{
   WALBERLA_ASSERT_NOT_NULLPTR( block );
   WALBERLA_ASSERT_NOT_NULLPTR( src );

   if( dstFromBlockData_ )
   {
      Field_T * dst = block->getData< Field_T >( dst_ );
      WALBERLA_ASSERT_NOT_NULLPTR( dst );
      return dst;
   }

   auto it = dstFields_.find( src );
   if( it != dstFields_.end() )
   {
#ifndef NDEBUG
      std::fill( (*it)->beginWithGhostLayer(), (*it)->end(), std::numeric_limits< Field_T::value_type >::quiet_NaN() );
#endif
      WALBERLA_ASSERT_NOT_NULLPTR( *it );
      return *it;
   }

   Field_T * dst = src->cloneUninitialized();
   WALBERLA_ASSERT_NOT_NULLPTR( dst );
   
   // take care of proper thread<->memory assignment (first-touch allocation policy !)
   WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ( dst,
      dst->get(x,y,z) = std::numeric_limits< Field_T::value_type >::quiet_NaN();
   )
   dstFields_.insert( dst );

   return dst;
}



} // namespace pde
} // namespace walberla
