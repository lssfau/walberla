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
//! \file SweepBase.h
//! \ingroup pde
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/debug/Debug.h"
#include "domain_decomposition/IBlock.h"
#include "field/GhostLayerField.h"
#include "field/SwapableCompare.h"

#include <set>



namespace walberla {
namespace pde {



class SweepBase
{
public:

   using Field_T = GhostLayerField<real_t, 1>;

   // block has NO dst u field
   SweepBase( const BlockDataID & uFieldId, const BlockDataID & fFieldId ) :
      src_( uFieldId ), f_( fFieldId ), dstFromBlockData_( false ) {}

   // every block has a dedicated dst u field
   SweepBase( const BlockDataID & src, const BlockDataID & dst, const BlockDataID & fFieldId ) :
      src_( src ), f_( fFieldId ), dstFromBlockData_( true ), dst_( dst ) {}

   virtual ~SweepBase() { for( auto field = dstFields_.begin(); field != dstFields_.end(); ++field ) delete *field; }

protected:

   inline Field_T * getSrcField( IBlock * const block ) const;
          Field_T * getDstField( IBlock * const block, Field_T * const src );

   inline Field_T * getUField( IBlock * const block ) const { return getSrcField( block ); }
   inline Field_T * getFField( IBlock * const block ) const;

   inline void getFields( IBlock * const block, Field_T * & u,                    Field_T * & f );
   inline void getFields( IBlock * const block, Field_T * & src, Field_T * & dst, Field_T * & f );



   const BlockDataID src_ {}; // "u"
   const BlockDataID f_ {};

   const bool dstFromBlockData_;
   const BlockDataID dst_ {};
   std::set< Field_T *, field::SwapableCompare< Field_T * > > dstFields_;
};



inline SweepBase::Field_T * SweepBase::getSrcField( IBlock * const block ) const
{
   WALBERLA_ASSERT_NOT_NULLPTR( block );

   Field_T * src = block->getData< Field_T >( src_ );
   
   WALBERLA_ASSERT_NOT_NULLPTR( src );
   
   return src;
}



inline SweepBase::Field_T * SweepBase::getFField( IBlock * const block ) const
{
   WALBERLA_ASSERT_NOT_NULLPTR( block );
   
   Field_T * f =  block->getData< Field_T >( f_ );

   WALBERLA_ASSERT_NOT_NULLPTR( f );
   
   return f;
}



inline void SweepBase::getFields( IBlock * const block, Field_T * & u, Field_T * & f )
{
   WALBERLA_ASSERT_NOT_NULLPTR( block );

   u = getSrcField( block );
   f = getFField( block );

   WALBERLA_ASSERT_NOT_NULLPTR( u );
   WALBERLA_ASSERT_NOT_NULLPTR( f );

   WALBERLA_ASSERT_EQUAL( u->xyzSize(), f->xyzSize() );
}



inline void SweepBase::getFields( IBlock * const block, Field_T * & src, Field_T * & dst, Field_T * & f )
{
   WALBERLA_ASSERT_NOT_NULLPTR( block );

   src = getSrcField( block );
   dst = getDstField( block, src );
   f   = getFField( block );

   WALBERLA_ASSERT_NOT_NULLPTR( src );
   WALBERLA_ASSERT_NOT_NULLPTR( dst );
   WALBERLA_ASSERT_NOT_NULLPTR( f );

   WALBERLA_ASSERT_EQUAL( src->xyzSize(), dst->xyzSize() );
   WALBERLA_ASSERT_EQUAL( src->xyzSize(), f->xyzSize() );
}



} // namespace pde
} // namespace walberla
