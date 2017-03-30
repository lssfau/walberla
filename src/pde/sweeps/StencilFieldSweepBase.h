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
//! \file StencilFieldSweepBase.h
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
class StencilFieldSweepBase : public SweepBase
{
public:

   typedef SweepBase::Field_T                                 Field_T;
   typedef GhostLayerField< real_t, Stencil_T::Size >  StencilField_T;

   // block has NO dst u field
   StencilFieldSweepBase( const BlockDataID & uFieldId, const BlockDataID & fFieldId, const BlockDataID & stencilFieldId ) :
      SweepBase( uFieldId, fFieldId ), stencil_( stencilFieldId ) {}

   // every block has a dedicated dst u field
   StencilFieldSweepBase( const BlockDataID & src, const BlockDataID & dst, const BlockDataID & fFieldId, const BlockDataID & stencilFieldId ) :
      SweepBase( src, dst, fFieldId ), stencil_( stencilFieldId ) {}

protected:

   inline StencilField_T * getStencilField( IBlock * const block ) const;

   using SweepBase::getFields;
   inline void getFields( IBlock * const block, Field_T * & u,                    Field_T * & f, StencilField_T * & stencil );
   inline void getFields( IBlock * const block, Field_T * & src, Field_T * & dst, Field_T * & f, StencilField_T * & stencil );



   const BlockDataID stencil_;
};



template< typename Stencil_T >
inline typename StencilFieldSweepBase<Stencil_T>::StencilField_T * StencilFieldSweepBase<Stencil_T>::getStencilField( IBlock * const block ) const
{
   WALBERLA_ASSERT_NOT_NULLPTR( block );
   
   StencilField_T * stencil = block->getData< StencilField_T >( stencil_ );

   WALBERLA_ASSERT_NOT_NULLPTR( stencil );
   
   return stencil;
}



template< typename Stencil_T >
inline void StencilFieldSweepBase<Stencil_T>::getFields( IBlock * const block, Field_T * & u, Field_T * & f, StencilField_T * & stencil )
{
   WALBERLA_ASSERT_NOT_NULLPTR( block );

   this->getFields( block, u, f );
   stencil = getStencilField( block );

   WALBERLA_ASSERT_NOT_NULLPTR( stencil );

   WALBERLA_ASSERT_EQUAL( u->xyzSize(), stencil->xyzSize() );
}



template< typename Stencil_T >
inline void StencilFieldSweepBase<Stencil_T>::getFields( IBlock * const block, Field_T * & src, Field_T * & dst, Field_T * & f,
                                                         StencilField_T * & stencil )
{
   WALBERLA_ASSERT_NOT_NULLPTR( block );

   this->getFields( block, src, dst, f );
   stencil = getStencilField( block );

   WALBERLA_ASSERT_NOT_NULLPTR( stencil );

   WALBERLA_ASSERT_EQUAL( src->xyzSize(), stencil->xyzSize() );
}



} // namespace pde
} // namespace walberla
