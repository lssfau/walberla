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
//! \file FlagFieldSweepBase.h
//! \ingroup lbm
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "SweepBase.h"
#include "field/FlagField.h"


namespace walberla {
namespace lbm {



template< typename LatticeModel_T, typename FlagField_T >
class FlagFieldSweepBase : public SweepBase<LatticeModel_T>
{
public:

   using PdfField_T = typename SweepBase<LatticeModel_T>::PdfField_T;
   using flag_t = typename FlagField_T::flag_t;

   // block has NO dst pdf field, lbm mask consists of multiple flags
   FlagFieldSweepBase( const BlockDataID & pdfField, const ConstBlockDataID & flagField, const Set< FlagUID > & lbmMask ) :
      SweepBase<LatticeModel_T>( pdfField ), flagField_( flagField ), lbmMask_( lbmMask ) {}

   // every block has a dedicated dst pdf field, lbm mask consists of multiple flags
   FlagFieldSweepBase( const BlockDataID & src, const BlockDataID & dst, const ConstBlockDataID & flagField, const Set< FlagUID > & lbmMask ) :
      SweepBase<LatticeModel_T>( src, dst ), flagField_( flagField ), lbmMask_( lbmMask ) {}

protected:

   using SweepBase<LatticeModel_T>::getFields;
   inline void getFields( IBlock * const block, PdfField_T * & src, PdfField_T * & dst, const FlagField_T * & flags );
   inline void getFields( IBlock * const block, PdfField_T * & src,                     const FlagField_T * & flags );

   flag_t getLbmMaskAndFields( IBlock * const block, PdfField_T * & src, PdfField_T * & dst, const FlagField_T * & flags )
   {
      getFields( block, src, dst, flags );
      return flags->getMask( lbmMask_ );
   }

   flag_t getLbmMaskAndFields( IBlock * const block, PdfField_T * & src, const FlagField_T * & flags )
   {
      getFields( block, src, flags );
      return flags->getMask( lbmMask_ );
   }



   const ConstBlockDataID flagField_;

   const Set< FlagUID > lbmMask_;
};



template< typename LatticeModel_T, typename FlagField_T >
inline void FlagFieldSweepBase< LatticeModel_T, FlagField_T >::getFields( IBlock * const block, PdfField_T * & src, PdfField_T * & dst,
                                                                          const FlagField_T * & flags )
{
   WALBERLA_ASSERT_NOT_NULLPTR( block );

   this->getFields( block, src, dst );
   flags = block->getData<FlagField_T>( flagField_ );

   WALBERLA_ASSERT_NOT_NULLPTR( flags );

   WALBERLA_ASSERT_EQUAL( src->xyzSize(), flags->xyzSize() );
}



template< typename LatticeModel_T, typename FlagField_T >
inline void FlagFieldSweepBase< LatticeModel_T, FlagField_T >::getFields( IBlock * const block, PdfField_T * & src,
                                                                          const FlagField_T * & flags )
{
   WALBERLA_ASSERT_NOT_NULLPTR( block );

   src = this->getSrcField( block );
   flags = block->getData<FlagField_T>( flagField_ );

   WALBERLA_ASSERT_NOT_NULLPTR( src );
   WALBERLA_ASSERT_NOT_NULLPTR( flags );

   WALBERLA_ASSERT_EQUAL( src->xyzSize(), flags->xyzSize() );
}



} // namespace lbm
} // namespace walberla
