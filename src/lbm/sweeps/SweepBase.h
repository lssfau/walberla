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
//! \ingroup lbm
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "lbm/field/DensityVelocityCallback.h"
#include "lbm/field/PdfField.h"
#include "core/debug/Debug.h"
#include "domain_decomposition/IBlock.h"
#include "field/EvaluationFilter.h"
#include "field/Field.h"
#include "field/SwapableCompare.h"
#include "field/iterators/IteratorMacros.h"

#include <set>
#include <vector>


namespace walberla {
namespace lbm {



template< typename LatticeModel_T,
          typename Filter_T = walberla::field::DefaultEvaluationFilter,
          typename DensityVelocityIn_T = DefaultDensityEquilibriumVelocityCalculation,
          typename DensityVelocityOut_T = DefaultDensityVelocityCallback >
class SweepBase
{
public:

   using PdfField_T = PdfField<LatticeModel_T>;

   // block has NO dst pdf field
   SweepBase( const BlockDataID & pdfField,
              const Filter_T & _filter = walberla::field::DefaultEvaluationFilter(),
              const DensityVelocityIn_T & _densityVelocityIn = DefaultDensityEquilibriumVelocityCalculation(),
              const DensityVelocityOut_T & _densityVelocityOut = DefaultDensityVelocityCallback() ) :
      src_( pdfField ), dstFromBlockData_( false ), filter_( _filter ),
      densityVelocityIn_( _densityVelocityIn ), densityVelocityOut_( _densityVelocityOut ) {}

   // every block has a dedicated dst pdf field
   SweepBase( const BlockDataID & src, const BlockDataID & dst,
              const Filter_T & _filter = walberla::field::DefaultEvaluationFilter(),
              const DensityVelocityIn_T & _densityVelocityIn = DefaultDensityEquilibriumVelocityCalculation(),
              const DensityVelocityOut_T & _densityVelocityOut = DefaultDensityVelocityCallback() ) :
      src_( src ), dstFromBlockData_( true ), dst_( dst ), filter_( _filter ),
      densityVelocityIn_( _densityVelocityIn ), densityVelocityOut_( _densityVelocityOut ) {}

   virtual ~SweepBase() { for( auto field = dstFields_.begin(); field != dstFields_.end(); ++field ) delete *field; }

   void filter( IBlock & block ) { filter_( block ); }
   bool filter( const cell_idx_t x, const cell_idx_t y, const cell_idx_t z ) const { return filter_(x,y,z); }

   void   densityVelocityIn( IBlock & block ) { densityVelocityIn_( block ); }
   real_t densityVelocityIn( Vector3<real_t> & velocity, const PdfField_T * const field,
                           const cell_idx_t x, const cell_idx_t y, const cell_idx_t z )
   {
      return densityVelocityIn_( velocity, field, x, y, z );
   }

   void densityVelocityOut( IBlock & block ) { densityVelocityOut_( block ); }
   void densityVelocityOut( const cell_idx_t x, const cell_idx_t y, const cell_idx_t z, const LatticeModel_T & lm,
                            const Vector3<real_t> & velocity, const real_t rho )
   {
      densityVelocityOut_( x, y, z, lm, velocity, rho );
   }

protected:

   inline PdfField_T * getSrcField( IBlock * const block ) const;
          PdfField_T * getDstField( IBlock * const block, PdfField_T * const src );

   inline void getFields( IBlock * const block, PdfField_T * & src, PdfField_T * & dst );



   const BlockDataID src_ {};

   const bool dstFromBlockData_;
   const BlockDataID dst_ {};
   std::set< PdfField_T *, field::SwapableCompare< PdfField_T * > > dstFields_;

   Filter_T filter_;
   DensityVelocityIn_T densityVelocityIn_;
   DensityVelocityOut_T densityVelocityOut_;
};



template< typename LatticeModel_T, typename Filter_T, typename DensityVelocityIn_T, typename DensityVelocityOut_T >
inline typename SweepBase< LatticeModel_T, Filter_T, DensityVelocityIn_T, DensityVelocityOut_T >::PdfField_T *
SweepBase< LatticeModel_T, Filter_T, DensityVelocityIn_T, DensityVelocityOut_T >::getSrcField( IBlock * const block ) const
{
   WALBERLA_ASSERT_NOT_NULLPTR( block );

   PdfField_T * src = block->getData<PdfField_T>( src_ );
   
   WALBERLA_ASSERT_NOT_NULLPTR( src );
   
   return src;
}



template< typename LatticeModel_T, typename Filter_T, typename DensityVelocityIn_T, typename DensityVelocityOut_T >
typename SweepBase< LatticeModel_T, Filter_T, DensityVelocityIn_T, DensityVelocityOut_T >::PdfField_T *
SweepBase< LatticeModel_T, Filter_T, DensityVelocityIn_T, DensityVelocityOut_T >::getDstField( IBlock * const block, PdfField_T * const src )
{
   WALBERLA_ASSERT_NOT_NULLPTR( block );
   WALBERLA_ASSERT_NOT_NULLPTR( src );

   if( dstFromBlockData_ )
   {
      PdfField_T * dst = block->getData<PdfField_T>( dst_ );
      WALBERLA_ASSERT_NOT_NULLPTR( dst );
      return dst;
   }

   auto it = dstFields_.find( src );
   if( it != dstFields_.end() )
   {
#ifndef NDEBUG
      std::fill( (*it)->beginWithGhostLayer(), (*it)->end(), std::numeric_limits< typename PdfField_T::value_type >::quiet_NaN() );
#endif
      WALBERLA_ASSERT_NOT_NULLPTR( *it );
      return *it;
   }

   PdfField_T * dst = src->cloneUninitialized();
   WALBERLA_ASSERT_NOT_NULLPTR( dst );
   
   // take care of proper thread<->memory assignment (first-touch allocation policy !)
   WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ( dst,
      for( uint_t f = uint_t(0); f < LatticeModel_T::Stencil::Size; ++f )
         dst->get(x,y,z,f) = std::numeric_limits< typename PdfField_T::value_type >::quiet_NaN();
   )
   dstFields_.insert( dst );

   return dst;
}



template< typename LatticeModel_T, typename Filter_T, typename DensityVelocityIn_T, typename DensityVelocityOut_T >
inline void SweepBase< LatticeModel_T, Filter_T, DensityVelocityIn_T, DensityVelocityOut_T >::getFields( IBlock * const block,
                                                                                                         PdfField_T * & src, PdfField_T * & dst )
{
   WALBERLA_ASSERT_NOT_NULLPTR( block );

   src = getSrcField( block );
   dst = getDstField( block, src );

   WALBERLA_ASSERT_NOT_NULLPTR( src );
   WALBERLA_ASSERT_NOT_NULLPTR( dst );

   WALBERLA_ASSERT_EQUAL( src->xyzSize(), dst->xyzSize() );
}



} // namespace lbm
} // namespace walberla
