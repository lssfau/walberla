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
//! \file ActiveCellSweep.h
//! \ingroup lbm
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "lbm/sweeps/FlagFieldSweepBase.h"
#include "core/OpenMP.h"


namespace walberla {
namespace lbm {



template< typename LatticeModel_T, typename FlagField_T, typename CellOperation >
class ActiveCellSweep : FlagFieldSweepBase< LatticeModel_T, FlagField_T >
{
public:

   using PdfField_T = typename FlagFieldSweepBase<LatticeModel_T, FlagField_T>::PdfField_T;

   // block has NO dst pdf field, lbm mask consists of multiple flags
   ActiveCellSweep( const CellOperation & op, const BlockDataID & pdfField, const ConstBlockDataID & flagField,
                    const Set< FlagUID > & lbmMask, const bool _useIterators = false ) :
      FlagFieldSweepBase<LatticeModel_T,FlagField_T>( pdfField, flagField, lbmMask ), cellOperation_( op ), useIterators_( _useIterators ) {}

   // every block has a dedicated dst pdf field, lbm mask consists of multiple flags
   ActiveCellSweep( const CellOperation & op, const BlockDataID & src, const BlockDataID & dst, const ConstBlockDataID & flagField,
                    const Set< FlagUID > & lbmMask, const bool _useIterators = false ) :
      FlagFieldSweepBase<LatticeModel_T,FlagField_T>( src, dst, flagField, lbmMask ), cellOperation_( op ), useIterators_( _useIterators ) {}

   virtual ~ActiveCellSweep() = default;

   const CellOperation & getCellOperation() const { return cellOperation_; }
         CellOperation & getCellOperation()       { return cellOperation_; }

   bool usesIterators() const { return useIterators_; }
   void  useIterators( const bool _useIterators ) { useIterators_ = _useIterators; }

   void operator()( IBlock * const block );

private:

   void itLoop( typename PdfField_T::iterator & srcIter, typename PdfField_T::iterator & dstIter,
                typename FlagField_T::const_iterator & flgIter, const typename FlagField_T::flag_t & lbm )
   {
      while( srcIter != PdfField_T::staticEnd )
      {
         if( field::isPartOfMaskSet( flgIter, lbm ) )
            cellOperation_( srcIter, dstIter );

         ++srcIter;
         ++dstIter;
         ++flgIter;
      }
   }



   CellOperation cellOperation_;
   bool useIterators_;
};



template< typename LatticeModel_T, typename FlagField_T, typename CellOperation >
void ActiveCellSweep< LatticeModel_T, FlagField_T, CellOperation >::operator()( IBlock * const block )
{
   PdfField_T * src( NULL );
   PdfField_T * dst( NULL );
   const FlagField_T * flagField( NULL );

   auto lbm = this->getLbmMaskAndFields( block, src, dst, flagField );

   WALBERLA_ASSERT_NOT_NULLPTR( src );
   WALBERLA_ASSERT_NOT_NULLPTR( dst );
   WALBERLA_ASSERT_NOT_NULLPTR( flagField );

   WALBERLA_ASSERT_GREATER_EQUAL( src->nrOfGhostLayers(), 1 );

   // execute stream & collide kernels
   
   const auto & lm = src->latticeModel();
   dst->resetLatticeModel( lm ); /* required so that member functions for getting density and equilibrium velocity can be called for dst! */

   cellOperation_.configure( lm );

   if( useIterators_ )
   {
#ifdef _OPENMP
      const cell_idx_t xSize = cell_idx_c( src->xSize() );
      const cell_idx_t ySize = cell_idx_c( src->ySize() );
      const int        zSize =      int_c( src->zSize() );

      #pragma omp parallel for schedule(static)
      for( int iz = 0; iz < zSize; ++iz ) {
         cell_idx_t z = cell_idx_c( iz );
         CellInterval interval( cell_idx_t(0), cell_idx_t(0), z, xSize - cell_idx_t(1), ySize - cell_idx_t(1), z );

         auto srcIter = src->beginSliceXYZ( interval );
         auto dstIter = dst->beginSliceXYZ( interval );
         auto flgIter = flagField->beginSliceXYZ( interval );

         itLoop( srcIter, dstIter, flgIter, lbm );
      }
#else
      auto srcIter = src->beginXYZ();
      auto dstIter = dst->beginXYZ();
      auto flgIter = flagField->beginXYZ();

      itLoop( srcIter, dstIter, flgIter, lbm );
#endif
   }
   else
   {
      const cell_idx_t xSize = cell_idx_c( src->xSize() );
      const cell_idx_t ySize = cell_idx_c( src->ySize() );
      const cell_idx_t zSize = cell_idx_c( src->zSize() );

#ifdef _OPENMP
      const int izSize = int_c( zSize );
      #pragma omp parallel for schedule(static)
      for( int iz = 0; iz < izSize; ++iz ) {
         cell_idx_t z =  cell_idx_c( iz );
#else
      for( cell_idx_t z = 0; z < zSize; ++z ) {
#endif
         for( cell_idx_t y = 0; y < ySize; ++y ) {
            for( cell_idx_t x = 0; x < xSize; ++x )
            {
               if( flagField->isPartOfMaskSet( x, y, z, lbm ) )
                  cellOperation_( src, dst, x, y, z );
            }
         }
      }
   }

   src->swapDataPointers( dst );
}



} // namespace lbm
} // namespace walberla
