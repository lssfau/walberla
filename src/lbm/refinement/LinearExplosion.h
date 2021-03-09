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
//! \file LinearExplosion.h
//! \ingroup lbm
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "lbm/field/PdfField.h"
#include "lbm/lattice_model/NeighborsStencil.h"

#include "blockforest/Block.h"
#include "blockforest/BlockNeighborhoodSection.h"

#include "core/OpenMP.h"
#include "core/cell/CellInterval.h"
#include "core/debug/Debug.h"

#include "stencil/D3Q27.h"

#include <map>



namespace walberla {
namespace lbm {
namespace refinement {



namespace internal {

// These two functions have been (static) member functions of class LinearExplosion. Sadly, there have been strange compile errors with
// the IBM compiler -> possible bug in the compiler. As a workaround, these two functions are now two free functions ...

template< typename PdfField_T, typename BoundaryHandling_T, typename CoarseField >
void fillTemporaryCoarseField( const cell_idx_t y, const cell_idx_t z, const CellInterval & coarse,
                               const PdfField_T * const pdfField, const BoundaryHandling_T * const boundaryHandling,
                               const shared_ptr< CoarseField > & tmpField );

template< typename PdfField_T, typename BoundaryHandling_T, typename CoarseField, typename BoolField >
void linearInterpolation( const cell_idx_t y, const cell_idx_t z, const CellInterval & fine,
                          PdfField_T * const pdfField, const BoundaryHandling_T * const boundaryHandling,
                          const shared_ptr< CoarseField > & tmpField, const shared_ptr< BoolField > & boolField );

} // namespace internal



template< typename LatticeModel_T, typename BoundaryHandling_T >
class LinearExplosion
{
private:

   using Stencil = typename LatticeModel_T::Stencil;

   using PdfField_T = PdfField<LatticeModel_T>;

   using CoarseField = field::GhostLayerField<typename PdfField_T::value_type, PdfField_T::F_SIZE>;
   using BoolField = field::GhostLayerField<bool, 1>;

public:

   LinearExplosion( const BlockDataID & pdfFieldId, const ConstBlockDataID & boundaryHandlingId ) :
      pdfFieldId_( pdfFieldId ), boundaryHandlingId_( boundaryHandlingId ) {}

   void operator()( Block * block );

protected:



   BlockDataID pdfFieldId_;
   ConstBlockDataID boundaryHandlingId_;

   std::map< Cell, std::pair< shared_ptr< CoarseField >, shared_ptr< BoolField > > > tmpFields_;
};



template< typename LatticeModel_T, typename BoundaryHandling_T >
void LinearExplosion< LatticeModel_T, BoundaryHandling_T >::operator()( Block * block )
{
   PdfField_T * pdfField = block->template getData< PdfField_T >( pdfFieldId_ );
   const BoundaryHandling_T * boundaryHandling = block->template getData< BoundaryHandling_T >( boundaryHandlingId_ );

   CellInterval fineInterval = pdfField->xyzSize();
   CellInterval coarseInterval = fineInterval;

   fineInterval.expand( cell_idx_t(4) );
   coarseInterval.xMax() = ( coarseInterval.xMax() >> 1 );
   coarseInterval.yMax() = ( coarseInterval.yMax() >> 1 );
   coarseInterval.zMax() = ( Stencil::D == uint_t(3) ) ? ( coarseInterval.zMax() >> 1 ) : cell_idx_t(0);
   coarseInterval.expand( cell_idx_t(2) );

   std::vector< CellInterval > coarseIntervals;
   std::vector< CellInterval >   fineIntervals;

   // determine regions that need to be interpolated

   for( auto d = NeighborsStencil< LatticeModel_T >::type::beginNoCenter(); d != NeighborsStencil< LatticeModel_T >::type::end(); ++d )
   {
      if( block->neighborhoodSectionHasLargerBlock( blockforest::getBlockNeighborhoodSectionIndex(*d) ) )
      {
         CellInterval coarse = coarseInterval;
         CellInterval fine = fineInterval;

         for( uint_t i = 0; i < 3; ++i )
         {
            if( stencil::c[i][*d] == -1 )
            {
               coarse.max()[i] = coarse.min()[i] + cell_idx_t(1);
               fine.max()[i] = fine.min()[i] + cell_idx_t(3);
            }
            else if( stencil::c[i][*d] == 1 )
            {
               coarse.min()[i] = coarse.max()[i] - cell_idx_t(1);
               fine.min()[i] = fine.max()[i] - cell_idx_t(3);
            }
            else
            {
               WALBERLA_ASSERT_EQUAL( stencil::c[i][*d], 0 );
               coarse.min()[i] += cell_idx_t(2);
               coarse.max()[i] -= cell_idx_t(2);
               fine.min()[i] += cell_idx_t(4);
               fine.max()[i] -= cell_idx_t(4);
            }
         }

         coarseIntervals.push_back( coarse );
         fineIntervals.push_back( fine );
      }
   }

   // it there are no coarse neighbors -> return

   if( coarseIntervals.empty() )
   {
      WALBERLA_ASSERT( fineIntervals.empty() );
      return;
   }

   // fetch temporary fields from temporary field storage "tmpFields_"

   if( tmpFields_.find( coarseInterval.max() ) == tmpFields_.end() )
   {
      tmpFields_[ coarseInterval.max() ] = std::make_pair( make_shared< CoarseField >( coarseInterval.xSize() - uint_t(4), coarseInterval.ySize() - uint_t(4),
                                                                                       coarseInterval.zSize() - uint_t(4), uint_t(2), pdfField->layout() ),
                                                           make_shared< BoolField >( coarseInterval.xSize() - uint_t(4), coarseInterval.ySize() - uint_t(4),
                                                                                     coarseInterval.zSize() - uint_t(4), uint_t(3), pdfField->layout() ) );
   }
   auto fields = tmpFields_[ coarseInterval.max() ];
   auto tmpField = fields.first;
   auto boolField = fields.second;
   
#ifndef NDEBUG
   std::fill( tmpField->beginWithGhostLayer(), tmpField->end(), std::numeric_limits< typename PdfField_T::value_type >::quiet_NaN() );
   std::fill( boolField->beginWithGhostLayer(), boolField->end(), false );
#endif

   // fill temporary coarse field with values (only regions that are required for the interpolation!)

   for( auto coarse = coarseIntervals.begin(); coarse != coarseIntervals.end(); ++coarse )
   {
#ifdef _OPENMP

      if( coarse->zSize() >= coarse->ySize() )
      {
         int zMin = int_c( coarse->zMin() );
         int zMax = int_c( coarse->zMax() ) + 1;

         #pragma omp parallel for schedule(static)
         for( int zi = zMin; zi < zMax; ++zi )
         {
            const cell_idx_t z = cell_idx_c(zi);
            for( cell_idx_t y = coarse->yMin(); y <= coarse->yMax(); ++y )
               internal::fillTemporaryCoarseField< PdfField_T, BoundaryHandling_T, CoarseField >( y, z, *coarse, pdfField, boundaryHandling, tmpField );
         }
      }
      else
      {
         int yMin = int_c( coarse->yMin() );
         int yMax = int_c( coarse->yMax() ) + 1;

         #pragma omp parallel for schedule(static)
         for( int yi = yMin; yi < yMax; ++yi )
         {
            const cell_idx_t y = cell_idx_c(yi);
            for( cell_idx_t z = coarse->zMin(); z <= coarse->zMax(); ++z )
               internal::fillTemporaryCoarseField< PdfField_T, BoundaryHandling_T, CoarseField >( y, z, *coarse, pdfField, boundaryHandling, tmpField );
         }
      }
#else
      for( cell_idx_t z = coarse->zMin(); z <= coarse->zMax(); ++z ) {
         for( cell_idx_t y = coarse->yMin(); y <= coarse->yMax(); ++y )
         {
            internal::fillTemporaryCoarseField< PdfField_T, BoundaryHandling_T, CoarseField >( y, z, *coarse, pdfField, boundaryHandling, tmpField );
         }
      }
#endif

      CellInterval expanded( *coarse );
      expanded.expand( cell_idx_t(1) );

      for( auto cell = boolField->beginSliceXYZ( expanded ); cell != boolField->end(); ++cell )
         *cell = false;
   }

   for( auto coarse = coarseIntervals.begin(); coarse != coarseIntervals.end(); ++coarse ) {
      for( auto z = coarse->zMin(); z <= coarse->zMax(); ++z ) {
         for( auto y = coarse->yMin(); y <= coarse->yMax(); ++y ) {
            for( auto x = coarse->xMin(); x <= coarse->xMax(); ++x )
            {
               if( boundaryHandling->isDomain( cell_idx_t(2) * x, cell_idx_t(2) * y, cell_idx_t(2) * z ) )
               {
                  WALBERLA_ASSERT( boundaryHandling->isDomain( cell_idx_t(2) * x + cell_idx_t(1), cell_idx_t(2) * y,                 cell_idx_t(2) * z                 ) );
                  WALBERLA_ASSERT( boundaryHandling->isDomain( cell_idx_t(2) * x,                 cell_idx_t(2) * y + cell_idx_t(1), cell_idx_t(2) * z                 ) );
                  WALBERLA_ASSERT( boundaryHandling->isDomain( cell_idx_t(2) * x + cell_idx_t(1), cell_idx_t(2) * y + cell_idx_t(1), cell_idx_t(2) * z                 ) );
                  if( Stencil::D == uint_t(3) )
                  {
                     WALBERLA_ASSERT( boundaryHandling->isDomain( cell_idx_t(2) * x,                 cell_idx_t(2) * y,                 cell_idx_t(2) * z + cell_idx_t(1) ) );
                     WALBERLA_ASSERT( boundaryHandling->isDomain( cell_idx_t(2) * x + cell_idx_t(1), cell_idx_t(2) * y,                 cell_idx_t(2) * z + cell_idx_t(1) ) );
                     WALBERLA_ASSERT( boundaryHandling->isDomain( cell_idx_t(2) * x,                 cell_idx_t(2) * y + cell_idx_t(1), cell_idx_t(2) * z + cell_idx_t(1) ) );
                     WALBERLA_ASSERT( boundaryHandling->isDomain( cell_idx_t(2) * x + cell_idx_t(1), cell_idx_t(2) * y + cell_idx_t(1), cell_idx_t(2) * z + cell_idx_t(1) ) );
                  }
                  
                  boolField->get(x,y,z) = true;
               }
               else
               {
                  WALBERLA_ASSERT( !boundaryHandling->isDomain( cell_idx_t(2) * x + cell_idx_t(1), cell_idx_t(2) * y,                 cell_idx_t(2) * z                 ) );
                  WALBERLA_ASSERT( !boundaryHandling->isDomain( cell_idx_t(2) * x,                 cell_idx_t(2) * y + cell_idx_t(1), cell_idx_t(2) * z                 ) );
                  WALBERLA_ASSERT( !boundaryHandling->isDomain( cell_idx_t(2) * x + cell_idx_t(1), cell_idx_t(2) * y + cell_idx_t(1), cell_idx_t(2) * z                 ) );
                  if( Stencil::D == uint_t(3) )
                  {
                     WALBERLA_ASSERT( !boundaryHandling->isDomain( cell_idx_t(2) * x,                 cell_idx_t(2) * y,                 cell_idx_t(2) * z + cell_idx_t(1) ) );
                     WALBERLA_ASSERT( !boundaryHandling->isDomain( cell_idx_t(2) * x + cell_idx_t(1), cell_idx_t(2) * y,                 cell_idx_t(2) * z + cell_idx_t(1) ) );
                     WALBERLA_ASSERT( !boundaryHandling->isDomain( cell_idx_t(2) * x,                 cell_idx_t(2) * y + cell_idx_t(1), cell_idx_t(2) * z + cell_idx_t(1) ) );
                     WALBERLA_ASSERT( !boundaryHandling->isDomain( cell_idx_t(2) * x + cell_idx_t(1), cell_idx_t(2) * y + cell_idx_t(1), cell_idx_t(2) * z + cell_idx_t(1) ) );
                  }
               }
            }
         }
      }
   }

   // linear interpolation

   for( auto fine = fineIntervals.begin(); fine != fineIntervals.end(); ++fine )
   {
      WALBERLA_ASSERT( (fine->xSize() & uint_t(1)) == uint_t(0) );
      WALBERLA_ASSERT( (fine->ySize() & uint_t(1)) == uint_t(0) );
      WALBERLA_ASSERT( ( Stencil::D == uint_t(2) && fine->zSize() == uint_t(1) ) || ( Stencil::D == uint_t(3) && (fine->zSize() & uint_t(1)) == uint_t(0) ) );
      
#ifdef _OPENMP

      if( fine->zSize() >= fine->ySize() )
      {
         int zSize = (Stencil::D == uint_t(3)) ? ( int_c( fine->zSize() ) / 2 ) : 1;

         #pragma omp parallel for schedule(static)
         for( int zi = 0; zi < zSize; ++zi )
         {
            const cell_idx_t z = fine->zMin() + cell_idx_c(zi) * cell_idx_t(2);
            for( cell_idx_t y = fine->yMin(); y <= fine->yMax(); y += cell_idx_t(2) )
               internal::linearInterpolation< PdfField_T, BoundaryHandling_T, CoarseField, BoolField >( y, z, *fine, pdfField, boundaryHandling, tmpField, boolField );
         }
      }
      else
      {
         int ySize = int_c( fine->ySize() ) / 2;

         #pragma omp parallel for schedule(static)
         for( int yi = 0; yi < ySize; ++yi )
         {
            const cell_idx_t y = fine->yMin() + cell_idx_c(yi) * cell_idx_t(2);
            for( cell_idx_t z = fine->zMin(); z <= fine->zMax(); z += cell_idx_t(2) )
               internal::linearInterpolation< PdfField_T, BoundaryHandling_T, CoarseField, BoolField >( y, z, *fine, pdfField, boundaryHandling, tmpField, boolField );
         }
      }
#else
      for( cell_idx_t z = fine->zMin(); z <= fine->zMax(); z += cell_idx_t(2) ) {
         for( cell_idx_t y = fine->yMin(); y <= fine->yMax(); y += cell_idx_t(2) )
         {
            internal::linearInterpolation< PdfField_T, BoundaryHandling_T, CoarseField, BoolField >( y, z, *fine, pdfField, boundaryHandling, tmpField, boolField );
         }
      }
#endif
   }
}



namespace internal {

template< typename PdfField_T, typename BoundaryHandling_T, typename CoarseField >
void fillTemporaryCoarseField( const cell_idx_t y, const cell_idx_t z, const CellInterval & coarse,
                               const PdfField_T * const pdfField, const BoundaryHandling_T * const boundaryHandling,
                               const shared_ptr< CoarseField > & tmpField )
{
   for( cell_idx_t x = coarse.xMin(); x <= coarse.xMax(); ++x ) {
      for( uint_t f = 0; f < PdfField_T::F_SIZE; ++f )
      {
         const auto fx = cell_idx_t(2) * x;
         const auto fy = cell_idx_t(2) * y;
         const auto fz = cell_idx_t(2) * z;

         const auto value = pdfField->get( fx, fy, fz, f );

         if( boundaryHandling->isDomain(fx,fy,fz) )
         {
            WALBERLA_ASSERT( !math::isnan( value ) );
         }

         tmpField->get(x,y,z,f) = value;

         /*
         WALBERLA_ASSERT( !math::isnan( pdfField->get( fx                , fy                , fz                , f ) ) );
         WALBERLA_ASSERT( !math::isnan( pdfField->get( fx                , fy                , fz + cell_idx_t(1), f ) ) );
         WALBERLA_ASSERT( !math::isnan( pdfField->get( fx                , fy + cell_idx_t(1), fz                , f ) ) );
         WALBERLA_ASSERT( !math::isnan( pdfField->get( fx                , fy + cell_idx_t(1), fz + cell_idx_t(1), f ) ) );
         WALBERLA_ASSERT( !math::isnan( pdfField->get( fx + cell_idx_t(1), fy                , fz                , f ) ) );
         WALBERLA_ASSERT( !math::isnan( pdfField->get( fx + cell_idx_t(1), fy                , fz + cell_idx_t(1), f ) ) );
         WALBERLA_ASSERT( !math::isnan( pdfField->get( fx + cell_idx_t(1), fy + cell_idx_t(1), fz                , f ) ) );
         WALBERLA_ASSERT( !math::isnan( pdfField->get( fx + cell_idx_t(1), fy + cell_idx_t(1), fz + cell_idx_t(1), f ) ) );

         auto value  = pdfField->get( fx                , fy                , fz                , f );
              value += pdfField->get( fx                , fy                , fz + cell_idx_t(1), f );
              value += pdfField->get( fx                , fy + cell_idx_t(1), fz                , f );
              value += pdfField->get( fx                , fy + cell_idx_t(1), fz + cell_idx_t(1), f );
              value += pdfField->get( fx + cell_idx_t(1), fy                , fz                , f );
              value += pdfField->get( fx + cell_idx_t(1), fy                , fz + cell_idx_t(1), f );
              value += pdfField->get( fx + cell_idx_t(1), fy + cell_idx_t(1), fz                , f );
              value += pdfField->get( fx + cell_idx_t(1), fy + cell_idx_t(1), fz + cell_idx_t(1), f );

         tmpField->get(x,y,z,f) = real_t(0.125) * value;
         */
      }
   }
}



template< typename PdfField_T, typename BoundaryHandling_T, typename CoarseField, typename BoolField >
void linearInterpolation( const cell_idx_t y, const cell_idx_t z, const CellInterval & fine,
                          PdfField_T * const pdfField, const BoundaryHandling_T * const boundaryHandling,
                          const shared_ptr< CoarseField > & tmpField, const shared_ptr< BoolField > & boolField )
{
   // see: Chen et al., "Grid refinement in lattice Boltzmann methods based on volumetric formulation", 2006

   static const real_t weights[8][3] = { { -0.25, -0.25, -0.25 },
                                         { -0.25, -0.25,  0.25 },
                                         { -0.25,  0.25, -0.25 },
                                         { -0.25,  0.25,  0.25 },
                                         {  0.25, -0.25, -0.25 },
                                         {  0.25, -0.25,  0.25 },
                                         {  0.25,  0.25, -0.25 },
                                         {  0.25,  0.25,  0.25 } };

   for( cell_idx_t x = fine.xMin(); x < fine.xMax(); x += cell_idx_t(2) )
   {
      if( boundaryHandling->isDomain(x,y,z) )
      {
         WALBERLA_ASSERT( boundaryHandling->isDomain( x + cell_idx_t(1), y,                 z                 ) );
         WALBERLA_ASSERT( boundaryHandling->isDomain( x,                 y + cell_idx_t(1), z                 ) );
         WALBERLA_ASSERT( boundaryHandling->isDomain( x + cell_idx_t(1), y + cell_idx_t(1), z                 ) );
         if( PdfField_T::Stencil::D == uint_t(3) )
         {
            WALBERLA_ASSERT( boundaryHandling->isDomain( x,                 y,                 z + cell_idx_t(1) ) );
            WALBERLA_ASSERT( boundaryHandling->isDomain( x + cell_idx_t(1), y,                 z + cell_idx_t(1) ) );
            WALBERLA_ASSERT( boundaryHandling->isDomain( x,                 y + cell_idx_t(1), z + cell_idx_t(1) ) );
            WALBERLA_ASSERT( boundaryHandling->isDomain( x + cell_idx_t(1), y + cell_idx_t(1), z + cell_idx_t(1) ) );
         }

         Cell cell( ( ( x + cell_idx_t(4) ) >> 1 ) - cell_idx_t(2),
                    ( ( y + cell_idx_t(4) ) >> 1 ) - cell_idx_t(2),
                    (PdfField_T::Stencil::D == uint_t(3)) ? ( ( ( z + cell_idx_t(4) ) >> 1 ) - cell_idx_t(2) ) : cell_idx_t(0) );

         Cell min[3], max[3];

         min[0][0] = cell[0] - cell_idx_t(1); min[0][1] = cell[1]; min[0][2] = cell[2];
         max[0][0] = cell[0] + cell_idx_t(1); max[0][1] = cell[1]; max[0][2] = cell[2];
         min[1][0] = cell[0]; min[1][1] = cell[1] - cell_idx_t(1); min[1][2] = cell[2];
         max[1][0] = cell[0]; max[1][1] = cell[1] + cell_idx_t(1); max[1][2] = cell[2];
         min[2][0] = cell[0]; min[2][1] = cell[1]; min[2][2] = cell[2] - cell_idx_t(1);
         max[2][0] = cell[0]; max[2][1] = cell[1]; max[2][2] = cell[2] + cell_idx_t(1);

         for( uint_t f = 0; f < PdfField_T::F_SIZE; ++f )
         {
            WALBERLA_ASSERT( !math::isnan( tmpField->get( cell, f ) ) );

            const auto v = tmpField->get( cell, f );

            Vector3< real_t > grad( real_t(0) );

            for( uint_t i = 0; i < PdfField_T::Stencil::D; ++i )
            {

#define WALBERLA_LBM_REFINEMENT_EXPLOSION_EXCLUDE_EXTRAPOLATION
#ifdef WALBERLA_LBM_REFINEMENT_EXPLOSION_EXCLUDE_EXTRAPOLATION

               if( boolField->get( max[i] ) && boolField->get( min[i] ) )
               {
                  WALBERLA_ASSERT( !math::isnan( tmpField->get( max[i], f ) ) );
                  WALBERLA_ASSERT( !math::isnan( tmpField->get( min[i], f ) ) );

                  grad[i] = real_t(0.5) * ( tmpField->get( max[i], f ) - tmpField->get( min[i], f ) );
               }

#else
               if( boolField->get( max[i] ) )
               {
                  WALBERLA_ASSERT( !math::isnan( tmpField->get( max[i], f ) ) );
                  if( boolField->get( min[i] ) )
                  {
                     WALBERLA_ASSERT( !math::isnan( tmpField->get( min[i], f ) ) );
                     grad[i] = real_t(0.5) * ( tmpField->get( max[i], f ) - tmpField->get( min[i], f ) );
                  }
                  else
                  {
                     grad[i] = tmpField->get( max[i], f ) - v;
                  }

               }
               else if( boolField->get( min[i] ) )
               {
                  WALBERLA_ASSERT( !math::isnan( tmpField->get( min[i], f ) ) );
                  grad[i] = v - tmpField->get( min[i], f );
               }
#endif
#undef WALBERLA_LBM_REFINEMENT_EXPLOSION_EXCLUDE_EXTRAPOLATION
            }

#define WALBERLA_LBM_REFINEMENT_EXPLOSION_CHEN_CORRECTION
#ifdef WALBERLA_LBM_REFINEMENT_EXPLOSION_CHEN_CORRECTION

            const auto direction = PdfField_T::Stencil::dir[f];
            Vector3< real_t > cNorm( stencil::cNorm[0][direction], stencil::cNorm[1][direction], stencil::cNorm[2][direction] );

            grad = grad - cNorm * ( cNorm * grad );
#endif
#undef WALBERLA_LBM_REFINEMENT_EXPLOSION_CHEN_CORRECTION

            const auto xx = x + cell_idx_t(1);
            const auto yy = y + cell_idx_t(1);
            pdfField->get( x , y , z , f ) = v + ( weights[0][0] * grad[0] + weights[0][1] * grad[1] + weights[0][2] * grad[2] );
            pdfField->get( x , yy, z , f ) = v + ( weights[2][0] * grad[0] + weights[2][1] * grad[1] + weights[2][2] * grad[2] );
            pdfField->get( xx, y , z , f ) = v + ( weights[4][0] * grad[0] + weights[4][1] * grad[1] + weights[4][2] * grad[2] );
            pdfField->get( xx, yy, z , f ) = v + ( weights[6][0] * grad[0] + weights[6][1] * grad[1] + weights[6][2] * grad[2] );

            if( PdfField_T::Stencil::D == uint_t(3) )
            {
               const auto zz = z + cell_idx_t(1);
               pdfField->get( x , y , zz, f ) = v + ( weights[1][0] * grad[0] + weights[1][1] * grad[1] + weights[1][2] * grad[2] );
               pdfField->get( x , yy, zz, f ) = v + ( weights[3][0] * grad[0] + weights[3][1] * grad[1] + weights[3][2] * grad[2] );
               pdfField->get( xx, y , zz, f ) = v + ( weights[5][0] * grad[0] + weights[5][1] * grad[1] + weights[5][2] * grad[2] );
               pdfField->get( xx, yy, zz, f ) = v + ( weights[7][0] * grad[0] + weights[7][1] * grad[1] + weights[7][2] * grad[2] );
            }
         }
      }
      else
      {
         WALBERLA_ASSERT( !boundaryHandling->isDomain( x + cell_idx_t(1), y,                 z                 ) );
         WALBERLA_ASSERT( !boundaryHandling->isDomain( x,                 y + cell_idx_t(1), z                 ) );
         WALBERLA_ASSERT( !boundaryHandling->isDomain( x + cell_idx_t(1), y + cell_idx_t(1), z                 ) );
         if( PdfField_T::Stencil::D == uint_t(3) )
         {
            WALBERLA_ASSERT( !boundaryHandling->isDomain( x,                 y,                 z + cell_idx_t(1) ) );
            WALBERLA_ASSERT( !boundaryHandling->isDomain( x + cell_idx_t(1), y,                 z + cell_idx_t(1) ) );
            WALBERLA_ASSERT( !boundaryHandling->isDomain( x,                 y + cell_idx_t(1), z + cell_idx_t(1) ) );
            WALBERLA_ASSERT( !boundaryHandling->isDomain( x + cell_idx_t(1), y + cell_idx_t(1), z + cell_idx_t(1) ) );
         }
      }
   }
}

} // namespace internal



} // namespace refinement
} // namespace lbm
} // namespace walberla
