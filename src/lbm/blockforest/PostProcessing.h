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
//! \file PostProcessing.h
//! \ingroup lbm
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "blockforest/BlockDataHandling.h"
#include "core/cell/CellInterval.h"
#include "lbm/field/PdfField.h"
#include "lbm/lattice_model/NeighborsStencil.h"
#include "field/FlagField.h"



namespace walberla {
namespace lbm {



namespace internal
{
   using MarkerField_T = field::Field<uint8_t, 1>;
}



// To be used after dynamic changes to the block structure (= after dynamic refinement)
template< typename LatticeModel_T, typename Filter_T >
class PostProcessing
{
public:

   using PdfField_T = PdfField<LatticeModel_T>;
   using Stencil_T = typename LatticeModel_T::Stencil;
   using NeighborsStencil_T = typename NeighborsStencil<LatticeModel_T>::type;

   PostProcessing( const BlockDataID & pdfFieldId, const BlockDataID & markerFieldId, const Filter_T & filter ) :
      pdfFieldId_( pdfFieldId ), markerFieldId_( markerFieldId ), filter_( filter )
   {}

   void operator()( BlockForest & forest, const PhantomBlockForest & phantomForest );

private:

   BlockDataID pdfFieldId_;
   BlockDataID markerFieldId_;
   Filter_T filter_;
};



template< typename LatticeModel_T, typename Filter_T >
void PostProcessing< LatticeModel_T, Filter_T >::operator()( BlockForest & forest, const PhantomBlockForest & )
{
   /*
   static const real_t weights[8][3] = { { -0.25, -0.25, -0.25 },
                                         { -0.25, -0.25,  0.25 },
                                         { -0.25,  0.25, -0.25 },
                                         { -0.25,  0.25,  0.25 },
                                         {  0.25, -0.25, -0.25 },
                                         {  0.25, -0.25,  0.25 },
                                         {  0.25,  0.25, -0.25 },
                                         {  0.25,  0.25,  0.25 } };
   */

   for( auto block = forest.begin(); block != forest.end(); ++block )
   {
      internal::MarkerField_T * markerField = block->template getData< internal::MarkerField_T >( markerFieldId_ );

      if( markerField )
      {
         if( ! markerField->xyzSize().empty() )
         {
            PdfField_T * pdfField = block->template getData< PdfField_T >( pdfFieldId_ );
            WALBERLA_ASSERT_NOT_NULLPTR( pdfField );

            filter_( *block );

            const CellInterval interval = pdfField->xyzSize();
            WALBERLA_ASSERT_EQUAL( pdfField->xyzSize(), markerField->xyzSize() );

            /*
            // trilinear interpolation with extrapolation on block boundaries and chen correction
            //
            auto phantom = phantomForest.getBlock( dynamic_cast< const BlockID & >( block->getId() ) );
            WALBERLA_ASSERT_NOT_NULLPTR( phantom );
            if( phantom->sourceBlockIsLarger() )
            {
               WALBERLA_ASSERT_EQUAL( interval.xSize() & uint_t(1), uint_t(0) );
               WALBERLA_ASSERT_EQUAL( interval.ySize() & uint_t(1), uint_t(0) );
               WALBERLA_ASSERT( Stencil_T::D == uint_t(2) || ( interval.zSize() & uint_t(1) ) == uint_t(0) );

               const uint_t cxs = interval.xSize() / uint_t(2);
               const uint_t cys = interval.ySize() / uint_t(2);
               const uint_t czs = ( Stencil_T::D == uint_t(2) ) ? uint_t(1) : ( interval.zSize() / uint_t(2) );

               field::Field< typename PdfField_T::value_type, PdfField_T::F_SIZE > coarse( cxs, cys, czs, pdfField->layout() );

               const CellInterval cInterval = coarse.xyzSize();

               WALBERLA_FOR_ALL_CELLS_IN_INTERVAL_XYZ( cInterval,

                  const cell_idx_t fx = x * cell_idx_t(2);
                  const cell_idx_t fy = y * cell_idx_t(2);
                  const cell_idx_t fz = z * cell_idx_t(2);

                  for( uint_t f = uint_t(0); f != PdfField_T::F_SIZE; ++f )
                     coarse(x,y,z,f) = pdfField->get(fx,fy,fz,f);
               )

               WALBERLA_FOR_ALL_CELLS_IN_INTERVAL_XYZ( cInterval,

                  const cell_idx_t fx = x * cell_idx_t(2);
                  const cell_idx_t fy = y * cell_idx_t(2);
                  const cell_idx_t fz = z * cell_idx_t(2);

                  const cell_idx_t one = cell_idx_t(1);

                  if( markerField->get(fx,fy,fz) != uint8_t(0) )
                  {
                     if( filter_(fx,fy,fz) && filter_(fx+one,fy,fz) && filter_(fx,fy+one,fz) && filter_(fx+one,fy+one,fz) &&
                         ( Stencil_T::D == uint_t(2) ||
                           ( filter_(fx,fy,fz+one) && filter_(fx+one,fy,fz+one) && filter_(fx,fy+one,fz+one) && filter_(fx+one,fy+one,fz+one) ) ) )
                     {
                        Cell min[3];
                        Cell max[3];

                        min[0][0] = x - one; min[0][1] = y; min[0][2] = z;
                        max[0][0] = x + one; max[0][1] = y; max[0][2] = z;
                        min[1][0] = x; min[1][1] = y - one; min[1][2] = z;
                        max[1][0] = x; max[1][1] = y + one; max[1][2] = z;
                        min[2][0] = x; min[2][1] = y; min[2][2] = z - one;
                        max[2][0] = x; max[2][1] = y; max[2][2] = z + one;

                        Cell fmin[3];
                        Cell fmax[3];

                        for( int i = 0; i != 3; ++i )
                        {
                           fmin[i] = min[i];
                           fmax[i] = max[i];
                           for( int j = 0; j != 3; ++j )
                           {
                              fmin[i][j] *= cell_idx_t(2);
                              fmax[i][j] *= cell_idx_t(2);
                           }
                        }

                        for( uint_t f = 0; f < PdfField_T::F_SIZE; ++f )
                        {
                           WALBERLA_ASSERT( !math::isnan( coarse(x,y,z,f) ) );
                           const auto v = coarse(x,y,z,f);

                           Vector3< real_t > grad( real_t(0) );

                           for( uint_t i = 0; i < Stencil_T::D; ++i )
                           {
                              if( interval.contains( fmax[i] ) && markerField->get( fmax[i] ) != uint8_t(0) )
                              {
                                 if( interval.contains( fmin[i] ) && markerField->get( fmin[i] ) != uint8_t(0) )
                                 {
                                    grad[i] = real_t(0.5) * ( coarse( max[i], f ) - coarse( min[i], f ) );
                                 }
                                 else
                                 {
                                    grad[i] = coarse( max[i], f ) - v;
                                 }
                              }
                              else if( interval.contains( fmin[i] ) && markerField->get( fmin[i] ) != uint8_t(0) )
                              {
                                 grad[i] = v - coarse( min[i], f );
                              }
                           }

                           const auto direction = Stencil_T::dir[f];
                           Vector3< real_t > cNorm( stencil::cNorm[0][direction], stencil::cNorm[1][direction], stencil::cNorm[2][direction] );

                           grad = grad - cNorm * ( cNorm * grad );

                           const auto xx = fx + one;
                           const auto yy = fy + one;
                           pdfField->get( fx, fy, fz , f ) = v + ( weights[0][0] * grad[0] + weights[0][1] * grad[1] + weights[0][2] * grad[2] );
                           pdfField->get( fx, yy, fz , f ) = v + ( weights[2][0] * grad[0] + weights[2][1] * grad[1] + weights[2][2] * grad[2] );
                           pdfField->get( xx, fy, fz , f ) = v + ( weights[4][0] * grad[0] + weights[4][1] * grad[1] + weights[4][2] * grad[2] );
                           pdfField->get( xx, yy, fz , f ) = v + ( weights[6][0] * grad[0] + weights[6][1] * grad[1] + weights[6][2] * grad[2] );

                           if( Stencil_T::D == uint_t(3) )
                           {
                              const auto zz = fz + one;
                              pdfField->get( fx, fy, zz, f ) = v + ( weights[1][0] * grad[0] + weights[1][1] * grad[1] + weights[1][2] * grad[2] );
                              pdfField->get( fx, yy, zz, f ) = v + ( weights[3][0] * grad[0] + weights[3][1] * grad[1] + weights[3][2] * grad[2] );
                              pdfField->get( xx, fy, zz, f ) = v + ( weights[5][0] * grad[0] + weights[5][1] * grad[1] + weights[5][2] * grad[2] );
                              pdfField->get( xx, yy, zz, f ) = v + ( weights[7][0] * grad[0] + weights[7][1] * grad[1] + weights[7][2] * grad[2] );
                           }
                        }

                     }
                  }
                  else
                  {
                     WALBERLA_ASSERT_UNEQUAL( markerField->get(fx+one,fy,    fz    ), uint8_t(0) );
                     WALBERLA_ASSERT_UNEQUAL( markerField->get(fx,    fy+one,fz    ), uint8_t(0) );
                     WALBERLA_ASSERT_UNEQUAL( markerField->get(fx+one,fy+one,fz    ), uint8_t(0) );
                     if( Stencil_T::D != uint_t(2) )
                     {
                        WALBERLA_ASSERT_UNEQUAL( markerField->get(fx,    fy,    fz+one), uint8_t(0) );
                        WALBERLA_ASSERT_UNEQUAL( markerField->get(fx+one,fy,    fz+one), uint8_t(0) );
                        WALBERLA_ASSERT_UNEQUAL( markerField->get(fx,    fy+one,fz+one), uint8_t(0) );
                        WALBERLA_ASSERT_UNEQUAL( markerField->get(fx+one,fy+one,fz+one), uint8_t(0) );
                     }
                  }
               )
            }
            */

            // treatment of boundary->fluid cell conversions
            //
            WALBERLA_FOR_ALL_CELLS_XYZ( pdfField,

               if( filter_(x,y,z) )
               {
                  if( markerField->get(x,y,z) == uint8_t(0) )
                  {
                     bool reconstructed( false );

                     for( auto it = NeighborsStencil_T::begin(); it != NeighborsStencil_T::end() && !reconstructed; ++it )
                     {
                        const cell_idx_t nx = x + it.cx();
                        const cell_idx_t ny = y + it.cy();
                        const cell_idx_t nz = z + it.cz();

                        if( interval.contains( nx, ny, nz ) )
                        {
                           if( filter_( nx, ny, nz ) && (markerField->get( nx, ny, nz ) != uint8_t(0)) )
                           {
                              for( uint_t f = uint_t(0); f != pdfField->fSize(); ++f )
                                 pdfField->get(x,y,z,f) = pdfField->get(nx,ny,nz,f);
                              reconstructed = true;
                           }
                        }
                     }
                     
                     if( !reconstructed )
                        pdfField->setToEquilibrium(x,y,z);

                     /*
                     Vector3<real_t> velocity;
                     real_t density( real_t(0) );
                     real_t count( real_t(0) );

                     for( auto it = NeighborsStencil_T::begin(); it != NeighborsStencil_T::end(); ++it )
                     {
                        const cell_idx_t nx = x + it.cx();
                        const cell_idx_t ny = y + it.cy();
                        const cell_idx_t nz = z + it.cz();

                        if( interval.contains( nx, ny, nz ) )
                        {
                           if( filter_( nx, ny, nz ) && (markerField->get( nx, ny, nz ) != uint8_t(0)) )
                           {
                              Vector3<real_t> vel;
                              density += pdfField->getDensityAndVelocity(vel,x,y,z);
                              velocity += vel;
                              count += real_t(1);
                           }
                        }
                     }
                     
                     if( count > real_t(0) )
                     {
                        const real_t factor = real_t(1) / count;
                        velocity *= factor;
                        density *= factor;
                     }
                     else
                     {
                        density = real_t(1);
                     }

                     pdfField->setDensityAndVelocity(x,y,z,velocity,density);
                     */
                  }
               }
            )

            markerField->resize( uint_t(0), uint_t(0), uint_t(0) );
         }
      }
   }
}






template< typename LatticeModel_T, typename Filter_T >
class MarkerFieldGenerator
{
public:

   using PdfField_T = PdfField<LatticeModel_T>;

   MarkerFieldGenerator( const BlockDataID & pdfFieldId, const BlockDataID & markerFieldId, const Filter_T & filter ) :
      pdfFieldId_( pdfFieldId ), markerFieldId_( markerFieldId ), filter_( filter )
   {}

   void operator()( BlockForest & forest, const PhantomBlockForest & phantomForest );

private:

   BlockDataID pdfFieldId_;
   BlockDataID markerFieldId_;
   Filter_T filter_;
};



template< typename LatticeModel_T, typename Filter_T >
void MarkerFieldGenerator< LatticeModel_T, Filter_T >::operator()( BlockForest & forest, const PhantomBlockForest & )
{
   for( auto block = forest.begin(); block != forest.end(); ++block )
   {
      PdfField_T * pdfField = block->template getData< PdfField_T >( pdfFieldId_ );
      if( pdfField )
      {
         internal::MarkerField_T * markerField = block->template getData< internal::MarkerField_T >( markerFieldId_ );
         WALBERLA_ASSERT_NOT_NULLPTR( markerField );
         if( markerField->xyzSize().empty() )
         {
            const CellInterval interval = pdfField->xyzSize();
            markerField->resize( interval.xSize(), interval.ySize(), interval.zSize() );
            filter_( *block );
            for( cell_idx_t z = interval.zMin(); z <= interval.zMax(); ++z ) {
               for( cell_idx_t y = interval.yMin(); y <= interval.yMax(); ++y ) {
                  for( cell_idx_t x = interval.xMin(); x <= interval.xMax(); ++x )
                  {
                     markerField->get(x,y,z) = filter_(x,y,z) ? uint8_t(1) : uint8_t(0);
                  }
               }
            }
         }
      }
   }
}






template< typename LatticeModel_T, typename Filter_T >
class MarkerData : public blockforest::BlockDataHandling< internal::MarkerField_T >
{
public:

   using PdfField_T = PdfField<LatticeModel_T>;

   MarkerData( const BlockDataID & pdfFieldId, const Filter_T & filter ) :
      pdfFieldId_( pdfFieldId ), filter_( filter )
   {}

   ~MarkerData() override = default;

   internal::MarkerField_T * initialize( IBlock * const ) override { return allocate(); }

   void serialize( IBlock * const block, const BlockDataID & id, mpi::SendBuffer & buffer ) override;

   void serializeCoarseToFine( Block * const block, const BlockDataID & id, mpi::SendBuffer & buffer, const uint_t child ) override;
   void serializeFineToCoarse( Block * const block, const BlockDataID & id, mpi::SendBuffer & buffer ) override;

   internal::MarkerField_T * deserialize( IBlock * const ) override { return allocate(); }

   internal::MarkerField_T * deserializeCoarseToFine( Block * const ) override { return allocate(); }
   internal::MarkerField_T * deserializeFineToCoarse( Block * const ) override { return allocate(); }
   
   void deserialize( IBlock * const block, const BlockDataID & id, mpi::RecvBuffer & buffer ) override;

   void deserializeCoarseToFine( Block * const block, const BlockDataID & id, mpi::RecvBuffer & buffer ) override;
   void deserializeFineToCoarse( Block * const block, const BlockDataID & id, mpi::RecvBuffer & buffer, const uint_t child ) override;   

protected:

   internal::MarkerField_T * allocate() { return new internal::MarkerField_T( uint_t(0), uint_t(0), uint_t(0) ); }

   void sizeCheck( const uint_t xSize, const uint_t ySize, const uint_t zSize )
   {
      WALBERLA_CHECK( (xSize & uint_t(1)) == uint_t(0), "The x-size of your field must be divisible by 2." );
      WALBERLA_CHECK( (ySize & uint_t(1)) == uint_t(0), "The y-size of your field must be divisible by 2." );
      if( LatticeModel_T::Stencil::D == uint_t(2) )
      { WALBERLA_CHECK( zSize == uint_t(1), "The z-size of your field must be equal to 1 (pseudo 2D mode)." ); }
      else
      { WALBERLA_CHECK( (zSize & uint_t(1)) == uint_t(0), "The z-size of your field must be divisible by 2." ); }
   }

   BlockDataID pdfFieldId_;
   Filter_T filter_; // attention: the same filter object cannot be used in parallel by multiple threads working on different blocks

}; // class MarkerData



template< typename LatticeModel_T, typename Filter_T >
void MarkerData< LatticeModel_T, Filter_T >::serialize( IBlock * const block, const BlockDataID &, mpi::SendBuffer & buffer )
{
   PdfField_T * field = block->template getData<PdfField_T >( pdfFieldId_ );
   WALBERLA_ASSERT_NOT_NULLPTR( field );

#ifndef NDEBUG
   buffer << field->xSize() << field->ySize() << field->zSize();
#endif

   const CellInterval interval = field->xyzSize();

   std::vector<bool> filter;
   filter.reserve( interval.numCells() );

   Filter_T localFilter( filter_ );
   localFilter( *block );

   for( cell_idx_t z = interval.zMin(); z <= interval.zMax(); ++z ) {
      for( cell_idx_t y = interval.yMin(); y <= interval.yMax(); ++y ) {
         for( cell_idx_t x = interval.xMin(); x <= interval.xMax(); ++x )
         {
            filter.push_back( localFilter(x,y,z) );
         }
      }
   }

   buffer << filter;
}



template< typename LatticeModel_T, typename Filter_T >
void MarkerData< LatticeModel_T, Filter_T >::serializeCoarseToFine( Block * const block, const BlockDataID &, mpi::SendBuffer & buffer, const uint_t child )
{
   PdfField_T * field = block->template getData<PdfField_T >( pdfFieldId_ );
   WALBERLA_ASSERT_NOT_NULLPTR( field );

   const uint_t xSize = field->xSize();
   const uint_t ySize = field->ySize();
   const uint_t zSize = field->zSize();
   sizeCheck( xSize, ySize, zSize );
   
#ifndef NDEBUG
   buffer << child << ( xSize / uint_t(2) ) << ( ySize / uint_t(2) ) << ( (LatticeModel_T::Stencil::D == uint_t(3)) ? ( zSize / uint_t(2) ) : zSize );
#endif   

   std::vector<bool> filter( (LatticeModel_T::Stencil::D == uint_t(3)) ? ((xSize * ySize * zSize) / uint_t(8)) : ((xSize * ySize * zSize) / uint_t(4)) );

   Filter_T localFilter( filter_ );
   localFilter( *block );

   uint_t i( uint_t(0) );
   const cell_idx_t zBegin = (LatticeModel_T::Stencil::D == uint_t(3)) ? ( (child & uint_t(4)) ? ( cell_idx_c( zSize ) / cell_idx_t(2) ) : cell_idx_t(0) ) : cell_idx_t(0);
   const cell_idx_t zEnd = (LatticeModel_T::Stencil::D == uint_t(3)) ? ( (child & uint_t(4)) ? cell_idx_c( zSize ) : ( cell_idx_c( zSize ) / cell_idx_t(2) ) ) : cell_idx_t(1);
   for( cell_idx_t z = zBegin; z < zEnd; ++z )
   {
      const cell_idx_t yEnd = (child & uint_t(2)) ? cell_idx_c( ySize ) : ( cell_idx_c( ySize ) / cell_idx_t(2) );
      for( cell_idx_t y = (child & uint_t(2)) ? ( cell_idx_c( ySize ) / cell_idx_t(2) ) : cell_idx_t(0); y < yEnd; ++y )
      {
         const cell_idx_t xEnd = (child & uint_t(1)) ? cell_idx_c( xSize ) : ( cell_idx_c( xSize ) / cell_idx_t(2) );
         for( cell_idx_t x = (child & uint_t(1)) ? ( cell_idx_c( xSize ) / cell_idx_t(2) ) : cell_idx_t(0); x < xEnd; ++x )
         {
            WALBERLA_ASSERT_LESS( i, filter.size() );
            filter[i] = localFilter(x,y,z);
            ++i;
         }
      }
   }
   WALBERLA_ASSERT_EQUAL( i, filter.size() );

   buffer << filter;
}



template< typename LatticeModel_T, typename Filter_T >
void MarkerData< LatticeModel_T, Filter_T >::serializeFineToCoarse( Block * const block, const BlockDataID &, mpi::SendBuffer & buffer )
{
   PdfField_T * field = block->template getData<PdfField_T >( pdfFieldId_ );
   WALBERLA_ASSERT_NOT_NULLPTR( field );

   const uint_t xSize = field->xSize();
   const uint_t ySize = field->ySize();
   const uint_t zSize = field->zSize();
   sizeCheck( xSize, ySize, zSize );

#ifndef NDEBUG
   buffer << block->getId().getBranchId() << ( xSize / uint_t(2) ) << ( ySize / uint_t(2) ) << ( (LatticeModel_T::Stencil::D == uint_t(3)) ? ( zSize / uint_t(2) ) : zSize );
#endif

   std::vector<bool> filter( (LatticeModel_T::Stencil::D == uint_t(3)) ? ((xSize * ySize * zSize) / uint_t(8)) : ((xSize * ySize * zSize) / uint_t(4)) );

   Filter_T localFilter( filter_ );
   localFilter( *block );

   uint_t i( uint_t(0) );
   for( cell_idx_t z = cell_idx_t(0); z < cell_idx_c( zSize ); z += cell_idx_t(2) ) {
      for( cell_idx_t y = cell_idx_t(0); y < cell_idx_c( ySize ); y += cell_idx_t(2) ) {
         for( cell_idx_t x = cell_idx_t(0); x < cell_idx_c( xSize ); x += cell_idx_t(2) )
         {
            bool result = localFilter( x,                 y,                 z                 ) &&
                          localFilter( x + cell_idx_t(1), y                , z                 ) &&
                          localFilter( x                , y + cell_idx_t(1), z                 ) &&
                          localFilter( x + cell_idx_t(1), y + cell_idx_t(1), z                 );
            if( LatticeModel_T::Stencil::D == uint_t(3) )
            {
               result = result && localFilter( x                , y                , z + cell_idx_t(1) ) &&
                                  localFilter( x + cell_idx_t(1), y                , z + cell_idx_t(1) ) &&
                                  localFilter( x                , y + cell_idx_t(1), z + cell_idx_t(1) ) &&
                                  localFilter( x + cell_idx_t(1), y + cell_idx_t(1), z + cell_idx_t(1) );
            }

            WALBERLA_ASSERT_LESS( i, filter.size() );
            filter[i] = result;
            ++i;
         }
      }
   }
   WALBERLA_ASSERT_EQUAL( i, filter.size() );

   buffer << filter;
}



template< typename LatticeModel_T, typename Filter_T >
void MarkerData< LatticeModel_T, Filter_T >::deserialize( IBlock * const block, const BlockDataID & id, mpi::RecvBuffer & buffer )
{
   PdfField_T * pdfField = block->template getData<PdfField_T >( pdfFieldId_ );
   WALBERLA_ASSERT_NOT_NULLPTR( pdfField );

#ifndef NDEBUG
   uint_t xSender( uint_t(0) );
   uint_t ySender( uint_t(0) );
   uint_t zSender( uint_t(0) );
   buffer >> xSender >> ySender >> zSender;
   WALBERLA_ASSERT_EQUAL( xSender, pdfField->xSize() );
   WALBERLA_ASSERT_EQUAL( ySender, pdfField->ySize() );
   WALBERLA_ASSERT_EQUAL( zSender, pdfField->zSize() );
#endif

   const CellInterval interval = pdfField->xyzSize();

   std::vector<bool> filter;
   buffer >> filter;
   WALBERLA_ASSERT_EQUAL( filter.size(), interval.numCells() );

   internal::MarkerField_T * markerField = block->template getData< internal::MarkerField_T >( id );
   if( markerField->xyzSize().empty() )
      markerField->resize( interval.xSize(), interval.ySize(), interval.zSize() );

   uint_t i( uint_t(0) );
   for( cell_idx_t z = interval.zMin(); z <= interval.zMax(); ++z ) {
      for( cell_idx_t y = interval.yMin(); y <= interval.yMax(); ++y ) {
         for( cell_idx_t x = interval.xMin(); x <= interval.xMax(); ++x )
         {
            markerField->get(x,y,z) = filter[i] ? uint8_t(1) : uint8_t(0);
            ++i;
         }
      }
   }
}



template< typename LatticeModel_T, typename Filter_T >
void MarkerData< LatticeModel_T, Filter_T >::deserializeCoarseToFine( Block * const block, const BlockDataID & id, mpi::RecvBuffer & buffer )
{
   PdfField_T * pdfField = block->template getData<PdfField_T >( pdfFieldId_ );
   WALBERLA_ASSERT_NOT_NULLPTR( pdfField );

   const uint_t xSize = pdfField->xSize();
   const uint_t ySize = pdfField->ySize();
   const uint_t zSize = pdfField->zSize();
   sizeCheck( xSize, ySize, zSize );

#ifndef NDEBUG
   uint_t branchId( uint_t(0) );
   uint_t xSender( uint_t(0) );
   uint_t ySender( uint_t(0) );
   uint_t zSender( uint_t(0) );
   buffer >> branchId >> xSender >> ySender >> zSender;
   WALBERLA_ASSERT_EQUAL( branchId, block->getId().getBranchId() );
   WALBERLA_ASSERT_EQUAL( xSender, xSize / uint_t(2) );
   WALBERLA_ASSERT_EQUAL( ySender, ySize / uint_t(2) );
   if( LatticeModel_T::Stencil::D == uint_t(3) )
   { WALBERLA_ASSERT_EQUAL( zSender, zSize / uint_t(2) ); }
   else
   { WALBERLA_ASSERT_EQUAL( zSender, zSize ); }
#endif

   std::vector<bool> filter;
   buffer >> filter;
   WALBERLA_ASSERT_EQUAL( filter.size(), (LatticeModel_T::Stencil::D == uint_t(3)) ? (pdfField->xyzSize().numCells() / uint_t(8)) : (pdfField->xyzSize().numCells() / uint_t(4)) );

   internal::MarkerField_T * markerField = block->template getData< internal::MarkerField_T >( id );
   if( markerField->xyzSize().empty() )
      markerField->resize( xSize, ySize, zSize );   

   uint_t i( uint_t(0) );
   for( cell_idx_t z = cell_idx_t(0); z < cell_idx_c( zSize ); z += cell_idx_t(2) ) {
      for( cell_idx_t y = cell_idx_t(0); y < cell_idx_c( ySize ); y += cell_idx_t(2) ) {
         for( cell_idx_t x = cell_idx_t(0); x < cell_idx_c( xSize ); x += cell_idx_t(2) )
         {
            WALBERLA_ASSERT_LESS( i, filter.size() );
            const uint8_t marker = filter[i] ? uint8_t(1) : uint8_t(0);
            markerField->get( x,                 y,                 z                 ) = marker;
            markerField->get( x + cell_idx_t(1), y                , z                 ) = marker;
            markerField->get( x                , y + cell_idx_t(1), z                 ) = marker;
            markerField->get( x + cell_idx_t(1), y + cell_idx_t(1), z                 ) = marker;
            if( LatticeModel_T::Stencil::D == uint_t(3) )
            {
               markerField->get( x                , y                , z + cell_idx_t(1) ) = marker;
               markerField->get( x + cell_idx_t(1), y                , z + cell_idx_t(1) ) = marker;
               markerField->get( x                , y + cell_idx_t(1), z + cell_idx_t(1) ) = marker;
               markerField->get( x + cell_idx_t(1), y + cell_idx_t(1), z + cell_idx_t(1) ) = marker;
            }
            ++i;
         }
      }
   }
   WALBERLA_ASSERT_EQUAL( i, filter.size() );
}



template< typename LatticeModel_T, typename Filter_T >
void MarkerData< LatticeModel_T, Filter_T >::deserializeFineToCoarse( Block * const block, const BlockDataID & id, mpi::RecvBuffer & buffer, const uint_t child )
{
   PdfField_T * pdfField = block->template getData<PdfField_T >( pdfFieldId_ );
   WALBERLA_ASSERT_NOT_NULLPTR( pdfField );

   const uint_t xSize = pdfField->xSize();
   const uint_t ySize = pdfField->ySize();
   const uint_t zSize = pdfField->zSize();
   sizeCheck( xSize, ySize, zSize );
   
#ifndef NDEBUG
   uint_t branchId( uint_t(0) );
   uint_t xSender( uint_t(0) );
   uint_t ySender( uint_t(0) );
   uint_t zSender( uint_t(0) );
   buffer >> branchId >> xSender >> ySender >> zSender;
   WALBERLA_ASSERT_EQUAL( branchId, child );
   WALBERLA_ASSERT_EQUAL( xSender, xSize / uint_t(2) );
   WALBERLA_ASSERT_EQUAL( ySender, ySize / uint_t(2) );
   if( LatticeModel_T::Stencil::D == uint_t(3) )
   { WALBERLA_ASSERT_EQUAL( zSender, zSize / uint_t(2) ); }
   else
   { WALBERLA_ASSERT_EQUAL( zSender, zSize ); }
#endif

   std::vector<bool> filter;
   buffer >> filter;
   WALBERLA_ASSERT_EQUAL( filter.size(), (LatticeModel_T::Stencil::D == uint_t(3)) ? (pdfField->xyzSize().numCells() / uint_t(8)) : (pdfField->xyzSize().numCells() / uint_t(4)) );
   
   internal::MarkerField_T * markerField = block->template getData< internal::MarkerField_T >( id );
   if( markerField->xyzSize().empty() )
      markerField->resize( xSize, ySize, zSize );  

   uint_t i( uint_t(0) );
   const cell_idx_t zBegin = (LatticeModel_T::Stencil::D == uint_t(3)) ? ((child & uint_t(4)) ? ( cell_idx_c( zSize ) / cell_idx_t(2) ) : cell_idx_t(0)) : cell_idx_t(0);
   const cell_idx_t zEnd = (LatticeModel_T::Stencil::D == uint_t(3)) ? ((child & uint_t(4)) ? cell_idx_c( zSize ) : ( cell_idx_c( zSize ) / cell_idx_t(2) )) : cell_idx_t(1);
   for( cell_idx_t z = zBegin; z < zEnd; ++z )
   {
      const cell_idx_t yEnd = (child & uint_t(2)) ? cell_idx_c( ySize ) : ( cell_idx_c( ySize ) / cell_idx_t(2) );
      for( cell_idx_t y = (child & uint_t(2)) ? ( cell_idx_c( ySize ) / cell_idx_t(2) ) : cell_idx_t(0); y < yEnd; ++y )
      {
         const cell_idx_t xEnd = (child & uint_t(1)) ? cell_idx_c( xSize ) : ( cell_idx_c( xSize ) / cell_idx_t(2) );
         for( cell_idx_t x = (child & uint_t(1)) ? ( cell_idx_c( xSize ) / cell_idx_t(2) ) : cell_idx_t(0); x < xEnd; ++x )
         {
            WALBERLA_ASSERT_LESS( i, filter.size() );
            markerField->get(x,y,z) = filter[i] ? uint8_t(1) : uint8_t(0);
            ++i;
         }
      }
   }
   WALBERLA_ASSERT_EQUAL( i, filter.size() );
}



} // namespace lbm
} // namespace walberla
