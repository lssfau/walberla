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
//! \file Streaming.h
//! \ingroup lbm
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "lbm/field/PdfField.h"
#include "lbm/lattice_model/LatticeModelBase.h"

#include "core/debug/Debug.h"

#include "field/iterators/IteratorMacros.h"

#include "stencil/D3Q19.h"
#include "stencil/Directions.h"

#include <type_traits>


namespace walberla {
namespace lbm {



template< typename LatticeModel_T, typename FlagField_T, class Enable = void >
struct Stream; // streaming (stream pull) performed only for cells marked as 'LBM' cells

template< typename LatticeModel_T, class Enable = void >
struct StreamEverything; // streaming performed for all cells



////////////
// STREAM //
////////////

template< typename LatticeModel_T, typename FlagField_T >
struct Stream< LatticeModel_T, FlagField_T, typename std::enable_if< ! std::is_same< typename LatticeModel_T::Stencil,
                                                                                                         stencil::D3Q19>::value >::type >
{
   static_assert( (std::is_same< typename LatticeModel_T::Stencil, stencil::D3Q19 >::value == false), "There is a specialization for D3Q19!" );

   using PdfField_T = PdfField<LatticeModel_T>;
   using Stencil = typename LatticeModel_T::Stencil;

   static void execute( PdfField_T * src, PdfField_T * dst, const FlagField_T * flagField, const typename FlagField_T::flag_t lbm,
                        const uint_t numberOfGhostLayersToInclude = uint_t(0) );
};

template< typename LatticeModel_T, typename FlagField_T >
void Stream< LatticeModel_T, FlagField_T, typename std::enable_if< ! std::is_same< typename LatticeModel_T::Stencil,
                                                                                                       stencil::D3Q19 >::value >::type
   >::execute( PdfField_T * src, PdfField_T * dst, const FlagField_T * flagField, const typename FlagField_T::flag_t lbm,
               const uint_t numberOfGhostLayersToInclude )
{
   WALBERLA_ASSERT_NOT_NULLPTR( src );
   WALBERLA_ASSERT_NOT_NULLPTR( dst );
   WALBERLA_ASSERT_NOT_NULLPTR( flagField );

   WALBERLA_ASSERT_GREATER( src->nrOfGhostLayers(), numberOfGhostLayersToInclude );
   WALBERLA_ASSERT_GREATER_EQUAL( dst->nrOfGhostLayers(), numberOfGhostLayersToInclude );
   WALBERLA_ASSERT_GREATER_EQUAL( flagField->nrOfGhostLayers(), numberOfGhostLayersToInclude );

   WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ( src, numberOfGhostLayersToInclude,

      if( flagField->isPartOfMaskSet( x, y, z, lbm ) )
      {
         for( auto d = Stencil::begin(); d != Stencil::end(); ++d )
            dst->get( x, y, z, d.toIdx() ) = src->get( x-d.cx(), y-d.cy(), z-d.cz(), d.toIdx() );
      }

   ) // WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ

   src->swapDataPointers( dst );
}



template< typename LatticeModel_T, typename FlagField_T >
struct Stream< LatticeModel_T, FlagField_T, typename std::enable_if< std::is_same< typename LatticeModel_T::Stencil,
                                                                                       stencil::D3Q19 >::value >::type >
{
   static_assert( (std::is_same< typename LatticeModel_T::Stencil, stencil::D3Q19 >::value), "Only works with D3Q19!" );

   using PdfField_T = PdfField<LatticeModel_T>;
   using Stencil = typename LatticeModel_T::Stencil;

   static void execute( PdfField_T * src, PdfField_T * dst, const FlagField_T * flagField, const typename FlagField_T::flag_t lbm,
                        const uint_t numberOfGhostLayersToInclude = uint_t(0) );
};

template< typename LatticeModel_T, typename FlagField_T >
void Stream< LatticeModel_T, FlagField_T, typename std::enable_if< std::is_same< typename LatticeModel_T::Stencil,
                                                                                     stencil::D3Q19 >::value >::type
   >::execute( PdfField_T * src, PdfField_T * dst, const FlagField_T * flagField, const typename FlagField_T::flag_t lbm,
               const uint_t numberOfGhostLayersToInclude )
{
   WALBERLA_ASSERT_NOT_NULLPTR( src );
   WALBERLA_ASSERT_NOT_NULLPTR( dst );
   WALBERLA_ASSERT_NOT_NULLPTR( flagField );

   WALBERLA_ASSERT_GREATER( src->nrOfGhostLayers(), numberOfGhostLayersToInclude );
   WALBERLA_ASSERT_GREATER_EQUAL( dst->nrOfGhostLayers(), numberOfGhostLayersToInclude );
   WALBERLA_ASSERT_GREATER_EQUAL( flagField->nrOfGhostLayers(), numberOfGhostLayersToInclude );

   WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ( src, numberOfGhostLayersToInclude,

      using namespace stencil;

      if( flagField->isPartOfMaskSet( x, y, z, lbm ) )
      {
         dst->get(x,y,z, Stencil::idx[C])  = src->get(x  , y  , z  , Stencil::idx[C]);
         dst->get(x,y,z, Stencil::idx[N])  = src->get(x  , y-1, z  , Stencil::idx[N]);
         dst->get(x,y,z, Stencil::idx[S])  = src->get(x  , y+1, z  , Stencil::idx[S]);
         dst->get(x,y,z, Stencil::idx[W])  = src->get(x+1, y  , z  , Stencil::idx[W]);
         dst->get(x,y,z, Stencil::idx[E])  = src->get(x-1, y  , z  , Stencil::idx[E]);
         dst->get(x,y,z, Stencil::idx[T])  = src->get(x  , y  , z-1, Stencil::idx[T]);
         dst->get(x,y,z, Stencil::idx[B])  = src->get(x  , y  , z+1, Stencil::idx[B]);
         dst->get(x,y,z, Stencil::idx[NW]) = src->get(x+1, y-1, z  , Stencil::idx[NW]);
         dst->get(x,y,z, Stencil::idx[NE]) = src->get(x-1, y-1, z  , Stencil::idx[NE]);
         dst->get(x,y,z, Stencil::idx[SW]) = src->get(x+1, y+1, z  , Stencil::idx[SW]);
         dst->get(x,y,z, Stencil::idx[SE]) = src->get(x-1, y+1, z  , Stencil::idx[SE]);
         dst->get(x,y,z, Stencil::idx[TN]) = src->get(x  , y-1, z-1, Stencil::idx[TN]);
         dst->get(x,y,z, Stencil::idx[TS]) = src->get(x  , y+1, z-1, Stencil::idx[TS]);
         dst->get(x,y,z, Stencil::idx[TW]) = src->get(x+1, y  , z-1, Stencil::idx[TW]);
         dst->get(x,y,z, Stencil::idx[TE]) = src->get(x-1, y  , z-1, Stencil::idx[TE]);
         dst->get(x,y,z, Stencil::idx[BN]) = src->get(x  , y-1, z+1, Stencil::idx[BN]);
         dst->get(x,y,z, Stencil::idx[BS]) = src->get(x  , y+1, z+1, Stencil::idx[BS]);
         dst->get(x,y,z, Stencil::idx[BW]) = src->get(x+1, y  , z+1, Stencil::idx[BW]);
         dst->get(x,y,z, Stencil::idx[BE]) = src->get(x-1, y  , z+1, Stencil::idx[BE]);
      }

   ) // WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ

   src->swapDataPointers( dst );
}



///////////////////////
// STREAM EVERYTHING //
///////////////////////

template< typename LatticeModel_T >
struct StreamEverything< LatticeModel_T, typename std::enable_if< ! std::is_same< typename LatticeModel_T::Stencil,
                                                                                                      stencil::D3Q19 >::value >::type >
{
   static_assert( (std::is_same< typename LatticeModel_T::Stencil, stencil::D3Q19 >::value == false), "There is a specialization for D3Q19!" );

   using PdfField_T = PdfField<LatticeModel_T>;
   using Stencil = typename LatticeModel_T::Stencil;

   static void execute( PdfField_T * src, PdfField_T * dst, const uint_t numberOfGhostLayersToInclude = uint_t(0) );
};

template< typename LatticeModel_T >
void StreamEverything< LatticeModel_T, typename std::enable_if< ! std::is_same< typename LatticeModel_T::Stencil,
                                                                                                    stencil::D3Q19 >::value >::type
   >::execute( PdfField_T * src, PdfField_T * dst, const uint_t numberOfGhostLayersToInclude )
{
   WALBERLA_ASSERT_NOT_NULLPTR( src );
   WALBERLA_ASSERT_NOT_NULLPTR( dst );

   WALBERLA_ASSERT_GREATER( src->nrOfGhostLayers(), numberOfGhostLayersToInclude );
   WALBERLA_ASSERT_GREATER_EQUAL( dst->nrOfGhostLayers(), numberOfGhostLayersToInclude );

   if( src->layout() == field::fzyx && dst->layout() == field::fzyx )
   {
      WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_YZ( src, numberOfGhostLayersToInclude,

         const cell_idx_t xSize = cell_idx_c( src->xSize() );
         for( auto d = Stencil::begin(); d != Stencil::end(); ++d )
         {
            real_t * srcPtr = &src->get( -cell_idx_c( numberOfGhostLayersToInclude )-d.cx(), y-d.cy(), z-d.cz(), d.toIdx() );
            real_t * dstPtr = &dst->get( -cell_idx_c( numberOfGhostLayersToInclude ), y, z, d.toIdx() );
            std::copy( srcPtr, srcPtr + xSize + cell_idx_t(2) * cell_idx_c( numberOfGhostLayersToInclude ), dstPtr );
         }

      ) // WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_YZ
   }
   else // ==> src->layout() == field::zyxf || dst->layout() == field::zyxf
   {
      WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ( src, numberOfGhostLayersToInclude,

         for( auto d = Stencil::begin(); d != Stencil::end(); ++d )
            dst->get( x, y, z, d.toIdx() ) = src->get( x-d.cx(), y-d.cy(), z-d.cz(), d.toIdx() );

      ) // WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ
   }

   src->swapDataPointers( dst );
}



template< typename LatticeModel_T >
struct StreamEverything< LatticeModel_T, typename std::enable_if< std::is_same< typename LatticeModel_T::Stencil,
                                                                                    stencil::D3Q19 >::value >::type >
{
   static_assert( (std::is_same< typename LatticeModel_T::Stencil, stencil::D3Q19 >::value), "Only works with D3Q19!" );

   using PdfField_T = PdfField<LatticeModel_T>;
   using Stencil = typename LatticeModel_T::Stencil;

   static void execute( PdfField_T * src, PdfField_T * dst, const uint_t numberOfGhostLayersToInclude = uint_t(0) );
};

template< typename LatticeModel_T >
void StreamEverything< LatticeModel_T, typename std::enable_if< std::is_same< typename LatticeModel_T::Stencil,
                                                                                  stencil::D3Q19 >::value >::type
   >::execute( PdfField_T * src, PdfField_T * dst, const uint_t numberOfGhostLayersToInclude )
{
   WALBERLA_ASSERT_NOT_NULLPTR( src );
   WALBERLA_ASSERT_NOT_NULLPTR( dst );

   WALBERLA_ASSERT_GREATER( src->nrOfGhostLayers(), numberOfGhostLayersToInclude );
   WALBERLA_ASSERT_GREATER_EQUAL( dst->nrOfGhostLayers(), numberOfGhostLayersToInclude );

   if( src->layout() == field::fzyx && dst->layout() == field::fzyx )
   {
      WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_YZ( src, numberOfGhostLayersToInclude,

         const cell_idx_t xSize = cell_idx_c( src->xSize() );
         for( auto d = Stencil::begin(); d != Stencil::end(); ++d )
         {
            real_t * srcPtr = &src->get( -cell_idx_c( numberOfGhostLayersToInclude )-d.cx(), y-d.cy(), z-d.cz(), d.toIdx() );
            real_t * dstPtr = &dst->get( -cell_idx_c( numberOfGhostLayersToInclude ), y, z, d.toIdx() );
            std::copy( srcPtr, srcPtr + xSize + cell_idx_t(2) * cell_idx_c( numberOfGhostLayersToInclude ), dstPtr );
         }

      ) // WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_YZ
   }
   else // ==> src->layout() == field::zyxf || dst->layout() == field::zyxf
   {
      WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ( src, numberOfGhostLayersToInclude,

         using namespace stencil;

         dst->get(x,y,z, Stencil::idx[C])  = src->get(x  , y  , z  , Stencil::idx[C]);
         dst->get(x,y,z, Stencil::idx[N])  = src->get(x  , y-1, z  , Stencil::idx[N]);
         dst->get(x,y,z, Stencil::idx[S])  = src->get(x  , y+1, z  , Stencil::idx[S]);
         dst->get(x,y,z, Stencil::idx[W])  = src->get(x+1, y  , z  , Stencil::idx[W]);
         dst->get(x,y,z, Stencil::idx[E])  = src->get(x-1, y  , z  , Stencil::idx[E]);
         dst->get(x,y,z, Stencil::idx[T])  = src->get(x  , y  , z-1, Stencil::idx[T]);
         dst->get(x,y,z, Stencil::idx[B])  = src->get(x  , y  , z+1, Stencil::idx[B]);
         dst->get(x,y,z, Stencil::idx[NW]) = src->get(x+1, y-1, z  , Stencil::idx[NW]);
         dst->get(x,y,z, Stencil::idx[NE]) = src->get(x-1, y-1, z  , Stencil::idx[NE]);
         dst->get(x,y,z, Stencil::idx[SW]) = src->get(x+1, y+1, z  , Stencil::idx[SW]);
         dst->get(x,y,z, Stencil::idx[SE]) = src->get(x-1, y+1, z  , Stencil::idx[SE]);
         dst->get(x,y,z, Stencil::idx[TN]) = src->get(x  , y-1, z-1, Stencil::idx[TN]);
         dst->get(x,y,z, Stencil::idx[TS]) = src->get(x  , y+1, z-1, Stencil::idx[TS]);
         dst->get(x,y,z, Stencil::idx[TW]) = src->get(x+1, y  , z-1, Stencil::idx[TW]);
         dst->get(x,y,z, Stencil::idx[TE]) = src->get(x-1, y  , z-1, Stencil::idx[TE]);
         dst->get(x,y,z, Stencil::idx[BN]) = src->get(x  , y-1, z+1, Stencil::idx[BN]);
         dst->get(x,y,z, Stencil::idx[BS]) = src->get(x  , y+1, z+1, Stencil::idx[BS]);
         dst->get(x,y,z, Stencil::idx[BW]) = src->get(x+1, y  , z+1, Stencil::idx[BW]);
         dst->get(x,y,z, Stencil::idx[BE]) = src->get(x-1, y  , z+1, Stencil::idx[BE]);

      ) // WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ
   }

   src->swapDataPointers( dst );
}



} // namespace lbm
} // namespace walberla
