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
//! \file StreamPull.h
//! \ingroup lbm
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "lbm/field/PdfField.h"
#include "lbm/lattice_model/LatticeModelBase.h"

#include "core/debug/Debug.h"

#include "domain_decomposition/IBlock.h"

#include "field/EvaluationFilter.h"
#include "field/iterators/IteratorMacros.h"

#include "stencil/D2Q9.h"
#include "stencil/D3Q19.h"
#include "stencil/D3Q27.h"
#include "stencil/Directions.h"

#include <type_traits>

namespace walberla {
namespace lbm {



// template< typename LatticeModel_T, class Enable = void >
// struct StreamPullEverything; // streaming performed for all cells

// template< typename LatticeModel_T, class Enable = void >
// struct StreamPull; // streaming (stream pull) performed only for cells selected by an evaluation filter



///////////////////////
// STREAM EVERYTHING //
///////////////////////

template< typename LatticeModel_T, class Enable = void >
struct StreamPullEverything
{
   static_assert( (std::is_same< typename LatticeModel_T::Stencil, stencil::D2Q9  >::value == false), "There is a specialization for D2Q9!" );
   static_assert( (std::is_same< typename LatticeModel_T::Stencil, stencil::D3Q19 >::value == false), "There is a specialization for D3Q19!" );
   static_assert( (std::is_same< typename LatticeModel_T::Stencil, stencil::D3Q27 >::value == false), "There is a specialization for D3Q27!" );

   using PdfField_T = PdfField<LatticeModel_T>;
   using Stencil_T = typename LatticeModel_T::Stencil;

   static void execute( PdfField_T * src, PdfField_T * dst, const uint_t numberOfGhostLayersToInclude = uint_t(0) );
};

template< typename LatticeModel_T, class Enable >
void StreamPullEverything< LatticeModel_T, Enable >::execute( PdfField_T * src, PdfField_T * dst, const uint_t numberOfGhostLayersToInclude )
{
   WALBERLA_ASSERT_NOT_NULLPTR( src );
   WALBERLA_ASSERT_NOT_NULLPTR( dst );

   WALBERLA_ASSERT_GREATER( src->nrOfGhostLayers(), numberOfGhostLayersToInclude );
   WALBERLA_ASSERT_GREATER_EQUAL( dst->nrOfGhostLayers(), numberOfGhostLayersToInclude );

   if( src->layout() == field::fzyx && dst->layout() == field::fzyx )
   {
      WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_YZ( src, numberOfGhostLayersToInclude,

         const cell_idx_t xSize = cell_idx_c( src->xSize() );
         for( auto d = Stencil_T::begin(); d != Stencil_T::end(); ++d )
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

         for( auto d = Stencil_T::begin(); d != Stencil_T::end(); ++d )
            dst->get( x, y, z, d.toIdx() ) = src->get( x-d.cx(), y-d.cy(), z-d.cz(), d.toIdx() );

      ) // WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ
   }

   src->swapDataPointers( dst );
}



// D2Q9 //

template< typename LatticeModel_T >
struct StreamPullEverything< LatticeModel_T, typename std::enable_if< std::is_same< typename LatticeModel_T::Stencil,
                                                                                          stencil::D2Q9 >::value >::type >
{
   static_assert( (std::is_same< typename LatticeModel_T::Stencil, stencil::D2Q9 >::value), "Only works with D2Q9!" );

   using PdfField_T = PdfField<LatticeModel_T>;
   using Stencil_T = typename LatticeModel_T::Stencil;

   static void execute( PdfField_T * src, PdfField_T * dst, const uint_t numberOfGhostLayersToInclude = uint_t(0) );
};

template< typename LatticeModel_T >
void StreamPullEverything< LatticeModel_T, typename std::enable_if< std::is_same< typename LatticeModel_T::Stencil,
                                                                                        stencil::D2Q9 >::value >::type
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
         for( auto d = Stencil_T::begin(); d != Stencil_T::end(); ++d )
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

         dst->get(x,y,z, Stencil_T::idx[C])  = src->get(x  , y  , z, Stencil_T::idx[C]);
         dst->get(x,y,z, Stencil_T::idx[N])  = src->get(x  , y-1, z, Stencil_T::idx[N]);
         dst->get(x,y,z, Stencil_T::idx[S])  = src->get(x  , y+1, z, Stencil_T::idx[S]);
         dst->get(x,y,z, Stencil_T::idx[W])  = src->get(x+1, y  , z, Stencil_T::idx[W]);
         dst->get(x,y,z, Stencil_T::idx[E])  = src->get(x-1, y  , z, Stencil_T::idx[E]);
         dst->get(x,y,z, Stencil_T::idx[NW]) = src->get(x+1, y-1, z, Stencil_T::idx[NW]);
         dst->get(x,y,z, Stencil_T::idx[NE]) = src->get(x-1, y-1, z, Stencil_T::idx[NE]);
         dst->get(x,y,z, Stencil_T::idx[SW]) = src->get(x+1, y+1, z, Stencil_T::idx[SW]);
         dst->get(x,y,z, Stencil_T::idx[SE]) = src->get(x-1, y+1, z, Stencil_T::idx[SE]);

      ) // WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ
   }

   src->swapDataPointers( dst );
}



// D3Q19 //

template< typename LatticeModel_T >
struct StreamPullEverything< LatticeModel_T, typename std::enable_if< std::is_same< typename LatticeModel_T::Stencil,
                                                                                          stencil::D3Q19 >::value >::type >
{
   static_assert( (std::is_same< typename LatticeModel_T::Stencil, stencil::D3Q19 >::value), "Only works with D3Q19!" );

   using PdfField_T = PdfField<LatticeModel_T>;
   using Stencil_T = typename LatticeModel_T::Stencil;

   static void execute( PdfField_T * src, PdfField_T * dst, const uint_t numberOfGhostLayersToInclude = uint_t(0) );
};

template< typename LatticeModel_T >
void StreamPullEverything< LatticeModel_T, typename std::enable_if< std::is_same< typename LatticeModel_T::Stencil,
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
         for( auto d = Stencil_T::begin(); d != Stencil_T::end(); ++d )
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

         dst->get(x,y,z, Stencil_T::idx[C])  = src->get(x  , y  , z  , Stencil_T::idx[C]);
         dst->get(x,y,z, Stencil_T::idx[N])  = src->get(x  , y-1, z  , Stencil_T::idx[N]);
         dst->get(x,y,z, Stencil_T::idx[S])  = src->get(x  , y+1, z  , Stencil_T::idx[S]);
         dst->get(x,y,z, Stencil_T::idx[W])  = src->get(x+1, y  , z  , Stencil_T::idx[W]);
         dst->get(x,y,z, Stencil_T::idx[E])  = src->get(x-1, y  , z  , Stencil_T::idx[E]);
         dst->get(x,y,z, Stencil_T::idx[T])  = src->get(x  , y  , z-1, Stencil_T::idx[T]);
         dst->get(x,y,z, Stencil_T::idx[B])  = src->get(x  , y  , z+1, Stencil_T::idx[B]);
         dst->get(x,y,z, Stencil_T::idx[NW]) = src->get(x+1, y-1, z  , Stencil_T::idx[NW]);
         dst->get(x,y,z, Stencil_T::idx[NE]) = src->get(x-1, y-1, z  , Stencil_T::idx[NE]);
         dst->get(x,y,z, Stencil_T::idx[SW]) = src->get(x+1, y+1, z  , Stencil_T::idx[SW]);
         dst->get(x,y,z, Stencil_T::idx[SE]) = src->get(x-1, y+1, z  , Stencil_T::idx[SE]);
         dst->get(x,y,z, Stencil_T::idx[TN]) = src->get(x  , y-1, z-1, Stencil_T::idx[TN]);
         dst->get(x,y,z, Stencil_T::idx[TS]) = src->get(x  , y+1, z-1, Stencil_T::idx[TS]);
         dst->get(x,y,z, Stencil_T::idx[TW]) = src->get(x+1, y  , z-1, Stencil_T::idx[TW]);
         dst->get(x,y,z, Stencil_T::idx[TE]) = src->get(x-1, y  , z-1, Stencil_T::idx[TE]);
         dst->get(x,y,z, Stencil_T::idx[BN]) = src->get(x  , y-1, z+1, Stencil_T::idx[BN]);
         dst->get(x,y,z, Stencil_T::idx[BS]) = src->get(x  , y+1, z+1, Stencil_T::idx[BS]);
         dst->get(x,y,z, Stencil_T::idx[BW]) = src->get(x+1, y  , z+1, Stencil_T::idx[BW]);
         dst->get(x,y,z, Stencil_T::idx[BE]) = src->get(x-1, y  , z+1, Stencil_T::idx[BE]);

      ) // WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ
   }

   src->swapDataPointers( dst );
}



// D3Q27 //

template< typename LatticeModel_T >
struct StreamPullEverything< LatticeModel_T, typename std::enable_if< std::is_same< typename LatticeModel_T::Stencil,
                                                                                          stencil::D3Q27 >::value >::type >
{
   static_assert( (std::is_same< typename LatticeModel_T::Stencil, stencil::D3Q27 >::value), "Only works with D3Q27!" );

   using PdfField_T = PdfField<LatticeModel_T>;
   using Stencil_T = typename LatticeModel_T::Stencil;

   static void execute( PdfField_T * src, PdfField_T * dst, const uint_t numberOfGhostLayersToInclude = uint_t(0) );
};

template< typename LatticeModel_T >
void StreamPullEverything< LatticeModel_T, typename std::enable_if< std::is_same< typename LatticeModel_T::Stencil,
                                                                                        stencil::D3Q27 >::value >::type
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
         for( auto d = Stencil_T::begin(); d != Stencil_T::end(); ++d )
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

         dst->get(x,y,z, Stencil_T::idx[C])   = src->get(x  , y  , z  , Stencil_T::idx[C]);
         dst->get(x,y,z, Stencil_T::idx[N])   = src->get(x  , y-1, z  , Stencil_T::idx[N]);
         dst->get(x,y,z, Stencil_T::idx[S])   = src->get(x  , y+1, z  , Stencil_T::idx[S]);
         dst->get(x,y,z, Stencil_T::idx[W])   = src->get(x+1, y  , z  , Stencil_T::idx[W]);
         dst->get(x,y,z, Stencil_T::idx[E])   = src->get(x-1, y  , z  , Stencil_T::idx[E]);
         dst->get(x,y,z, Stencil_T::idx[T])   = src->get(x  , y  , z-1, Stencil_T::idx[T]);
         dst->get(x,y,z, Stencil_T::idx[B])   = src->get(x  , y  , z+1, Stencil_T::idx[B]);
         dst->get(x,y,z, Stencil_T::idx[NW])  = src->get(x+1, y-1, z  , Stencil_T::idx[NW]);
         dst->get(x,y,z, Stencil_T::idx[NE])  = src->get(x-1, y-1, z  , Stencil_T::idx[NE]);
         dst->get(x,y,z, Stencil_T::idx[SW])  = src->get(x+1, y+1, z  , Stencil_T::idx[SW]);
         dst->get(x,y,z, Stencil_T::idx[SE])  = src->get(x-1, y+1, z  , Stencil_T::idx[SE]);
         dst->get(x,y,z, Stencil_T::idx[TN])  = src->get(x  , y-1, z-1, Stencil_T::idx[TN]);
         dst->get(x,y,z, Stencil_T::idx[TS])  = src->get(x  , y+1, z-1, Stencil_T::idx[TS]);
         dst->get(x,y,z, Stencil_T::idx[TW])  = src->get(x+1, y  , z-1, Stencil_T::idx[TW]);
         dst->get(x,y,z, Stencil_T::idx[TE])  = src->get(x-1, y  , z-1, Stencil_T::idx[TE]);
         dst->get(x,y,z, Stencil_T::idx[BN])  = src->get(x  , y-1, z+1, Stencil_T::idx[BN]);
         dst->get(x,y,z, Stencil_T::idx[BS])  = src->get(x  , y+1, z+1, Stencil_T::idx[BS]);
         dst->get(x,y,z, Stencil_T::idx[BW])  = src->get(x+1, y  , z+1, Stencil_T::idx[BW]);
         dst->get(x,y,z, Stencil_T::idx[BE])  = src->get(x-1, y  , z+1, Stencil_T::idx[BE]);
         dst->get(x,y,z, Stencil_T::idx[TNE]) = src->get(x-1, y-1, z-1, Stencil_T::idx[TNE]);
         dst->get(x,y,z, Stencil_T::idx[TNW]) = src->get(x+1, y-1, z-1, Stencil_T::idx[TNW]);
         dst->get(x,y,z, Stencil_T::idx[TSE]) = src->get(x-1, y+1, z-1, Stencil_T::idx[TSE]);
         dst->get(x,y,z, Stencil_T::idx[TSW]) = src->get(x+1, y+1, z-1, Stencil_T::idx[TSW]);
         dst->get(x,y,z, Stencil_T::idx[BNE]) = src->get(x-1, y-1, z+1, Stencil_T::idx[BNE]);
         dst->get(x,y,z, Stencil_T::idx[BNW]) = src->get(x+1, y-1, z+1, Stencil_T::idx[BNW]);
         dst->get(x,y,z, Stencil_T::idx[BSE]) = src->get(x-1, y+1, z+1, Stencil_T::idx[BSE]);
         dst->get(x,y,z, Stencil_T::idx[BSW]) = src->get(x+1, y+1, z+1, Stencil_T::idx[BSW]);

      ) // WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ
   }

   src->swapDataPointers( dst );
}



////////////
// STREAM //
////////////

template< typename LatticeModel_T, class Enable = void >
struct StreamPull
{
   static_assert( (std::is_same< typename LatticeModel_T::Stencil, stencil::D2Q9  >::value == false), "There is a specialization for D2Q9!" );
   static_assert( (std::is_same< typename LatticeModel_T::Stencil, stencil::D3Q19 >::value == false), "There is a specialization for D3Q19!" );
   static_assert( (std::is_same< typename LatticeModel_T::Stencil, stencil::D3Q27 >::value == false), "There is a specialization for D3Q27!" );

   using PdfField_T = PdfField<LatticeModel_T>;
   using Stencil_T = typename LatticeModel_T::Stencil;

   template< typename Filter_T >
   static inline void execute( PdfField_T * src, PdfField_T * dst, IBlock * block, Filter_T & filter, const uint_t numberOfGhostLayersToInclude = uint_t(0) );

   static inline void execute( PdfField_T * src, PdfField_T * dst, IBlock * block, walberla::field::DefaultEvaluationFilter & filter,
                               const uint_t numberOfGhostLayersToInclude = uint_t(0) );

   static inline void execute( PdfField_T * src, PdfField_T * dst, cell_idx_t x, cell_idx_t y, cell_idx_t z );
};

template< typename LatticeModel_T, class Enable >
template< typename Filter_T >
void StreamPull< LatticeModel_T, Enable >::execute( PdfField_T * src, PdfField_T * dst, IBlock * block, Filter_T & filter, const uint_t numberOfGhostLayersToInclude )
{
   WALBERLA_ASSERT_NOT_NULLPTR( src );
   WALBERLA_ASSERT_NOT_NULLPTR( dst );

   WALBERLA_ASSERT_GREATER( src->nrOfGhostLayers(), numberOfGhostLayersToInclude );
   WALBERLA_ASSERT_GREATER_EQUAL( dst->nrOfGhostLayers(), numberOfGhostLayersToInclude );

   filter( *block );
   
   WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ( src, numberOfGhostLayersToInclude,

      if( filter(x,y,z) )
         execute( src, dst, x,y,z );

   ) // WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ

   src->swapDataPointers( dst );
}

template< typename LatticeModel_T, class Enable >
void StreamPull< LatticeModel_T, Enable >::execute( PdfField_T * src, PdfField_T * dst, IBlock * /*block*/, walberla::field::DefaultEvaluationFilter & /*filter*/, const uint_t numberOfGhostLayersToInclude )
{
   StreamPullEverything< LatticeModel_T >::execute( src, dst, numberOfGhostLayersToInclude );
}


template< typename LatticeModel_T, class Enable >
void StreamPull< LatticeModel_T, Enable >::execute( PdfField_T * src, PdfField_T * dst, cell_idx_t x, cell_idx_t y, cell_idx_t z )
{
   for( auto d = Stencil_T::begin(); d != Stencil_T::end(); ++d )
      dst->get( x, y, z, d.toIdx() ) = src->get( x-d.cx(), y-d.cy(), z-d.cz(), d.toIdx() );
}



// D2Q9 //

template< typename LatticeModel_T >
struct StreamPull< LatticeModel_T, typename std::enable_if< std::is_same< typename LatticeModel_T::Stencil,
                                                                                stencil::D2Q9 >::value >::type >
{
   static_assert( (std::is_same< typename LatticeModel_T::Stencil, stencil::D2Q9 >::value), "Only works with D2Q9!" );

   using PdfField_T = PdfField<LatticeModel_T>;
   using Stencil_T = typename LatticeModel_T::Stencil;

   template< typename Filter_T >
   static void execute( PdfField_T * src, PdfField_T * dst, IBlock * block, Filter_T & filter, const uint_t numberOfGhostLayersToInclude = uint_t(0) );

   static void execute( PdfField_T * src, PdfField_T * dst, IBlock * block, walberla::field::DefaultEvaluationFilter & filter,
                        const uint_t numberOfGhostLayersToInclude = uint_t(0) );

   static inline void execute( PdfField_T * src, PdfField_T * dst, cell_idx_t x, cell_idx_t y, cell_idx_t z );
};

template< typename LatticeModel_T >
template< typename Filter_T >
void StreamPull< LatticeModel_T, typename std::enable_if< std::is_same< typename LatticeModel_T::Stencil,
                                                                              stencil::D2Q9 >::value >::type
   >::execute( PdfField_T * src, PdfField_T * dst, IBlock * block, Filter_T & filter, const uint_t numberOfGhostLayersToInclude )
{
   WALBERLA_ASSERT_NOT_NULLPTR( src );
   WALBERLA_ASSERT_NOT_NULLPTR( dst );

   WALBERLA_ASSERT_GREATER( src->nrOfGhostLayers(), numberOfGhostLayersToInclude );
   WALBERLA_ASSERT_GREATER_EQUAL( dst->nrOfGhostLayers(), numberOfGhostLayersToInclude );

   filter( *block );

   WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ( src, numberOfGhostLayersToInclude,

      if( filter(x,y,z) )
         execute(src,dst,x,y,z);

   ) // WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ

   src->swapDataPointers( dst );
}

template< typename LatticeModel_T >
void StreamPull< LatticeModel_T, typename std::enable_if< std::is_same< typename LatticeModel_T::Stencil,
                                                                              stencil::D2Q9 >::value >::type
   >::execute( PdfField_T * src, PdfField_T * dst, IBlock * /*block*/, walberla::field::DefaultEvaluationFilter & /*filter*/, const uint_t numberOfGhostLayersToInclude )
{
   StreamPullEverything< LatticeModel_T >::execute( src, dst, numberOfGhostLayersToInclude );
}

template< typename LatticeModel_T >
void StreamPull< LatticeModel_T, typename std::enable_if< std::is_same< typename LatticeModel_T::Stencil,
                                                                              stencil::D2Q9 >::value >::type
   >::execute( PdfField_T * src, PdfField_T * dst, cell_idx_t x, cell_idx_t y, cell_idx_t z )
{
   using namespace stencil;

   dst->get(x,y,z, Stencil_T::idx[C])  = src->get(x  , y  , z, Stencil_T::idx[C]);
   dst->get(x,y,z, Stencil_T::idx[N])  = src->get(x  , y-1, z, Stencil_T::idx[N]);
   dst->get(x,y,z, Stencil_T::idx[S])  = src->get(x  , y+1, z, Stencil_T::idx[S]);
   dst->get(x,y,z, Stencil_T::idx[W])  = src->get(x+1, y  , z, Stencil_T::idx[W]);
   dst->get(x,y,z, Stencil_T::idx[E])  = src->get(x-1, y  , z, Stencil_T::idx[E]);
   dst->get(x,y,z, Stencil_T::idx[NW]) = src->get(x+1, y-1, z, Stencil_T::idx[NW]);
   dst->get(x,y,z, Stencil_T::idx[NE]) = src->get(x-1, y-1, z, Stencil_T::idx[NE]);
   dst->get(x,y,z, Stencil_T::idx[SW]) = src->get(x+1, y+1, z, Stencil_T::idx[SW]);
   dst->get(x,y,z, Stencil_T::idx[SE]) = src->get(x-1, y+1, z, Stencil_T::idx[SE]);
}



// D3Q19 //

template< typename LatticeModel_T >
struct StreamPull< LatticeModel_T, typename std::enable_if< std::is_same< typename LatticeModel_T::Stencil,
                                                                                stencil::D3Q19 >::value >::type >
{
   static_assert( (std::is_same< typename LatticeModel_T::Stencil, stencil::D3Q19 >::value), "Only works with D3Q19!" );

   using PdfField_T = PdfField<LatticeModel_T>;
   using Stencil_T = typename LatticeModel_T::Stencil;

   template< typename Filter_T >
   static void execute( PdfField_T * src, PdfField_T * dst, IBlock * block, Filter_T & filter, const uint_t numberOfGhostLayersToInclude = uint_t(0) );

   static void execute( PdfField_T * src, PdfField_T * dst, IBlock * block, walberla::field::DefaultEvaluationFilter & filter,
                        const uint_t numberOfGhostLayersToInclude = uint_t(0) );

   static inline void execute( PdfField_T * src, PdfField_T * dst, cell_idx_t x, cell_idx_t y, cell_idx_t z );
};

template< typename LatticeModel_T >
template< typename Filter_T >
void StreamPull< LatticeModel_T, typename std::enable_if< std::is_same< typename LatticeModel_T::Stencil,
                                                                              stencil::D3Q19 >::value >::type
   >::execute( PdfField_T * src, PdfField_T * dst, IBlock * block, Filter_T & filter, const uint_t numberOfGhostLayersToInclude )
{
   WALBERLA_ASSERT_NOT_NULLPTR( src );
   WALBERLA_ASSERT_NOT_NULLPTR( dst );

   WALBERLA_ASSERT_GREATER( src->nrOfGhostLayers(), numberOfGhostLayersToInclude );
   WALBERLA_ASSERT_GREATER_EQUAL( dst->nrOfGhostLayers(), numberOfGhostLayersToInclude );
   
   filter( *block );

   WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ( src, numberOfGhostLayersToInclude,

      if( filter(x,y,z) )
         execute(src,dst,x,y,z);

   ) // WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ

   src->swapDataPointers( dst );
}

template< typename LatticeModel_T >
void StreamPull< LatticeModel_T, typename std::enable_if< std::is_same< typename LatticeModel_T::Stencil,
                                                                              stencil::D3Q19 >::value >::type
   >::execute( PdfField_T * src, PdfField_T * dst, IBlock * /*block*/, walberla::field::DefaultEvaluationFilter & /*filter*/, const uint_t numberOfGhostLayersToInclude )
{
   StreamPullEverything< LatticeModel_T >::execute( src, dst, numberOfGhostLayersToInclude );
}

template< typename LatticeModel_T >
void StreamPull< LatticeModel_T, typename std::enable_if< std::is_same< typename LatticeModel_T::Stencil,
                                                                              stencil::D3Q19 >::value >::type
   >::execute( PdfField_T * src, PdfField_T * dst, cell_idx_t x, cell_idx_t y, cell_idx_t z )
{
   using namespace stencil;

   dst->get(x,y,z, Stencil_T::idx[C])  = src->get(x  , y  , z  , Stencil_T::idx[C]);
   dst->get(x,y,z, Stencil_T::idx[N])  = src->get(x  , y-1, z  , Stencil_T::idx[N]);
   dst->get(x,y,z, Stencil_T::idx[S])  = src->get(x  , y+1, z  , Stencil_T::idx[S]);
   dst->get(x,y,z, Stencil_T::idx[W])  = src->get(x+1, y  , z  , Stencil_T::idx[W]);
   dst->get(x,y,z, Stencil_T::idx[E])  = src->get(x-1, y  , z  , Stencil_T::idx[E]);
   dst->get(x,y,z, Stencil_T::idx[T])  = src->get(x  , y  , z-1, Stencil_T::idx[T]);
   dst->get(x,y,z, Stencil_T::idx[B])  = src->get(x  , y  , z+1, Stencil_T::idx[B]);
   dst->get(x,y,z, Stencil_T::idx[NW]) = src->get(x+1, y-1, z  , Stencil_T::idx[NW]);
   dst->get(x,y,z, Stencil_T::idx[NE]) = src->get(x-1, y-1, z  , Stencil_T::idx[NE]);
   dst->get(x,y,z, Stencil_T::idx[SW]) = src->get(x+1, y+1, z  , Stencil_T::idx[SW]);
   dst->get(x,y,z, Stencil_T::idx[SE]) = src->get(x-1, y+1, z  , Stencil_T::idx[SE]);
   dst->get(x,y,z, Stencil_T::idx[TN]) = src->get(x  , y-1, z-1, Stencil_T::idx[TN]);
   dst->get(x,y,z, Stencil_T::idx[TS]) = src->get(x  , y+1, z-1, Stencil_T::idx[TS]);
   dst->get(x,y,z, Stencil_T::idx[TW]) = src->get(x+1, y  , z-1, Stencil_T::idx[TW]);
   dst->get(x,y,z, Stencil_T::idx[TE]) = src->get(x-1, y  , z-1, Stencil_T::idx[TE]);
   dst->get(x,y,z, Stencil_T::idx[BN]) = src->get(x  , y-1, z+1, Stencil_T::idx[BN]);
   dst->get(x,y,z, Stencil_T::idx[BS]) = src->get(x  , y+1, z+1, Stencil_T::idx[BS]);
   dst->get(x,y,z, Stencil_T::idx[BW]) = src->get(x+1, y  , z+1, Stencil_T::idx[BW]);
   dst->get(x,y,z, Stencil_T::idx[BE]) = src->get(x-1, y  , z+1, Stencil_T::idx[BE]);
}



// D3Q27 //

template< typename LatticeModel_T >
struct StreamPull< LatticeModel_T, typename std::enable_if< std::is_same< typename LatticeModel_T::Stencil,
                                                                                stencil::D3Q27 >::value >::type >
{
   static_assert( (std::is_same< typename LatticeModel_T::Stencil, stencil::D3Q27 >::value), "Only works with D3Q27!" );

   using PdfField_T = PdfField<LatticeModel_T>;
   using Stencil_T = typename LatticeModel_T::Stencil;

   template< typename Filter_T >
   static void execute( PdfField_T * src, PdfField_T * dst, IBlock * block, Filter_T & filter, const uint_t numberOfGhostLayersToInclude = uint_t(0) );

   static void execute( PdfField_T * src, PdfField_T * dst, IBlock * block, walberla::field::DefaultEvaluationFilter & filter,
                        const uint_t numberOfGhostLayersToInclude = uint_t(0) );

   static inline void execute( PdfField_T * src, PdfField_T * dst, cell_idx_t x, cell_idx_t y, cell_idx_t z );
};

template< typename LatticeModel_T >
template< typename Filter_T >
void StreamPull< LatticeModel_T, typename std::enable_if< std::is_same< typename LatticeModel_T::Stencil,
                                                                              stencil::D3Q27 >::value >::type
   >::execute( PdfField_T * src, PdfField_T * dst, IBlock * block, Filter_T & filter, const uint_t numberOfGhostLayersToInclude )
{
   WALBERLA_ASSERT_NOT_NULLPTR( src );
   WALBERLA_ASSERT_NOT_NULLPTR( dst );

   WALBERLA_ASSERT_GREATER( src->nrOfGhostLayers(), numberOfGhostLayersToInclude );
   WALBERLA_ASSERT_GREATER_EQUAL( dst->nrOfGhostLayers(), numberOfGhostLayersToInclude );

   filter( *block );

   WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ( src, numberOfGhostLayersToInclude,

      if( filter(x,y,z) )
         execute(src,dst,x,y,z);

   ) // WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ

   src->swapDataPointers( dst );
}

template< typename LatticeModel_T >
void StreamPull< LatticeModel_T, typename std::enable_if< std::is_same< typename LatticeModel_T::Stencil,
                                                                              stencil::D3Q27 >::value >::type
   >::execute( PdfField_T * src, PdfField_T * dst, IBlock * /*block*/, walberla::field::DefaultEvaluationFilter & /*filter*/, const uint_t numberOfGhostLayersToInclude )
{
   StreamPullEverything< LatticeModel_T >::execute( src, dst, numberOfGhostLayersToInclude );
}

template< typename LatticeModel_T >
void StreamPull< LatticeModel_T, typename std::enable_if< std::is_same< typename LatticeModel_T::Stencil,
                                                                              stencil::D3Q27 >::value >::type
   >::execute( PdfField_T * src, PdfField_T * dst, cell_idx_t x, cell_idx_t y, cell_idx_t z )
{
   using namespace stencil;

   dst->get(x,y,z, Stencil_T::idx[C])   = src->get(x  , y  , z  , Stencil_T::idx[C]);
   dst->get(x,y,z, Stencil_T::idx[N])   = src->get(x  , y-1, z  , Stencil_T::idx[N]);
   dst->get(x,y,z, Stencil_T::idx[S])   = src->get(x  , y+1, z  , Stencil_T::idx[S]);
   dst->get(x,y,z, Stencil_T::idx[W])   = src->get(x+1, y  , z  , Stencil_T::idx[W]);
   dst->get(x,y,z, Stencil_T::idx[E])   = src->get(x-1, y  , z  , Stencil_T::idx[E]);
   dst->get(x,y,z, Stencil_T::idx[T])   = src->get(x  , y  , z-1, Stencil_T::idx[T]);
   dst->get(x,y,z, Stencil_T::idx[B])   = src->get(x  , y  , z+1, Stencil_T::idx[B]);
   dst->get(x,y,z, Stencil_T::idx[NW])  = src->get(x+1, y-1, z  , Stencil_T::idx[NW]);
   dst->get(x,y,z, Stencil_T::idx[NE])  = src->get(x-1, y-1, z  , Stencil_T::idx[NE]);
   dst->get(x,y,z, Stencil_T::idx[SW])  = src->get(x+1, y+1, z  , Stencil_T::idx[SW]);
   dst->get(x,y,z, Stencil_T::idx[SE])  = src->get(x-1, y+1, z  , Stencil_T::idx[SE]);
   dst->get(x,y,z, Stencil_T::idx[TN])  = src->get(x  , y-1, z-1, Stencil_T::idx[TN]);
   dst->get(x,y,z, Stencil_T::idx[TS])  = src->get(x  , y+1, z-1, Stencil_T::idx[TS]);
   dst->get(x,y,z, Stencil_T::idx[TW])  = src->get(x+1, y  , z-1, Stencil_T::idx[TW]);
   dst->get(x,y,z, Stencil_T::idx[TE])  = src->get(x-1, y  , z-1, Stencil_T::idx[TE]);
   dst->get(x,y,z, Stencil_T::idx[BN])  = src->get(x  , y-1, z+1, Stencil_T::idx[BN]);
   dst->get(x,y,z, Stencil_T::idx[BS])  = src->get(x  , y+1, z+1, Stencil_T::idx[BS]);
   dst->get(x,y,z, Stencil_T::idx[BW])  = src->get(x+1, y  , z+1, Stencil_T::idx[BW]);
   dst->get(x,y,z, Stencil_T::idx[BE])  = src->get(x-1, y  , z+1, Stencil_T::idx[BE]);
   dst->get(x,y,z, Stencil_T::idx[TNE]) = src->get(x-1, y-1, z-1, Stencil_T::idx[TNE]);
   dst->get(x,y,z, Stencil_T::idx[TNW]) = src->get(x+1, y-1, z-1, Stencil_T::idx[TNW]);
   dst->get(x,y,z, Stencil_T::idx[TSE]) = src->get(x-1, y+1, z-1, Stencil_T::idx[TSE]);
   dst->get(x,y,z, Stencil_T::idx[TSW]) = src->get(x+1, y+1, z-1, Stencil_T::idx[TSW]);
   dst->get(x,y,z, Stencil_T::idx[BNE]) = src->get(x-1, y-1, z+1, Stencil_T::idx[BNE]);
   dst->get(x,y,z, Stencil_T::idx[BNW]) = src->get(x+1, y-1, z+1, Stencil_T::idx[BNW]);
   dst->get(x,y,z, Stencil_T::idx[BSE]) = src->get(x-1, y+1, z+1, Stencil_T::idx[BSE]);
   dst->get(x,y,z, Stencil_T::idx[BSW]) = src->get(x+1, y+1, z+1, Stencil_T::idx[BSW]);
}



} // namespace lbm
} // namespace walberla
