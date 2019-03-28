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
//! \file EqualLevelBorderStreamCorrection.h
//! \ingroup lbm
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "lbm/field/PdfField.h"

#include "blockforest/Block.h"
#include "blockforest/BlockNeighborhoodSection.h"

#include "core/cell/CellInterval.h"

#include "stencil/D2CornerStencil.h"
#include "stencil/D3EdgeCornerStencil.h"

#include <type_traits>



namespace walberla {
namespace lbm {
namespace refinement {



namespace internal {
template< typename LatticeModel_T, class Enable = void  >
struct EdgeCornerStencil;
template< typename LatticeModel_T >
struct EdgeCornerStencil< LatticeModel_T, typename std::enable_if< LatticeModel_T::Stencil::D == 2 >::type >
{
   typedef stencil::D2CornerStencil type;
};
template< typename LatticeModel_T >
struct EdgeCornerStencil< LatticeModel_T, typename std::enable_if< LatticeModel_T::Stencil::D == 3 >::type >
{
   typedef stencil::D3EdgeCornerStencil type;
};
}



template< typename LatticeModel_T >
class EqualLevelBorderStreamCorrection
{
public:

   typedef PdfField< LatticeModel_T > PdfField_T;
   typedef typename internal::EdgeCornerStencil< LatticeModel_T >::type EdgeCornerStencil_T;

   EqualLevelBorderStreamCorrection( const BlockDataID & pdfFieldId ) : pdfFieldId_( pdfFieldId ) {}

   void prepare( Block * block ) const { helper( block, true  ); }
   void   apply( Block * block ) const { helper( block, false ); }

protected:

   bool finerNeighborExistsInVicinity( Block * block, stencil::Direction dir ) const;
   void helper( Block * block, const bool prep ) const;

   BlockDataID pdfFieldId_;
};



template< typename LatticeModel_T >
bool EqualLevelBorderStreamCorrection< LatticeModel_T >::finerNeighborExistsInVicinity( Block * block, stencil::Direction dir ) const
{
   Vector3<int> min( -1 );
   Vector3<int> max(  1 );

   for( uint_t i = 0; i != 3; ++i )
   {
      if( stencil::c[i][dir] == -1 ) max[i] = 0;
      if( stencil::c[i][dir] ==  1 ) min[i] = 0;
   }
   if( LatticeModel_T::Stencil::D == uint_t(2) )
      min[2] = max[2] = 0;

   for( int z = min[2]; z <= max[2]; ++z ) {
      for( int y = min[1]; y <= max[1]; ++y ) {
         for( int x = min[0]; x <= max[0]; ++x )
         {
            if( x == 0 && y == 0 && z == 0 )
               continue;
            if( block->neighborhoodSectionHasSmallerBlocks( blockforest::getBlockNeighborhoodSectionIndex(x,y,z) ) )
               return true;
         }
      }
   }

   return false;
}



template< typename LatticeModel_T >
void EqualLevelBorderStreamCorrection< LatticeModel_T >::helper( Block * block, const bool prep ) const
{
   PdfField_T * pdfField = block->template getData< PdfField_T >( pdfFieldId_ );
   const CellInterval xyzSize = pdfField->xyzSize();

   for( auto dir = EdgeCornerStencil_T::beginNoCenter(); dir != EdgeCornerStencil_T::end(); ++dir )
   {
      if( block->neighborhoodSectionHasEquallySizedBlock( blockforest::getBlockNeighborhoodSectionIndex(*dir) ) &&
          finerNeighborExistsInVicinity( block, *dir ) )
      {
         CellInterval region( xyzSize );
         region.expand( cell_idx_t(1) );
         for( uint_t i = 0; i < 3; ++i )
         {
            if( stencil::c[i][*dir] == -1 ) region.max()[i] = region.min()[i];
            else if( stencil::c[i][*dir] == 1 ) region.min()[i] = region.max()[i];
            else
            {
               region.min()[i] += cell_idx_t(1);
               region.max()[i] -= cell_idx_t(1);
            }
         }

         auto invDir = stencil::inverseDir[ *dir ];

         if( prep )
         {
            for( auto cell = pdfField->beginSliceXYZ( region ); cell != pdfField->end(); ++cell ) {
               for( uint_t n = 0; n < PdfField_T::Stencil::d_per_d_length[invDir]; ++n )
               {
                  auto d = PdfField_T::Stencil::d_per_d[invDir][n];
                  auto idx = PdfField_T::Stencil::idx[ d ];
                  cell[idx] = cell.neighbor( d, idx );
               }
            }
         }
         else
         {
            for( auto cell = pdfField->beginSliceXYZ( region ); cell != pdfField->end(); ++cell ) {
               for( uint_t n = 0; n < PdfField_T::Stencil::d_per_d_length[invDir]; ++n )
               {
                  auto d = PdfField_T::Stencil::d_per_d[invDir][n];
                  auto idx = PdfField_T::Stencil::idx[ d ];
                  cell.neighbor( d, idx ) = cell[idx];
               }
            }
         }
      }
   }
}



} // namespace refinement
} // namespace lbm
} // namespace walberla
