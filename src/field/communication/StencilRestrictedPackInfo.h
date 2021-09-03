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
//! \file StencilRestrictedPackInfo.h
//! \ingroup lbm
//! \author Martin Bauer <martin.bauer@fau.de>
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//! \brief Packs only certain components of a field
//
//======================================================================================================================

#pragma once

#include "communication/UniformPackInfo.h"
#include "core/debug/Debug.h"
#include "stencil/Directions.h"


namespace walberla {
namespace field {
namespace communication {


/**
 * \brief PackInfo that packs only the components that point to the neighbor.
 *
 * see also documentation for FieldPackInfo
 *
 * \warning For fields with nrOfGhostLayers > 1:  only the components pointing towards
 *          the boundary are communicated, which may not be the desired behavior
 *          for the 'inner' ghost layers
 */
template< typename GhostLayerField_T, typename Stencil_T >
class StencilRestrictedPackInfo : public walberla::communication::UniformPackInfo
{
public:
   StencilRestrictedPackInfo( const BlockDataID & fieldId ) : fieldId_( fieldId ) {}
   ~StencilRestrictedPackInfo() override = default;

   bool constantDataExchange() const override { return true; }
   bool threadsafeReceiving()  const override { return true; }

   void unpackData( IBlock * receiver, stencil::Direction dir, mpi::RecvBuffer & buffer ) override;

   void communicateLocal( const IBlock * sender, IBlock * receiver, stencil::Direction dir ) override;

protected:

   void packDataImpl( const IBlock * sender, stencil::Direction dir, mpi::SendBuffer & outBuffer ) const override;

   const BlockDataID fieldId_;

   static_assert(GhostLayerField_T::F_SIZE == Stencil_T::Size, "Size of stencil and f size of field have to be equal");
};


template< typename GhostLayerField_T, typename Stencil >
void StencilRestrictedPackInfo<GhostLayerField_T, Stencil>::unpackData( IBlock * receiver, stencil::Direction dir, mpi::RecvBuffer & buffer )
{
   if( Stencil::idx[ stencil::inverseDir[dir] ] >= Stencil::Size )
      return;

   GhostLayerField_T * pdfField = receiver->getData< GhostLayerField_T >( fieldId_ );
   WALBERLA_ASSERT_NOT_NULLPTR( pdfField );
   WALBERLA_ASSERT_EQUAL( pdfField->nrOfGhostLayers(), 1 );

   stencil::Direction packerDirection = stencil::inverseDir[dir];

   for(auto i = pdfField->beginGhostLayerOnlyXYZ(dir); i != pdfField->end(); ++i )
      for(uint_t f = 0; f < Stencil::d_per_d_length[packerDirection]; ++f)
         buffer >> i.getF( Stencil::idx[ Stencil::d_per_d[packerDirection][f] ] );
}



template< typename GhostLayerField_T, typename Stencil >
void StencilRestrictedPackInfo<GhostLayerField_T, Stencil>::communicateLocal( const IBlock * sender, IBlock * receiver, stencil::Direction dir )
{
   if( Stencil::idx[dir] >= Stencil::Size )
      return;

   const GhostLayerField_T * sf = sender  ->getData< GhostLayerField_T >( fieldId_ );
         GhostLayerField_T * rf = receiver->getData< GhostLayerField_T >( fieldId_ );

   WALBERLA_ASSERT_EQUAL( sf->xyzSize(), rf->xyzSize() );

   typename GhostLayerField_T::const_iterator srcIter = sf->beginSliceBeforeGhostLayerXYZ(dir);
   typename GhostLayerField_T::iterator       dstIter = rf->beginGhostLayerOnlyXYZ(stencil::inverseDir[dir]);

   while( srcIter != sf->end() )
   {
      for( uint_t f = 0; f < Stencil::d_per_d_length[dir]; ++f )
         dstIter.getF( Stencil::idx[ Stencil::d_per_d[dir][f] ] ) = srcIter.getF( Stencil::idx[ Stencil::d_per_d[dir][f] ] );

      ++srcIter;
      ++dstIter;
   }
   WALBERLA_ASSERT( srcIter == sf->end() );
   WALBERLA_ASSERT( dstIter == rf->end() );
}



template< typename GhostLayerField_T, typename Stencil >
void StencilRestrictedPackInfo<GhostLayerField_T, Stencil>::packDataImpl( const IBlock * sender, stencil::Direction dir, mpi::SendBuffer & outBuffer ) const
{
   if( Stencil::idx[dir] >= Stencil::Size )
      return;

   const GhostLayerField_T * pdfField = sender->getData< GhostLayerField_T >( fieldId_ );
   WALBERLA_ASSERT_NOT_NULLPTR( pdfField );
   WALBERLA_ASSERT_EQUAL( pdfField->nrOfGhostLayers(), 1 );

   for( auto i = pdfField->beginSliceBeforeGhostLayerXYZ(dir); i != pdfField->end(); ++i )
      for(uint_t f = 0; f < Stencil::d_per_d_length[dir]; ++f)
         outBuffer << i.getF( Stencil::idx[ Stencil::d_per_d[dir][f] ] );
}



} // namespace communication
} // namespace field
} // namespace walberla
