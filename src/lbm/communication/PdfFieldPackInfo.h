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
//! \file PdfFieldPackInfo.h
//! \ingroup lbm
//! \author Martin Bauer <martin.bauer@fau.de>
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//! \brief Optimized packing for PDF fields, by communicating only required components.
//
//======================================================================================================================

#pragma once

#include "lbm/field/PdfField.h"
#include "communication/UniformPackInfo.h"
#include "core/debug/Debug.h"
#include "stencil/Directions.h"


namespace walberla {
namespace lbm {



/**
 * \brief Optimized PackInfo for PDF fields (communicates only components pointing to neighbor)
 *
 * In principle a PDF field can be communicated using a FieldPackInfo which
 * copies all entries (=values for f) of the ghost layer.
 *
 * For PDF fields, however, it is sufficient to communicate only those discrete velocities
 * which point into the direction of the neighboring block
 *
 * This PackInfo sends only these required components and therefore the sizes of the
 * messages decrease.
 *
 * see also documentation for FieldPackInfo
 *
 * \warning For fields with nrOfGhostLayers > 1:  only the components pointing towards
 *          the boundary are communicated, which may not be the desired behavior
 *          for the 'inner' ghost layers
 *
 * \ingroup lbm
 */
template< typename LatticeModel_T >
class PdfFieldPackInfo : public walberla::communication::UniformPackInfo
{
public:

   typedef PdfField<LatticeModel_T>          PdfField_T;
   typedef typename LatticeModel_T::Stencil  Stencil;

   PdfFieldPackInfo( const BlockDataID & pdfFieldId ) : pdfFieldId_( pdfFieldId ) {}
   virtual ~PdfFieldPackInfo() {}

   bool constantDataExchange() const { return true; }
   bool threadsafeReceiving()  const { return true; }

   void unpackData( IBlock * receiver, stencil::Direction dir, mpi::RecvBuffer & buffer );

   void communicateLocal( const IBlock * sender, IBlock * receiver, stencil::Direction dir );

protected:

   void packDataImpl( const IBlock * sender, stencil::Direction dir, mpi::SendBuffer & outBuffer ) const;



   const BlockDataID pdfFieldId_;
};



template< typename LatticeModel_T >
void PdfFieldPackInfo< LatticeModel_T >::unpackData( IBlock * receiver, stencil::Direction dir, mpi::RecvBuffer & buffer )
{
   if( Stencil::idx[ stencil::inverseDir[dir] ] >= Stencil::Size )
      return;

   PdfField_T * pdfField = receiver->getData< PdfField_T >( pdfFieldId_ );
   WALBERLA_ASSERT_NOT_NULLPTR( pdfField );
   WALBERLA_ASSERT_EQUAL( pdfField->nrOfGhostLayers(), 1 );

   stencil::Direction packerDirection = stencil::inverseDir[dir];

   for(auto i = pdfField->beginGhostLayerOnlyXYZ(dir); i != pdfField->end(); ++i )
      for(uint_t f = 0; f < Stencil::d_per_d_length[packerDirection]; ++f)
         buffer >> i.getF( Stencil::idx[ Stencil::d_per_d[packerDirection][f] ] );
}



template< typename LatticeModel_T >
void PdfFieldPackInfo< LatticeModel_T >::communicateLocal( const IBlock * sender, IBlock * receiver, stencil::Direction dir )
{
   if( Stencil::idx[dir] >= Stencil::Size )
      return;

   const PdfField_T * sf = sender  ->getData< PdfField_T >( pdfFieldId_ );
         PdfField_T * rf = receiver->getData< PdfField_T >( pdfFieldId_ );

   WALBERLA_ASSERT_EQUAL( sf->xyzSize(), rf->xyzSize() );

   typename PdfField_T::const_iterator srcIter = sf->beginSliceBeforeGhostLayerXYZ(dir);
   typename PdfField_T::iterator       dstIter = rf->beginGhostLayerOnlyXYZ(stencil::inverseDir[dir]);

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



template< typename LatticeModel_T >
void PdfFieldPackInfo< LatticeModel_T >::packDataImpl( const IBlock * sender, stencil::Direction dir, mpi::SendBuffer & outBuffer ) const
{
   if( Stencil::idx[dir] >= Stencil::Size )
      return;

   const PdfField_T * pdfField = sender->getData< PdfField_T >( pdfFieldId_ );
   WALBERLA_ASSERT_NOT_NULLPTR( pdfField );
   WALBERLA_ASSERT_EQUAL( pdfField->nrOfGhostLayers(), 1 );

   for( auto i = pdfField->beginSliceBeforeGhostLayerXYZ(dir); i != pdfField->end(); ++i )
      for(uint_t f = 0; f < Stencil::d_per_d_length[dir]; ++f)
         outBuffer << i.getF( Stencil::idx[ Stencil::d_per_d[dir][f] ] );
}



} // namespace lbm
} // namespace walberla
