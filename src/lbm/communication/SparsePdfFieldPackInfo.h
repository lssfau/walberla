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
 * \brief Optimized PackInfo for sparse PDF fields
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
 * Additionally this PackInfo only communicates PDFs in fluid cells.
 *
 * see also documentation for FieldPackInfo
 *
 * \warning For fields with nrOfGhostLayers > 1:  only the components pointing towards
 *          the boundary are communicated, which may not be the desired behavior
 *          for the 'inner' ghost layers
 *
 * \ingroup lbm
 */
template< typename LatticeModel_T, typename FlagField_T >
class SparsePdfFieldPackInfo : public communication::UniformPackInfo
{
public:

   typedef lbm::PdfField<LatticeModel_T>     PdfField_T;
   typedef typename FlagField_T::flag_t      flag_t;
   typedef typename LatticeModel_T::Stencil  Stencil;

   SparsePdfFieldPackInfo( const BlockDataID & pdfFieldId, const BlockDataID & flagFieldId, FlagUID flag, bool flagFieldConstant )
      : pdfFieldId_( pdfFieldId ), flagFieldId_( flagFieldId ), flag_( flag ), flagFieldConstant_( flagFieldConstant ) {}

   virtual ~SparsePdfFieldPackInfo() {}

   bool constantDataExchange() const { return flagFieldConstant_; }
   bool threadsafeReceiving()  const { return true; }

   void unpackData( IBlock * receiver, stencil::Direction dir, mpi::RecvBuffer & buffer );

   void communicateLocal( const IBlock * sender, IBlock * receiver, stencil::Direction dir );

protected:

   void packDataImpl( const IBlock * sender, stencil::Direction dir, mpi::SendBuffer & outBuffer ) const;


   BlockDataID pdfFieldId_;
   BlockDataID flagFieldId_;
   FlagUID     flag_;
   bool        flagFieldConstant_;
};



template< typename LatticeModel_T, typename FlagField_T >
void SparsePdfFieldPackInfo< LatticeModel_T, FlagField_T >::unpackData( IBlock * receiver, stencil::Direction dir, mpi::RecvBuffer & buffer )
{
   if( Stencil::idx[ stencil::inverseDir[dir] ] >= Stencil::Size )
      return;

   PdfField_T * pdfField = receiver->getData< PdfField_T >( pdfFieldId_ );
   WALBERLA_ASSERT_NOT_NULLPTR( pdfField );
   WALBERLA_ASSERT_EQUAL( pdfField->nrOfGhostLayers(), 1 );

   const FlagField_T * flagField = receiver->getData< FlagField_T >( flagFieldId_ );
   WALBERLA_ASSERT_NOT_NULLPTR( flagField );
   WALBERLA_ASSERT_EQUAL( flagField->nrOfGhostLayers(), 1 );

   WALBERLA_ASSERT_EQUAL( flagField->xyzSize(), pdfField->xyzSize() );

   const flag_t mask = flagField->getFlag( flag_ );

   stencil::Direction packerDirection = stencil::inverseDir[dir];
   const auto packerIndex = Stencil::idx[ packerDirection ];

   WALBERLA_ASSERT( packerIndex < Stencil::Size );

   WALBERLA_DEBUG_SECTION()
   {
      uint_t ctr = 0;
      for( auto flagIt = flagField->beginGhostLayerOnlyXYZ(dir); flagIt != flagField->end(); ++flagIt )
      {
         if( *flagIt & mask )
            ctr += Stencil::d_per_d_length[packerIndex];
      }
      uint_t recvCtr = 0;
      buffer >> recvCtr;

      WALBERLA_ASSERT_EQUAL( ctr, recvCtr, "The number of packed values and values to be unpacked do not match!\n" );
   }

   auto pdfIt  =  pdfField->beginGhostLayerOnlyXYZ(dir);
   auto flagIt = flagField->beginGhostLayerOnlyXYZ(dir);

   while( pdfIt != pdfField->end() )
   {
      WALBERLA_ASSERT_UNEQUAL( flagIt, flagField->end() );

      if( *flagIt & mask )
         for(uint_t f = 0; f < Stencil::d_per_d_length[packerIndex]; ++f)
            buffer >> pdfIt.getF( Stencil::idx[ Stencil::d_per_d[packerIndex][f] ] );

      ++pdfIt;
      ++flagIt;
   }

   WALBERLA_ASSERT_EQUAL( flagIt, flagField->end() );
}



template< typename LatticeModel_T, typename FlagField_T >
void SparsePdfFieldPackInfo< LatticeModel_T, FlagField_T >::communicateLocal( const IBlock * sender, IBlock * receiver, stencil::Direction dir )
{
   if( Stencil::idx[dir] >= Stencil::Size )
      return;

   const PdfField_T * senderPdf   = sender  ->getData< PdfField_T >( pdfFieldId_ );
         PdfField_T * receiverPdf = receiver->getData< PdfField_T >( pdfFieldId_ );

   WALBERLA_ASSERT_NOT_NULLPTR( senderPdf );
   WALBERLA_ASSERT_NOT_NULLPTR( receiverPdf );
   WALBERLA_ASSERT_EQUAL( senderPdf->xyzSize(), receiverPdf->xyzSize() );
   WALBERLA_ASSERT_EQUAL( senderPdf->nrOfGhostLayers(), 1 );
   WALBERLA_ASSERT_EQUAL( receiverPdf->nrOfGhostLayers(), 1 );

   const FlagField_T * senderFlagField = sender->getData< FlagField_T >( flagFieldId_ );
   WALBERLA_ASSERT_NOT_NULLPTR( senderFlagField );
   WALBERLA_ASSERT_EQUAL( senderFlagField->xyzSize(), senderPdf->xyzSize() );
   WALBERLA_ASSERT_EQUAL( senderFlagField->nrOfGhostLayers(), 1 );

   const flag_t mask = senderFlagField->getFlag( flag_ );

   typename PdfField_T::const_iterator sendPdfIter = senderPdf->beginSliceBeforeGhostLayerXYZ(dir);
   typename PdfField_T::iterator       recvPdfIter = receiverPdf->beginGhostLayerOnlyXYZ(stencil::inverseDir[dir]);

   typename FlagField_T::const_iterator sendFlagIter = senderFlagField->beginSliceBeforeGhostLayerXYZ(dir);

   const auto dirIndex = Stencil::idx[ dir ];

   WALBERLA_ASSERT( dirIndex < Stencil::Size );

   while( sendPdfIter != senderPdf->end() )
   {
      WALBERLA_ASSERT_UNEQUAL( recvPdfIter, receiverPdf->end() );
      WALBERLA_ASSERT_UNEQUAL( sendFlagIter, senderFlagField->end() );

      if( *sendFlagIter & mask )
         for( uint_t f = 0; f < Stencil::d_per_d_length[dirIndex]; ++f )
            recvPdfIter.getF( Stencil::idx[ Stencil::d_per_d[dirIndex][f] ] ) = sendPdfIter.getF( Stencil::idx[ Stencil::d_per_d[dirIndex][f] ] );

      ++sendPdfIter;
      ++recvPdfIter;
      ++sendFlagIter;
   }

   WALBERLA_ASSERT_EQUAL( sendPdfIter, senderPdf->end() );
   WALBERLA_ASSERT_EQUAL( recvPdfIter, receiverPdf->end() );
   WALBERLA_ASSERT_EQUAL( sendFlagIter, senderFlagField->end() );
}



template< typename LatticeModel_T, typename FlagField_T >
void SparsePdfFieldPackInfo< LatticeModel_T, FlagField_T >::packDataImpl( const IBlock * sender, stencil::Direction dir, mpi::SendBuffer & outBuffer ) const
{
   if( Stencil::idx[dir] >= Stencil::Size )
      return;

   const PdfField_T * pdfField = sender->getData< PdfField_T >( pdfFieldId_ );
   WALBERLA_ASSERT_NOT_NULLPTR( pdfField );
   WALBERLA_ASSERT_EQUAL( pdfField->nrOfGhostLayers(), 1 );

   const FlagField_T * flagField = sender->getData< FlagField_T >( flagFieldId_ );
   WALBERLA_ASSERT_NOT_NULLPTR( flagField );
   WALBERLA_ASSERT_EQUAL( flagField->nrOfGhostLayers(), 1 );

   WALBERLA_ASSERT_EQUAL( flagField->xyzSize(), pdfField->xyzSize() );

   const flag_t mask = flagField->getFlag( flag_ );

   const auto dirIndex = Stencil::idx[ dir ];

   WALBERLA_ASSERT( dirIndex < Stencil::Size );

   WALBERLA_DEBUG_SECTION()
   {
      uint_t ctr = 0;
      for( auto flagIt = flagField->beginSliceBeforeGhostLayerXYZ(dir); flagIt != flagField->end(); ++flagIt )
      {
         if( *flagIt & mask )
            ctr += Stencil::d_per_d_length[dirIndex];
      }
      outBuffer << ctr;
   }

   auto pdfIt  =  pdfField->beginSliceBeforeGhostLayerXYZ(dir);
   auto flagIt = flagField->beginSliceBeforeGhostLayerXYZ(dir);
   while( pdfIt != pdfField->end() )
   {
      WALBERLA_ASSERT_UNEQUAL( flagIt, flagField->end() );

      if( *flagIt & mask )
         for(uint_t f = 0; f < Stencil::d_per_d_length[dirIndex]; ++f)
            outBuffer << pdfIt.getF( Stencil::idx[ Stencil::d_per_d[dirIndex][f] ] );

      ++pdfIt;
      ++flagIt;
   }


   WALBERLA_ASSERT_EQUAL( flagIt, flagField->end() );
}





} // namespace lbm
} // namespace walberla
