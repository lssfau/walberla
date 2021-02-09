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
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "TimeStepPdfPackInfo.h"

#include "lbm/field/PdfField.h"
#include "blockforest/BlockNeighborhoodSection.h"
#include "core/cell/CellInterval.h"
#include "stencil/Directions.h"


namespace walberla {
namespace lbm {
namespace refinement {



#ifdef NDEBUG
template< typename LatticeModel_T >
#else
template< typename LatticeModel_T, typename BoundaryHandling_T >
#endif
class PdfFieldPackInfo : public TimeStepPdfPackInfo
{
public:

   typedef PdfField<LatticeModel_T>          PdfField_T;
   typedef typename LatticeModel_T::Stencil  Stencil;

#ifdef NDEBUG   
   PdfFieldPackInfo( const BlockDataID & pdfFieldId, const bool _optimizedEqualLevelCommunication = true,
                     const bool _optimizedForLinearExplosion = true ) :
      pdfFieldId_( pdfFieldId ), optimizedEqualLevelCommunication_( _optimizedEqualLevelCommunication ),
      optimizedForLinearExplosion_( _optimizedForLinearExplosion ), equalLevelCells_( equalLevelCells() ) {}
#else
   PdfFieldPackInfo( const BlockDataID & pdfFieldId, const ConstBlockDataID & boundaryHandlingId,
                     const bool _optimizedEqualLevelCommunication = true, const bool _optimizedForLinearExplosion = true ) :
      pdfFieldId_( pdfFieldId ), boundaryHandlingId_( boundaryHandlingId ),
      optimizedEqualLevelCommunication_( _optimizedEqualLevelCommunication ), optimizedForLinearExplosion_( _optimizedForLinearExplosion ),
      equalLevelCells_( equalLevelCells() ) {}
#endif

   ~PdfFieldPackInfo() override = default;

   bool optimizedEqualLevelCommunication() const override { return optimizedEqualLevelCommunication_; }
   void optimizeEqualLevelCommunication( const bool value = true ) override { optimizedEqualLevelCommunication_ = value; }
   
   bool optimizedForLinearExplosion() const override { return optimizedForLinearExplosion_; }
   void optimizeForLinearExplosion( const bool value = true ) override { optimizedForLinearExplosion_ = value; }
   
   bool constantDataExchange() const override { return true; }
   bool threadsafeReceiving()  const override { return true; }

   void       unpackDataEqualLevel( Block * receiver, stencil::Direction dir, mpi::RecvBuffer & buffer ) override;
   void communicateLocalEqualLevel( const Block * sender, Block * receiver, stencil::Direction dir ) override;

   void       unpackDataCoarseToFine( Block * fineReceiver, const BlockID & coarseSender, stencil::Direction dir, mpi::RecvBuffer & buffer ) override;
   void communicateLocalCoarseToFine( const Block * coarseSender, Block * fineReceiver, stencil::Direction dir ) override;

   void       unpackDataFineToCoarse( Block * coarseReceiver, const BlockID & fineSender, stencil::Direction dir, mpi::RecvBuffer & buffer ) override;
   void communicateLocalFineToCoarse( const Block * fineSender, Block * coarseReceiver, stencil::Direction dir ) override;

protected:

   void packDataEqualLevelImpl( const Block * sender, stencil::Direction dir, mpi::SendBuffer & buffer ) const override;
   void packDataCoarseToFineImpl( const Block * coarseSender, const BlockID &   fineReceiver, stencil::Direction dir, mpi::SendBuffer & buffer ) const override;
   void packDataFineToCoarseImpl( const Block *   fineSender, const BlockID & coarseReceiver, stencil::Direction dir, mpi::SendBuffer & buffer ) const override;

   ///////////////////////////////////////////////////////////////////////
   // Helper functions for determining packing/unpacking cell intervals //
   ///////////////////////////////////////////////////////////////////////

   static inline CellInterval equalLevelPackInterval  ( stencil::Direction dir, const CellInterval & cellBB, const uint_t numberOfLayers );
   static        CellInterval equalLevelUnpackInterval( stencil::Direction dir, const CellInterval & cellBB, const uint_t numberOfLayers );

   inline bool equalLevelFaceIntervalSplitable( const CellInterval & interval, stencil::Direction dir ) const;
   inline std::vector< CellInterval > splitEqualLevelFaceInterval( const CellInterval & interval, stencil::Direction dir ) const; // dir: from sender to receiver

   static inline Vector3< cell_idx_t > getNeighborShift( const BlockID & smallBlock, stencil::Direction dir ); // dir: direction from big to small block

   static CellInterval coarseToFinePackInterval  ( stencil::Direction dir, const CellInterval & cellBB, const BlockID & smallBlock );
   static CellInterval coarseToFineUnpackInterval( stencil::Direction dir, const CellInterval & cellBB, const BlockID & smallBlock );

   static inline  CellInterval fineToCoarsePackInterval  ( stencil::Direction dir, const CellInterval & cellBB );
   static         CellInterval fineToCoarseUnpackInterval( stencil::Direction dir, const CellInterval & cellBB, const BlockID & smallBlock );

   //////////////////////////////
   // General helper functions //
   //////////////////////////////

   static inline uint_t equalLevelCells();
   
   static inline bool isFaceDirection( stencil::Direction dir );
   static inline bool isEdgeDirection( stencil::Direction dir );
   static inline bool isCornerDirection( stencil::Direction dir );

   static inline bool blocksConnectedByFaces( const Block * block, const BlockID & neighbor );
   static inline bool blocksConnectedByEdges( const Block * block, const BlockID & neighbor );

   static inline bool divisibleByTwo( const CellInterval & cellBB );

   static bool coarserNeighborExistsInVicinity( const Block * block, stencil::Direction dir );



   BlockDataID pdfFieldId_;
#ifndef NDEBUG
   ConstBlockDataID boundaryHandlingId_;
#endif
   
   bool optimizedEqualLevelCommunication_;
   bool optimizedForLinearExplosion_;
   
   uint_t equalLevelCells_;
};



/////////////////
// Equal level //
/////////////////

#ifdef NDEBUG
template< typename LatticeModel_T >
void PdfFieldPackInfo< LatticeModel_T >::packDataEqualLevelImpl( const Block * sender, stencil::Direction dir, mpi::SendBuffer & buffer ) const
#else
template< typename LatticeModel_T, typename BoundaryHandling_T >
void PdfFieldPackInfo< LatticeModel_T, BoundaryHandling_T >::packDataEqualLevelImpl( const Block * sender, stencil::Direction dir,
                                                                                     mpi::SendBuffer & buffer ) const
#endif
{
#ifndef NDEBUG
   if( Stencil::D == uint_t(2) )
      WALBERLA_ASSERT_EQUAL( stencil::cz[dir], 0 );
#endif

   const bool optimizedCommunication = ( optimizedEqualLevelCommunication_ && !coarserNeighborExistsInVicinity( sender, dir ) );

   if( optimizedCommunication && Stencil::d_per_d_length[dir] == uint_t(0) )
      return;

   const PdfField_T * field = sender->getData< PdfField_T >( pdfFieldId_ );

   if( optimizedCommunication )
   {
      CellInterval packingInterval = equalLevelPackInterval( dir, field->xyzSize(), uint_t(1) );

      for( auto cell = field->beginSliceXYZ( packingInterval ); cell != field->end(); ++cell )
         for( uint_t d = 0; d < Stencil::d_per_d_length[dir]; ++d )
            buffer << cell.getF( Stencil::idx[ Stencil::d_per_d[dir][d] ] );

      /*
       * --> fzyx (sadly, seems to be slower ... ?)
       *
      for( auto z = packingInterval.zMin(); z <= packingInterval.zMax(); ++z ) {
         for( auto y = packingInterval.yMin(); y <= packingInterval.yMax(); ++y ) {
            for( uint_t d = 0; d < Stencil::d_per_d_length[dir]; ++d ) {

               auto * s = &( field->get( packingInterval.xMin(), y, z, Stencil::d_per_d[dir][d] ) );

               for( uint_t x = 0; x < packingInterval.xSize(); ++x )
                  buffer << s[x];
            }
         }
      }
      */
   }
   else
   {
      CellInterval packingInterval = equalLevelPackInterval( dir, field->xyzSize(), equalLevelCells_ );

      if( optimizedEqualLevelCommunication_ && isFaceDirection( dir ) && Stencil::D == uint_t(3) )
      {
         WALBERLA_ASSERT( equalLevelFaceIntervalSplitable( packingInterval, dir ) );

         std::vector< CellInterval > intervals = splitEqualLevelFaceInterval( packingInterval, dir );
         WALBERLA_ASSERT_EQUAL( intervals.size(), 5 );

         for( uint_t i = 0; i < 4; ++i )
         {
            for( auto cell = field->beginSliceXYZ( intervals[i] ); cell != field->end(); ++cell )
               for( uint_t d = 0; d < Stencil::Size; ++d )
                  buffer << cell.getF( d );
         }

         if( Stencil::d_per_d_length[dir] > uint_t(0) )
            for( auto cell = field->beginSliceXYZ( intervals[4] ); cell != field->end(); ++cell )
               for( uint_t d = 0; d < Stencil::d_per_d_length[dir]; ++d )
                  buffer << cell.getF( Stencil::idx[ Stencil::d_per_d[dir][d] ] );
      }
      else
      {
         for( auto cell = field->beginSliceXYZ( packingInterval ); cell != field->end(); ++cell )
            for( uint_t d = 0; d < Stencil::Size; ++d )
               buffer << cell.getF( d );

         /*
          * --> fzyx (sadly, seems to be slower ... ?)
          *
         for( auto z = packingInterval.zMin(); z <= packingInterval.zMax(); ++z ) {
            for( auto y = packingInterval.yMin(); y <= packingInterval.yMax(); ++y ) {
               for( uint_t d = 0; d < Stencil::Size; ++d ) {

                  auto * s = &( field->get( packingInterval.xMin(), y, z, d ) );

                  for( uint_t x = 0; x < packingInterval.xSize(); ++x )
                     buffer << s[x];
               }
            }
         }
         */
      }
   }
}



#ifdef NDEBUG
template< typename LatticeModel_T >
void PdfFieldPackInfo< LatticeModel_T >::unpackDataEqualLevel( Block * receiver, stencil::Direction dir, mpi::RecvBuffer & buffer )
#else
template< typename LatticeModel_T, typename BoundaryHandling_T >
void PdfFieldPackInfo< LatticeModel_T, BoundaryHandling_T >::unpackDataEqualLevel( Block * receiver, stencil::Direction dir,
                                                                                   mpi::RecvBuffer & buffer )
#endif
{
#ifndef NDEBUG
   if( Stencil::D == uint_t(2) )
      WALBERLA_ASSERT_EQUAL( stencil::cz[dir], 0 );
#endif

   const bool optimizedCommunication = ( optimizedEqualLevelCommunication_ && !coarserNeighborExistsInVicinity( receiver, dir ) );

   const auto invDir = stencil::inverseDir[dir];

   if( optimizedCommunication && Stencil::d_per_d_length[invDir] == uint_t(0) )
      return;

   PdfField_T * field = receiver->getData< PdfField_T >( pdfFieldId_ );

   if( optimizedCommunication )
   {
      CellInterval unpackingInterval = equalLevelUnpackInterval( dir, field->xyzSize(), uint_t(1) );

      for( auto cell = field->beginSliceXYZ( unpackingInterval ); cell != field->end(); ++cell )
         for( uint_t d = 0; d < Stencil::d_per_d_length[invDir]; ++d )
            buffer >> cell.getF( Stencil::idx[ Stencil::d_per_d[invDir][d] ] );

      /*
       * --> fzyx (sadly, seems to be slower ... ?)
       *
      for( auto z = unpackingInterval.zMin(); z <= unpackingInterval.zMax(); ++z ) {
         for( auto y = unpackingInterval.yMin(); y <= unpackingInterval.yMax(); ++y ) {
            for( uint_t d = 0; d < Stencil::d_per_d_length[invDir]; ++d ) {

               auto * r = &( field->get( unpackingInterval.xMin(), y, z, Stencil::d_per_d[invDir][d] ) );

               for( uint_t x = 0; x < unpackingInterval.xSize(); ++x )
                  buffer >> r[x];
            }
         }
      }
      */
   }
   else
   {
      CellInterval unpackingInterval = equalLevelUnpackInterval( dir, field->xyzSize(), equalLevelCells_ );

      if( optimizedEqualLevelCommunication_ && isFaceDirection( dir ) && Stencil::D == uint_t(3) )
      {
         WALBERLA_ASSERT( equalLevelFaceIntervalSplitable( unpackingInterval, invDir ) );

         std::vector< CellInterval > intervals = splitEqualLevelFaceInterval( unpackingInterval, invDir );
         WALBERLA_ASSERT_EQUAL( intervals.size(), 5 );

         for( uint_t i = 0; i < 4; ++i )
         {
            for( auto cell = field->beginSliceXYZ( intervals[i] ); cell != field->end(); ++cell )
               for( uint_t d = 0; d < Stencil::Size; ++d )
                  buffer >> cell.getF( d );
         }

         if( Stencil::d_per_d_length[invDir] > uint_t(0) )
            for( auto cell = field->beginSliceXYZ( intervals[4] ); cell != field->end(); ++cell )
               for( uint_t d = 0; d < Stencil::d_per_d_length[invDir]; ++d )
                  buffer >> cell.getF( Stencil::idx[ Stencil::d_per_d[invDir][d] ] );
      }
      else
      {
         for( auto cell = field->beginSliceXYZ( unpackingInterval ); cell != field->end(); ++cell )
            for( uint_t d = 0; d < Stencil::Size; ++d )
               buffer >> cell.getF( d );

         /*
          * --> fzyx (sadly, seems to be slower ... ?)
          *
         for( auto z = unpackingInterval.zMin(); z <= unpackingInterval.zMax(); ++z ) {
            for( auto y = unpackingInterval.yMin(); y <= unpackingInterval.yMax(); ++y ) {
               for( uint_t d = 0; d < Stencil::Size; ++d ) {

                  auto * r = &( field->get( unpackingInterval.xMin(), y, z, d ) );

                  for( uint_t x = 0; x < unpackingInterval.xSize(); ++x )
                     buffer >> r[x];
               }
            }
         }
         */
      }
   }
}



#ifdef NDEBUG
template< typename LatticeModel_T >
void PdfFieldPackInfo< LatticeModel_T >::communicateLocalEqualLevel( const Block * sender, Block * receiver, stencil::Direction dir )
#else
template< typename LatticeModel_T, typename BoundaryHandling_T >
void PdfFieldPackInfo< LatticeModel_T, BoundaryHandling_T >::communicateLocalEqualLevel( const Block * sender, Block * receiver, stencil::Direction dir )
#endif
{
#ifndef NDEBUG
   if( Stencil::D == uint_t(2) )
      WALBERLA_ASSERT_EQUAL( stencil::cz[dir], 0 );
#endif

   const bool optimizedCommunication = ( optimizedEqualLevelCommunication_ && !coarserNeighborExistsInVicinity( sender, dir ) );

   if( optimizedCommunication && Stencil::d_per_d_length[dir] == uint_t(0) )
      return;

   const PdfField_T * sf =   sender->getData< PdfField_T >( pdfFieldId_ );
         PdfField_T * rf = receiver->getData< PdfField_T >( pdfFieldId_ );

   WALBERLA_ASSERT_EQUAL( sf->xyzSize(), rf->xyzSize() );

   if( optimizedCommunication )
   {
      CellInterval   packingInterval = equalLevelPackInterval( dir, sf->xyzSize(), uint_t(1) );
      CellInterval unpackingInterval = equalLevelUnpackInterval( stencil::inverseDir[dir], rf->xyzSize(), uint_t(1) );

      WALBERLA_ASSERT_EQUAL( packingInterval.xSize(), unpackingInterval.xSize() );
      WALBERLA_ASSERT_EQUAL( packingInterval.ySize(), unpackingInterval.ySize() );
      WALBERLA_ASSERT_EQUAL( packingInterval.zSize(), unpackingInterval.zSize() );

      auto sCell = sf->beginSliceXYZ(   packingInterval );
      auto rCell = rf->beginSliceXYZ( unpackingInterval );
      while( sCell != sf->end() )
      {
         WALBERLA_ASSERT( rCell != rf->end() );

         for( uint_t d = 0; d < Stencil::d_per_d_length[dir]; ++d )
         {
            const auto idx = Stencil::idx[ Stencil::d_per_d[dir][d] ];
            rCell.getF( idx ) = sCell.getF( idx );
         }

         ++sCell;
         ++rCell;
      }
      WALBERLA_ASSERT( rCell == rf->end() );

      /*
       * --> fzyx (sadly, seems to be slower ... ?)
       *
      for( auto zp = packingInterval.zMin(), zu = unpackingInterval.zMin(); zp <= packingInterval.zMax(); ++zp, ++zu ) {
         for( auto yp = packingInterval.yMin(), yu = unpackingInterval.yMin(); yp <= packingInterval.yMax(); ++yp, ++yu ) {
            for( uint_t d = 0; d < Stencil::d_per_d_length[dir]; ++d ) {

               auto * s = &( sf->get(   packingInterval.xMin(), yp, zp, Stencil::d_per_d[dir][d] ) );
               auto * r = &( rf->get( unpackingInterval.xMin(), yu, zu, Stencil::d_per_d[dir][d] ) );

               std::copy( s, s+packingInterval.xSize(), r );
            }
         }
      }
      */
   }
   else
   {
      CellInterval   packingInterval = equalLevelPackInterval( dir, sf->xyzSize(), equalLevelCells_ );
      CellInterval unpackingInterval = equalLevelUnpackInterval( stencil::inverseDir[dir], rf->xyzSize(), equalLevelCells_ );

      WALBERLA_ASSERT_EQUAL( packingInterval.xSize(), unpackingInterval.xSize() );
      WALBERLA_ASSERT_EQUAL( packingInterval.ySize(), unpackingInterval.ySize() );
      WALBERLA_ASSERT_EQUAL( packingInterval.zSize(), unpackingInterval.zSize() );

      if( optimizedEqualLevelCommunication_ && isFaceDirection( dir ) && Stencil::D == uint_t(3) )
      {
         WALBERLA_ASSERT( equalLevelFaceIntervalSplitable(   packingInterval, dir ) );
         WALBERLA_ASSERT( equalLevelFaceIntervalSplitable( unpackingInterval, dir ) );

         std::vector< CellInterval > pIntervals = splitEqualLevelFaceInterval(   packingInterval, dir );
         std::vector< CellInterval > uIntervals = splitEqualLevelFaceInterval( unpackingInterval, dir );
         WALBERLA_ASSERT_EQUAL( pIntervals.size(), 5 );
         WALBERLA_ASSERT_EQUAL( uIntervals.size(), 5 );

         for( uint_t i = 0; i < 4; ++i )
         {
            auto sCell = sf->beginSliceXYZ( pIntervals[i] );
            auto rCell = rf->beginSliceXYZ( uIntervals[i] );
            while( sCell != sf->end() )
            {
               WALBERLA_ASSERT( rCell != rf->end() );

               for( uint_t d = 0; d < Stencil::Size; ++d )
                  rCell.getF( d ) = sCell.getF( d );

               ++sCell;
               ++rCell;
            }
            WALBERLA_ASSERT( rCell == rf->end() );
         }

         if( Stencil::d_per_d_length[dir] > uint_t(0) )
         {
            auto sCell = sf->beginSliceXYZ( pIntervals[4] );
            auto rCell = rf->beginSliceXYZ( uIntervals[4] );
            while( sCell != sf->end() )
            {
               WALBERLA_ASSERT( rCell != rf->end() );

               for( uint_t d = 0; d < Stencil::d_per_d_length[dir]; ++d )
               {
                  const auto idx = Stencil::idx[ Stencil::d_per_d[dir][d] ];
                  rCell.getF( idx ) = sCell.getF( idx );
               }

               ++sCell;
               ++rCell;
            }
            WALBERLA_ASSERT( rCell == rf->end() );
         }
      }
      else
      {
         auto sCell = sf->beginSliceXYZ(   packingInterval );
         auto rCell = rf->beginSliceXYZ( unpackingInterval );
         while( sCell != sf->end() )
         {
            WALBERLA_ASSERT( rCell != rf->end() );

            for( uint_t d = 0; d < Stencil::Size; ++d )
               rCell.getF( d ) = sCell.getF( d );

            ++sCell;
            ++rCell;
         }
         WALBERLA_ASSERT( rCell == rf->end() );

         /*
          * --> fzyx (sadly, seems to be slower ... ?)
          *
         for( auto zp = packingInterval.zMin(), zu = unpackingInterval.zMin(); zp <= packingInterval.zMax(); ++zp, ++zu ) {
            for( auto yp = packingInterval.yMin(), yu = unpackingInterval.yMin(); yp <= packingInterval.yMax(); ++yp, ++yu ) {
               for( uint_t d = 0; d < Stencil::Size; ++d ) {

                  auto * s = &( sf->get(   packingInterval.xMin(), yp, zp, d ) );
                  auto * r = &( rf->get( unpackingInterval.xMin(), yu, zu, d ) );

                  std::copy( s, s+packingInterval.xSize(), r );
               }
            }
         }
         */
      }
   }
}



////////////////////
// Coarse to fine //
////////////////////

#ifdef NDEBUG
template< typename LatticeModel_T >
void PdfFieldPackInfo< LatticeModel_T >::packDataCoarseToFineImpl( const Block * coarseSender, const BlockID & fineReceiver,
                                                                   stencil::Direction dir, mpi::SendBuffer & buffer ) const
#else
template< typename LatticeModel_T, typename BoundaryHandling_T >
void PdfFieldPackInfo< LatticeModel_T, BoundaryHandling_T >::packDataCoarseToFineImpl( const Block * coarseSender, const BlockID & fineReceiver,
                                                                                       stencil::Direction dir, mpi::SendBuffer & buffer ) const
#endif
{
#ifndef NDEBUG
   if( Stencil::D == uint_t(2) )
      WALBERLA_ASSERT_EQUAL( stencil::cz[dir], 0 );
#endif

   const PdfField_T * field = coarseSender->getData< PdfField_T >( pdfFieldId_ );

#ifndef NDEBUG
   const BoundaryHandling_T * boundaryHandling = coarseSender->template getData< BoundaryHandling_T >( boundaryHandlingId_ );
#endif

   CellInterval packingInterval = coarseToFinePackInterval( dir, field->xyzSize(), fineReceiver );

   for( cell_idx_t z = packingInterval.zMin(); z <= packingInterval.zMax(); ++z ) {
      for( cell_idx_t y = packingInterval.yMin(); y <= packingInterval.yMax(); ++y ) {
         for( cell_idx_t x = packingInterval.xMin(); x <= packingInterval.xMax(); ++x ) {
            for( uint_t d = 0; d < Stencil::Size; ++d )
          //for( uint_t d = 0; d < Stencil::d_per_d_length[dir]; ++d )
            {
               buffer << field->get( x, y, z, d );
             //buffer << field->get( x, y, z, Stencil::idx[ Stencil::d_per_d[dir][d] ] );

#ifndef NDEBUG
               if( boundaryHandling->isDomain(x,y,z) )
                  WALBERLA_ASSERT( !math::isnan( field->get( x, y, z, d ) ) );
#endif
            }
         }
      }
   }
}



#ifdef NDEBUG
template< typename LatticeModel_T >
void PdfFieldPackInfo< LatticeModel_T >::unpackDataCoarseToFine( Block * fineReceiver, const BlockID & /*coarseSender*/,
                                                                 stencil::Direction dir, mpi::RecvBuffer & buffer )
#else
template< typename LatticeModel_T, typename BoundaryHandling_T >
void PdfFieldPackInfo< LatticeModel_T, BoundaryHandling_T >::unpackDataCoarseToFine( Block * fineReceiver, const BlockID & /*coarseSender*/,
                                                                                     stencil::Direction dir, mpi::RecvBuffer & buffer )
#endif
{
#ifndef NDEBUG
   if( Stencil::D == uint_t(2) )
      WALBERLA_ASSERT_EQUAL( stencil::cz[dir], 0 );
#endif

   PdfField_T * field = fineReceiver->getData< PdfField_T >( pdfFieldId_ );

   CellInterval unpackingInterval = coarseToFineUnpackInterval( dir, field->xyzSize(), fineReceiver->getId() );

   if( optimizedForLinearExplosion_ )
   {
      for( cell_idx_t z = unpackingInterval.zMin(); z <= unpackingInterval.zMax(); z += cell_idx_t(2) ) {
         for( cell_idx_t y = unpackingInterval.yMin(); y <= unpackingInterval.yMax(); y += cell_idx_t(2) ) {
            for( cell_idx_t x = unpackingInterval.xMin(); x <= unpackingInterval.xMax(); x += cell_idx_t(2) ) {
               for( uint_t idx = 0; idx < Stencil::Size; ++idx )
               {
                  typename PdfField_T::value_type value;
                  buffer >> value;
                  field->get( x, y, z, idx ) = value;
               }
            }
         }
      }
   }
   else
   {
      for( cell_idx_t z = unpackingInterval.zMin(); z <= unpackingInterval.zMax(); z += cell_idx_t(2) ) {
         for( cell_idx_t y = unpackingInterval.yMin(); y <= unpackingInterval.yMax(); y += cell_idx_t(2) ) {
            for( cell_idx_t x = unpackingInterval.xMin(); x <= unpackingInterval.xMax(); x += cell_idx_t(2) ) {
               for( uint_t idx = 0; idx < Stencil::Size; ++idx )
             //for( uint_t d = 0; d < Stencil::d_per_d_length[ stencil::inverseDir[dir] ]; ++d )
               {
                  typename PdfField_T::value_type value;
                  buffer >> value;

                //const auto idx = Stencil::idx[ Stencil::d_per_d[ stencil::inverseDir[dir] ][d] ];

                  field->get( x,                 y,                 z,                 idx ) = value;
                  field->get( x + cell_idx_t(1), y,                 z,                 idx ) = value;
                  field->get( x,                 y + cell_idx_t(1), z,                 idx ) = value;
                  field->get( x + cell_idx_t(1), y + cell_idx_t(1), z,                 idx ) = value;
                  if( Stencil::D == uint_t(3) )
                  {
                     field->get( x,                 y,                 z + cell_idx_t(1), idx ) = value;
                     field->get( x + cell_idx_t(1), y,                 z + cell_idx_t(1), idx ) = value;
                     field->get( x,                 y + cell_idx_t(1), z + cell_idx_t(1), idx ) = value;
                     field->get( x + cell_idx_t(1), y + cell_idx_t(1), z + cell_idx_t(1), idx ) = value;
                  }
               }
            }
         }
      }
   }
}



#ifdef NDEBUG
template< typename LatticeModel_T >
void PdfFieldPackInfo< LatticeModel_T >::communicateLocalCoarseToFine( const Block * coarseSender, Block * fineReceiver, stencil::Direction dir )
#else
template< typename LatticeModel_T, typename BoundaryHandling_T >
void PdfFieldPackInfo< LatticeModel_T, BoundaryHandling_T >::communicateLocalCoarseToFine( const Block * coarseSender, Block * fineReceiver,
                                                                                           stencil::Direction dir )
#endif
{
#ifndef NDEBUG
   if( Stencil::D == uint_t(2) )
      WALBERLA_ASSERT_EQUAL( stencil::cz[dir], 0 );
#endif

   const PdfField_T * sf = coarseSender->getData< PdfField_T >( pdfFieldId_ );
         PdfField_T * rf = fineReceiver->getData< PdfField_T >( pdfFieldId_ );

#ifndef NDEBUG
   const BoundaryHandling_T * boundaryHandling = coarseSender->template getData< BoundaryHandling_T >( boundaryHandlingId_ );
#endif

   WALBERLA_ASSERT_EQUAL( sf->xyzSize(), rf->xyzSize() );

   CellInterval   packingInterval = coarseToFinePackInterval( dir, sf->xyzSize(), fineReceiver->getId() );
   CellInterval unpackingInterval = coarseToFineUnpackInterval( stencil::inverseDir[dir], rf->xyzSize(), fineReceiver->getId() );

   WALBERLA_ASSERT_EQUAL( packingInterval.xSize() * uint_t(2), unpackingInterval.xSize() );
   WALBERLA_ASSERT_EQUAL( packingInterval.ySize() * uint_t(2), unpackingInterval.ySize() );
   if( Stencil::D == uint_t(3) )
   {
      WALBERLA_ASSERT_EQUAL( packingInterval.zSize() * uint_t(2), unpackingInterval.zSize() );
   }
   else
   {
      WALBERLA_ASSERT_EQUAL( packingInterval.zSize(), unpackingInterval.zSize() );
   }

   if( optimizedForLinearExplosion_ )
   {
      cell_idx_t rz = unpackingInterval.zMin();
      for( cell_idx_t sz = packingInterval.zMin(); sz <= packingInterval.zMax(); ++sz )
      {
         cell_idx_t ry = unpackingInterval.yMin();
         for( cell_idx_t sy = packingInterval.yMin(); sy <= packingInterval.yMax(); ++sy )
         {
            cell_idx_t rx = unpackingInterval.xMin();
            for( cell_idx_t sx = packingInterval.xMin(); sx <= packingInterval.xMax(); ++sx ) {
               for( uint_t idx = 0; idx < Stencil::Size; ++idx )
               {
                  rf->get( rx, ry, rz, idx ) = sf->get( sx, sy, sz, idx );

#ifndef NDEBUG
               if( boundaryHandling->isDomain(sx,sy,sz) )
                  WALBERLA_ASSERT( !math::isnan( sf->get( sx, sy, sz, idx ) ), sx << ", " << sy << ", " << sz << ", " << idx << " coarse sender block = " << coarseSender->getId() << " in dir <" << stencil::cx[dir] << "," << stencil::cy[dir] << "," << stencil::cz[dir] << "> " );
#endif
               }
               rx += cell_idx_t(2);
            }
            ry += cell_idx_t(2);
            WALBERLA_ASSERT_EQUAL( rx, unpackingInterval.xMax() + cell_idx_t(1) );
         }
         rz += cell_idx_t(2);
         WALBERLA_ASSERT_EQUAL( ry, unpackingInterval.yMax() + cell_idx_t(1) );
      }
      if( Stencil::D == uint_t(3) )
      {
         WALBERLA_ASSERT_EQUAL( rz, unpackingInterval.zMax() + cell_idx_t(1) );
      }
      else
      {
         WALBERLA_ASSERT_EQUAL( rz, unpackingInterval.zMax() + cell_idx_t(2) );
      }
   }
   else
   {
      cell_idx_t rz = unpackingInterval.zMin();
      for( cell_idx_t sz = packingInterval.zMin(); sz <= packingInterval.zMax(); ++sz )
      {
         cell_idx_t ry = unpackingInterval.yMin();
         for( cell_idx_t sy = packingInterval.yMin(); sy <= packingInterval.yMax(); ++sy )
         {
            cell_idx_t rx = unpackingInterval.xMin();
            for( cell_idx_t sx = packingInterval.xMin(); sx <= packingInterval.xMax(); ++sx ) {
               for( uint_t idx = 0; idx < Stencil::Size; ++idx )
             //for( uint_t d = 0; d < Stencil::d_per_d_length[dir]; ++d )
               {
                //const auto idx = Stencil::idx[ Stencil::d_per_d[dir][d] ];

                  const auto & value = sf->get(sx,sy,sz,idx);

#ifndef NDEBUG
               if( boundaryHandling->isDomain(sx,sy,sz) )
                  WALBERLA_ASSERT( !math::isnan( value ), "value at " << sx << ", " << sy << ", " << sz << ", " << idx << " coarse sender block = " << coarseSender->getId() << " in dir <" << stencil::cx[dir] << "," << stencil::cy[dir] << "," << stencil::cz[dir] << "> " );
#endif

                  rf->get( rx,                 ry,                 rz,                 idx ) = value;
                  rf->get( rx + cell_idx_t(1), ry,                 rz,                 idx ) = value;
                  rf->get( rx,                 ry + cell_idx_t(1), rz,                 idx ) = value;
                  rf->get( rx + cell_idx_t(1), ry + cell_idx_t(1), rz,                 idx ) = value;
                  if( Stencil::D == uint_t(3) )
                  {
                     rf->get( rx,                 ry,                 rz + cell_idx_t(1), idx ) = value;
                     rf->get( rx + cell_idx_t(1), ry,                 rz + cell_idx_t(1), idx ) = value;
                     rf->get( rx,                 ry + cell_idx_t(1), rz + cell_idx_t(1), idx ) = value;
                     rf->get( rx + cell_idx_t(1), ry + cell_idx_t(1), rz + cell_idx_t(1), idx ) = value;
                  }
               }
               rx += cell_idx_t(2);
            }
            ry += cell_idx_t(2);
            WALBERLA_ASSERT_EQUAL( rx, unpackingInterval.xMax() + cell_idx_t(1) );
         }
         rz += cell_idx_t(2);
         WALBERLA_ASSERT_EQUAL( ry, unpackingInterval.yMax() + cell_idx_t(1) );
      }
      if( Stencil::D == uint_t(3) )
      {
         WALBERLA_ASSERT_EQUAL( rz, unpackingInterval.zMax() + cell_idx_t(1) );
      }
      else
      {
         WALBERLA_ASSERT_EQUAL( rz, unpackingInterval.zMax() + cell_idx_t(2) );
      }
   }
}



////////////////////
// Fine to coarse //
////////////////////

#ifdef NDEBUG
template< typename LatticeModel_T >
void PdfFieldPackInfo< LatticeModel_T >::packDataFineToCoarseImpl( const Block * fineSender, const BlockID & coarseReceiver,
                                                                   stencil::Direction dir, mpi::SendBuffer & buffer ) const
#else
template< typename LatticeModel_T, typename BoundaryHandling_T >
void PdfFieldPackInfo< LatticeModel_T, BoundaryHandling_T >::packDataFineToCoarseImpl( const Block * fineSender, const BlockID & coarseReceiver,
                                                                                       stencil::Direction dir, mpi::SendBuffer & buffer ) const
#endif
{
#ifndef NDEBUG
   if( Stencil::D == uint_t(2) )
      WALBERLA_ASSERT_EQUAL( stencil::cz[dir], 0 );
#endif

   if( Stencil::d_per_d_length[dir] == uint_t(0) )
      return;

   if( ( ( isEdgeDirection(dir) || isCornerDirection(dir) ) && blocksConnectedByFaces( fineSender, coarseReceiver ) ) ||
       ( isCornerDirection(dir) && blocksConnectedByEdges( fineSender, coarseReceiver ) ) )
      return;

   const PdfField_T * field = fineSender->getData< PdfField_T >( pdfFieldId_ );
   
#ifndef NDEBUG
   const BoundaryHandling_T * boundaryHandling = fineSender->template getData< BoundaryHandling_T >( boundaryHandlingId_ );
#endif   

   CellInterval packingInterval = fineToCoarsePackInterval( dir, field->xyzSize() );

   const real_t factor = ( Stencil::D == uint_t(3) ) ? real_t( 0.125 ) : real_t( 0.25 );

   for( cell_idx_t z = packingInterval.zMin(); z <= packingInterval.zMax(); z += cell_idx_t(2) ) {
      for( cell_idx_t y = packingInterval.yMin(); y <= packingInterval.yMax(); y += cell_idx_t(2) ) {
         for( cell_idx_t x = packingInterval.xMin(); x <= packingInterval.xMax(); x += cell_idx_t(2) ) {
            for( uint_t d = 0; d < Stencil::d_per_d_length[dir]; ++d )
            {
               const auto idx = Stencil::idx[ Stencil::d_per_d[dir][d] ];

#ifndef NDEBUG
               if( boundaryHandling->isDomain(x,y,z) )
               {
                  WALBERLA_ASSERT( !math::isnan( field->get( x,                 y,                 z,                 idx ) ), x << ", " << y << ", " << z << ", " << idx << " fine sender block = " << fineSender->getId() << " in dir <" << stencil::cx[dir] << "," << stencil::cy[dir] << "," << stencil::cz[dir] << "> ");
                  WALBERLA_ASSERT( !math::isnan( field->get( x + cell_idx_t(1), y,                 z,                 idx ) ), x + cell_idx_t(1) << ", " << y << ", " << z << ", " << idx <<" fine sender block = " << fineSender->getId() << " in dir <" << stencil::cx[dir] << "," << stencil::cy[dir] << "," << stencil::cz[dir] << "> " );
                  WALBERLA_ASSERT( !math::isnan( field->get( x,                 y + cell_idx_t(1), z,                 idx ) ), x << ", " << y + cell_idx_t(1) << ", " << z << ", " << idx << " fine sender block = " << fineSender->getId() << " in dir <" << stencil::cx[dir] << "," << stencil::cy[dir] << "," << stencil::cz[dir] << "> " );
                  WALBERLA_ASSERT( !math::isnan( field->get( x + cell_idx_t(1), y + cell_idx_t(1), z,                 idx ) ), x + cell_idx_t(1) << ", " << y + cell_idx_t(1) << ", " << z << ", " << idx << " fine sender block = " << fineSender->getId() << " in dir <" << stencil::cx[dir] << "," << stencil::cy[dir] << "," << stencil::cz[dir] << "> " );
                  if( Stencil::D == uint_t(3) )
                  {
                     WALBERLA_ASSERT( !math::isnan( field->get( x,                 y,                 z + cell_idx_t(1), idx ) ), x << ", " << y << ", " << z + cell_idx_t(1) << ", " << idx << " fine sender block = " << fineSender->getId() << " in dir <" << stencil::cx[dir] << "," << stencil::cy[dir] << "," << stencil::cz[dir] << "> "  );
                     WALBERLA_ASSERT( !math::isnan( field->get( x + cell_idx_t(1), y,                 z + cell_idx_t(1), idx ) ), x + cell_idx_t(1) << ", " << y << ", " << z + cell_idx_t(1) << ", " << idx << " fine sender block = " << fineSender->getId() << " in dir <" << stencil::cx[dir] << "," << stencil::cy[dir] << "," << stencil::cz[dir] << "> " );
                     WALBERLA_ASSERT( !math::isnan( field->get( x,                 y + cell_idx_t(1), z + cell_idx_t(1), idx ) ), x << ", " << y + cell_idx_t(1) << ", " << z + cell_idx_t(1) << ", " << idx << " fine sender block = " << fineSender->getId() << " in dir <" << stencil::cx[dir] << "," << stencil::cy[dir] << "," << stencil::cz[dir] << "> " );
                     WALBERLA_ASSERT( !math::isnan( field->get( x + cell_idx_t(1), y + cell_idx_t(1), z + cell_idx_t(1), idx ) ), x + cell_idx_t(1) << ", " << y + cell_idx_t(1) << ", " << z + cell_idx_t(1) << ", " << idx << " fine sender block = " << fineSender->getId() << " in dir <" << stencil::cx[dir] << "," << stencil::cy[dir] << "," << stencil::cz[dir] << "> " );
                  }
               }
#endif

               typename PdfField_T::value_type value  = field->get( x,                 y,                 z,                 idx );
                                               value += field->get( x + cell_idx_t(1), y,                 z,                 idx );
                                               value += field->get( x,                 y + cell_idx_t(1), z,                 idx );
                                               value += field->get( x + cell_idx_t(1), y + cell_idx_t(1), z,                 idx );
               if( Stencil::D == uint_t(3) )
               {
                                               value += field->get( x,                 y,                 z + cell_idx_t(1), idx );
                                               value += field->get( x + cell_idx_t(1), y,                 z + cell_idx_t(1), idx );
                                               value += field->get( x,                 y + cell_idx_t(1), z + cell_idx_t(1), idx );
                                               value += field->get( x + cell_idx_t(1), y + cell_idx_t(1), z + cell_idx_t(1), idx );
               }

               buffer << ( factor * value );
            }
         }
      }
   }
}



#ifdef NDEBUG
template< typename LatticeModel_T >
void PdfFieldPackInfo< LatticeModel_T >::unpackDataFineToCoarse( Block * coarseReceiver, const BlockID & fineSender,
                                                                 stencil::Direction dir, mpi::RecvBuffer & buffer )
#else
template< typename LatticeModel_T, typename BoundaryHandling_T >
void PdfFieldPackInfo< LatticeModel_T, BoundaryHandling_T >::unpackDataFineToCoarse( Block * coarseReceiver, const BlockID & fineSender,
                                                                                     stencil::Direction dir, mpi::RecvBuffer & buffer )
#endif
{
#ifndef NDEBUG
   if( Stencil::D == uint_t(2) )
      WALBERLA_ASSERT_EQUAL( stencil::cz[dir], 0 );
#endif

   auto invDir = stencil::inverseDir[dir];
   if( Stencil::d_per_d_length[invDir] == uint_t(0) )
      return;

   if( ( ( isEdgeDirection(dir) || isCornerDirection(dir) ) && blocksConnectedByFaces( coarseReceiver, fineSender ) ) ||
       ( isCornerDirection(dir) && blocksConnectedByEdges( coarseReceiver, fineSender ) ) )
      return;

   PdfField_T * field = coarseReceiver->getData< PdfField_T >( pdfFieldId_ );

#ifndef NDEBUG
   const BoundaryHandling_T * boundaryHandling = coarseReceiver->template getData< BoundaryHandling_T >( boundaryHandlingId_ );
#endif

   CellInterval unpackingInterval = fineToCoarseUnpackInterval( dir, field->xyzSize(), fineSender );

   for( cell_idx_t z = unpackingInterval.zMin(); z <= unpackingInterval.zMax(); ++z ) {
      for( cell_idx_t y = unpackingInterval.yMin(); y <= unpackingInterval.yMax(); ++y ) {
         for( cell_idx_t x = unpackingInterval.xMin(); x <= unpackingInterval.xMax(); ++x ) {
            for( uint_t d = 0; d < Stencil::d_per_d_length[invDir]; ++d )
            {
               buffer >> field->get( x, y, z, Stencil::idx[ Stencil::d_per_d[invDir][d] ] );

#ifndef NDEBUG
               if( boundaryHandling->isDomain(x,y,z) )
                  WALBERLA_ASSERT( !math::isnan( field->get( x, y, z, Stencil::idx[ Stencil::d_per_d[invDir][d] ] ) ) );
#endif
            }
         }
      }
   }
}



#ifdef NDEBUG
template< typename LatticeModel_T >
void PdfFieldPackInfo< LatticeModel_T >::communicateLocalFineToCoarse( const Block * fineSender, Block * coarseReceiver, stencil::Direction dir )
#else
template< typename LatticeModel_T, typename BoundaryHandling_T >
void PdfFieldPackInfo< LatticeModel_T, BoundaryHandling_T >::communicateLocalFineToCoarse( const Block * fineSender, Block * coarseReceiver,
                                                                                           stencil::Direction dir )
#endif
{
#ifndef NDEBUG
   if( Stencil::D == uint_t(2) )
      WALBERLA_ASSERT_EQUAL( stencil::cz[dir], 0 );
#endif

   if( Stencil::d_per_d_length[dir] == uint_t(0) )
      return;

   if( ( ( isEdgeDirection(dir) || isCornerDirection(dir) ) && blocksConnectedByFaces( fineSender, coarseReceiver->getId() ) ) ||
       ( isCornerDirection(dir) && blocksConnectedByEdges( fineSender, coarseReceiver->getId() ) ) )
      return;

   const PdfField_T * sf =     fineSender->getData< PdfField_T >( pdfFieldId_ );
         PdfField_T * rf = coarseReceiver->getData< PdfField_T >( pdfFieldId_ );
         
#ifndef NDEBUG
   const BoundaryHandling_T * boundaryHandling = fineSender->template getData< BoundaryHandling_T >( boundaryHandlingId_ );
#endif

   WALBERLA_ASSERT_EQUAL( sf->xyzSize(), rf->xyzSize() );

   CellInterval   packingInterval = fineToCoarsePackInterval( dir, sf->xyzSize() );
   CellInterval unpackingInterval = fineToCoarseUnpackInterval( stencil::inverseDir[dir], rf->xyzSize(), fineSender->getId() );

   WALBERLA_ASSERT_EQUAL( packingInterval.xSize(), unpackingInterval.xSize() * uint_t(2) );
   WALBERLA_ASSERT_EQUAL( packingInterval.ySize(), unpackingInterval.ySize() * uint_t(2) );
   if( Stencil::D == uint_t(3) )
   {
      WALBERLA_ASSERT_EQUAL( packingInterval.zSize(), unpackingInterval.zSize() * uint_t(2) );
   }
   else
   {
      WALBERLA_ASSERT_EQUAL( packingInterval.zSize(), unpackingInterval.zSize() );
   }

   const real_t factor = ( Stencil::D == uint_t(3) ) ? real_t( 0.125 ) : real_t( 0.25 );

   cell_idx_t sz = packingInterval.zMin();
   for( cell_idx_t rz = unpackingInterval.zMin(); rz <= unpackingInterval.zMax(); ++rz )
   {
      cell_idx_t sy = packingInterval.yMin();
      for( cell_idx_t ry = unpackingInterval.yMin(); ry <= unpackingInterval.yMax(); ++ry )
      {
         cell_idx_t sx = packingInterval.xMin();
         for( cell_idx_t rx = unpackingInterval.xMin(); rx <= unpackingInterval.xMax(); ++rx ) {
            for( uint_t d = 0; d < Stencil::d_per_d_length[dir]; ++d )
            {
               const auto idx = Stencil::idx[ Stencil::d_per_d[dir][d] ];

#ifndef NDEBUG
               if( boundaryHandling->isDomain(sx,sy,sz) )
               {
                  WALBERLA_ASSERT( !math::isnan( sf->get( sx,                 sy,                 sz,                 idx ) ), sx << ", " << sy << ", " << sz << ", " << idx << " fine sender block = " << fineSender->getId() << " in dir <" << stencil::cx[dir] << "," << stencil::cy[dir] << "," << stencil::cz[dir] << "> ");
                  WALBERLA_ASSERT( !math::isnan( sf->get( sx + cell_idx_t(1), sy,                 sz,                 idx ) ), sx + cell_idx_t(1) << ", " << sy << ", " << sz << ", " << idx << " fine sender block = " << fineSender->getId() << " in dir <" << stencil::cx[dir] << "," << stencil::cy[dir] << "," << stencil::cz[dir] << "> ");
                  WALBERLA_ASSERT( !math::isnan( sf->get( sx,                 sy + cell_idx_t(1), sz,                 idx ) ), sx << ", " << sy + cell_idx_t(1) << ", " << sz << ", " << idx << " fine sender block = " << fineSender->getId() << " in dir <" << stencil::cx[dir] << "," << stencil::cy[dir] << "," << stencil::cz[dir] << "> ");
                  WALBERLA_ASSERT( !math::isnan( sf->get( sx + cell_idx_t(1), sy + cell_idx_t(1), sz,                 idx ) ), sx + cell_idx_t(1) << ", " << sy + cell_idx_t(1) << ", " << sz << ", " << idx << " fine sender block = " << fineSender->getId() << " in dir <" << stencil::cx[dir] << "," << stencil::cy[dir] << "," << stencil::cz[dir] << "> ");
                  if( Stencil::D == uint_t(3) )
                  {
                     WALBERLA_ASSERT( !math::isnan( sf->get( sx,                 sy,                 sz + cell_idx_t(1), idx ) ), sx << ", " << sy << ", " << sz + cell_idx_t(1) << ", " << idx << " fine sender block = " << fineSender->getId() << " in dir <" << stencil::cx[dir] << "," << stencil::cy[dir] << "," << stencil::cz[dir] << "> ");
                     WALBERLA_ASSERT( !math::isnan( sf->get( sx + cell_idx_t(1), sy,                 sz + cell_idx_t(1), idx ) ), sx + cell_idx_t(1) << ", " << sy << ", " << sz + cell_idx_t(1) << ", " << idx << " fine sender block = " << fineSender->getId() << " in dir <" << stencil::cx[dir] << "," << stencil::cy[dir] << "," << stencil::cz[dir] << "> ");
                     WALBERLA_ASSERT( !math::isnan( sf->get( sx,                 sy + cell_idx_t(1), sz + cell_idx_t(1), idx ) ), sx << ", " << sy + cell_idx_t(1) << ", " << sz + cell_idx_t(1) << ", " << idx << " fine sender block = " << fineSender->getId() << " in dir <" << stencil::cx[dir] << "," << stencil::cy[dir] << "," << stencil::cz[dir] << "> ");
                     WALBERLA_ASSERT( !math::isnan( sf->get( sx + cell_idx_t(1), sy + cell_idx_t(1), sz + cell_idx_t(1), idx ) ), sx + cell_idx_t(1) << ", " << sy + cell_idx_t(1) << ", " << sz + cell_idx_t(1) << ", " << idx << " fine sender block = " << fineSender->getId() << " in dir <" << stencil::cx[dir] << "," << stencil::cy[dir] << "," << stencil::cz[dir] << "> ");
                  }
               }
#endif

               typename PdfField_T::value_type value  = sf->get( sx,                 sy,                 sz,                 idx );
                                               value += sf->get( sx + cell_idx_t(1), sy,                 sz,                 idx );
                                               value += sf->get( sx,                 sy + cell_idx_t(1), sz,                 idx );
                                               value += sf->get( sx + cell_idx_t(1), sy + cell_idx_t(1), sz,                 idx );
               if( Stencil::D == uint_t(3) )
               {
                                               value += sf->get( sx,                 sy,                 sz + cell_idx_t(1), idx );
                                               value += sf->get( sx + cell_idx_t(1), sy,                 sz + cell_idx_t(1), idx );
                                               value += sf->get( sx,                 sy + cell_idx_t(1), sz + cell_idx_t(1), idx );
                                               value += sf->get( sx + cell_idx_t(1), sy + cell_idx_t(1), sz + cell_idx_t(1), idx );
               }

               rf->get( rx, ry, rz, idx ) = factor * value;
            }
            sx += cell_idx_t(2);
         }
         sy += cell_idx_t(2);
         WALBERLA_ASSERT_EQUAL( sx, packingInterval.xMax() + cell_idx_t(1) );
      }
      sz += cell_idx_t(2);
      WALBERLA_ASSERT_EQUAL( sy, packingInterval.yMax() + cell_idx_t(1) );
   }
   if( Stencil::D == uint_t(3) )
   {
      WALBERLA_ASSERT_EQUAL( sz, packingInterval.zMax() + cell_idx_t(1) );
   }
   else
   {
      WALBERLA_ASSERT_EQUAL( sz, packingInterval.zMax() + cell_idx_t(2) );
   }
}



///////////////////////////////////////////////////////////////////////
// Helper functions for determining packing/unpacking cell intervals //
///////////////////////////////////////////////////////////////////////

#ifdef NDEBUG
template< typename LatticeModel_T >
inline CellInterval PdfFieldPackInfo< LatticeModel_T >::equalLevelPackInterval( stencil::Direction dir, const CellInterval & cellBB,
                                                                                const uint_t numberOfLayers )
#else
template< typename LatticeModel_T, typename BoundaryHandling_T >
inline CellInterval PdfFieldPackInfo< LatticeModel_T, BoundaryHandling_T >::equalLevelPackInterval( stencil::Direction dir, const CellInterval & cellBB,
                                                                                                    const uint_t numberOfLayers )
#endif
{
   CellInterval interval = equalLevelUnpackInterval( dir, cellBB, numberOfLayers );

   for( uint_t i = 0; i != Stencil::D; ++i )
   {
      const auto offset = cell_idx_c( stencil::c[i][dir] ) * cell_idx_c( numberOfLayers );
      interval.min()[i] -= offset;
      interval.max()[i] -= offset;
   }

   return interval;
}



#ifdef NDEBUG
template< typename LatticeModel_T >
CellInterval PdfFieldPackInfo< LatticeModel_T >::equalLevelUnpackInterval( stencil::Direction dir, const CellInterval & cellBB,
                                                                           const uint_t numberOfLayers )
#else
template< typename LatticeModel_T, typename BoundaryHandling_T >
CellInterval PdfFieldPackInfo< LatticeModel_T, BoundaryHandling_T >::equalLevelUnpackInterval( stencil::Direction dir, const CellInterval & cellBB,
                                                                                               const uint_t numberOfLayers )
#endif
{
   CellInterval interval( cellBB );
   interval.expand( cell_idx_c(numberOfLayers) );

   for( uint_t i = 0; i != 3; ++i )
   {
      const auto c = cell_idx_c( stencil::c[i][dir] );

      if( c == -1 ) interval.max()[i] = interval.min()[i] + cell_idx_c( numberOfLayers - 1 );
      else if( c == 1 ) interval.min()[i] = interval.max()[i] - cell_idx_c( numberOfLayers - 1 );
      else
      {
         WALBERLA_ASSERT_EQUAL( c, cell_idx_t(0) );
         interval.min()[i] += cell_idx_c( numberOfLayers );
         interval.max()[i] -= cell_idx_c( numberOfLayers );
      }
   }

   return interval;
}



#ifdef NDEBUG
template< typename LatticeModel_T >
inline bool PdfFieldPackInfo< LatticeModel_T >::equalLevelFaceIntervalSplitable( const CellInterval & interval, stencil::Direction dir ) const
#else
template< typename LatticeModel_T, typename BoundaryHandling_T >
inline bool PdfFieldPackInfo< LatticeModel_T, BoundaryHandling_T >::equalLevelFaceIntervalSplitable( const CellInterval & interval, stencil::Direction dir ) const
#endif
{
   if( stencil::cx[dir] != 0 )
   {
      WALBERLA_ASSERT_EQUAL( stencil::cy[dir], 0 );
      WALBERLA_ASSERT_EQUAL( stencil::cz[dir], 0 );

      WALBERLA_ASSERT_EQUAL( interval.xSize(), equalLevelCells_ );

      return interval.ySize() > uint_t(2) && interval.zSize() > uint_t(2);
   }
   else if( stencil::cy[dir] != 0 )
   {
      WALBERLA_ASSERT_EQUAL( stencil::cx[dir], 0 );
      WALBERLA_ASSERT_EQUAL( stencil::cz[dir], 0 );

      WALBERLA_ASSERT_EQUAL( interval.ySize(), equalLevelCells_ );

      return interval.xSize() > uint_t(2) && interval.zSize() > uint_t(2);
   }

   WALBERLA_ASSERT_UNEQUAL( stencil::cz[dir], 0 );
   WALBERLA_ASSERT_EQUAL( stencil::cx[dir], 0 );
   WALBERLA_ASSERT_EQUAL( stencil::cy[dir], 0 );

   WALBERLA_ASSERT_EQUAL( interval.zSize(), equalLevelCells_ );

   return interval.xSize() > uint_t(2) && interval.ySize() > uint_t(2);
}



#ifdef NDEBUG
template< typename LatticeModel_T >
inline std::vector< CellInterval > PdfFieldPackInfo< LatticeModel_T >::splitEqualLevelFaceInterval( const CellInterval & interval,
                                                                                                    stencil::Direction dir ) const // dir: from sender to receiver
#else
template< typename LatticeModel_T, typename BoundaryHandling_T >
inline std::vector< CellInterval > PdfFieldPackInfo< LatticeModel_T, BoundaryHandling_T >::splitEqualLevelFaceInterval( const CellInterval & interval,
                                                                                                                        stencil::Direction dir ) const // dir: from sender to receiver
#endif
{
   std::vector< CellInterval > intervals;

   if( stencil::cx[dir] != 0 )
   {
      WALBERLA_ASSERT_EQUAL( stencil::cy[dir], 0 );
      WALBERLA_ASSERT_EQUAL( stencil::cz[dir], 0 );

      WALBERLA_ASSERT_EQUAL( interval.xSize(), equalLevelCells_ );

      intervals.emplace_back( interval.xMin(), interval.yMin(), interval.zMax(),
                                         interval.xMax(), interval.yMax(), interval.zMax() );
      intervals.emplace_back( interval.xMin(), interval.yMin(), interval.zMin() + cell_idx_t(1),
                                         interval.xMax(), interval.yMin(), interval.zMax() - cell_idx_t(1) );
      intervals.emplace_back( interval.xMin(), interval.yMax(), interval.zMin() + cell_idx_t(1),
                                         interval.xMax(), interval.yMax(), interval.zMax() - cell_idx_t(1) );
      intervals.emplace_back( interval.xMin(), interval.yMin(), interval.zMin(),
                                         interval.xMax(), interval.yMax(), interval.zMin() );

      const cell_idx_t x = ( stencil::cx[dir] > 0 ) ? interval.xMax() : interval.xMin();
      intervals.emplace_back( x, interval.yMin() + cell_idx_t(1), interval.zMin() + cell_idx_t(1),
                                         x, interval.yMax() - cell_idx_t(1), interval.zMax() - cell_idx_t(1) );
   }
   else if( stencil::cy[dir] != 0 )
   {
      WALBERLA_ASSERT_EQUAL( stencil::cx[dir], 0 );
      WALBERLA_ASSERT_EQUAL( stencil::cz[dir], 0 );

      WALBERLA_ASSERT_EQUAL( interval.ySize(), equalLevelCells_ );

      intervals.emplace_back( interval.xMin(), interval.yMin(), interval.zMax(),
                                         interval.xMax(), interval.yMax(), interval.zMax() );
      intervals.emplace_back( interval.xMin(), interval.yMin(), interval.zMin() + cell_idx_t(1),
                                         interval.xMin(), interval.yMax(), interval.zMax() - cell_idx_t(1) );
      intervals.emplace_back( interval.xMax(), interval.yMin(), interval.zMin() + cell_idx_t(1),
                                         interval.xMax(), interval.yMax(), interval.zMax() - cell_idx_t(1) );
      intervals.emplace_back( interval.xMin(), interval.yMin(), interval.zMin(),
                                         interval.xMax(), interval.yMax(), interval.zMin() );

      const cell_idx_t y = ( stencil::cy[dir] > 0 ) ? interval.yMax() : interval.yMin();
      intervals.emplace_back( interval.xMin() + cell_idx_t(1), y, interval.zMin() + cell_idx_t(1),
                                         interval.xMax() - cell_idx_t(1), y, interval.zMax() - cell_idx_t(1) );
   }
   else
   {
      WALBERLA_ASSERT_UNEQUAL( stencil::cz[dir], 0 );
      WALBERLA_ASSERT_EQUAL( stencil::cx[dir], 0 );
      WALBERLA_ASSERT_EQUAL( stencil::cy[dir], 0 );

      WALBERLA_ASSERT_EQUAL( interval.zSize(), equalLevelCells_ );

      intervals.emplace_back( interval.xMin(), interval.yMax(),                 interval.zMin(),
                                         interval.xMax(), interval.yMax(),                 interval.zMax() );
      intervals.emplace_back( interval.xMin(), interval.yMin() + cell_idx_t(1), interval.zMin(),
                                         interval.xMin(), interval.yMax() - cell_idx_t(1), interval.zMax() );
      intervals.emplace_back( interval.xMax(), interval.yMin() + cell_idx_t(1), interval.zMin(),
                                         interval.xMax(), interval.yMax() - cell_idx_t(1), interval.zMax() );
      intervals.emplace_back( interval.xMin(), interval.yMin(),                 interval.zMin(),
                                         interval.xMax(), interval.yMin(),                 interval.zMax() );

      const cell_idx_t z = ( stencil::cz[dir] > 0 ) ? interval.zMax() : interval.zMin();
      intervals.emplace_back( interval.xMin() + cell_idx_t(1), interval.yMin() + cell_idx_t(1), z,
                                         interval.xMax() - cell_idx_t(1), interval.yMax() - cell_idx_t(1), z );
   }

   return intervals;
}



#ifdef NDEBUG
template< typename LatticeModel_T >
inline Vector3< cell_idx_t > PdfFieldPackInfo< LatticeModel_T >::getNeighborShift( const BlockID & smallBlock, stencil::Direction dir ) // dir: direction from big to small block
#else
template< typename LatticeModel_T, typename BoundaryHandling_T >
inline Vector3< cell_idx_t > PdfFieldPackInfo< LatticeModel_T, BoundaryHandling_T >::getNeighborShift( const BlockID & smallBlock,
                                                                                                       stencil::Direction dir ) // dir: direction from big to small block
#endif
{
   Vector3< cell_idx_t > shift;

   uint_t branchId = smallBlock.getBranchId();

   shift[0] = ( stencil::cx[dir] == 0 ) ? ( ( ( branchId & uint_t(1) ) == uint_t(0) ) ? cell_idx_t(-1) : cell_idx_t(1) ) : cell_idx_t(0);
   shift[1] = ( stencil::cy[dir] == 0 ) ? ( ( ( branchId & uint_t(2) ) == uint_t(0) ) ? cell_idx_t(-1) : cell_idx_t(1) ) : cell_idx_t(0);
   shift[2] = ( Stencil::D == uint_t(3) ) ?
            ( ( stencil::cz[dir] == 0 ) ? ( ( ( branchId & uint_t(4) ) == uint_t(0) ) ? cell_idx_t(-1) : cell_idx_t(1) ) : cell_idx_t(0) ) : cell_idx_t(0);

   return shift;
}



#ifdef NDEBUG
template< typename LatticeModel_T >
CellInterval PdfFieldPackInfo< LatticeModel_T >::coarseToFinePackInterval( stencil::Direction dir, const CellInterval & cellBB,
                                                                           const BlockID & smallBlock )
#else
template< typename LatticeModel_T, typename BoundaryHandling_T >
CellInterval PdfFieldPackInfo< LatticeModel_T, BoundaryHandling_T >::coarseToFinePackInterval( stencil::Direction dir, const CellInterval & cellBB,
                                                                                               const BlockID & smallBlock )
#endif
{
   CellInterval interval = equalLevelPackInterval( dir, cellBB, uint_t(2) );
 //CellInterval interval = equalLevelPackInterval( dir, cellBB, uint_t(1) );

   Vector3< cell_idx_t > shift = getNeighborShift( smallBlock, dir );

   WALBERLA_ASSERT( divisibleByTwo( cellBB ) );

   for( uint_t i = 0; i != Stencil::D; ++i )
   {
      if( shift[i] == cell_idx_t(-1) )
         interval.max()[i] = interval.min()[i] + cell_idx_c( cellBB.size(i) / uint_t(2) ) + cell_idx_t(1);
       //interval.max()[i] = interval.min()[i] + cell_idx_c( cellBB.size(i) / uint_t(2) );
      if( shift[i] == cell_idx_t( 1) )
         interval.min()[i] = interval.max()[i] - cell_idx_c( cellBB.size(i) / uint_t(2) ) - cell_idx_t(1);
       //interval.min()[i] = interval.max()[i] - cell_idx_c( cellBB.size(i) / uint_t(2) );
   }

   WALBERLA_ASSERT( cellBB.contains( interval ) );

   return interval;
}



#ifdef NDEBUG
template< typename LatticeModel_T >
CellInterval PdfFieldPackInfo< LatticeModel_T >::coarseToFineUnpackInterval( stencil::Direction dir, const CellInterval & cellBB,
                                                                             const BlockID & smallBlock )
#else
template< typename LatticeModel_T, typename BoundaryHandling_T >
CellInterval PdfFieldPackInfo< LatticeModel_T, BoundaryHandling_T >::coarseToFineUnpackInterval( stencil::Direction dir, const CellInterval & cellBB,
                                                                                                 const BlockID & smallBlock )
#endif
{
   CellInterval interval = equalLevelUnpackInterval( dir, cellBB, uint_t(4) );
 //CellInterval interval = equalLevelUnpackInterval( dir, cellBB, uint_t(2) );

   Vector3< cell_idx_t > shift = getNeighborShift( smallBlock, dir );

   for( uint_t i = 0; i != Stencil::D; ++i )
   {
      if( shift[i] == cell_idx_t(-1) )
         interval.max()[i] += cell_idx_t(4);
       //interval.max()[i] += cell_idx_t(2);
      if( shift[i] == cell_idx_t( 1) )
         interval.min()[i] -= cell_idx_t(4);
       //interval.min()[i] -= cell_idx_t(2);
   }

#ifndef NDEBUG
   CellInterval expandedCellBB( cellBB );
   expandedCellBB.expand( cell_idx_t(4) );
   WALBERLA_ASSERT( expandedCellBB.contains( interval ) );
#endif

   return interval;
}



#ifdef NDEBUG
template< typename LatticeModel_T >
inline CellInterval PdfFieldPackInfo< LatticeModel_T >::fineToCoarsePackInterval( stencil::Direction dir, const CellInterval & cellBB )
#else
template< typename LatticeModel_T, typename BoundaryHandling_T >
inline CellInterval PdfFieldPackInfo< LatticeModel_T, BoundaryHandling_T >::fineToCoarsePackInterval( stencil::Direction dir, const CellInterval & cellBB )
#endif
{
   return equalLevelUnpackInterval( dir, cellBB, uint_t(2) );
}



#ifdef NDEBUG
template< typename LatticeModel_T >
CellInterval PdfFieldPackInfo< LatticeModel_T >::fineToCoarseUnpackInterval( stencil::Direction dir, const CellInterval & cellBB,
                                                                             const BlockID & smallBlock )
#else
template< typename LatticeModel_T, typename BoundaryHandling_T >
CellInterval PdfFieldPackInfo< LatticeModel_T, BoundaryHandling_T >::fineToCoarseUnpackInterval( stencil::Direction dir, const CellInterval & cellBB,
                                                                                                 const BlockID & smallBlock )
#endif
{
   CellInterval interval = equalLevelPackInterval( dir, cellBB, uint_t(1) );
   Vector3< cell_idx_t > shift = getNeighborShift( smallBlock, dir );

   WALBERLA_ASSERT( divisibleByTwo( cellBB ) );

   for( uint_t i = 0; i != Stencil::D; ++i )
   {
      if( shift[i] == cell_idx_t(-1) )
         interval.max()[i] = interval.min()[i] + cell_idx_c( cellBB.size(i) / uint_t(2) ) - cell_idx_t(1);
      if( shift[i] == cell_idx_t( 1) )
         interval.min()[i] = interval.max()[i] - cell_idx_c( cellBB.size(i) / uint_t(2) ) + cell_idx_t(1);
   }

   WALBERLA_ASSERT( cellBB.contains( interval ) );

   return interval;
}



//////////////////////////////
// General helper functions //
//////////////////////////////

#ifdef NDEBUG
template< typename LatticeModel_T >
inline uint_t PdfFieldPackInfo< LatticeModel_T >::equalLevelCells()
#else
template< typename LatticeModel_T, typename BoundaryHandling_T >
inline uint_t PdfFieldPackInfo< LatticeModel_T, BoundaryHandling_T >::equalLevelCells()
#endif
{
   if( Stencil::containsDir( stencil::TNE ) || Stencil::containsDir( stencil::TNW ) || Stencil::containsDir( stencil::TSE ) ||
       Stencil::containsDir( stencil::TSW ) || Stencil::containsDir( stencil::BNE ) || Stencil::containsDir( stencil::BNW ) ||
       Stencil::containsDir( stencil::BSE ) || Stencil::containsDir( stencil::BSW ) )
      return uint_t(3);
      
   return uint_t(2);
}



#ifdef NDEBUG
template< typename LatticeModel_T >
inline bool PdfFieldPackInfo< LatticeModel_T >::isFaceDirection( stencil::Direction dir )
#else
template< typename LatticeModel_T, typename BoundaryHandling_T >
inline bool PdfFieldPackInfo< LatticeModel_T, BoundaryHandling_T >::isFaceDirection( stencil::Direction dir )
#endif
{
   return ( dir == stencil::N ) || ( dir == stencil::S ) || ( dir == stencil::W ) ||
          ( dir == stencil::E ) || ( dir == stencil::T ) || ( dir == stencil::B );
}



#ifdef NDEBUG
template< typename LatticeModel_T >
inline bool PdfFieldPackInfo< LatticeModel_T >::isEdgeDirection( stencil::Direction dir )
#else
template< typename LatticeModel_T, typename BoundaryHandling_T >
inline bool PdfFieldPackInfo< LatticeModel_T, BoundaryHandling_T >::isEdgeDirection( stencil::Direction dir )
#endif
{
   return ( dir == stencil::NW ) || ( dir == stencil::NE ) || ( dir == stencil::SW ) || ( dir == stencil::SE ) ||
          ( dir == stencil::TN ) || ( dir == stencil::TS ) || ( dir == stencil::TW ) || ( dir == stencil::TE ) ||
          ( dir == stencil::BN ) || ( dir == stencil::BS ) || ( dir == stencil::BW ) || ( dir == stencil::BE );
}



#ifdef NDEBUG
template< typename LatticeModel_T >
inline bool PdfFieldPackInfo< LatticeModel_T >::isCornerDirection( stencil::Direction dir )
#else
template< typename LatticeModel_T, typename BoundaryHandling_T >
inline bool PdfFieldPackInfo< LatticeModel_T, BoundaryHandling_T >::isCornerDirection( stencil::Direction dir )
#endif
{
   return ( dir == stencil::TNE ) || ( dir == stencil::TNW ) || ( dir == stencil::TSE ) || ( dir == stencil::TSW ) ||
          ( dir == stencil::BNE ) || ( dir == stencil::BNW ) || ( dir == stencil::BSE ) || ( dir == stencil::BSW );
}



#ifdef NDEBUG
template< typename LatticeModel_T >
inline bool PdfFieldPackInfo< LatticeModel_T >::blocksConnectedByFaces( const Block * block, const BlockID & neighbor )
#else
template< typename LatticeModel_T, typename BoundaryHandling_T >
inline bool PdfFieldPackInfo< LatticeModel_T, BoundaryHandling_T >::blocksConnectedByFaces( const Block * block, const BlockID & neighbor )
#endif
{
   const uint_t face[] = { uint_t(4), uint_t(10), uint_t(12), uint_t(13), uint_t(15), uint_t(21) };
   for( int i = 0; i != 6; ++i )
   {
      for( uint_t n = 0; n != block->getNeighborhoodSectionSize( face[i] ); ++n )
         if( block->getNeighborId( face[i], n ) == neighbor )
            return true;
   }
   return false;
}



#ifdef NDEBUG
template< typename LatticeModel_T >
inline bool PdfFieldPackInfo< LatticeModel_T >::blocksConnectedByEdges( const Block * block, const BlockID & neighbor )
#else
template< typename LatticeModel_T, typename BoundaryHandling_T >
inline bool PdfFieldPackInfo< LatticeModel_T, BoundaryHandling_T >::blocksConnectedByEdges( const Block * block, const BlockID & neighbor )
#endif
{
   const uint_t face[] = { uint_t( 1), uint_t( 3), uint_t( 5), uint_t( 7), uint_t( 9), uint_t(11),
                           uint_t(14), uint_t(16), uint_t(18), uint_t(20), uint_t(22), uint_t(24) };

   for( int i = 0; i != 12; ++i )
   {
      for( uint_t n = 0; n != block->getNeighborhoodSectionSize( face[i] ); ++n )
         if( block->getNeighborId( face[i], n ) == neighbor )
            return true;
   }
   return false;
}



#ifdef NDEBUG
template< typename LatticeModel_T >
inline bool PdfFieldPackInfo< LatticeModel_T >::divisibleByTwo( const CellInterval & cellBB )
#else
template< typename LatticeModel_T, typename BoundaryHandling_T >
inline bool PdfFieldPackInfo< LatticeModel_T, BoundaryHandling_T >::divisibleByTwo( const CellInterval & cellBB )
#endif
{
   return ( ( cellBB.xSize() & uint_t(1) ) == uint_t(0) ) &&
          ( ( cellBB.ySize() & uint_t(1) ) == uint_t(0) ) &&
          ( ( Stencil::D == uint_t(2) && cellBB.zSize() == uint_t(1) ) || ( Stencil::D == uint_t(3) && ( cellBB.zSize() & uint_t(1) ) == uint_t(0) ) );
}



#ifdef NDEBUG
template< typename LatticeModel_T >
bool PdfFieldPackInfo< LatticeModel_T >::coarserNeighborExistsInVicinity( const Block * block, stencil::Direction dir )
#else
template< typename LatticeModel_T, typename BoundaryHandling_T >
bool PdfFieldPackInfo< LatticeModel_T, BoundaryHandling_T >::coarserNeighborExistsInVicinity( const Block * block, stencil::Direction dir )
#endif
{
   if( block->getLevel() == uint_t(0) )
      return false;

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
            if( block->neighborhoodSectionHasLargerBlock( blockforest::getBlockNeighborhoodSectionIndex(x,y,z) ) )
               return true;
         }
      }
   }

   return false;
}



} // namespace refinement
} // namespace lbm
} // namespace walberla
