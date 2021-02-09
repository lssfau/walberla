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
//! \file PackInfo.h
//! \ingroup lbm
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "blockforest/BlockNeighborhoodSection.h"
#include "blockforest/communication/NonUniformPackInfo.h"
#include "core/cell/CellInterval.h"
#include "stencil/Directions.h"


namespace walberla {
namespace field {
namespace refinement {



template< typename Field_T, typename Stencil >
class PackInfo : public blockforest::communication::NonUniformPackInfo
{
public:

   PackInfo( const BlockDataID & fieldId ) : fieldId_( fieldId ) {}
   ~PackInfo() override = default;

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

   static inline Vector3< cell_idx_t > getNeighborShift( const BlockID & smallBlock, stencil::Direction dir ); // dir: direction from big to small block

   static CellInterval coarseToFinePackInterval  ( stencil::Direction dir, const CellInterval & cellBB, const BlockID & smallBlock );
   static CellInterval coarseToFineUnpackInterval( stencil::Direction dir, const CellInterval & cellBB, const BlockID & smallBlock );

   static inline  CellInterval fineToCoarsePackInterval  ( stencil::Direction dir, const CellInterval & cellBB );
   static         CellInterval fineToCoarseUnpackInterval( stencil::Direction dir, const CellInterval & cellBB, const BlockID & smallBlock );

   //////////////////////////////
   // General helper functions //
   //////////////////////////////

   static inline bool isFaceDirection( stencil::Direction dir );
   static inline bool isEdgeDirection( stencil::Direction dir );
   static inline bool isCornerDirection( stencil::Direction dir );

   static inline bool blocksConnectedByFaces( const Block * block, const BlockID & neighbor );
   static inline bool blocksConnectedByEdges( const Block * block, const BlockID & neighbor );

   static inline bool divisibleByTwo( const CellInterval & cellBB );

   static bool coarserNeighborExists( const Block * block, stencil::Direction dir );



   const BlockDataID fieldId_;
};



/////////////////
// Equal level //
/////////////////

template< typename Field_T, typename Stencil >
void PackInfo< Field_T, Stencil >::packDataEqualLevelImpl( const Block * sender, stencil::Direction dir, mpi::SendBuffer & buffer ) const
{
#ifndef NDEBUG
   if( Stencil::D == uint_t(2) )
      WALBERLA_ASSERT_EQUAL( stencil::cz[dir], 0 );
#endif

   const Field_T * field = sender->getData< Field_T >( fieldId_ );

   CellInterval packingInterval = equalLevelPackInterval( dir, field->xyzSize(), uint_t(1) );

   for( auto cell = field->beginSliceXYZ( packingInterval ); cell != field->end(); ++cell )
      for( uint_t idx = 0; idx < Field_T::F_SIZE; ++idx )
         buffer << cell.getF( idx );
}



template< typename Field_T, typename Stencil >
void PackInfo< Field_T, Stencil >::unpackDataEqualLevel( Block * receiver, stencil::Direction dir, mpi::RecvBuffer & buffer )
{
#ifndef NDEBUG
   if( Stencil::D == uint_t(2) )
      WALBERLA_ASSERT_EQUAL( stencil::cz[dir], 0 );
#endif

   Field_T * field = receiver->getData< Field_T >( fieldId_ );

   CellInterval unpackingInterval = equalLevelUnpackInterval( dir, field->xyzSize(), uint_t(1) );

   for( auto cell = field->beginSliceXYZ( unpackingInterval ); cell != field->end(); ++cell )
      for( uint_t idx = 0; idx < Field_T::F_SIZE; ++idx )
         buffer >> cell.getF( idx );
}



template< typename Field_T, typename Stencil >
void PackInfo< Field_T, Stencil >::communicateLocalEqualLevel( const Block * sender, Block * receiver, stencil::Direction dir )
{
#ifndef NDEBUG
   if( Stencil::D == uint_t(2) )
      WALBERLA_ASSERT_EQUAL( stencil::cz[dir], 0 );
#endif

   const Field_T * sf =   sender->getData< Field_T >( fieldId_ );
         Field_T * rf = receiver->getData< Field_T >( fieldId_ );

   WALBERLA_ASSERT_EQUAL( sf->xyzSize(), rf->xyzSize() );

   CellInterval   packingInterval = equalLevelPackInterval( dir, sf->xyzSize(), uint_t(1) );
   CellInterval unpackingInterval = equalLevelUnpackInterval( stencil::inverseDir[dir], rf->xyzSize(), uint_t(1) );

   auto sCell = sf->beginSliceXYZ(   packingInterval );
   auto rCell = rf->beginSliceXYZ( unpackingInterval );
   while( sCell != sf->end() )
   {
      WALBERLA_ASSERT( rCell != rf->end() );

      for( uint_t idx = 0; idx < Field_T::F_SIZE; ++idx )
         rCell.getF( idx ) = sCell.getF( idx );

      ++sCell;
      ++rCell;
   }
   WALBERLA_ASSERT( rCell == rf->end() );
}



////////////////////
// Coarse to fine //
////////////////////

template< typename Field_T, typename Stencil >
void PackInfo< Field_T, Stencil >::packDataCoarseToFineImpl( const Block * coarseSender, const BlockID & fineReceiver,
                                                                       stencil::Direction dir, mpi::SendBuffer & buffer ) const
{
#ifndef NDEBUG
   if( Stencil::D == uint_t(2) )
      WALBERLA_ASSERT_EQUAL( stencil::cz[dir], 0 );
#endif

   const Field_T * field = coarseSender->getData< Field_T >( fieldId_ );

   CellInterval packingInterval = coarseToFinePackInterval( dir, field->xyzSize(), fineReceiver );

   for( cell_idx_t z = packingInterval.zMin(); z <= packingInterval.zMax(); ++z )
      for( cell_idx_t y = packingInterval.yMin(); y <= packingInterval.yMax(); ++y )
         for( cell_idx_t x = packingInterval.xMin(); x <= packingInterval.xMax(); ++x )
            for( uint_t idx = 0; idx < Field_T::F_SIZE; ++idx )
               buffer << field->get( x, y, z, idx );
}



template< typename Field_T, typename Stencil >
void PackInfo< Field_T, Stencil >::unpackDataCoarseToFine( Block * fineReceiver, const BlockID & /*coarseSender*/,
                                                                     stencil::Direction dir, mpi::RecvBuffer & buffer )
{
#ifndef NDEBUG
   if( Stencil::D == uint_t(2) )
      WALBERLA_ASSERT_EQUAL( stencil::cz[dir], 0 );
#endif

   Field_T * field = fineReceiver->getData< Field_T >( fieldId_ );

   CellInterval unpackingInterval = coarseToFineUnpackInterval( dir, field->xyzSize(), fineReceiver->getId() );

   for( cell_idx_t z = unpackingInterval.zMin(); z <= unpackingInterval.zMax(); z += cell_idx_t(2) ) {
      for( cell_idx_t y = unpackingInterval.yMin(); y <= unpackingInterval.yMax(); y += cell_idx_t(2) ) {
         for( cell_idx_t x = unpackingInterval.xMin(); x <= unpackingInterval.xMax(); x += cell_idx_t(2) ) {
            for( uint_t idx = 0; idx < Field_T::F_SIZE; ++idx )
            {
               typename Field_T::value_type value;
               buffer >> value;

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



template< typename Field_T, typename Stencil >
void PackInfo< Field_T, Stencil >::communicateLocalCoarseToFine( const Block * coarseSender, Block * fineReceiver, stencil::Direction dir )
{
#ifndef NDEBUG
   if( Stencil::D == uint_t(2) )
      WALBERLA_ASSERT_EQUAL( stencil::cz[dir], 0 );
#endif

   const Field_T * sf = coarseSender->getData< Field_T >( fieldId_ );
         Field_T * rf = fineReceiver->getData< Field_T >( fieldId_ );

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

   cell_idx_t rz = unpackingInterval.zMin();
   for( cell_idx_t sz = packingInterval.zMin(); sz <= packingInterval.zMax(); ++sz )
   {
      cell_idx_t ry = unpackingInterval.yMin();
      for( cell_idx_t sy = packingInterval.yMin(); sy <= packingInterval.yMax(); ++sy )
      {
         cell_idx_t rx = unpackingInterval.xMin();
         for( cell_idx_t sx = packingInterval.xMin(); sx <= packingInterval.xMax(); ++sx ) {
            for( uint_t idx = 0; idx < Field_T::F_SIZE; ++idx )
            {
               const auto & value = sf->get(sx,sy,sz,idx);

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



////////////////////
// Fine to coarse //
////////////////////

template< typename Field_T, typename Stencil >
void PackInfo< Field_T, Stencil >::packDataFineToCoarseImpl( const Block * fineSender, const BlockID & coarseReceiver,
                                                                       stencil::Direction dir, mpi::SendBuffer & buffer ) const
{
#ifndef NDEBUG
   if( Stencil::D == uint_t(2) )
      WALBERLA_ASSERT_EQUAL( stencil::cz[dir], 0 );
#endif

   if( ( ( isEdgeDirection(dir) || isCornerDirection(dir) ) && blocksConnectedByFaces( fineSender, coarseReceiver ) ) ||
       ( isCornerDirection(dir) && blocksConnectedByEdges( fineSender, coarseReceiver ) ) )
      return;

   const Field_T * field = fineSender->getData< Field_T >( fieldId_ );

   CellInterval packingInterval = fineToCoarsePackInterval( dir, field->xyzSize() );

   const real_t factor = ( Stencil::D == uint_t(3) ) ? real_t( 0.125 ) : real_t( 0.25 );

   for( cell_idx_t z = packingInterval.zMin(); z <= packingInterval.zMax(); z += cell_idx_t(2) ) {
      for( cell_idx_t y = packingInterval.yMin(); y <= packingInterval.yMax(); y += cell_idx_t(2) ) {
         for( cell_idx_t x = packingInterval.xMin(); x <= packingInterval.xMax(); x += cell_idx_t(2) ) {
            for( uint_t idx = 0; idx < Field_T::F_SIZE; ++idx )
            {
               typename Field_T::value_type value  = field->get( x,                 y,                 z,                 idx );
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



template< typename Field_T, typename Stencil >
void PackInfo< Field_T, Stencil >::unpackDataFineToCoarse( Block * coarseReceiver, const BlockID & fineSender,
                                                                     stencil::Direction dir, mpi::RecvBuffer & buffer )
{
#ifndef NDEBUG
   if( Stencil::D == uint_t(2) )
      WALBERLA_ASSERT_EQUAL( stencil::cz[dir], 0 );
#endif

   if( ( ( isEdgeDirection(dir) || isCornerDirection(dir) ) && blocksConnectedByFaces( coarseReceiver, fineSender ) ) ||
       ( isCornerDirection(dir) && blocksConnectedByEdges( coarseReceiver, fineSender ) ) )
      return;

   Field_T * field = coarseReceiver->getData< Field_T >( fieldId_ );

   CellInterval unpackingInterval = fineToCoarseUnpackInterval( dir, field->xyzSize(), fineSender );

   for( cell_idx_t z = unpackingInterval.zMin(); z <= unpackingInterval.zMax(); ++z )
      for( cell_idx_t y = unpackingInterval.yMin(); y <= unpackingInterval.yMax(); ++y )
         for( cell_idx_t x = unpackingInterval.xMin(); x <= unpackingInterval.xMax(); ++x )
            for( uint_t idx = 0; idx < Field_T::F_SIZE; ++idx )
               buffer >> field->get( x, y, z, idx );
}



template< typename Field_T, typename Stencil >
void PackInfo< Field_T, Stencil >::communicateLocalFineToCoarse( const Block * fineSender, Block * coarseReceiver, stencil::Direction dir )
{
#ifndef NDEBUG
   if( Stencil::D == uint_t(2) )
      WALBERLA_ASSERT_EQUAL( stencil::cz[dir], 0 );
#endif

   if( ( ( isEdgeDirection(dir) || isCornerDirection(dir) ) && blocksConnectedByFaces( fineSender, coarseReceiver->getId() ) ) ||
       ( isCornerDirection(dir) && blocksConnectedByEdges( fineSender, coarseReceiver->getId() ) ) )
      return;

   const Field_T * sf =     fineSender->getData< Field_T >( fieldId_ );
         Field_T * rf = coarseReceiver->getData< Field_T >( fieldId_ );

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
            for( uint_t idx = 0; idx < Field_T::F_SIZE; ++idx )
            {
               typename Field_T::value_type value  = sf->get( sx,                 sy,                 sz,                 idx );
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

template< typename Field_T, typename Stencil >
inline CellInterval PackInfo< Field_T, Stencil >::equalLevelPackInterval( stencil::Direction dir, const CellInterval & cellBB,
                                                                                    const uint_t numberOfLayers )
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



template< typename Field_T, typename Stencil >
CellInterval PackInfo< Field_T, Stencil >::equalLevelUnpackInterval( stencil::Direction dir, const CellInterval & cellBB,
                                                                               const uint_t numberOfLayers )
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



template< typename Field_T, typename Stencil >
inline Vector3< cell_idx_t > PackInfo< Field_T, Stencil >::getNeighborShift( const BlockID & smallBlock, stencil::Direction dir ) // dir: direction from big to small block
{
   Vector3< cell_idx_t > shift;

   uint_t branchId = smallBlock.getBranchId();

   shift[0] = ( stencil::cx[dir] == 0 ) ? ( ( ( branchId & uint_t(1) ) == uint_t(0) ) ? cell_idx_t(-1) : cell_idx_t(1) ) : cell_idx_t(0);
   shift[1] = ( stencil::cy[dir] == 0 ) ? ( ( ( branchId & uint_t(2) ) == uint_t(0) ) ? cell_idx_t(-1) : cell_idx_t(1) ) : cell_idx_t(0);
   shift[2] = ( Stencil::D == uint_t(3) ) ?
            ( ( stencil::cz[dir] == 0 ) ? ( ( ( branchId & uint_t(4) ) == uint_t(0) ) ? cell_idx_t(-1) : cell_idx_t(1) ) : cell_idx_t(0) ) : cell_idx_t(0);

   return shift;
}



template< typename Field_T, typename Stencil >
CellInterval PackInfo< Field_T, Stencil >::coarseToFinePackInterval( stencil::Direction dir, const CellInterval & cellBB,
                                                                               const BlockID & smallBlock )
{
   CellInterval interval = equalLevelPackInterval( dir, cellBB, uint_t(1) );
   Vector3< cell_idx_t > shift = getNeighborShift( smallBlock, dir );

   WALBERLA_ASSERT( divisibleByTwo( cellBB ) );

   for( uint_t i = 0; i != Stencil::D; ++i )
   {
      if( shift[i] == cell_idx_t(-1) )
         interval.max()[i] = interval.min()[i] + cell_idx_c( cellBB.size(i) / uint_t(2) );
      if( shift[i] == cell_idx_t( 1) )
         interval.min()[i] = interval.max()[i] - cell_idx_c( cellBB.size(i) / uint_t(2) );
   }

   WALBERLA_ASSERT( cellBB.contains( interval ) );

   return interval;
}



template< typename Field_T, typename Stencil >
CellInterval PackInfo< Field_T, Stencil >::coarseToFineUnpackInterval( stencil::Direction dir, const CellInterval & cellBB,
                                                                                 const BlockID & smallBlock )
{
   CellInterval interval = equalLevelUnpackInterval( dir, cellBB, uint_t(2) );
   Vector3< cell_idx_t > shift = getNeighborShift( smallBlock, dir );

   for( uint_t i = 0; i != Stencil::D; ++i )
   {
      if( shift[i] == cell_idx_t(-1) )
         interval.max()[i] += cell_idx_t(2);
      if( shift[i] == cell_idx_t( 1) )
         interval.min()[i] -= cell_idx_t(2);
   }

#ifndef NDEBUG
   CellInterval expandedCellBB( cellBB );
   expandedCellBB.expand( cell_idx_t(2) );
   WALBERLA_ASSERT( expandedCellBB.contains( interval ) );
#endif

   return interval;
}



template< typename Field_T, typename Stencil >
inline CellInterval PackInfo< Field_T, Stencil >::fineToCoarsePackInterval( stencil::Direction dir, const CellInterval & cellBB )
{
   return equalLevelPackInterval( dir, cellBB, uint_t(2) );
}



template< typename Field_T, typename Stencil >
CellInterval PackInfo< Field_T, Stencil >::fineToCoarseUnpackInterval( stencil::Direction dir, const CellInterval & cellBB,
                                                                                 const BlockID & smallBlock )
{
   CellInterval interval = equalLevelUnpackInterval( dir, cellBB, uint_t(1) );
   Vector3< cell_idx_t > shift = getNeighborShift( smallBlock, dir );

   WALBERLA_ASSERT( divisibleByTwo( cellBB ) );

   for( uint_t i = 0; i != Stencil::D; ++i )
   {
      if( shift[i] == cell_idx_t(-1) )
         interval.max()[i] = interval.min()[i] + cell_idx_c( cellBB.size(i) / uint_t(2) ) - cell_idx_t(1);
      if( shift[i] == cell_idx_t( 1) )
         interval.min()[i] = interval.max()[i] - cell_idx_c( cellBB.size(i) / uint_t(2) ) + cell_idx_t(1);
   }

#ifndef NDEBUG
   CellInterval expandedCellBB( cellBB );
   expandedCellBB.expand( cell_idx_t(1) );
   WALBERLA_ASSERT( expandedCellBB.contains( interval ) );
#endif

   return interval;
}



//////////////////////////////
// General helper functions //
//////////////////////////////

template< typename Field_T, typename Stencil >
inline bool PackInfo< Field_T, Stencil >::isFaceDirection( stencil::Direction dir )
{
   return ( dir == stencil::N ) || ( dir == stencil::S ) || ( dir == stencil::W ) ||
          ( dir == stencil::E ) || ( dir == stencil::T ) || ( dir == stencil::B );
}



template< typename Field_T, typename Stencil >
inline bool PackInfo< Field_T, Stencil >::isEdgeDirection( stencil::Direction dir )
{
   return ( dir == stencil::NW ) || ( dir == stencil::NE ) || ( dir == stencil::SW ) || ( dir == stencil::SE ) ||
          ( dir == stencil::TN ) || ( dir == stencil::TS ) || ( dir == stencil::TW ) || ( dir == stencil::TE ) ||
          ( dir == stencil::BN ) || ( dir == stencil::BS ) || ( dir == stencil::BW ) || ( dir == stencil::BE );
}



template< typename Field_T, typename Stencil >
inline bool PackInfo< Field_T, Stencil >::isCornerDirection( stencil::Direction dir )
{
   return ( dir == stencil::TNE ) || ( dir == stencil::TNW ) || ( dir == stencil::TSE ) || ( dir == stencil::TSW ) ||
          ( dir == stencil::BNE ) || ( dir == stencil::BNW ) || ( dir == stencil::BSE ) || ( dir == stencil::BSW );
}



template< typename Field_T, typename Stencil >
inline bool PackInfo< Field_T, Stencil >::blocksConnectedByFaces( const Block * block, const BlockID & neighbor )
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



template< typename Field_T, typename Stencil >
inline bool PackInfo< Field_T, Stencil >::blocksConnectedByEdges( const Block * block, const BlockID & neighbor )
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



template< typename Field_T, typename Stencil >
inline bool PackInfo< Field_T, Stencil >::divisibleByTwo( const CellInterval & cellBB )
{
   return ( ( cellBB.xSize() & uint_t(1) ) == uint_t(0) ) &&
          ( ( cellBB.ySize() & uint_t(1) ) == uint_t(0) ) &&
          ( ( Stencil::D == uint_t(2) && cellBB.zSize() == uint_t(1) ) || ( Stencil::D == uint_t(3) && ( cellBB.zSize() & uint_t(1) ) == uint_t(0) ) );
}



template< typename Field_T, typename Stencil >
bool PackInfo< Field_T, Stencil >::coarserNeighborExists( const Block * block, stencil::Direction dir )
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
   if( Stencil::D == uint_t(2) )
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
} // namespace field
} // namespace walberla
