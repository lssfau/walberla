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
//! \file NonUniformPackInfo.h
//! \ingroup blockforest
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//! \author David Staubach <david.staubach@fau.de>
//
//======================================================================================================================

#pragma once

#include "blockforest/Block.h"
#include "blockforest/BlockID.h"

#include "core/debug/Debug.h"
#include "core/mpi/RecvBuffer.h"
#include "core/mpi/SendBuffer.h"

#include "stencil/Directions.h"


namespace walberla {
namespace blockforest {
namespace communication {


class NonUniformPackInfo
{
public:

   //**Construction & Destruction************************************************************
   /*! \name Construction & Destruction */
   //@{
            NonUniformPackInfo() = default;
   virtual ~NonUniformPackInfo() = default;
   //@}
   //*******************************************************************************************************************

   /**
    * Should return true if the amount of data that is packed for a given block in direction
    * "dir" is guaranteed to remain constant over time. False otherwise.
    * If you are not sure what to return, return false! Returning false is always safe.
    * Falsely return true will lead to errors! However, if the data can be guaranteed to remain
    * constant over time, returning true enables performance optimizations during the communication.
    */
   virtual bool constantDataExchange() const = 0;

   /**
    * Must return false if calling `unpackData*()` and/or `communicateLocal*()` methods is not thread-safe.
    * True otherwise.
    * If you are not sure what to return, return false! Returning false is always safe.
    * Falsely return true will most likely lead to errors! However, if both `unpackData*()` AND
    * `communicateLocal*()` are thread-safe, returning true can lead to performance improvements.
    */
   virtual bool threadsafeReceiving() const = 0;

   /// Must be thread-safe! Calls \ref packDataEqualLevelImpl.
   inline void packDataEqualLevel( const Block * sender, stencil::Direction dir, mpi::SendBuffer & buffer ) const;

   /// If NOT thread-safe, \ref threadsafeReceiving must return false!
   virtual void unpackDataEqualLevel( Block * receiver, stencil::Direction dir, mpi::RecvBuffer & buffer ) = 0;

   /// If NOT thread-safe, \ref threadsafeReceiving must return false!
   virtual void communicateLocalEqualLevel( const Block * sender, Block * receiver, stencil::Direction dir ) = 0;

   inline  void packDataCoarseToFine        ( const Block * coarseSender, const BlockID & fineReceiver, stencil::Direction dir, mpi::SendBuffer & buffer ) const;
   virtual void unpackDataCoarseToFine      (       Block * fineReceiver, const BlockID & coarseSender, stencil::Direction dir, mpi::RecvBuffer & buffer ) = 0;
   virtual void communicateLocalCoarseToFine( const Block * coarseSender, Block * fineReceiver, stencil::Direction dir ) = 0;

   inline  void packDataFineToCoarse        ( const Block * fineSender,     const BlockID & coarseReceiver, stencil::Direction dir, mpi::SendBuffer & buffer ) const;
   virtual void unpackDataFineToCoarse      (       Block * coarseReceiver, const BlockID & fineSender,     stencil::Direction dir, mpi::RecvBuffer & buffer ) = 0;
   virtual void communicateLocalFineToCoarse( const Block * fineSender, Block * coarseReceiver, stencil::Direction dir ) = 0;

#ifndef NDEBUG
   void clearBufferSizeCheckMap() { bufferSize_.clear(); }
#endif
   
protected:

   /// Must be thread-safe!
   virtual void packDataEqualLevelImpl( const Block * sender, stencil::Direction dir, mpi::SendBuffer & buffer ) const = 0;

   virtual void packDataCoarseToFineImpl( const Block * coarseSender, const BlockID &   fineReceiver, stencil::Direction dir, mpi::SendBuffer & buffer ) const = 0;
   virtual void packDataFineToCoarseImpl( const Block *   fineSender, const BlockID & coarseReceiver, stencil::Direction dir, mpi::SendBuffer & buffer ) const = 0;

#ifndef NDEBUG
   mutable std::map< const Block *, std::map< stencil::Direction, std::map< uint_t, size_t > > > bufferSize_;
#endif
};



inline void NonUniformPackInfo::packDataEqualLevel( const Block * sender, stencil::Direction dir, mpi::SendBuffer & buffer ) const
{
#ifndef NDEBUG
   size_t const sizeBefore = buffer.size();
#endif

   packDataEqualLevelImpl( sender, dir, buffer );

#ifndef NDEBUG
   size_t const sizeAfter = buffer.size();
   if( constantDataExchange() )
   {
#ifdef _OPENMP
      #pragma omp critical (NonUniformPackInfo_packData) // -> packData must be thread-safe!
      {
#endif
      auto & blockMap = bufferSize_[ sender ];
      auto & sizeMap  = blockMap[ dir ];
      auto dirEntry = sizeMap.find( uint_t(0) );
      if( dirEntry == sizeMap.end() )
         sizeMap[ uint_t(0) ] = sizeAfter - sizeBefore;
      else
         WALBERLA_ASSERT_EQUAL( sizeMap[ uint_t(0) ], (sizeAfter - sizeBefore) )
#ifdef _OPENMP
      }
#endif
   }
#endif
}



inline void NonUniformPackInfo::packDataCoarseToFine( const Block * coarseSender, const BlockID & fineReceiver, stencil::Direction dir, mpi::SendBuffer & buffer ) const
{
#ifndef NDEBUG
   size_t const sizeBefore = buffer.size();
#endif

   packDataCoarseToFineImpl( coarseSender, fineReceiver, dir, buffer );

#ifndef NDEBUG
   size_t const sizeAfter = buffer.size();
   if( constantDataExchange() )
   {
#ifdef _OPENMP
      #pragma omp critical (NonUniformPackInfo_packData) // -> packData must be thread-safe!
      {
#endif
      auto & blockMap = bufferSize_[ coarseSender ];
      auto & sizeMap  = blockMap[ dir ];
      auto dirEntry = sizeMap.find( fineReceiver.getBranchId() );
      if( dirEntry == sizeMap.end() )
         sizeMap[ fineReceiver.getBranchId() ] = sizeAfter - sizeBefore;
      else
         WALBERLA_ASSERT_EQUAL( sizeMap[ fineReceiver.getBranchId() ], (sizeAfter - sizeBefore) )
#ifdef _OPENMP
      }
#endif
   }
#endif
}



inline void NonUniformPackInfo::packDataFineToCoarse( const Block * fineSender, const BlockID & coarseReceiver, stencil::Direction dir, mpi::SendBuffer & buffer ) const
{
#ifndef NDEBUG
   size_t const sizeBefore = buffer.size();
#endif

   packDataFineToCoarseImpl( fineSender, coarseReceiver, dir, buffer );

#ifndef NDEBUG
   size_t const sizeAfter = buffer.size();
   if( constantDataExchange() )
   {
#ifdef _OPENMP
      #pragma omp critical (NonUniformPackInfo_packData) // -> packData must be thread-safe!
      {
#endif
      auto & blockMap = bufferSize_[ fineSender ];
      auto & sizeMap  = blockMap[ dir ];
      auto dirEntry = sizeMap.find( uint_t(0) );
      if( dirEntry == sizeMap.end() )
         sizeMap[ uint_t(0) ] = sizeAfter - sizeBefore;
      else
         WALBERLA_ASSERT_EQUAL( sizeMap[ uint_t(0) ], (sizeAfter - sizeBefore) )
#ifdef _OPENMP
      }
#endif
   }
#endif
}


} // namespace communication
} // namespace blockforest
} // namespace walberla
