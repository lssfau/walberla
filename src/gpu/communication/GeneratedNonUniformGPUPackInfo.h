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
//! \file GeneratedNonUniformGPUPackInfo.h
//! \ingroup gpu
//! \author Markus Holzer <markus.holzer@fau.de>
//
//======================================================================================================================

#pragma once

#include "blockforest/Block.h"
#include "blockforest/BlockID.h"

#include "gpu/GPUWrapper.h"
#include "gpu/communication/CustomMemoryBuffer.h"

#include "stencil/Directions.h"

using GpuBuffer_T = walberla::gpu::communication::GPUMemoryBuffer;


namespace walberla::gpu {


class GeneratedNonUniformGPUPackInfo
{
 public:
   using VoidFunction                  = std::function< void( gpuStream_t) >;
   GeneratedNonUniformGPUPackInfo() = default;
   virtual ~GeneratedNonUniformGPUPackInfo() = default;

   virtual bool constantDataExchange() const = 0;
   virtual bool threadsafeReceiving() const = 0;

   inline void packDataEqualLevel( const Block * sender, stencil::Direction dir, GpuBuffer_T & buffer) const;
   virtual void unpackDataEqualLevel( Block * receiver, stencil::Direction dir, GpuBuffer_T & buffer) = 0;
   virtual void communicateLocalEqualLevel( const Block * sender, Block * receiver, stencil::Direction dir, gpuStream_t stream) = 0;
   virtual void getLocalEqualLevelCommFunction( std::vector< VoidFunction >& commFunctions, const Block * sender, Block * receiver, stencil::Direction dir) = 0;

   inline  void packDataCoarseToFine        ( const Block * coarseSender, const BlockID & fineReceiver, stencil::Direction dir, GpuBuffer_T & buffer) const;
   virtual void unpackDataCoarseToFine      (       Block * fineReceiver, const BlockID & coarseSender, stencil::Direction dir, GpuBuffer_T & buffer) = 0;
   virtual void communicateLocalCoarseToFine( const Block * coarseSender, Block * fineReceiver, stencil::Direction dir ) = 0;
   virtual void communicateLocalCoarseToFine( const Block * coarseSender, Block * fineReceiver, stencil::Direction dir, GpuBuffer_T & buffer, gpuStream_t stream) = 0;
   virtual void getLocalCoarseToFineCommFunction( std::vector< VoidFunction >& commFunctions, const Block * coarseSender, Block * fineReceiver, stencil::Direction dir, GpuBuffer_T & buffer) = 0;

   inline  void packDataFineToCoarse        ( const Block * fineSender,     const BlockID & coarseReceiver, stencil::Direction dir, GpuBuffer_T & buffer) const;
   virtual void unpackDataFineToCoarse      (       Block * coarseReceiver, const BlockID & fineSender,     stencil::Direction dir, GpuBuffer_T & buffer) = 0;
   virtual void communicateLocalFineToCoarse( const Block * fineSender, Block * coarseReceiver, stencil::Direction dir) = 0;
   virtual void communicateLocalFineToCoarse( const Block * fineSender, Block * coarseReceiver, stencil::Direction dir, GpuBuffer_T & buffer, gpuStream_t stream) = 0;
   virtual void getLocalFineToCoarseCommFunction( std::vector< VoidFunction >& commFunctions, const Block * fineSender, Block * coarseReceiver, stencil::Direction dir, GpuBuffer_T & buffer) = 0;

   virtual uint_t sizeEqualLevelSend( const Block * sender, stencil::Direction dir) = 0;
   virtual uint_t sizeCoarseToFineSend ( const Block * coarseSender, const BlockID & fineReceiver, stencil::Direction dir) = 0;
   virtual uint_t sizeFineToCoarseSend ( const Block * fineSender, stencil::Direction dir) = 0;


#ifndef NDEBUG
   void clearBufferSizeCheckMap() { bufferSize_.clear(); }
#endif

 protected:
   virtual void packDataEqualLevelImpl(const Block* sender, stencil::Direction dir, GpuBuffer_T & buffer) const = 0;
   virtual void packDataCoarseToFineImpl(const Block* coarseSender, const BlockID& fineReceiver, stencil::Direction dir, GpuBuffer_T & buffer) const = 0;
   virtual void packDataFineToCoarseImpl(const Block* fineSender, const BlockID& coarseReceiver, stencil::Direction dir, GpuBuffer_T & buffer) const = 0;

#ifndef NDEBUG
   mutable std::map< const Block *, std::map< stencil::Direction, std::map< uint_t, size_t > > > bufferSize_;
#endif

};

inline void GeneratedNonUniformGPUPackInfo::packDataEqualLevel( const Block * sender, stencil::Direction dir, GpuBuffer_T & buffer ) const
{
#ifndef NDEBUG
   size_t const sizeBefore = buffer.size();
#endif

   packDataEqualLevelImpl( sender, dir, buffer );

#ifndef NDEBUG
size_t const sizeAfter = buffer.size();
if( constantDataExchange() )
{
      auto & blockMap = bufferSize_[ sender ];
      auto & sizeMap  = blockMap[ dir ];
      auto dirEntry = sizeMap.find( uint_t(0) );
      if( dirEntry == sizeMap.end() )
         sizeMap[ uint_t(0) ] = sizeAfter - sizeBefore;
      else
         WALBERLA_ASSERT_EQUAL( sizeMap[ uint_t(0) ], (sizeAfter - sizeBefore) )
}
#endif
}



inline void GeneratedNonUniformGPUPackInfo::packDataCoarseToFine( const Block * coarseSender, const BlockID & fineReceiver, stencil::Direction dir, GpuBuffer_T & buffer ) const
{
#ifndef NDEBUG
   size_t const sizeBefore = buffer.size();
#endif

   packDataCoarseToFineImpl( coarseSender, fineReceiver, dir, buffer );

#ifndef NDEBUG
size_t const sizeAfter = buffer.size();
if( constantDataExchange() )
{
      auto & blockMap = bufferSize_[ coarseSender ];
      auto & sizeMap  = blockMap[ dir ];
      auto dirEntry = sizeMap.find( fineReceiver.getBranchId() );
      if( dirEntry == sizeMap.end() )
         sizeMap[ fineReceiver.getBranchId() ] = sizeAfter - sizeBefore;
      else
         WALBERLA_ASSERT_EQUAL( sizeMap[ fineReceiver.getBranchId() ], (sizeAfter - sizeBefore) )
}
#endif
}



inline void GeneratedNonUniformGPUPackInfo::packDataFineToCoarse( const Block * fineSender, const BlockID & coarseReceiver, stencil::Direction dir, GpuBuffer_T & buffer ) const
{
#ifndef NDEBUG
   size_t const sizeBefore = buffer.size();
#endif

   packDataFineToCoarseImpl( fineSender, coarseReceiver, dir, buffer );

#ifndef NDEBUG
size_t const sizeAfter = buffer.size();
if( constantDataExchange() )
{
      auto & blockMap = bufferSize_[ fineSender ];
      auto & sizeMap  = blockMap[ dir ];
      auto dirEntry = sizeMap.find( uint_t(0) );
      if( dirEntry == sizeMap.end() )
         sizeMap[ uint_t(0) ] = sizeAfter - sizeBefore;
      else
         WALBERLA_ASSERT_EQUAL( sizeMap[ uint_t(0) ], (sizeAfter - sizeBefore) )
}
#endif
}


} //namespace walberla::gpu