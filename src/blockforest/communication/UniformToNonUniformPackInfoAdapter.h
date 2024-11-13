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
//! \file UniformToNonUniformPackInfoAdapter.h
//! \ingroup blockforest
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//
//======================================================================================================================

#pragma once

#include "NonUniformPackInfo.h"

#include "communication/UniformPackInfo.h"

namespace walberla {
namespace blockforest {
namespace communication {


//*******************************************************************************************************************
/*! 
 * Adapter to use a \ref communication::UniformPackInfo in a \ref NonUniformBufferedScheme. No communication between coarse <-> fine blocks
 * happens.
 */
//*******************************************************************************************************************
class UniformToNonUniformPackInfoAdapter : public NonUniformPackInfo
{
public:

   //**Construction & Destruction************************************************************************************
   /*! \name Construction & Destruction */
   //@{
   UniformToNonUniformPackInfoAdapter( const shared_ptr<walberla::communication::UniformPackInfo> & uniformPackInfo ) : uniformPackInfo_( uniformPackInfo ) { }
   ~UniformToNonUniformPackInfoAdapter() override = default;
   //@}
   //****************************************************************************************************************

   /**
   * Should return true if the amount of data that is packed for a given block in direction
   * "dir" is guaranteed to remain constant over time. False otherwise.
   * If you are not sure what to return, return false! Returning false is always safe.
   * Falsely return true will lead to errors! However, if the data can be guaranteed to remain
   * constant over time, returning true enables performance optimizations during the communication.
   */
   bool constantDataExchange() const override { return uniformPackInfo_->constantDataExchange(); }

   /**
   * Must return false if calling `unpackData*()` and/or `communicateLocal*()` methods is not thread-safe.
   * True otherwise.
   * If you are not sure what to return, return false! Returning false is always safe.
   * Falsely return true will most likely lead to errors! However, if both `unpackData*()` AND
   * `communicateLocal*()` are thread-safe, returning true can lead to performance improvements.
   */
   bool threadsafeReceiving() const override { return uniformPackInfo_->threadsafeReceiving(); }

   /// If NOT thread-safe, \ref threadsafeReceiving must return false!
   void unpackDataEqualLevel( Block * receiver, stencil::Direction dir, mpi::RecvBuffer & buffer ) override { uniformPackInfo_->unpackData( receiver, dir, buffer ); }

   /// If NOT thread-safe, \ref threadsafeReceiving must return false!
   void communicateLocalEqualLevel( const Block * sender, Block * receiver, stencil::Direction dir ) override { uniformPackInfo_->communicateLocal( sender, receiver, dir ); }

   void unpackDataCoarseToFine( Block * /*fineReceiver*/, const BlockID & /*coarseSender*/, stencil::Direction /*dir*/, mpi::RecvBuffer & /*buffer*/ ) override { }
   void communicateLocalCoarseToFine( const Block * /*coarseSender*/, Block * /*fineReceiver*/, stencil::Direction /*dir*/ ) override { }

   void unpackDataFineToCoarse( Block * /*coarseReceiver*/, const BlockID & /*fineSender*/, stencil::Direction /*dir*/, mpi::RecvBuffer & /*buffer*/ ) override { }
   void communicateLocalFineToCoarse( const Block * /*fineSender*/, Block * /*coarseReceiver*/, stencil::Direction /*dir*/ ) override { }

protected:

   shared_ptr<walberla::communication::UniformPackInfo> uniformPackInfo_;

   /// Must be thread-safe!
   void packDataEqualLevelImpl( const Block * sender, stencil::Direction dir, mpi::SendBuffer & buffer ) const override { uniformPackInfo_->packData( sender, dir, buffer ); }

   void packDataCoarseToFineImpl( const Block * /*coarseSender*/, const BlockID &   /*fineReceiver*/, stencil::Direction /*dir*/, mpi::SendBuffer & /*buffer*/ ) const override { }
   void packDataFineToCoarseImpl( const Block *   /*fineSender*/, const BlockID & /*coarseReceiver*/, stencil::Direction /*dir*/, mpi::SendBuffer & /*buffer*/ ) const override { }
};


} // namespace communication
} // namespace blockforest
} // namespace walberla
