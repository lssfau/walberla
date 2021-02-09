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
//! \file ReducePackInfo.h
//! \ingroup communication
//! \author Matthias Markl <matthias.markl@fau.de>
//! \brief Stores/reduces ghost layers to/from a communication buffer
//
//======================================================================================================================

#pragma once

#include "core/mpi/RecvBuffer.h"
#include "core/mpi/SendBuffer.h"
#include "domain_decomposition/IBlock.h"
#include "stencil/Directions.h"

#include <set>


namespace walberla {
namespace communication {

/**
 * Data packing/unpacking for reduce operations
 * This class can be used, when multiple data sets from different senders should be reduced at the receiver
 * \ingroup comm
 */
class ReducePackInfo
{
public:
   ReducePackInfo( ) : size_(0u) {}
   virtual ~ReducePackInfo() = default;

   size_t getSize() const { return size_; }

   void reset(){ initBlocks_.clear(); }

   void communicateLocal( const IBlock * sender, IBlock * receiver, stencil::Direction dir );
   virtual void packData        ( const IBlock * sender,   stencil::Direction dir, mpi::SendBuffer & outBuffer ) = 0;
   void unpackData      (       IBlock * receiver, stencil::Direction dir, mpi::RecvBuffer & buffer    );

protected:
   virtual size_t initData( IBlock * receiver, stencil::Direction dir ) = 0;

   virtual void safeCommunicateLocal( const IBlock * sender, IBlock * receiver, stencil::Direction dir ) = 0;
   virtual void safeUnpackData      (       IBlock * receiver, stencil::Direction dir, mpi::RecvBuffer & buffer ) = 0;

private:
   size_t            size_;
   std::set<IBlock*> initBlocks_;
};


inline void ReducePackInfo::communicateLocal( const IBlock * sender, IBlock * receiver, stencil::Direction dir )
{
   if( initBlocks_.find(receiver) == initBlocks_.end() ){
      initBlocks_.insert(receiver);
      size_ = initData(receiver,dir);
   }
   safeCommunicateLocal( sender, receiver, dir );
}


inline void ReducePackInfo::unpackData( IBlock* receiver, stencil::Direction dir, mpi::RecvBuffer & buffer )
{
   if( initBlocks_.find(receiver) == initBlocks_.end() ){
      initBlocks_.insert(receiver);
      size_ = initData(receiver,dir);
   }
   safeUnpackData( receiver, dir, buffer );
}

} // namespace communication
} // namespace walberla
