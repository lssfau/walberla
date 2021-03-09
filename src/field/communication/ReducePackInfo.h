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
//! \file ReduceFieldPackInfo.h
//! \ingroup field
//! \author Matthias Markl <matthias.markl@fau.de>
//! \brief Stores/reduces ghost layers to/from a communication buffer
//
//======================================================================================================================

#pragma once

#include "field/GhostLayerField.h"
#include "communication/ReducePackInfo.h"
#include "core/debug/Debug.h"
#include "domain_decomposition/BlockDataID.h"


namespace walberla {
namespace field {
namespace communication {


/**
 * Data packing/unpacking for reduce operations
 * This class can be used, when multiple data sets from different senders should be reduced at the receiver
 * \ingroup comm
 */
template< template<typename> class ReduceOperation, typename GhostLayerField_T >
class ReducePackInfo : public walberla::communication::ReducePackInfo
{
public:
   using T = typename GhostLayerField_T::value_type;

   ReducePackInfo( const BlockDataID & bdId, T init )
        : bdId_(bdId),
          init_(init),
          op_()
   {}

   ~ReducePackInfo() override = default;

   bool constantDataExchange() const { return false; }
   bool threadsafeReceiving()  const { return false; }

   void safeCommunicateLocal( const IBlock * sender, IBlock * receiver, stencil::Direction dir ) override;
   void packData            ( const IBlock * sender,   stencil::Direction dir, mpi::SendBuffer & outBuffer ) override;
   void safeUnpackData      (       IBlock * receiver, stencil::Direction dir, mpi::RecvBuffer & buffer    ) override;

   size_t getDataSize() const { return getSize() * sizeof(T); }

protected:
   size_t initData( IBlock * receiver, stencil::Direction dir ) override;

   const BlockDataID        bdId_;
   const T                  init_;
   const ReduceOperation<T> op_;
};



template< template<typename> class ReduceOperation, typename GhostLayerField_T >
size_t ReducePackInfo< ReduceOperation, GhostLayerField_T >::initData( IBlock * receiver, stencil::Direction dir )
{
   GhostLayerField_T * f = receiver->getData< GhostLayerField_T > (bdId_);
   WALBERLA_ASSERT_NOT_NULLPTR(f);

   size_t size( 0u );
   for( auto i = f->beginSliceBeforeGhostLayer(dir); i != f->end(); ++i ){
      *i = init_;
      ++size;
   }
   return size;
}



template< template<typename> class ReduceOperation, typename GhostLayerField_T >
void ReducePackInfo< ReduceOperation, GhostLayerField_T >::safeCommunicateLocal( const IBlock * sender, IBlock * receiver, stencil::Direction dir )
{
   const GhostLayerField_T * sf = sender  ->getData< GhostLayerField_T > (bdId_);
         GhostLayerField_T * rf = receiver->getData< GhostLayerField_T > (bdId_);

   WALBERLA_ASSERT_EQUAL(sf->xSize(), rf->xSize());
   WALBERLA_ASSERT_EQUAL(sf->ySize(), rf->ySize());
   WALBERLA_ASSERT_EQUAL(sf->zSize(), rf->zSize());

   auto srcIter = sf->beginSliceBeforeGhostLayer(dir);
   auto dstIter = rf->beginGhostLayerOnly(stencil::inverseDir[dir]);

   while( srcIter != sf->end() ) {
      *dstIter = op_(*srcIter, *dstIter);
      ++srcIter;
      ++dstIter;
   }
   WALBERLA_ASSERT(srcIter == sf->end() && dstIter == rf->end());
}



template< template<typename> class ReduceOperation, typename GhostLayerField_T >
void ReducePackInfo< ReduceOperation, GhostLayerField_T>::packData( const IBlock * sender, stencil::Direction dir,
                                                                    mpi::GenericSendBuffer<unsigned char> & outBuffer )
{
   const GhostLayerField_T * f = sender->getData< GhostLayerField_T > (bdId_);
   WALBERLA_ASSERT_NOT_NULLPTR(f);

   for( auto i = f->beginSliceBeforeGhostLayer(dir); i != f->end(); ++i )
      outBuffer << *i;
}



template< template<typename> class ReduceOperation, typename GhostLayerField_T >
void ReducePackInfo< ReduceOperation, GhostLayerField_T >::safeUnpackData( IBlock* receiver, stencil::Direction dir,
                                                                           mpi::GenericRecvBuffer<unsigned char> & buffer )
{
   GhostLayerField_T * f = receiver->getData< GhostLayerField_T > (bdId_);
   WALBERLA_ASSERT_NOT_NULLPTR(f);

   T buf(0);
   for(auto i = f->beginGhostLayerOnly(dir); i != f->end(); ++i ){
      buffer >> buf;
      *i = op_( *i, buf );
   }
}


} // namespace communication
} // namespace field
} // namespace walberla
