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
//! \file UniformPullReductionPackInfo.h
//! \ingroup field
//! \author Tobias Schruff <schruff@iww.rwth-aachen.de>
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//
//======================================================================================================================

#pragma once

#include "communication/UniformPackInfo.h"

#include "core/debug/Debug.h"

#include "field/GhostLayerField.h"

#include "stencil/Directions.h"


namespace walberla {
namespace field {
namespace communication {

/**
 * Data packing/unpacking for ghost layer based communication of a single walberla::field::Field
 *
 * \ingroup comm
 *
 * template ReduceOperation is e.g. std::plus
 *
 * This pack info is used to apply a given reduce operation (e.g. +) to the values in the interior of a ghost layer field
 * together with the values coming from the sender's ghost layer.
 */
template< template<typename> class ReduceOperation, typename GhostLayerField_T >
class UniformPullReductionPackInfo : public walberla::communication::UniformPackInfo
{
public:
   typedef typename GhostLayerField_T::value_type T;

   UniformPullReductionPackInfo( const BlockDataID & bdID ) : bdID_( bdID ), communicateAllGhostLayers_( true ),
                                 numberOfGhostLayers_( 0 ) {}

   UniformPullReductionPackInfo( const BlockDataID & bdID, const uint_t numberOfGhostLayers ) : bdID_( bdID ),
                                 communicateAllGhostLayers_( false ), numberOfGhostLayers_(  numberOfGhostLayers ) {}

   ~UniformPullReductionPackInfo() override = default;

   bool constantDataExchange() const override { return mpi::BufferSizeTrait<T>::constantSize; }
   bool threadsafeReceiving()  const override { return true; }

   void unpackData(IBlock * receiver, stencil::Direction dir, mpi::RecvBuffer & buffer) override;

   void communicateLocal(const IBlock * sender, IBlock * receiver, stencil::Direction dir) override;

protected:
   void packDataImpl(const IBlock * sender, stencil::Direction dir, mpi::SendBuffer & outBuffer) const override;

   uint_t numberOfGhostLayersToCommunicate( const GhostLayerField_T * const field ) const;

   const  BlockDataID bdID_;
   bool   communicateAllGhostLayers_;
   uint_t numberOfGhostLayers_;
   ReduceOperation<T> op_;
};



template< template<typename> class ReduceOperation, typename GhostLayerField_T >
void UniformPullReductionPackInfo<ReduceOperation, GhostLayerField_T>::unpackData( IBlock * receiver, stencil::Direction dir, mpi::RecvBuffer & buffer )
{
   GhostLayerField_T * f = receiver->getData< GhostLayerField_T >( bdID_ );
   WALBERLA_ASSERT_NOT_NULLPTR(f);

   cell_idx_t nrOfGhostLayers = cell_idx_c( numberOfGhostLayersToCommunicate( f ) );

   T buf( 0 );
   for( auto i = f->beginSliceBeforeGhostLayer( dir, nrOfGhostLayers ); i != f->end(); ++i ) {
      buffer >> buf;
      *i = op_( *i, buf );
   }
}



template< template<typename> class ReduceOperation, typename GhostLayerField_T>
void UniformPullReductionPackInfo<ReduceOperation, GhostLayerField_T>::communicateLocal( const IBlock * sender, IBlock * receiver, stencil::Direction dir )
{
   const GhostLayerField_T * sf = sender  ->getData< GhostLayerField_T >( bdID_ );
         GhostLayerField_T * rf = receiver->getData< GhostLayerField_T >( bdID_ );

   WALBERLA_ASSERT_EQUAL(sf->xSize(), rf->xSize());
   WALBERLA_ASSERT_EQUAL(sf->ySize(), rf->ySize());
   WALBERLA_ASSERT_EQUAL(sf->zSize(), rf->zSize());

   uint_t nrOfGhostLayers = numberOfGhostLayersToCommunicate( sf );
   auto srcIter = sf->beginGhostLayerOnly( nrOfGhostLayers, dir );
   auto dstIter = rf->beginSliceBeforeGhostLayer( stencil::inverseDir[dir], cell_idx_c( nrOfGhostLayers ) );

   while( srcIter != sf->end() ) {
      *dstIter = op_( *srcIter, *dstIter );
      ++srcIter;
      ++dstIter;
   }

   WALBERLA_ASSERT(srcIter == sf->end() && dstIter == rf->end());
}



template< template<typename> class ReduceOperation, typename GhostLayerField_T>
void UniformPullReductionPackInfo<ReduceOperation, GhostLayerField_T>::packDataImpl( const IBlock * sender, stencil::Direction dir, mpi::SendBuffer & outBuffer ) const
{
   const GhostLayerField_T * f = sender->getData< GhostLayerField_T >( bdID_ );
   WALBERLA_ASSERT_NOT_NULLPTR(f);

   uint_t nrOfGhostLayers = numberOfGhostLayersToCommunicate( f );
   for( auto i = f->beginGhostLayerOnly( nrOfGhostLayers, dir ); i != f->end(); ++i )
      outBuffer << *i;
}



template< template<typename> class ReduceOperation, typename GhostLayerField_T>
uint_t UniformPullReductionPackInfo<ReduceOperation, GhostLayerField_T>::numberOfGhostLayersToCommunicate( const GhostLayerField_T * const field ) const
{
   if( communicateAllGhostLayers_ )
   {
      return field->nrOfGhostLayers();
   }
   else
   {
      WALBERLA_ASSERT_LESS_EQUAL( numberOfGhostLayers_, field->nrOfGhostLayers() );
      return numberOfGhostLayers_;
   }
}

} // namespace communication
} // namespace field
} // namespace walberla
