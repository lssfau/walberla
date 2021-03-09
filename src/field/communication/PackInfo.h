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
//! \ingroup field
//! \author Martin Bauer <martin.bauer@fau.de>
//! \brief Stores/loads ghost layers to/from a communication buffer
//
//======================================================================================================================

#pragma once

#include "field/GhostLayerField.h"
#include "communication/UniformPackInfo.h"
#include "core/debug/Debug.h"
#include "stencil/Directions.h"


namespace walberla {
namespace field {
namespace communication {


/**
 * Data packing/unpacking for ghost layer based communication of a single walberla::field::Field
 * \ingroup field
 * Template parameters are equivalent of the parameters of GhostLayerField that is communicated
 */
template<typename GhostLayerField_T>
class PackInfo : public walberla::communication::UniformPackInfo
{
public:
   using T = typename GhostLayerField_T::value_type;

   PackInfo( const BlockDataID & bdId ) : bdId_( bdId ), communicateAllGhostLayers_( true ),
                                          numberOfGhostLayers_( 0 ) {}

   PackInfo( const BlockDataID & bdId, const uint_t numberOfGhostLayers ) : bdId_( bdId ),
      communicateAllGhostLayers_( false ), numberOfGhostLayers_(  numberOfGhostLayers ) {}

   ~PackInfo() override = default;

   bool constantDataExchange() const override { return mpi::BufferSizeTrait<T>::constantSize; }
   bool threadsafeReceiving()  const override { return true; }

   void unpackData(IBlock * receiver, stencil::Direction dir, mpi::RecvBuffer & buffer) override;

   void communicateLocal(const IBlock * sender, IBlock * receiver, stencil::Direction dir) override;

protected:
   void packDataImpl(const IBlock * sender, stencil::Direction dir, mpi::SendBuffer & outBuffer) const override;

   uint_t numberOfGhostLayersToCommunicate( const GhostLayerField_T * const field ) const;

   const BlockDataID bdId_;
   bool   communicateAllGhostLayers_;
   uint_t numberOfGhostLayers_;
};



template<typename GhostLayerField_T>
void PackInfo<GhostLayerField_T>::unpackData(IBlock * receiver, stencil::Direction dir, mpi::RecvBuffer & buffer)
{
   GhostLayerField_T * f = receiver->getData< GhostLayerField_T >( bdId_ );
   WALBERLA_ASSERT_NOT_NULLPTR(f);

   uint_t nrOfGhostLayers = numberOfGhostLayersToCommunicate( f );

#ifndef NDEBUG
   uint_t xSize, ySize, zSize;
   buffer >> xSize >> ySize >> zSize;
   WALBERLA_ASSERT_EQUAL(xSize, f->xSize());
   WALBERLA_ASSERT_EQUAL(ySize, f->ySize());
   WALBERLA_ASSERT_EQUAL(zSize, f->zSize());
#endif

   for( auto i = f->beginGhostLayerOnly( nrOfGhostLayers, dir ); i != f->end(); ++i )
      buffer >> *i;
}



template<typename GhostLayerField_T>
void PackInfo<GhostLayerField_T>::communicateLocal(const IBlock * sender, IBlock * receiver, stencil::Direction dir)
{
   const GhostLayerField_T * sf = sender  ->getData< GhostLayerField_T >( bdId_ );
         GhostLayerField_T * rf = receiver->getData< GhostLayerField_T >( bdId_ );

   WALBERLA_ASSERT_EQUAL(sf->xSize(), rf->xSize());
   WALBERLA_ASSERT_EQUAL(sf->ySize(), rf->ySize());
   WALBERLA_ASSERT_EQUAL(sf->zSize(), rf->zSize());

   uint_t nrOfGhostLayers = numberOfGhostLayersToCommunicate( sf );
   auto srcIter = sf->beginSliceBeforeGhostLayer( dir, cell_idx_c( nrOfGhostLayers ) );
   auto dstIter = rf->beginGhostLayerOnly( nrOfGhostLayers, stencil::inverseDir[dir] );

   while( srcIter != sf->end() ) {
      *dstIter = *srcIter;
      ++srcIter;
      ++dstIter;
   }

   WALBERLA_ASSERT(srcIter == sf->end() && dstIter == rf->end());
}



template<typename GhostLayerField_T>
void PackInfo<GhostLayerField_T>::packDataImpl(const IBlock * sender, stencil::Direction dir, mpi::SendBuffer & outBuffer) const
{
   const GhostLayerField_T * f = sender->getData< GhostLayerField_T >( bdId_ );
   WALBERLA_ASSERT_NOT_NULLPTR(f);

#ifndef NDEBUG
   outBuffer << f->xSize() << f->ySize() << f->zSize();
#endif

   cell_idx_t nrOfGhostLayers = cell_idx_c( numberOfGhostLayersToCommunicate( f ) );
   for( auto i = f->beginSliceBeforeGhostLayer(dir, nrOfGhostLayers ); i != f->end(); ++i )
      outBuffer << *i;

}



template<typename GhostLayerField_T>
uint_t PackInfo<GhostLayerField_T>::numberOfGhostLayersToCommunicate( const GhostLayerField_T * const field ) const
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
