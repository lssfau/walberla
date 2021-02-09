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
//! \file UniformMPIDatatypeInfo.h
//! \ingroup field
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/debug/Debug.h"
#include "communication/UniformMPIDatatypeInfo.h"
#include "field/communication/MPIDatatypes.h"

namespace walberla {
namespace field {
namespace communication {

template<typename GhostLayerField_T>
class UniformMPIDatatypeInfo : public walberla::communication::UniformMPIDatatypeInfo
{
public:
   UniformMPIDatatypeInfo( BlockDataID blockDataID ) : blockDataID_ ( blockDataID ),
      communicateAllGhostLayers_( true ), numberOfGhostLayers_(0) {}

   UniformMPIDatatypeInfo( BlockDataID blockDataID, const uint_t numberOfGhostLayers ) : blockDataID_( blockDataID ),
      communicateAllGhostLayers_( false ), numberOfGhostLayers_( numberOfGhostLayers ) {}

   ~UniformMPIDatatypeInfo() override = default;

   shared_ptr<mpi::Datatype> getSendDatatype ( IBlock * block, const stencil::Direction dir ) override
   {
      auto numGl = numberOfGhostLayersToCommunicate( block );
      return make_shared<mpi::Datatype>( mpiDatatypeSliceBeforeGhostlayer( *getField( block ), dir, numGl ) );
   }

   shared_ptr<mpi::Datatype> getRecvDatatype ( IBlock * block, const stencil::Direction dir ) override
   {
      auto numGl = numberOfGhostLayersToCommunicate( block );
      return make_shared<mpi::Datatype>( mpiDatatypeGhostLayerOnly( *getField(block), numGl, dir ) );
   }

   void * getSendPointer( IBlock * block, const stencil::Direction ) override
   {
      return getField(block)->data();
   }

   void * getRecvPointer( IBlock * block, const stencil::Direction ) override
   {
      return getField(block)->data();
   }

private:

   GhostLayerField_T * getField( IBlock * block )
   {
      GhostLayerField_T * const f =  block->getData<GhostLayerField_T>( blockDataID_ );
      WALBERLA_ASSERT_NOT_NULLPTR( f );
      return f;
   }

   uint_t numberOfGhostLayersToCommunicate( IBlock * block )
   {
      if( communicateAllGhostLayers_ )
      {
         return getField( block )->nrOfGhostLayers();
      }
      else
      {
         WALBERLA_ASSERT_LESS_EQUAL( numberOfGhostLayers_, getField( block )->nrOfGhostLayers() );
         return numberOfGhostLayers_;
      }
   }

   BlockDataID blockDataID_;
   bool communicateAllGhostLayers_;
   uint_t numberOfGhostLayers_;
};


} // namespace communication
} // namespace field
} // namespace walberla


