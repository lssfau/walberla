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
//! \file PdfFieldMPIDatatypeInfo.h
//! \ingroup lbm
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/debug/Debug.h"
#include "communication/UniformMPIDatatypeInfo.h"
#include "field/communication/MPIDatatypes.h"

#include <set>

namespace walberla {
namespace lbm {
namespace communication {

template<typename PdfField_T>
class PdfFieldMPIDatatypeInfo : public walberla::communication::UniformMPIDatatypeInfo
{
public:
   PdfFieldMPIDatatypeInfo( BlockDataID pdfFieldID ) : pdfFieldID_( pdfFieldID ) {}

   virtual ~PdfFieldMPIDatatypeInfo() {}

   virtual shared_ptr<mpi::Datatype> getSendDatatype ( IBlock * block, const stencil::Direction dir )
   {
      return make_shared<mpi::Datatype>( field::communication::mpiDatatypeSliceBeforeGhostlayerXYZ(
         *getField( block ), dir, uint_t( 1 ), getOptimizedCommunicationIndices( dir ), false ) );
   }

   virtual shared_ptr<mpi::Datatype> getRecvDatatype ( IBlock * block, const stencil::Direction dir )
   {
      return make_shared<mpi::Datatype>( field::communication::mpiDatatypeGhostLayerOnlyXYZ(
         *getField( block ), dir, false, getOptimizedCommunicationIndices( stencil::inverseDir[dir] ) ) );
   }

   virtual void * getSendPointer( IBlock * block, const stencil::Direction )
   {
      return getField(block)->data();
   }

   virtual void * getRecvPointer( IBlock * block, const stencil::Direction )
   {
      return getField(block)->data();
   }

private:

   inline static std::set< cell_idx_t > getOptimizedCommunicationIndices( const stencil::Direction dir )
   {
      std::set< cell_idx_t > result;
      for( uint_t i = 0; i < PdfField_T::Stencil::d_per_d_length[dir]; ++i )
      {
         result.insert( cell_idx_c( PdfField_T::Stencil::idx[PdfField_T::Stencil::d_per_d[dir][i]] ) );
      }

      return result;
   }

   PdfField_T * getField( IBlock * block )
   {
      PdfField_T * const f = block->getData<PdfField_T>( pdfFieldID_ );
      WALBERLA_ASSERT_NOT_NULLPTR( f );
      return f;
   }

   BlockDataID pdfFieldID_;
};


} // namespace communication
} // namespace lbm
} // namespace walberla


