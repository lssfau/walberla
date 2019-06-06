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
//! \author Martin Bauer <martin.bauer@fau.de>

//
//======================================================================================================================

#pragma once

#include "core/debug/Debug.h"
#include "communication/UniformMPIDatatypeInfo.h"
#include "field/communication/MPIDatatypes.h"

#include <set>

namespace walberla {
namespace field {
namespace communication {

template< typename GhostLayerField_T, typename Stencil_T >
class StencilRestrictedMPIDatatypeInfo : public walberla::communication::UniformMPIDatatypeInfo
{
public:
    StencilRestrictedMPIDatatypeInfo( BlockDataID fieldID ) : fieldID_( fieldID ) {}

    virtual ~StencilRestrictedMPIDatatypeInfo() {}

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
        for( uint_t i = 0; i < Stencil_T::d_per_d_length[dir]; ++i )
        {
            result.insert( cell_idx_c( Stencil_T::idx[Stencil_T::d_per_d[dir][i]] ) );
        }

        return result;
    }

    GhostLayerField_T * getField( IBlock * block )
    {
        GhostLayerField_T * const f = block->getData<GhostLayerField_T>( fieldID_ );
        WALBERLA_ASSERT_NOT_NULLPTR( f );
        return f;
    }

    BlockDataID fieldID_;
};


} // namespace communication
} // namespace field
} // namespace walberla


