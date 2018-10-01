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
//! \file BlockDataID.h
//! \ingroup domain_decomposition
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/DataTypes.h"
#include "core/mpi/BufferSizeTrait.h"
#include "core/mpi/RecvBuffer.h"
#include "core/mpi/SendBuffer.h"

#include <type_traits>


namespace walberla {
namespace domain_decomposition {



class ConstBlockDataID; // forward declaration;



class BlockDataID
{
   friend class ConstBlockDataID;

public:

            BlockDataID()                          = default;
   explicit BlockDataID( const uint_t id )         : id_( id ) {}
            BlockDataID( const BlockDataID& bdid ) = default;

   void pack( mpi::SendBuffer & buffer ) const { buffer << id_; }
   void unpack( mpi::RecvBuffer & buffer ) { buffer >> id_; }

   BlockDataID& operator=( const BlockDataID& bdid ) = default;

   bool operator==( const BlockDataID& bdid ) const { return id_ == bdid.id_; }
   bool operator!=( const BlockDataID& bdid ) const { return id_ != bdid.id_; }
   bool operator< ( const BlockDataID& bdid ) const { return id_ <  bdid.id_; }

   operator uint_t() const { return id_; }

private:

   uint_t id_ = 0;

}; // class BlockDataID
static_assert( std::is_trivially_copyable<BlockDataID>::value, "BlockDataID has to be trivially copyable!");


class ConstBlockDataID
{
public:

            ConstBlockDataID()                               = default;
   explicit ConstBlockDataID( const uint_t id )              : id_( id ) {}
            ConstBlockDataID( const BlockDataID&      bdid ) : id_( bdid.id_ ) {}
            ConstBlockDataID( const ConstBlockDataID& bdid ) = default;

   void pack( mpi::SendBuffer & buffer ) const { buffer << id_; }
   void unpack( mpi::RecvBuffer & buffer ) { buffer >> id_; }

   ConstBlockDataID& operator=( const BlockDataID&      bdid ) { id_ = bdid.id_; return *this; }
   ConstBlockDataID& operator=( const ConstBlockDataID& bdid ) = default;

   bool operator==( const ConstBlockDataID& bdid ) const { return id_ == bdid.id_; }
   bool operator!=( const ConstBlockDataID& bdid ) const { return id_ != bdid.id_; }
   bool operator< ( const ConstBlockDataID& bdid ) const { return id_ <  bdid.id_; }

   operator uint_t() const { return id_; }

private:

   uint_t id_ = 0;

}; // class ConstBlockDataID
static_assert( std::is_trivially_copyable<ConstBlockDataID>::value, "ConstBlockDataID has to be trivially copyable!");


} // namespace domain_decomposition

using domain_decomposition::BlockDataID;
using domain_decomposition::ConstBlockDataID;

} // namespace walberla



//======================================================================================================================
//
//  Send/Recv Buffer Serialization Specialization
//
//======================================================================================================================

namespace walberla {
namespace mpi {

template< typename T,  // Element type of SendBuffer
          typename G > // Growth policy of SendBuffer
mpi::GenericSendBuffer<T,G> & operator<<( mpi::GenericSendBuffer<T,G> & buffer, const BlockDataID & id )
{
   buffer.addDebugMarker( "id" );
   id.pack( buffer );
   return buffer;
}

template< typename T > // Element type  of RecvBuffer
mpi::GenericRecvBuffer<T> & operator>>( mpi::GenericRecvBuffer<T> & buffer, BlockDataID & id )
{
   buffer.readDebugMarker( "id" );
   id.unpack( buffer );
   return buffer;
}

template<>
struct BufferSizeTrait< BlockDataID > {
   static const bool constantSize = true;
   static const uint_t size = BufferSizeTrait<uint_t>::size;
};

template< typename T,  // Element type of SendBuffer
          typename G > // Growth policy of SendBuffer
mpi::GenericSendBuffer<T,G> & operator<<( mpi::GenericSendBuffer<T,G> & buffer, const ConstBlockDataID & id )
{
   buffer.addDebugMarker( "id" );
   id.pack( buffer );
   return buffer;
}

template< typename T > // Element type  of RecvBuffer
mpi::GenericRecvBuffer<T> & operator>>( mpi::GenericRecvBuffer<T> & buffer, ConstBlockDataID & id )
{
   buffer.readDebugMarker( "id" );
   id.unpack( buffer );
   return buffer;
}

template<>
struct BufferSizeTrait< ConstBlockDataID > {
   static const bool constantSize = true;
   static const uint_t size = BufferSizeTrait<uint_t>::size;
};

}
}
