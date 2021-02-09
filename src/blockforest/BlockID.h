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
//! \file BlockID.h
//! \ingroup blockforest
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "Types.h"
#include "Utility.h"

#include "core/debug/Debug.h"
#include "core/math/Uint.h"
#include "core/mpi/RecvBuffer.h"
#include "core/mpi/SendBuffer.h"

#include "domain_decomposition/IBlockID.h"

#include <limits>
#include <ostream>
#include <vector>


// include cmake definition (-> WALBERLA_BLOCKFOREST_PRIMITIVE_BLOCKID)
#include "blockforest/CMakeDefs.h"


namespace walberla {
namespace blockforest {



#ifndef WALBERLA_BLOCKFOREST_PRIMITIVE_BLOCKID



class BlockID : public IBlockID {

public:

          BlockID() : usedBits_( 0 ) {}
   inline BlockID( const BlockID& id ) : usedBits_( id.usedBits_ ), blocks_( id.blocks_ ) {}
   inline BlockID( const uint_t treeIndex, const uint_t treeIdMarker );
   inline BlockID( const BlockID& id, const uint_t branchId );
          BlockID( const std::vector< uint8_t >& array, const uint_t offset, const uint_t bytes );

   inline void init( const uint_t treeIndex, const uint_t treeIdMarker );
   inline void clear() { usedBits_ = 0; blocks_.clear(); }

   inline uint_t getUsedBits()  const { return usedBits_; }
   inline uint_t getUsedBytes() const { return ( usedBits_ >> uint_c(3) ) + ( ( usedBits_ & uint_c(7) ) ? uint_c(1) : uint_c(0) ); }
   inline uint_t getTreeId()    const;
   inline uint_t getTreeIndex() const;

   BlockID getSuperId() const { BlockID id( *this ); id.removeBranchId(); return id; }
   BlockID getFatherId() const { return getSuperId(); }

   void   appendBranchId( const uint_t branchId );
   void   removeBranchId();
   uint_t    getBranchId() const { WALBERLA_ASSERT( !blocks_.empty() ); return blocks_[0] & uint_c(7); }

   inline bool operator< ( const IBlockID& rhs ) const;
   inline bool operator==( const IBlockID& rhs ) const;
   inline bool operator!=( const IBlockID& rhs ) const { return !(operator==( rhs )); }

   inline IDType getID() const;

   inline std::ostream& toStream( std::ostream& os ) const;

   void toByteArray( std::vector< uint8_t >& array, const uint_t offset, const uint_t bytes ) const;

   template< typename Buffer_T > void   toBuffer( Buffer_T& buffer ) const;
   template< typename Buffer_T > void fromBuffer( Buffer_T& buffer );

private:

   uint_t usedBits_;
   std::vector< uint_t >  blocks_;

   static const uint_t SHIFT = math::UINT_BITS - 3;

   WALBERLA_STATIC_ASSERT( math::UINT_BITS > 31 );
   WALBERLA_STATIC_ASSERT( !(math::UINT_BITS & (math::UINT_BITS - 1)) ); // power of two

}; // class BlockID



inline BlockID::BlockID( const uint_t treeIndex, const uint_t treeIdMarker )
{
   WALBERLA_ASSERT_GREATER( uintMSBPosition(treeIdMarker), uintMSBPosition(treeIndex) );

   const uint_t value = treeIndex | treeIdMarker;

   usedBits_ = uintMSBPosition( value );
   blocks_.push_back( value );
}



inline BlockID::BlockID( const BlockID& id, const uint_t branchId )
{
   WALBERLA_ASSERT_LESS( branchId, 8 );

   usedBits_ = id.usedBits_;
   blocks_   = id.blocks_;

   appendBranchId( branchId );
}



inline void BlockID::init( const uint_t treeIndex, const uint_t treeIdMarker )
{
   WALBERLA_ASSERT_EQUAL( usedBits_, 0 );
   WALBERLA_ASSERT( blocks_.empty() );
   WALBERLA_ASSERT_GREATER( uintMSBPosition(treeIdMarker), uintMSBPosition(treeIndex) );

   const uint_t value = treeIndex | treeIdMarker;

   usedBits_ = uintMSBPosition( value );
   blocks_.push_back( value );
}



inline uint_t BlockID::getTreeId() const
{
   WALBERLA_ASSERT_GREATER( usedBits_, 0 );
   WALBERLA_ASSERT_LESS_EQUAL( usedBits_, math::UINT_BITS );
   WALBERLA_CHECK_EQUAL( blocks_.size(), 1 );
   WALBERLA_ASSERT_EQUAL( usedBits_, uintMSBPosition( blocks_[0] ) );

   return blocks_[0];
}



inline uint_t BlockID::getTreeIndex() const
{
   WALBERLA_ASSERT_GREATER( usedBits_, 0 );
   WALBERLA_ASSERT_LESS_EQUAL( usedBits_, math::UINT_BITS );
   WALBERLA_CHECK_EQUAL( blocks_.size(), 1 );
   WALBERLA_ASSERT_EQUAL( usedBits_, uintMSBPosition( blocks_[0] ) );

   return blocks_[0] & ( uintPow2(usedBits_ - 1) - 1 );
}



inline bool BlockID::operator<( const IBlockID& rhs ) const
{
   WALBERLA_ASSERT_EQUAL( dynamic_cast< const BlockID* >( &rhs ), &rhs );
   const BlockID& _rhs = *static_cast< const BlockID* >( &rhs );

   if( blocks_.size() == _rhs.blocks_.size() ) {

      for( uint_t i = blocks_.size(); i-- != 0; )
         if( blocks_[i] != _rhs.blocks_[i] )
            return blocks_[i] < _rhs.blocks_[i];

      return false;
   }

   return blocks_.size() < _rhs.blocks_.size();
}



inline bool BlockID::operator==( const IBlockID& rhs ) const
{
   WALBERLA_ASSERT_EQUAL( dynamic_cast< const BlockID* >( &rhs ), &rhs );
   const BlockID& _rhs = *static_cast< const BlockID* >( &rhs );

   return ( usedBits_ == _rhs.usedBits_ && blocks_ == _rhs.blocks_ );
}



inline IBlockID::IDType BlockID::getID() const
{
   WALBERLA_CHECK_EQUAL( blocks_.size(), 1 );
   return numeric_cast<IBlockID::IDType> (blocks_[0]);
}



inline std::ostream& BlockID::toStream( std::ostream& os ) const {

   for( uint_t i = blocks_.size(); i-- != 0; ) {
      os << uintToBitString( blocks_[i] );
      if( i > 0 ) os << " ";
   }
   return os;
}



template< typename Buffer_T >
void BlockID::toBuffer( Buffer_T& buffer ) const
{
   WALBERLA_ASSERT_LESS_EQUAL( getUsedBytes(), 255 );

   const uint8_t bytes = uint8_t( getUsedBytes() );
   buffer << bytes;

   uint8_t b(0);
   for( uint_t i = 0; i != blocks_.size(); ++i )
      for( uint_t j = 0; j != math::UINT_BYTES && b != bytes; ++j, ++b )
         buffer << uint8_c( ( blocks_[i] >> ( uint_c(j) * uint_c(8) ) ) & uint_c(255) );
}



template< typename Buffer_T >
void BlockID::fromBuffer( Buffer_T& buffer ) {

   uint8_t bytes(0);
   buffer >> bytes;

   const uint_t blocks = ( bytes / math::UINT_BYTES ) + ( ( (bytes & ( math::UINT_BYTES - 1 )) == 0 ) ? 0 : 1 );

   uint8_t b(0);
   blocks_.clear();
   for( uint_t i = 0; i != blocks; ++i ) {
      blocks_.push_back(0);
      for( uint_t j = 0; j != math::UINT_BYTES && b != bytes; ++j, ++b ) {
         uint8_t byte(0);
         buffer >> byte;
         blocks_.back() |= uint_c( byte ) << ( uint_c(j) * uint_c(8) );
      }
   }

   for( uint_t i = 0; i != math::UINT_BITS; ++i )
      if( ( blocks_[ blocks-1 ] & ( uint_c(1) << uint_c(i) ) ) != uint_c(0) )
         usedBits_ = i + 1 + (blocks-1) * math::UINT_BITS;
}



#else



class BlockID : public IBlockID {

public:

   inline BlockID() = default;
   inline BlockID( const BlockID& id ) = default;
   inline BlockID( const uint_t id ) : id_( id ) {}
   inline BlockID( const uint_t treeIndex, const uint_t treeIdMarker );
   inline BlockID( const BlockID& id, const uint_t branchId );
   inline BlockID( const std::vector< uint8_t >& array, const uint_t offset, const uint_t bytes ) : id_( byteArrayToUint( array, offset, bytes ) ) {}

   inline void init( const uint_t treeIndex, const uint_t treeIdMarker );
          void clear() { id_ = uint_c(0); }

   uint_t getUsedBits()  const { return uintMSBPosition( id_ ); }
   uint_t getUsedBytes() const { const uint_t bits( getUsedBits() ); return ( bits >> uint_c(3) ) + ( ( bits & uint_c(7) ) ? uint_c(1) : uint_c(0) ); }
   uint_t getTreeId()    const { return id_; }
   uint_t getTreeIndex() const { return id_ - uintPow2( getUsedBits() - 1 ); }

   BlockID getSuperId()  const { WALBERLA_ASSERT_GREATER_EQUAL( getUsedBits(), uint_c(4) ); return BlockID( id_ >> 3 ); }
   BlockID getFatherId() const { return getSuperId(); }

   void   appendBranchId( const uint_t branchId ) { WALBERLA_ASSERT_LESS_EQUAL( getUsedBits() + 3, math::UINT_BITS ); WALBERLA_ASSERT_LESS( branchId, 8 ); id_ = (id_ << 3) + branchId; }
   void   removeBranchId() { WALBERLA_ASSERT_GREATER_EQUAL( getUsedBits(), uint_c(4) ); id_ >>= 3; }
   uint_t    getBranchId() const { WALBERLA_ASSERT_GREATER_EQUAL( getUsedBits(), uint_c(4) ); return id_ & uint_c(7); }

   bool operator< ( const IBlockID& rhs ) const override
      { WALBERLA_ASSERT_EQUAL( dynamic_cast< const BlockID* >( &rhs ), &rhs ); return id_ <  static_cast< const BlockID* >( &rhs )->id_; }
   bool operator> ( const IBlockID& rhs ) const
      { WALBERLA_ASSERT_EQUAL( dynamic_cast< const BlockID* >( &rhs ), &rhs ); return id_ >  static_cast< const BlockID* >( &rhs )->id_; }
   bool operator==( const IBlockID& rhs ) const override
      { WALBERLA_ASSERT_EQUAL( dynamic_cast< const BlockID* >( &rhs ), &rhs ); return id_ == static_cast< const BlockID* >( &rhs )->id_; }
   bool operator!=( const IBlockID& rhs ) const override
      { WALBERLA_ASSERT_EQUAL( dynamic_cast< const BlockID* >( &rhs ), &rhs ); return id_ != static_cast< const BlockID* >( &rhs )->id_; }

   inline IDType getID() const override;

   inline std::ostream& toStream( std::ostream& os ) const override;

   void toByteArray( std::vector< uint8_t >& array, const uint_t offset, const uint_t bytes ) const { uintToByteArray( id_, array, offset, bytes ); }

   template< typename Buffer_T > void   toBuffer( Buffer_T& buffer ) const;
   template< typename Buffer_T > void fromBuffer( Buffer_T& buffer );

private:

   uint_t id_ = uint_c(0);

}; // class BlockID



inline BlockID::BlockID( const uint_t treeIndex, const uint_t treeIdMarker ) : id_( treeIndex | treeIdMarker )
{
   WALBERLA_ASSERT_GREATER( uintMSBPosition(treeIdMarker), uintMSBPosition(treeIndex) );
}



inline BlockID::BlockID( const BlockID& id, const uint_t branchId ) : id_( (id.id_ << 3) + branchId )
{
   WALBERLA_ASSERT_LESS_EQUAL( id.getUsedBits() + 3, math::UINT_BITS );
   WALBERLA_ASSERT_LESS( branchId, 8 );
}



inline void BlockID::init( const uint_t treeIndex, const uint_t treeIdMarker )
{
   WALBERLA_ASSERT_GREATER( uintMSBPosition(treeIdMarker), uintMSBPosition(treeIndex) );

   id_ = treeIndex | treeIdMarker;
}




inline IBlockID::IDType BlockID::getID() const
{
   return numeric_cast<IBlockID::IDType> (id_);
}




inline std::ostream& BlockID::toStream( std::ostream& os ) const {

   os << uintToBitString( id_ ) << " (" << id_ << ")";
   return os;
}



template< typename Buffer_T >
void BlockID::toBuffer( Buffer_T& buffer ) const {

   const uint8_t bytes = uint8_t( getUsedBytes() );
   buffer << bytes;
   for( uint8_t i = 0; i != bytes; ++i )
      buffer << uint8_c( ( id_ >> ( uint_c(i) * uint_c(8) ) ) & uint_c(255) );
}



template< typename Buffer_T >
void BlockID::fromBuffer( Buffer_T& buffer ) {

   uint8_t bytes(0);
   buffer >> bytes;

   WALBERLA_ASSERT_LESS_EQUAL( bytes, math::UINT_BYTES );

   id_ = uint_c(0);
   for( uint8_t i = 0; i != bytes; ++i ) {
      uint8_t byte(0);
      buffer >> byte;
      id_ |= uint_c( byte ) << ( uint_c(i) * uint_c(8) );
   }
}



#endif



} // namespace blockforest

typedef blockforest::BlockID BlockID;

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
inline mpi::GenericSendBuffer<T,G> & operator<<( mpi::GenericSendBuffer<T,G> & buffer, const blockforest::BlockID & id )
{
   buffer.addDebugMarker( "bi" );
   id.toBuffer( buffer );
   return buffer;
}

template< typename T > // Element type  of RecvBuffer
inline mpi::GenericRecvBuffer<T>& operator>>( mpi::GenericRecvBuffer<T> & buffer, blockforest::BlockID & id )
{
   buffer.readDebugMarker( "bi" );
   id.fromBuffer( buffer );
   return buffer;
}

template<>
struct BufferSizeTrait< blockforest::BlockID > { static const bool constantSize = false; };

}
}
