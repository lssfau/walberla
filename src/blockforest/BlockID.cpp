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
//! \file BlockID.cpp
//! \ingroup blockforest
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#include "BlockID.h"


#ifndef NDEBUG
#include <map>
#endif

namespace walberla {
namespace blockforest {



#ifndef WALBERLA_BLOCKFOREST_PRIMITIVE_BLOCKID



BlockID::BlockID( const std::vector< uint8_t >& array, const uint_t offset, const uint_t bytes )
{
   WALBERLA_ASSERT( !array.empty() );
   WALBERLA_ASSERT_GREATER( bytes, uint_c(0) );
   WALBERLA_ASSERT_LESS_EQUAL( offset + bytes, array.size() );

   static const uint8_t mask[8] = { 1, 2, 4, 8, 16, 32, 64, 128 };

   usedBits_ = 0;

   for( uint_t i = offset + bytes; i-- != offset && usedBits_ == 0; )
      for( uint_t j = 8;  j-- != 0 && usedBits_ == 0; )
         if( array[i] & mask[j] ) usedBits_ = ( i - offset ) * 8 + j + 1;

   WALBERLA_ASSERT_UNEQUAL( usedBits_, 0 );

   const uint_t blocks = ( usedBits_ / math::UINT_BITS ) + ( ( (usedBits_ & ( math::UINT_BITS - 1 )) == 0 ) ? 0 : 1 );

   uint_t b = offset;
   for( uint_t i = 0; i != blocks; ++i ) {
      blocks_.push_back(0);
      for( uint_t j = 0; j != math::UINT_BYTES && b != offset + bytes; ++j, ++b ) {
         blocks_.back() |= uint_c( array[b] ) << ( j * 8 );
      }
   }
}



void BlockID::appendBranchId( const uint_t branchId )
{
   WALBERLA_ASSERT( !blocks_.empty() );
   WALBERLA_ASSERT_LESS( branchId, 8 );
   WALBERLA_ASSERT_EQUAL( usedBits_, uintMSBPosition( blocks_.back() ) + math::UINT_BITS * (blocks_.size()-1) );

   const uint_t unusedBits = math::UINT_BITS - ( usedBits_ & ( math::UINT_BITS - 1 ) );

   if( unusedBits < 3 || unusedBits == math::UINT_BITS )
      blocks_.push_back(0);

   for( uint_t i = static_cast< uint_t >( blocks_.size() ) - 1; i != 0; --i )
         blocks_[i] = ( blocks_[i] << 3 ) | ( blocks_[i-1] >> SHIFT );

   blocks_[0] = ( blocks_[0] << 3 ) | branchId;
   usedBits_ += 3;

   WALBERLA_ASSERT_EQUAL( usedBits_, uintMSBPosition( blocks_.back() ) + math::UINT_BITS * (blocks_.size()-1) );
   WALBERLA_ASSERT_GREATER_EQUAL( blocks_.size() * math::UINT_BITS, usedBits_ );
}



void BlockID::removeBranchId()
{
   WALBERLA_ASSERT( !blocks_.empty() );
   WALBERLA_ASSERT_GREATER( usedBits_, 3 );
   WALBERLA_ASSERT_EQUAL( usedBits_, uintMSBPosition( blocks_.back() ) + math::UINT_BITS * (blocks_.size()-1) );

   for( uint_t i = 0; i != static_cast< uint_t >( blocks_.size() ) - 1; ++i )
      blocks_[i] = ( blocks_[i] >> 3 ) | ( (blocks_[i+1] & uint_c(7)) << SHIFT );

   const uint_t bits = usedBits_ & ( math::UINT_BITS - 1 );

   if( 0 < bits && bits < 4 ) blocks_.pop_back();
   else blocks_.back() >>= 3;

   usedBits_ -= 3;

   WALBERLA_ASSERT_EQUAL( usedBits_, uintMSBPosition( blocks_.back() ) + math::UINT_BITS * (blocks_.size()-1) );
   WALBERLA_ASSERT( !blocks_.empty() );
}



void BlockID::toByteArray( std::vector< uint8_t >& array, const uint_t offset, const uint_t bytes ) const
{
   WALBERLA_ASSERT( !array.empty() );
   WALBERLA_ASSERT_LESS_EQUAL( offset + bytes, array.size() );
   WALBERLA_ASSERT_GREATER_EQUAL( bytes * uint_c(8), usedBits_ );

   for( uint_t i = offset; i != ( offset + bytes ); ++i )
      array[i] = uint8_c(0);

   uint_t b = offset;

   for( uint_t i = 0; i != blocks_.size(); ++i ) {

      uint_t block = blocks_[i];

      for( uint_t j = 0; j != math::UINT_BYTES && b != offset + bytes; ++j, ++b ) {

         array[b] = static_cast< uint8_t >( block & uint_c(255) );
         block >>= 8;
      }
   }
}


#else

#ifdef WALBERLA_CXX_COMPILER_IS_MSVC
namespace { char dummy; } // disable MSVC warning LNK4221: This object file does not define any previously
                          // undefined public symbols, so it will not be used by any link operation that
                          // consumes this library
#endif

#endif



} // namespace blockforest
} // namespace walberla
