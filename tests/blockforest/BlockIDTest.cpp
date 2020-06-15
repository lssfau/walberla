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
//! \file BlockIDTest.cpp
//! \ingroup blockforest
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#include "blockforest/BlockID.h"
#include "core/debug/TestSubsystem.h"
#include "core/math/Random.h"
#include "core/mpi/RecvBuffer.h"
#include "core/mpi/SendBuffer.h"

#include <cstdlib>
#include <cstring>
#include <map>
#include <vector>


namespace walberla {
namespace blockforest {



static void test() {

   uint_t bit  = 1;
   uint_t mask = 0;

   for( uint_t i = 0; i != math::UINT_BITS; ++i ) {

      for( uint_t j = 0; j !=10000; ++j ) {

         const uint_t treeIndex = walberla::math::intRandom<uint_t>() & mask;

#ifndef WALBERLA_BLOCKFOREST_PRIMITIVE_BLOCKID
         const uint_t branches = walberla::math::intRandom<uint8_t>();
#else
         const uint_t branches = walberla::math::intRandom( uint_t(0), uint_c( ( math::UINT_BITS - i - 1 ) / 3 ) );
#endif
         std::vector< uint_t > branch;
         for( uint_t b = 0; b != branches; ++b )
            branch.push_back( walberla::math::intRandom<uint_t>() & uint_c(7) );

         BlockID id( treeIndex, bit );
         for( uint_t b = 0; b != branches; ++b )
            id.appendBranchId( branch[b] );

         uint_t bits = branches * 3 + i + 1;
         uint_t mod  = bits % 8;
         std::vector< uint8_t > byteArray( ( bits / 8 ) + ( ( mod > 0 ) ? 1 : 0 ), 0 );

         id.toByteArray( byteArray, 0, byteArray.size() );
         BlockID oldId( id );

         walberla::mpi::GenericSendBuffer<uint8_t> sb;
         oldId.toBuffer( sb );
         WALBERLA_CHECK_EQUAL( (oldId.getUsedBytes() + oldId.getUsedBytes() * mpi::BUFFER_DEBUG_OVERHEAD) + (1 + mpi::BUFFER_DEBUG_OVERHEAD ), sb.size() );

         for( uint_t b = branches; b-- != 0; ) {
            WALBERLA_CHECK_EQUAL( id.getBranchId(), branch[b] );
            id.removeBranchId();
         }
         WALBERLA_CHECK_EQUAL( id.getTreeIndex(), treeIndex );

         BlockID newId( byteArray, 0, byteArray.size() );
         WALBERLA_CHECK_EQUAL( oldId, newId );

         for( uint_t b = branches; b-- != 0; ) {
            WALBERLA_CHECK_EQUAL( newId.getBranchId(), branch[b] );
            newId.removeBranchId();
         }
         WALBERLA_CHECK_EQUAL( newId.getTreeIndex(), treeIndex );

         walberla::mpi::GenericRecvBuffer<uint8_t> rb( sb );
         newId.fromBuffer( rb );
         WALBERLA_CHECK_EQUAL( oldId, newId );
      }

      std::vector< BlockID >      blockId;
      std::vector< uint_t >       value;
      std::map< BlockID, uint_t > idMap;

      for( uint_t j = 0; j != 1000; ++j ) {

         const uint_t treeIndex = walberla::math::intRandom<uint_t>() & mask;

#ifndef WALBERLA_BLOCKFOREST_PRIMITIVE_BLOCKID
         const uint_t branches = walberla::math::intRandom<uint8_t>();
#else
         const uint_t branches = walberla::math::intRandom( uint_t(0), uint_c( ( math::UINT_BITS - i - 1 ) / 3 ) );
#endif
         BlockID id( treeIndex, bit );
         for( uint_t b = 0; b != branches; ++b )
            id.appendBranchId( walberla::math::intRandom<uint_t>() & uint_c(7) );

         if( idMap.find( id ) == idMap.end() ) {
            blockId.push_back( id );
            value.push_back( walberla::math::intRandom<uint_t>() );
            idMap[id] = value.back();
         }
      }
      for( uint_t j = 0; j != blockId.size(); ++j )
         WALBERLA_CHECK_EQUAL( value[j], idMap[blockId[j]] );

      bit  <<= 1;
      mask <<= 1;
      mask |= 1;
   }
}



} // namespace blockforest
} // namespace walberla



int main( int /*argc*/, char** /*argv*/ ) {

   walberla::debug::enterTestMode();

   walberla::blockforest::test();

   return EXIT_SUCCESS;
}
