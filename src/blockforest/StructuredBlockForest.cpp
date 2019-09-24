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
//! \file StructuredBlockForest.cpp
//! \ingroup blockforest
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#include "StructuredBlockForest.h"
#include "core/Abort.h"


namespace walberla {
namespace blockforest {



bool StructuredBlockForest::blockExists( const Cell& cell, const uint_t level ) const {

   real_t x;
   real_t y;
   real_t z;
   getCellCenter( x, y, z, cell, level );

   if( getBlockInformation().active() ) {

      BlockID id;
      if( !getBlockInformation().getId( id, x, y, z ) )
         return false;

      return getLevelFromBlockId( id ) == level;
   }

   const IBlock* block = getBlock(x,y,z);
   if( block == nullptr )
      return false;

   WALBERLA_ASSERT_EQUAL( dynamic_cast< const Block* >( block ), block );

   return static_cast< const Block* >( block )->getLevel() == level;
}



void StructuredBlockForest::getBlockID( IBlockID& id, const Cell& cell, const uint_t level ) const {

   WALBERLA_ASSERT_EQUAL( dynamic_cast< BlockID* >( &id ), &id );

   real_t x;
   real_t y;
   real_t z;
   getCellCenter( x, y, z, cell, level );

   if( getBlockInformation().active() ) {
      if( !getBlockInformation().getId( *static_cast< BlockID* >( &id ), x, y, z ) ||
          getLevelFromBlockId( *static_cast< BlockID* >( &id ) ) != level) {
         WALBERLA_ABORT( "Getting block ID failed: There exists no block at global cell " << cell << " on level " << level << "!" );
      }
   }
   else {
      const IBlock* const block = getBlock(x,y,z);
      if( block == nullptr ) {
         WALBERLA_ABORT( "Getting block ID failed: Locally, there exists no block at global cell " << cell << " on level " << level << "!\n"
                         "                         (for simulation global information you have to explicitly construct the block forest to "
                         "contain global knowledge)");
      }
      WALBERLA_ASSERT_EQUAL( dynamic_cast< const Block* >( block ), block );
      if( static_cast< const Block* >( block )->getLevel() != level ) {
         WALBERLA_ABORT( "Getting block ID failed: Locally, there exists no block at global cell " << cell << " on level " << level << "!\n"
                         "                         (for simulation global information you have to explicitly construct the block forest to "
                         "contain global knowledge)");
      }
      id = block->getId();
   }
}



} // namespace blockforest
} // namespace walberla
