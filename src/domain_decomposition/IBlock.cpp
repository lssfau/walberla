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
//! \file IBlock.cpp
//! \ingroup domain_decomposition
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#include "BlockStorage.h"
#include "IBlock.h"


namespace walberla {
namespace domain_decomposition {



/// The following members are not used for checking if two IBlock objects are equal: storage_ & storageIndex_
bool IBlock::operator==( const IBlock& rhs ) const
{
   if( aabb_ != rhs.aabb_ || state_ != rhs.state_ || data_.size() != rhs.data_.size() )
      return false;

   for( uint_t i = 0; i != data_.size(); ++i )
      if( *(data_[i]) != *(rhs.data_[i]) )
         return false;

   return equal( &rhs );
}



//**********************************************************************************************************************
/*!
*   Every derived class must call this constructor and pass a reference to the governing block storage data structure as
*   well as a reference to an axis-aligned bounding box that defines the exact part of the simulation space that is
*   assigned to this block.
*   If the third parameter 'local' is set to false ('local' is set to true by default), the caller indicates that this
*   block is NOT a normal, local block; Normal, local blocks contain all the block data that is required by the
*   simulation. Setting 'local' to false might be useful in many situations, for example, if a temporary object of type
*   'IBlock' must be created, or if the governing block storage data structure also allocates remote blocks for whatever
*   reason (remote blocks are "empty shells": blocks that reside on other processes and hence store no local block data).
*/
//**********************************************************************************************************************
IBlock::IBlock( BlockStorage& storage, const AABB& aabb, const IBlockID::IDType& id, const bool local ) :

   aabb_( aabb ), storage_( storage ), storageIndex_( -1 ) {

   if( local )
      storage_.registerBlock( std::make_pair(id, this) );
}



IBlock::~IBlock()
{
   for( uint_t i = data_.size(); i-- != 0; )
      delete data_[i];

   storage_.removeBlock( this );
}



} // namespace domain_decomposition
} // namespace walberla
