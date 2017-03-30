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
//! \file BlockCounter.h
//! \ingroup domain_decomposition
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "StructuredBlockStorage.h"

#include "core/DataTypes.h"
#include "core/Set.h"
#include "core/debug/CheckFunctions.h"
#include "core/debug/Debug.h"
#include "core/mpi/Reduce.h"
#include "core/uid/SUID.h"

#include <algorithm>



namespace walberla {
namespace domain_decomposition {



//**********************************************************************************************************************
/*!
*   \brief Class/functor for counting the blocks stored in parallel on all processes in a structured block storage
*
*   Returns the total number of blocks stored in a 'StructuredBlockStorage'. Can also return the number of blocks that
*   are contained only on a specific level.
*/
//**********************************************************************************************************************

class BlockCounter
{
public:

   BlockCounter( const weak_ptr< StructuredBlockStorage > & blocks,
                 const Set<SUID> & requiredSelectors     = Set<SUID>::emptySet(),
                 const Set<SUID> & incompatibleSelectors = Set<SUID>::emptySet() ) :
      totalNumberOfBlocks_( uint_t(0) ), blocks_( blocks ),
      requiredSelectors_(requiredSelectors), incompatibleSelectors_( incompatibleSelectors ) {}

   uint_t numberOfBlocks() const
   {
      return totalNumberOfBlocks_;
   }

   uint_t numberOfBlocks( const uint_t level ) const
   {
      WALBERLA_ASSERT_LESS( level, numberOfBlocks_.size() );
      return numberOfBlocks_[ level ];
   }

   void operator()()
   {
      auto blocks = blocks_.lock();
      WALBERLA_CHECK_NOT_NULLPTR( blocks, "Trying to access 'BlockCounter' for a block storage object that doesn't exist anymore" );

      totalNumberOfBlocks_ = uint_t(0);
      numberOfBlocks_.assign( blocks->getNumberOfLevels(), uint_t(0) );

      for( auto block = blocks->begin( requiredSelectors_, incompatibleSelectors_ ); block != blocks->end(); ++block )
      {
         WALBERLA_ASSERT_LESS( blocks->getLevel( *block ), blocks->getNumberOfLevels() );
         numberOfBlocks_[ blocks->getLevel( *block ) ] += uint_t(1);
      }

      mpi::allReduceInplace( numberOfBlocks_, mpi::SUM );

      for( auto numberOfBlocksOnLevel = numberOfBlocks_.begin(); numberOfBlocksOnLevel != numberOfBlocks_.end(); ++numberOfBlocksOnLevel )
         totalNumberOfBlocks_ += *numberOfBlocksOnLevel;
   }

private:

   uint_t totalNumberOfBlocks_;
   std::vector< uint_t > numberOfBlocks_;

   weak_ptr< StructuredBlockStorage > blocks_;

   Set<SUID> requiredSelectors_;
   Set<SUID> incompatibleSelectors_;
};



} // namespace domain_decomposition
} // namespace walberla
