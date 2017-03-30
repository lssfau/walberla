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
//! \file CellCounter.h
//! \ingroup field
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "FlagField.h"

#include "core/DataTypes.h"
#include "core/Set.h"
#include "core/debug/CheckFunctions.h"
#include "core/debug/Debug.h"
#include "core/mpi/Reduce.h"
#include "core/uid/SUID.h"

#include "domain_decomposition/StructuredBlockStorage.h"

#include <algorithm>
#include <numeric>



namespace walberla {
namespace field {



//**********************************************************************************************************************
/*!
*   \brief Class/functor for counting cells of a specific type that are stored in parallel on all processes in a
*          structured block storage
*
*   Returns the total number of cells of a specific type that are stored in a 'StructuredBlockStorage'. Can also return
*   the number of cells that are contained only on a specific level.
*   For counting cells of a specific type, a flag field is used: Only cells marked with certain flags (these flags have
*   to be specified in the constructor) are counted. By default, if at least on flag is set for a certain cell, this
*   cell is counted. If you only want to count cells for which _all_ specified flags are set, call member function
*   'allFlagsMustBeSet(true)'. After calling this function, only cells for which all flags are set are counted.
*   If no flags are provided (= an empty set is passed to the constructor -> this is the default behavior), all cells
*   are counted - independent of their specific type.
*   Additionally, you can provide an axis-aligned bounding box via member function 'addAABBFilter'. If an AABB is
*   added as a filter, only cells intersected by this AABB are considered during counting.
*/
//**********************************************************************************************************************

template< typename FlagField_T >
class CellCounter
{
public:

   CellCounter( const weak_ptr< StructuredBlockStorage > & blocks, const ConstBlockDataID & flagFieldId,
                const Set< FlagUID > & cellsToCount = Set< FlagUID >::emptySet(),
                const Set<SUID> & requiredSelectors     = Set<SUID>::emptySet(),
                const Set<SUID> & incompatibleSelectors = Set<SUID>::emptySet() ) :
      totalNumberOfCells_( uint64_t(0) ), totalNumberOfBlocksContainingCell_( uint64_t(0) ),
      blocks_( blocks ), allFlagsMustBeSet_( false ), cellsToCount_( cellsToCount ), flagFieldId_( flagFieldId ),
      requiredSelectors_( requiredSelectors ), incompatibleSelectors_( incompatibleSelectors )
   {}

   void allFlagsMustBeSet( const bool b )
   {
      allFlagsMustBeSet_ = b;
   }

   bool allFlagsMustBeSet() const
   {
      return allFlagsMustBeSet_;
   }

   void addAABBFilter( const math::AABB & aabb )
   {
      aabb_ = aabb;
   }

   void clearAABBFilter()
   {
      aabb_.init();
   }

   uint64_t numberOfCells() const
   {
      return totalNumberOfCells_;
   }

   uint64_t numberOfCells( const uint_t level ) const
   {
      WALBERLA_ASSERT_LESS( level, numberOfCells_.size() );
      return numberOfCells_[ level ];
   }

   const std::vector< uint64_t > & numberOfCellsPerLevel() const
   {
      return numberOfCells_;
   }


   uint64_t numberOfBlocksContainingCell() const
   {
      return totalNumberOfBlocksContainingCell_;
   }

   uint64_t numberOfBlocksContainingCell( const uint_t level ) const
   {
      WALBERLA_ASSERT_LESS( level, numberOfBlocksContainingCell_.size() );
      return numberOfBlocksContainingCell_[ level ];
   }

   const std::vector< uint64_t > & numberOfBlocksContainingCellPerLevel() const
   {
      return numberOfBlocksContainingCell_;
   }


   void operator()();

private:

   uint64_t                totalNumberOfCells_;
   std::vector< uint64_t > numberOfCells_;

   uint64_t                totalNumberOfBlocksContainingCell_;
   std::vector< uint64_t > numberOfBlocksContainingCell_;

   weak_ptr< StructuredBlockStorage > blocks_;

   bool allFlagsMustBeSet_;
   const Set< FlagUID >   cellsToCount_;
   const ConstBlockDataID flagFieldId_;

   math::AABB aabb_;

   Set<SUID> requiredSelectors_;
   Set<SUID> incompatibleSelectors_;
};



template< typename FlagField_T >
void CellCounter< FlagField_T >::operator()()
{
   totalNumberOfCells_ = uint64_t(0);
   totalNumberOfBlocksContainingCell_ = uint64_t(0);

   auto blocks = blocks_.lock();
   WALBERLA_CHECK_NOT_NULLPTR( blocks, "Trying to access 'CellCounter' for a block storage object that doesn't exist anymore" );

   numberOfCells_.assign( blocks->getNumberOfLevels(), uint64_t(0) );
   numberOfBlocksContainingCell_.assign( blocks->getNumberOfLevels(), uint64_t(0) );

   for( auto block = blocks->begin( requiredSelectors_, incompatibleSelectors_ ); block != blocks->end(); ++block )
   {
      const FlagField_T * flagField = block->template getData< FlagField_T >( flagFieldId_ );

      const auto level = blocks->getLevel( *block );
      WALBERLA_ASSERT_LESS( level, blocks->getNumberOfLevels() );

      cell::CellInterval interval( flagField->xyzSize() );
      if( !aabb_.empty() )
      {
         math::AABB blockAABB = block->getAABB();
         blockAABB.intersect( aabb_ );
         if( !blockAABB.empty() )
         {
            cell::CellInterval cellInterval = blocks->getCellBBFromAABB( blockAABB, level );
            blocks->transformGlobalToBlockLocalCellInterval( cellInterval, *block );
            interval.intersect( cellInterval );
         }
         else
         {
            interval = CellInterval();
         }
      }

      if( !interval.empty() )
      {
         if( cellsToCount_.empty() )
         {
            numberOfCells_[ level ] += uint64_c( interval.numCells() );
            ++numberOfBlocksContainingCell_[ level ];
         }
         else
         {
            typename FlagField_T::flag_t mask = 0;
            for( auto flag = cellsToCount_.begin(); flag != cellsToCount_.end(); ++flag )
               mask = static_cast< typename FlagField_T::flag_t >( mask | flagField->getFlag( *flag ) );

            bool blockHasCell = false;
            if( allFlagsMustBeSet_ )
            {
               for( auto cell = flagField->beginSliceXYZ( interval ); cell != flagField->end(); ++cell )
               {
                  if( isMaskSet( cell, mask ) )
                  {
                     ++numberOfCells_[ level ];
                     blockHasCell = true;
                  }
               }
            }
            else
            {
               for( auto cell = flagField->beginSliceXYZ( interval ); cell != flagField->end(); ++cell )
               {
                  if( isPartOfMaskSet( cell, mask ) )
                  {
                     ++numberOfCells_[ level ];
                     blockHasCell = true;
                  }
               }
            }

            if( blockHasCell )
            {
               ++numberOfBlocksContainingCell_[ level ];
            }
         }
      }
   }

   mpi::allReduceInplace( numberOfCells_,                mpi::SUM );
   mpi::allReduceInplace( numberOfBlocksContainingCell_, mpi::SUM );

   totalNumberOfCells_                = std::accumulate( numberOfCells_.begin(), numberOfCells_.end(), uint64_t(0) );
   totalNumberOfBlocksContainingCell_ = std::accumulate( numberOfBlocksContainingCell_.begin(),
                                                         numberOfBlocksContainingCell_.end(), uint64_t(0) );
}



} // namespace field
} // namespace walberla
