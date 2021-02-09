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
//! \file Block.h
//! \ingroup blockforest
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "BlockID.h"
#include "BlockNeighborhoodSection.h"
#include "BlockReconstruction.h"
#include "Types.h"

#include "core/debug/Debug.h"
#include "core/math/AABB.h"
#include "core/mpi/RecvBuffer.h"
#include "core/mpi/SendBuffer.h"

#include "domain_decomposition/IBlock.h"

#include <vector>


namespace walberla {
namespace blockforest {



class BlockForest;
class PhantomBlock;
class SetupBlock;



class Block : public IBlock {

public:

   struct NeighborBlock {

      NeighborBlock( const BlockForest & forest, const BlockID & id, const uint_t process, const Set<SUID> & state = Set<SUID>::emptySet() );

      bool operator==( const NeighborBlock & rhs ) const { return ( id_ == rhs.id_ && process_ == rhs.process_ && state_ == rhs.state_ ); }
      bool operator!=( const NeighborBlock & rhs ) const { return !operator==( rhs ); }

      const BlockID &   getId()      const { return id_; }
            uint_t      getProcess() const { return process_; }
      const Set<SUID> & getState()   const { return state_; }
      const AABB &      getAABB()    const { return aabb_; }

      BlockID   id_;
      uint_t    process_;
      Set<SUID> state_;
      AABB      aabb_;
   };



   Block( BlockForest & forest, const SetupBlock * const block );
   Block( BlockForest & forest, const BlockID & id, const AABB & aabb, const Set<SUID> & state, const uint_t level,
          const BlockReconstruction::NeighborhoodReconstruction< Block > & neighborhoodReconstruction,
          const std::vector< BlockReconstruction::NeighborhoodReconstructionBlock > & neighbors );
   Block( BlockForest & forest, const PhantomBlock & phantom );
   Block( BlockForest & forest, const BlockID & id, const AABB & aabb, const uint_t level, mpi::RecvBuffer & buffer,
          const std::function< uint_t ( const uint_t ) > & processMapping = std::function< uint_t ( const uint_t ) >() );

   ~Block() override = default;

   void toBuffer( mpi::SendBuffer & buffer ) const;

   const BlockForest & getForest() const { return forest_; }
         BlockForest & getForest()       { return forest_; }

   const BlockID & getId()      const override { return id_; }
         uint_t    getProcess() const;
         uint_t    getLevel()   const { return level_; }

   inline       void                           clearNeighborhoodSection  ( const uint_t index );
   inline const std::vector< NeighborBlock* >& getNeighborhoodSection    ( const uint_t index ) const;
   inline       uint_t                         getNeighborhoodSectionSize( const uint_t index ) const;
   inline       void                           addNeighbor               ( const uint_t sectionIndex, const uint_t neighborIndex );
   inline const NeighborBlock &                getNeighbor               ( const uint_t sectionIndex, const uint_t neighborIndex ) const;
   inline const BlockID &                      getNeighborId             ( const uint_t sectionIndex, const uint_t neighborIndex ) const;
   inline       uint_t                         getNeighborProcess        ( const uint_t sectionIndex, const uint_t neighborIndex ) const;
   inline       bool                           neighborExistsLocally     ( const uint_t sectionIndex, const uint_t neighborIndex ) const;
   inline       bool                           neighborExistsRemotely    ( const uint_t sectionIndex, const uint_t neighborIndex ) const;
   inline const Set<SUID> &                    getNeighborState          ( const uint_t sectionIndex, const uint_t neighborIndex ) const;
   inline const AABB &                         getNeighborAABB           ( const uint_t sectionIndex, const uint_t neighborIndex ) const;

   inline bool neighborhoodSectionHasBlocks           ( const uint_t sectionIndex ) const;
   inline bool neighborhoodSectionHasSmallerBlocks    ( const uint_t sectionIndex ) const;
   inline bool neighborhoodSectionHasEquallySizedBlock( const uint_t sectionIndex ) const;
   inline bool neighborhoodSectionHasLargerBlock      ( const uint_t sectionIndex ) const;

                void                          clearNeighborhood     () { neighborhood_.clear(); }
          const std::vector< NeighborBlock >& getNeighborhood       () const { return neighborhood_; }
                uint_t                        getNeighborhoodSize   () const { return neighborhood_.size(); }
   inline       void                          addNeighbor           ( const BlockID& id, const uint_t process, const Set<SUID>& state );
   inline const NeighborBlock &               getNeighbor           ( const uint_t index ) const;
   inline const BlockID &                     getNeighborId         ( const uint_t index ) const;
   inline       uint_t                        getNeighborProcess    ( const uint_t index ) const;
   inline       bool                          neighborExistsLocally ( const uint_t index ) const;
                bool                          neighborExistsRemotely( const uint_t index ) const { return !neighborExistsLocally( index ); }
   inline const Set<SUID> &                   getNeighborState      ( const uint_t index ) const;
   inline const AABB &                        getNeighborAABB       ( const uint_t index ) const;

   void resetNeighborhood( const PhantomBlock & phantom );

   bool targetBlockHasTheSameSize() const { return targetLevel_ == level_; }
   bool targetBlockIsLarger() const { return targetLevel_ < level_; }
   bool targetBlockIsSmaller() const { return targetLevel_ > level_; }

   inline uint_t getTargetLevel() const { return targetLevel_; }
   inline void   setTargetLevel( const uint_t tl );
   
   const std::vector< uint_t > & getTargetProcess() const { return targetProcess_; }
   void clearTargetProcess() { targetProcess_.clear(); }
   void addTargetProcess( const uint_t process ) { targetProcess_.push_back( process ); }
   void setTargetProcess( const uint_t index, const uint_t process ) { WALBERLA_ASSERT_GREATER( targetProcess_.size(), index ); targetProcess_[index] = process; }

protected:

   bool equal( const IBlock* rhs ) const override;

private:

   BlockForest & forest_;

   BlockID id_;
   uint_t  level_; // 0=coarse (= block is located on the initial grid) -> 1 -> 2 -> finer

   std::vector< NeighborBlock* >  neighborhoodSection_[26]; // the 26 neighborhood sections
   std::vector< NeighborBlock  >  neighborhood_;            // all neighbor blocks

   uint_t targetLevel_; // | level_ - targetLevel_ | <= 1
   std::vector< uint_t > targetProcess_;
};



inline void Block::clearNeighborhoodSection( const uint_t index )
{
   WALBERLA_ASSERT_LESS( index, 26 );

   neighborhoodSection_[index].clear();
}



inline const std::vector< Block::NeighborBlock* >& Block::getNeighborhoodSection( const uint_t index ) const
{
   WALBERLA_ASSERT_LESS( index, 26 );

   return neighborhoodSection_[index];
}



inline uint_t Block::getNeighborhoodSectionSize( const uint_t index ) const
{
   WALBERLA_ASSERT_LESS( index, 26 );

   return neighborhoodSection_[index].size();
}



inline void Block::addNeighbor( const uint_t sectionIndex, const uint_t neighborIndex )
{
   WALBERLA_ASSERT_LESS( sectionIndex, 26 );
   WALBERLA_ASSERT_LESS( neighborIndex, neighborhood_.size() );
   WALBERLA_ASSERT_LESS( neighborhoodSection_[sectionIndex].size(), getBlockMaxNeighborhoodSectionSize( sectionIndex ) );

   // ATTENTION: if "neighborhood_" is changed afterwards, this pointer might get invalidated!
   neighborhoodSection_[sectionIndex].push_back( &(neighborhood_[ neighborIndex ]) );
}



inline const Block::NeighborBlock& Block::getNeighbor( const uint_t sectionIndex, const uint_t neighborIndex ) const
{
   WALBERLA_ASSERT_LESS( sectionIndex, 26 );
   WALBERLA_ASSERT_LESS( neighborIndex, neighborhoodSection_[sectionIndex].size() );

   return *(neighborhoodSection_[sectionIndex][neighborIndex]);
}



inline const BlockID& Block::getNeighborId( const uint_t sectionIndex, const uint_t neighborIndex ) const
{
   WALBERLA_ASSERT_LESS( sectionIndex, 26 );
   WALBERLA_ASSERT_LESS( neighborIndex, neighborhoodSection_[sectionIndex].size() );

   return neighborhoodSection_[sectionIndex][neighborIndex]->id_;
}



inline uint_t Block::getNeighborProcess( const uint_t sectionIndex, const uint_t neighborIndex ) const
{
   WALBERLA_ASSERT_LESS( sectionIndex, 26 );
   WALBERLA_ASSERT_LESS( neighborIndex, neighborhoodSection_[sectionIndex].size() );

   return neighborhoodSection_[sectionIndex][neighborIndex]->process_;
}



inline bool Block::neighborExistsLocally( const uint_t sectionIndex, const uint_t neighborIndex ) const
{
   WALBERLA_ASSERT_LESS( sectionIndex, 26 );
   WALBERLA_ASSERT_LESS( neighborIndex, neighborhoodSection_[sectionIndex].size() );

   return neighborhoodSection_[sectionIndex][neighborIndex]->process_ == getProcess();
}



inline bool Block::neighborExistsRemotely( const uint_t sectionIndex, const uint_t neighborIndex ) const {

   return !neighborExistsLocally( sectionIndex, neighborIndex );
}



inline const Set<SUID> & Block::getNeighborState( const uint_t sectionIndex, const uint_t neighborIndex ) const
{
   WALBERLA_ASSERT_LESS( sectionIndex, 26 );
   WALBERLA_ASSERT_LESS( neighborIndex, neighborhoodSection_[sectionIndex].size() );

   return neighborhoodSection_[sectionIndex][neighborIndex]->state_;
}



inline const AABB & Block::getNeighborAABB( const uint_t sectionIndex, const uint_t neighborIndex ) const
{
   WALBERLA_ASSERT_LESS( sectionIndex, 26 );
   WALBERLA_ASSERT_LESS( neighborIndex, neighborhoodSection_[sectionIndex].size() );

   return neighborhoodSection_[sectionIndex][neighborIndex]->aabb_;
}



inline bool Block::neighborhoodSectionHasBlocks( const uint_t sectionIndex ) const
{
   WALBERLA_ASSERT_LESS( sectionIndex, 26 );

   return !neighborhoodSection_[sectionIndex].empty();
}

inline bool Block::neighborhoodSectionHasSmallerBlocks( const uint_t sectionIndex ) const
{
   WALBERLA_ASSERT_LESS( sectionIndex, 26 );

   return !neighborhoodSection_[sectionIndex].empty() && neighborhoodSection_[sectionIndex][0]->id_.getUsedBits() > id_.getUsedBits();
}

inline bool Block::neighborhoodSectionHasEquallySizedBlock( const uint_t sectionIndex ) const
{
   WALBERLA_ASSERT_LESS( sectionIndex, 26 );

   return !neighborhoodSection_[sectionIndex].empty() && neighborhoodSection_[sectionIndex][0]->id_.getUsedBits() == id_.getUsedBits();
}

inline bool Block::neighborhoodSectionHasLargerBlock( const uint_t sectionIndex ) const
{
   WALBERLA_ASSERT_LESS( sectionIndex, 26 );

   return !neighborhoodSection_[sectionIndex].empty() && neighborhoodSection_[sectionIndex][0]->id_.getUsedBits() < id_.getUsedBits();
}



inline void Block::addNeighbor( const BlockID & id, const uint_t process, const Set<SUID> & state )
{
#ifndef NDEBUG
   for( uint_t i = 0; i != neighborhood_.size(); ++i )
      WALBERLA_ASSERT( neighborhood_[i].getId() < id || id < neighborhood_[i].getId() );
#endif

   neighborhood_.emplace_back( forest_, id, process, state );
}



inline const Block::NeighborBlock& Block::getNeighbor( const uint_t index ) const
{
   WALBERLA_ASSERT_LESS( index, neighborhood_.size() );

   return neighborhood_[index];
}



inline const BlockID& Block::getNeighborId( const uint_t index ) const
{
   WALBERLA_ASSERT_LESS( index, neighborhood_.size() );

   return neighborhood_[index].id_;
}



inline uint_t Block::getNeighborProcess( const uint_t index ) const
{
   WALBERLA_ASSERT_LESS( index, neighborhood_.size() );

   return neighborhood_[index].process_;
}



inline bool Block::neighborExistsLocally ( const uint_t index ) const
{
   WALBERLA_ASSERT_LESS( index, neighborhood_.size() );

   return neighborhood_[index].process_ == getProcess();
}



inline const Set<SUID> & Block::getNeighborState( const uint_t index ) const
{
   WALBERLA_ASSERT_LESS( index, neighborhood_.size() );

   return neighborhood_[index].state_;
}



inline const AABB & Block::getNeighborAABB( const uint_t index ) const
{
   WALBERLA_ASSERT_LESS( index, neighborhood_.size() );

   return neighborhood_[index].aabb_;
}



inline void Block::setTargetLevel( const uint_t tl )
{
   WALBERLA_ASSERT( tl == level_ || tl == (level_ - uint_t(1)) || tl == (level_ + uint_t(1)) );

   targetLevel_ = tl;
}



} // namespace blockforest

typedef blockforest::Block Block;

} // namespace walberla


