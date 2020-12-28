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
//! \file PhantomBlock.h
//! \ingroup blockforest
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "BlockID.h"
#include "BlockNeighborhoodSection.h"
#include "Types.h"

#include "core/debug/Debug.h"
#include "core/math/AABB.h"
#include "core/uid/SUID.h"
#include "core/Any.h"

#include <vector>



namespace walberla {
namespace blockforest {


class PhantomBlockForest;


class PhantomBlock
{
public:

   struct NeighborBlock
   {
      NeighborBlock( const PhantomBlockForest & phantomForest, const BlockID & id, const uint_t process, const Set<SUID> & state = Set<SUID>::emptySet() );

      bool operator==( const NeighborBlock & rhs ) const { return ( id_ == rhs.id_ && process_ == rhs.process_ && state_ == rhs.state_ ); }
      bool operator!=( const NeighborBlock& rhs ) const { return !operator==( rhs ); }

      const BlockID &   getId()      const { return id_; }
            uint_t      getProcess() const { return process_; }
            void        setProcess( const uint_t process ) { process_ = process; }
      const Set<SUID> & getState()   const { return state_; }
      const AABB &      getAABB()    const { return aabb_; }

      BlockID   id_;
      uint_t    process_;
      Set<SUID> state_;
      AABB      aabb_;
   };
   
   
   
   PhantomBlock( PhantomBlockForest & phantomForest, const BlockID & id, const Set<SUID> & state, const AABB & aabb, const uint_t level,
                 const uint_t sourceLevel, const std::vector< uint_t > sourceProcess, const uint_t targetProcess )
      : phantomForest_( phantomForest ), id_( id ), state_( state ), aabb_( aabb ), level_( level ),
        sourceLevel_( sourceLevel ), sourceProcess_( sourceProcess ), targetProcess_( targetProcess )
   {}
   
   const PhantomBlockForest & getPhantomForest() const { return phantomForest_; }
         PhantomBlockForest & getPhantomForest()       { return phantomForest_; }

   const BlockID &   getId()    const { return id_; }
   const Set<SUID> & getState() const { return state_; }
   const AABB &      getAABB()  const { return aabb_; }
         uint_t      getProcess() const;
         uint_t      getLevel() const { return level_; }
   
   template< typename T >
   void addData( const T & data ) { data_ = data; }
   
   template< typename T >
   T getData() const { return walberla::any_cast<T>( data_ ); }
   
   bool hasData() const {
#ifndef WALBERLA_USE_STD_EXPERIMENTAL_ANY
      return data_.has_value();
#else
      return !(data_.empty());
#endif
   }
   
   inline       void                             clearNeighborhoodSection  ( const uint_t index );
   inline const std::vector< NeighborBlock * > & getNeighborhoodSection    ( const uint_t index ) const;
   inline       uint_t                           getNeighborhoodSectionSize( const uint_t index ) const;
   inline       void                             addNeighbor               ( const uint_t sectionIndex, const uint_t neighborIndex );
   inline const NeighborBlock &                  getNeighbor               ( const uint_t sectionIndex, const uint_t neighborIndex ) const;
   inline const BlockID &                        getNeighborId             ( const uint_t sectionIndex, const uint_t neighborIndex ) const;
   inline       uint_t                           getNeighborProcess        ( const uint_t sectionIndex, const uint_t neighborIndex ) const;
   inline       bool                             neighborExistsLocally     ( const uint_t sectionIndex, const uint_t neighborIndex ) const;
   inline       bool                             neighborExistsRemotely    ( const uint_t sectionIndex, const uint_t neighborIndex ) const;
   inline const Set<SUID> &                      getNeighborState          ( const uint_t sectionIndex, const uint_t neighborIndex ) const;
   inline const AABB &                           getNeighborAABB           ( const uint_t sectionIndex, const uint_t neighborIndex ) const;

   inline bool neighborhoodSectionHasBlocks           ( const uint_t sectionIndex ) const;
   inline bool neighborhoodSectionHasSmallerBlocks    ( const uint_t sectionIndex ) const;
   inline bool neighborhoodSectionHasEquallySizedBlock( const uint_t sectionIndex ) const;
   inline bool neighborhoodSectionHasLargerBlock      ( const uint_t sectionIndex ) const;

                void                           clearNeighborhood     () { neighborhood_.clear(); }
          const std::vector< NeighborBlock > & getNeighborhood       () const { return neighborhood_; }
                uint_t                         getNeighborhoodSize   () const { return neighborhood_.size(); }
   inline       void                           addNeighbor           ( const BlockID & id, const uint_t process, const Set<SUID> & state );
   inline const NeighborBlock &                getNeighbor           ( const uint_t index ) const;
   inline const BlockID &                      getNeighborId         ( const uint_t index ) const;
   inline       uint_t                         getNeighborProcess    ( const uint_t index ) const;
   inline       void                           setNeighborProcess    ( const uint_t index, const uint_t process );
   inline       bool                           neighborExistsLocally ( const uint_t index ) const;
                bool                           neighborExistsRemotely( const uint_t index ) const { return !neighborExistsLocally( index ); }
   inline const Set<SUID> &                    getNeighborState      ( const uint_t index ) const;
   inline const AABB &                         getNeighborAABB       ( const uint_t index ) const;

   bool sourceBlockHasTheSameSize() const { return sourceLevel_ == level_; }
   bool sourceBlockIsLarger() const { return sourceLevel_ < level_; }
   bool sourceBlockIsSmaller() const { return sourceLevel_ > level_; }

   uint_t getSourceLevel() const { return sourceLevel_; }
   
   const std::vector< uint_t > & getSourceProcess() const { return sourceProcess_; }
   
   uint_t getTargetProcess() const { return targetProcess_; }
   void   setTargetProcess( const uint_t p ) { targetProcess_= p; }

private:

   PhantomBlockForest & phantomForest_;

   BlockID id_;
   Set<SUID> state_;

   // can both be reconstructed from 'id_'
   AABB aabb_;
   uint_t level_;

   // set by the user/application via callback
   walberla::any data_;

   std::vector< NeighborBlock* > neighborhoodSection_[26]; // the 26 neighborhood sections (can be restored from 'neighborhood_')
   std::vector< NeighborBlock  > neighborhood_;            // all neighbor blocks

   uint_t sourceLevel_; // | sourceLevel_ - level_ | == 1
   std::vector< uint_t > sourceProcess_; // PhantomBlock <-> Block Connection

   uint_t targetProcess_; // PhantomBlock process migration during dynamic load balancing
};



inline void PhantomBlock::clearNeighborhoodSection( const uint_t index )
{
   WALBERLA_ASSERT_LESS( index, 26 );

   neighborhoodSection_[index].clear();
}



inline const std::vector< PhantomBlock::NeighborBlock * > & PhantomBlock::getNeighborhoodSection( const uint_t index ) const
{
   WALBERLA_ASSERT_LESS( index, 26 );

   return neighborhoodSection_[index];
}



inline uint_t PhantomBlock::getNeighborhoodSectionSize( const uint_t index ) const
{
   WALBERLA_ASSERT_LESS( index, 26 );

   return neighborhoodSection_[index].size();
}



inline void PhantomBlock::addNeighbor( const uint_t sectionIndex, const uint_t neighborIndex )
{
   WALBERLA_ASSERT_LESS( sectionIndex, 26 );
   WALBERLA_ASSERT_LESS( neighborIndex, neighborhood_.size() );
   WALBERLA_ASSERT_LESS( neighborhoodSection_[sectionIndex].size(), getBlockMaxNeighborhoodSectionSize( sectionIndex ) );

   // ATTENTION: if "neighborhood_" is changed afterwards, this pointer might get invalidated!
   neighborhoodSection_[sectionIndex].push_back( &(neighborhood_[ neighborIndex ]) );
}



inline const PhantomBlock::NeighborBlock & PhantomBlock::getNeighbor( const uint_t sectionIndex, const uint_t neighborIndex ) const
{
   WALBERLA_ASSERT_LESS( sectionIndex, 26 );
   WALBERLA_ASSERT_LESS( neighborIndex, neighborhoodSection_[sectionIndex].size() );

   return *(neighborhoodSection_[sectionIndex][neighborIndex]);
}



inline const BlockID & PhantomBlock::getNeighborId( const uint_t sectionIndex, const uint_t neighborIndex ) const
{
   WALBERLA_ASSERT_LESS( sectionIndex, 26 );
   WALBERLA_ASSERT_LESS( neighborIndex, neighborhoodSection_[sectionIndex].size() );

   return neighborhoodSection_[sectionIndex][neighborIndex]->id_;
}



inline uint_t PhantomBlock::getNeighborProcess( const uint_t sectionIndex, const uint_t neighborIndex ) const
{
   WALBERLA_ASSERT_LESS( sectionIndex, 26 );
   WALBERLA_ASSERT_LESS( neighborIndex, neighborhoodSection_[sectionIndex].size() );

   return neighborhoodSection_[sectionIndex][neighborIndex]->process_;
}



inline bool PhantomBlock::neighborExistsLocally( const uint_t sectionIndex, const uint_t neighborIndex ) const
{
   WALBERLA_ASSERT_LESS( sectionIndex, 26 );
   WALBERLA_ASSERT_LESS( neighborIndex, neighborhoodSection_[sectionIndex].size() );

   return neighborhoodSection_[sectionIndex][neighborIndex]->process_ == getProcess();
}



inline bool PhantomBlock::neighborExistsRemotely( const uint_t sectionIndex, const uint_t neighborIndex ) const {

   return !neighborExistsLocally( sectionIndex, neighborIndex );
}



inline const Set<SUID> & PhantomBlock::getNeighborState( const uint_t sectionIndex, const uint_t neighborIndex ) const
{
   WALBERLA_ASSERT_LESS( sectionIndex, 26 );
   WALBERLA_ASSERT_LESS( neighborIndex, neighborhoodSection_[sectionIndex].size() );

   return neighborhoodSection_[sectionIndex][neighborIndex]->state_;
}



inline const AABB & PhantomBlock::getNeighborAABB( const uint_t sectionIndex, const uint_t neighborIndex ) const
{
   WALBERLA_ASSERT_LESS( sectionIndex, 26 );
   WALBERLA_ASSERT_LESS( neighborIndex, neighborhoodSection_[sectionIndex].size() );

   return neighborhoodSection_[sectionIndex][neighborIndex]->aabb_;
}



inline bool PhantomBlock::neighborhoodSectionHasBlocks( const uint_t sectionIndex ) const
{
   WALBERLA_ASSERT_LESS( sectionIndex, 26 );

   return !neighborhoodSection_[sectionIndex].empty();
}

inline bool PhantomBlock::neighborhoodSectionHasSmallerBlocks( const uint_t sectionIndex ) const
{
   WALBERLA_ASSERT_LESS( sectionIndex, 26 );

   return !neighborhoodSection_[sectionIndex].empty() && neighborhoodSection_[sectionIndex][0]->id_.getUsedBits() > id_.getUsedBits();
}

inline bool PhantomBlock::neighborhoodSectionHasEquallySizedBlock( const uint_t sectionIndex ) const
{
   WALBERLA_ASSERT_LESS( sectionIndex, 26 );

   return !neighborhoodSection_[sectionIndex].empty() && neighborhoodSection_[sectionIndex][0]->id_.getUsedBits() == id_.getUsedBits();
}

inline bool PhantomBlock::neighborhoodSectionHasLargerBlock( const uint_t sectionIndex ) const
{
   WALBERLA_ASSERT_LESS( sectionIndex, 26 );

   return !neighborhoodSection_[sectionIndex].empty() && neighborhoodSection_[sectionIndex][0]->id_.getUsedBits() < id_.getUsedBits();
}



inline void PhantomBlock::addNeighbor( const BlockID & id, const uint_t process, const Set<SUID> & state )
{
#ifndef NDEBUG
   for( uint_t i = 0; i != neighborhood_.size(); ++i )
      WALBERLA_ASSERT( neighborhood_[i].getId() < id || id < neighborhood_[i].getId() );
#endif

   neighborhood_.push_back( NeighborBlock( phantomForest_, id, process, state ) );
}



inline const PhantomBlock::NeighborBlock & PhantomBlock::getNeighbor( const uint_t index ) const
{
   WALBERLA_ASSERT_LESS( index, neighborhood_.size() );

   return neighborhood_[index];
}



inline const BlockID & PhantomBlock::getNeighborId( const uint_t index ) const
{
   WALBERLA_ASSERT_LESS( index, neighborhood_.size() );

   return neighborhood_[index].id_;
}



inline uint_t PhantomBlock::getNeighborProcess( const uint_t index ) const
{
   WALBERLA_ASSERT_LESS( index, neighborhood_.size() );

   return neighborhood_[index].process_;
}



inline void PhantomBlock::setNeighborProcess( const uint_t index, const uint_t process )
{
   WALBERLA_ASSERT_LESS( index, neighborhood_.size() );

   neighborhood_[index].process_ = process;
}



inline bool PhantomBlock::neighborExistsLocally( const uint_t index ) const
{
   WALBERLA_ASSERT_LESS( index, neighborhood_.size() );

   return neighborhood_[index].process_ == getProcess();
}



inline const Set<SUID> & PhantomBlock::getNeighborState( const uint_t index ) const
{
   WALBERLA_ASSERT_LESS( index, neighborhood_.size() );

   return neighborhood_[index].state_;
}



inline const AABB & PhantomBlock::getNeighborAABB( const uint_t index ) const
{
   WALBERLA_ASSERT_LESS( index, neighborhood_.size() );

   return neighborhood_[index].aabb_;
}



} // namespace blockforest

typedef blockforest::PhantomBlock PhantomBlock;

} // namespace walberla
