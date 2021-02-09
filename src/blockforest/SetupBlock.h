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
//! \file SetupBlock.h
//! \ingroup blockforest
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "BlockID.h"
#include "BlockNeighborhoodSection.h"
#include "BlockReconstruction.h"
#include "Types.h"

#include "core/NonCopyable.h"
#include "core/debug/Debug.h"
#include "core/math/AABB.h"
#include "core/uid/SUID.h"
#include "stencil/Directions.h"

#include <vector>


namespace walberla {
namespace blockforest {



class SetupBlock : private NonCopyable {

public:

   inline SetupBlock( SetupBlock* const father, const BlockID& id,
                      const real_t xmin, const real_t ymin, const real_t zmin, // incl.
                      const real_t xmax, const real_t ymax, const real_t zmax, // excl.
                      const uint_t level );

   ~SetupBlock() { for( uint_t i = 0; i != children_.size(); ++i ) delete children_[i]; }

   const BlockID& getId()            const { return Id_; }
         uint_t   getProcess()       const { return process_; }
         uint_t   getTargetProcess() const { return getProcess(); }
         void     assignProcess      ( const uint_t process ) { process_ = process; }
         void     assignTargetProcess( const uint_t process ) { assignProcess( process ); }

   const Set<SUID>& getState() const { return state_; }
   void setState( const Set<SUID>& state ) { state_  = state; }
   void addState( const Set<SUID>& state ) { state_ += state; }
   void addState( const SUID&      state ) { state_ += state; }
   void clearState() { state_.clear(); }

   const AABB & getAABB()  const { return aabb_; }
         uint_t getLevel() const { return level_; }

   workload_t getWorkload() const { return workload_; }
   void       setWorkload( const workload_t w ) { WALBERLA_ASSERT_GREATER_EQUAL( w, static_cast< workload_t >(0) ); workload_ = w; }

   memory_t getMemory() const { return memory_; }
   void     setMemory( const memory_t m ) { WALBERLA_ASSERT_GREATER_EQUAL( m, static_cast< memory_t >(0) ); memory_ = m; }

   bool isMarked() const { return marker_; }
   void setMarker( const bool marker ) { marker_ = marker; }

   const SetupBlock* getFather() const { return father_; }
         SetupBlock* getFather()       { return father_; }
         void        setFather( SetupBlock* const father) { father_ = father; }

   inline const SetupBlock* getChild( const uint_t index ) const { WALBERLA_ASSERT_LESS( index, children_.size() ); return children_[index]; }
   inline       SetupBlock* getChild( const uint_t index )       { WALBERLA_ASSERT_LESS( index, children_.size() ); return children_[index]; }
   inline       void        setChild( const uint_t index, SetupBlock* const child );

   bool hasFather()   const { return father_ != nullptr; }
   bool hasChildren() const { return !children_.empty(); }

   const std::vector< SetupBlock* >& getNeighborhoodSection( const uint_t index ) const { WALBERLA_ASSERT_LESS( index, 26 ); return neighborhoodSection_[index]; }

   inline       uint_t      getNeighborhoodSectionSize( const uint_t index ) const { WALBERLA_ASSERT_LESS( index, 26 ); return neighborhoodSection_[index].size(); }
   inline       void        clearNeighborhoodSection  ( const uint_t index )       { WALBERLA_ASSERT_LESS( index, 26 ); neighborhoodSection_[index].clear(); }
   inline const SetupBlock* getNeighbor                       ( const uint_t sectionIndex, const uint_t neighborIndex ) const;
   inline       SetupBlock* getNeighbor                       ( const uint_t sectionIndex, const uint_t neighborIndex );
   inline       void        addNeighbor                       ( const uint_t sectionIndex, SetupBlock* const block );
   inline const BlockID&    getNeighborId                     ( const uint_t sectionIndex, const uint_t neighborIndex ) const;
   inline       uint_t      getNeighborProcess                ( const uint_t sectionIndex, const uint_t neighborIndex ) const;
   inline       uint_t      getNeighborTargetProcess          ( const uint_t sectionIndex, const uint_t neighborIndex ) const;
   inline const Set<SUID>&  getNeighborState                  ( const uint_t sectionIndex, const uint_t neighborIndex ) const;
   inline       bool        neighborIsLocatedOnTheSameProcess ( const uint_t sectionIndex, const uint_t neighborIndex ) const;
   inline       bool        neighborIsAssignedToTheSameProcess( const uint_t sectionIndex, const uint_t neighborIndex ) const;

   inline bool neighborhoodSectionHasBlocks           ( const uint_t sectionIndex ) const;
   inline bool neighborhoodSectionHasSmallerBlocks    ( const uint_t sectionIndex ) const;
   inline bool neighborhoodSectionHasEquallySizedBlock( const uint_t sectionIndex ) const;
   inline bool neighborhoodSectionHasLargerBlock      ( const uint_t sectionIndex ) const;

   inline const std::vector< SetupBlock* >& getNeighborhood()     const { return neighborhood_; }
   inline       uint_t                      getNeighborhoodSize() const { return neighborhood_.size(); }
   inline const SetupBlock*                 getNeighbor                       ( const uint_t index ) const;
   inline const BlockID&                    getNeighborId                     ( const uint_t index ) const;
   inline       uint_t                      getNeighborProcess                ( const uint_t index ) const;
   inline       uint_t                      getNeighborTargetProcess          ( const uint_t index ) const;
   inline const Set<SUID>&                  getNeighborState                  ( const uint_t index ) const;
   inline       bool                        neighborIsLocatedOnTheSameProcess ( const uint_t index ) const;
   inline       bool                        neighborIsAssignedToTheSameProcess( const uint_t index ) const;

   void assembleNeighborhood(); ///< assemble 'neighborhood_' from 'neighborhoodSection_[i]'

   void split();

   uint_t getIndex() const { return index_; }
   void   setIndex( const uint_t index ) { index_ = index; }



private:

   BlockID Id_;
   uint_t  process_;

   Set<SUID> state_; ///< a set of SUIDs the describes the current state of this block

   AABB   aabb_;
   uint_t level_; // 0=coarse (= block is located on the initial grid) -> 1 -> 2 -> finer

   workload_t  workload_;
   memory_t    memory_;

   bool  marker_; ///< used during the initial setup phase

   SetupBlock*                 father_;
   std::vector< SetupBlock* >  children_;

   std::vector< SetupBlock* >  neighborhoodSection_[26]; // the 26 neighborhood sections
   std::vector< SetupBlock* >  neighborhood_;            // all neighbor blocks

   uint_t index_; ///< used during static load balancing with METIS
};



inline SetupBlock::SetupBlock( SetupBlock* const father, const BlockID& id,
                               const real_t xmin, const real_t ymin, const real_t zmin, // incl.
                               const real_t xmax, const real_t ymax, const real_t zmax, // excl.
                               const uint_t level ) :

   Id_( id ), process_( 0 ), level_( level ), workload_( 0 ), memory_( 0 ), marker_( false ), father_( father ), index_( 0 ) {

   aabb_.initMinMaxCorner( xmin, ymin, zmin, xmax, ymax, zmax );
}



inline void SetupBlock::setChild( const uint_t index, SetupBlock* const child )
{
   WALBERLA_ASSERT_LESS( index, 8 );
   WALBERLA_ASSERT( children_.empty() || children_.size() == 8 );

   if( children_.empty() )
      children_.resize( 8, nullptr );

   children_[index] = child;
}



inline const SetupBlock* SetupBlock::getNeighbor( const uint_t sectionIndex, const uint_t neighborIndex ) const
{
   WALBERLA_ASSERT_LESS( sectionIndex, uint_c(26) );
   WALBERLA_ASSERT_LESS( neighborIndex, neighborhoodSection_[sectionIndex].size() );

   return neighborhoodSection_[sectionIndex][neighborIndex];
}



inline SetupBlock* SetupBlock::getNeighbor( const uint_t sectionIndex, const uint_t neighborIndex )
{
   WALBERLA_ASSERT_LESS( sectionIndex, uint_c(26) );
   WALBERLA_ASSERT_LESS( neighborIndex, neighborhoodSection_[sectionIndex].size() );

   return neighborhoodSection_[sectionIndex][neighborIndex];
}



inline void SetupBlock::addNeighbor( const uint_t sectionIndex, SetupBlock* const block )
{
   WALBERLA_ASSERT_LESS( sectionIndex, uint_c(26) );
   WALBERLA_ASSERT_NOT_NULLPTR( block );
   WALBERLA_ASSERT_LESS( neighborhoodSection_[sectionIndex].size(), getBlockMaxNeighborhoodSectionSize( sectionIndex ) );

   return neighborhoodSection_[sectionIndex].push_back( block );
}



inline const BlockID& SetupBlock::getNeighborId( const uint_t sectionIndex, const uint_t neighborIndex ) const
{
   WALBERLA_ASSERT_LESS( sectionIndex, uint_c(26) );
   WALBERLA_ASSERT_LESS( neighborIndex, neighborhoodSection_[sectionIndex].size() );

   return neighborhoodSection_[sectionIndex][neighborIndex]->getId();
}



inline uint_t SetupBlock::getNeighborProcess( const uint_t sectionIndex, const uint_t neighborIndex ) const
{
   WALBERLA_ASSERT_LESS( sectionIndex, uint_c(26) );
   WALBERLA_ASSERT_LESS( neighborIndex, neighborhoodSection_[sectionIndex].size() );

   return neighborhoodSection_[sectionIndex][neighborIndex]->getProcess();
}

inline uint_t SetupBlock::getNeighborTargetProcess( const uint_t sectionIndex, const uint_t neighborIndex ) const
{
   return getNeighborProcess( sectionIndex, neighborIndex );
}



inline const Set<SUID>& SetupBlock::getNeighborState( const uint_t sectionIndex, const uint_t neighborIndex ) const
{
   WALBERLA_ASSERT_LESS( sectionIndex, uint_c(26) );
   WALBERLA_ASSERT_LESS( neighborIndex, neighborhoodSection_[sectionIndex].size() );

   return neighborhoodSection_[sectionIndex][neighborIndex]->getState();
}



inline bool SetupBlock::neighborIsLocatedOnTheSameProcess ( const uint_t sectionIndex, const uint_t neighborIndex ) const
{
   WALBERLA_ASSERT_LESS( sectionIndex, uint_c(26) );
   WALBERLA_ASSERT_LESS( neighborIndex, neighborhoodSection_[sectionIndex].size() );

   return neighborhoodSection_[sectionIndex][neighborIndex]->getProcess() == getProcess();
}



inline bool SetupBlock::neighborIsAssignedToTheSameProcess( const uint_t sectionIndex, const uint_t neighborIndex ) const
{
   WALBERLA_ASSERT_LESS( sectionIndex, uint_c(26) );
   WALBERLA_ASSERT_LESS( neighborIndex, neighborhoodSection_[sectionIndex].size() );

   return neighborhoodSection_[sectionIndex][neighborIndex]->getTargetProcess() == getTargetProcess();
}



inline bool SetupBlock::neighborhoodSectionHasBlocks( const uint_t sectionIndex ) const
{
   WALBERLA_ASSERT_LESS( sectionIndex, uint_c(26) );

   return !neighborhoodSection_[sectionIndex].empty();
}

inline bool SetupBlock::neighborhoodSectionHasSmallerBlocks( const uint_t sectionIndex ) const
{
   WALBERLA_ASSERT_LESS( sectionIndex, uint_c(26) );

   return !neighborhoodSection_[sectionIndex].empty() && neighborhoodSection_[sectionIndex][0]->level_ > level_;
}

inline bool SetupBlock::neighborhoodSectionHasEquallySizedBlock( const uint_t sectionIndex ) const
{
   WALBERLA_ASSERT_LESS( sectionIndex, uint_c(26) );

   return !neighborhoodSection_[sectionIndex].empty() && neighborhoodSection_[sectionIndex][0]->level_ == level_;
}

inline bool SetupBlock::neighborhoodSectionHasLargerBlock( const uint_t sectionIndex ) const
{
   WALBERLA_ASSERT_LESS( sectionIndex, uint_c(26) );

   return !neighborhoodSection_[sectionIndex].empty() && neighborhoodSection_[sectionIndex][0]->level_ < level_;
}



inline const SetupBlock* SetupBlock::getNeighbor( const uint_t index ) const
{
   WALBERLA_ASSERT_LESS( index, neighborhood_.size() );

   return neighborhood_[index];
}



inline const BlockID& SetupBlock::getNeighborId( const uint_t index ) const
{
   WALBERLA_ASSERT_LESS( index, neighborhood_.size() );

   return neighborhood_[index]->getId();
}



inline uint_t SetupBlock::getNeighborProcess( const uint_t index ) const
{
   WALBERLA_ASSERT_LESS( index, neighborhood_.size() );

   return neighborhood_[index]->getProcess();
}

inline uint_t SetupBlock::getNeighborTargetProcess( const uint_t index ) const {
   return getNeighborProcess( index );
}



inline const Set<SUID>& SetupBlock::getNeighborState( const uint_t index ) const
{
   WALBERLA_ASSERT_LESS( index, neighborhood_.size() );

   return neighborhood_[index]->getState();
}



inline bool SetupBlock::neighborIsLocatedOnTheSameProcess ( const uint_t index ) const
{
   WALBERLA_ASSERT_LESS( index, neighborhood_.size() );

   return neighborhood_[index]->getProcess() == getProcess();
}



inline bool SetupBlock::neighborIsAssignedToTheSameProcess( const uint_t index ) const
{
   WALBERLA_ASSERT_LESS( index, neighborhood_.size() );

   return neighborhood_[index]->getTargetProcess() == getTargetProcess();
}



} // namespace blockforest

using blockforest::SetupBlock;

} // namespace walberla
