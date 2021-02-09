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
//! \file BlockReconstruction.h
//! \ingroup blockforest
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "BlockID.h"
#include "BlockNeighborhoodSection.h"
#include "BlockNeighborhoodConstruction.h"
#include "Types.h"

#include "core/debug/Debug.h"
#include "core/math/AABB.h"
#include "core/uid/SUID.h"

#include <map>
#include <set>
#include <vector>


namespace walberla {
namespace blockforest {



class BlockReconstruction {

public:

   class AABBReconstruction {
   public:
      AABBReconstruction( const AABB& domain, const uint_t xSize, const uint_t ySize, const uint_t zSize, const uint_t treeIdDigits ) :
         domain_( domain ), xSize_( xSize ), ySize_( ySize ), zSize_( zSize ), treeIdDigits_( treeIdDigits ) {}
      uint_t operator()( AABB& aabb, const BlockID& blockId ) const { return BlockReconstruction::reconstructAABB( aabb, blockId, domain_,
                                                                                                                   xSize_, ySize_, zSize_,
                                                                                                                   treeIdDigits_ ); }
   private:
      const AABB   domain_;
      const uint_t xSize_, ySize_, zSize_;
      const uint_t treeIdDigits_;
   };



   class NeighborhoodReconstructionBlock {
   public:
      NeighborhoodReconstructionBlock( const BlockID& id, const uint_t process, const AABB& aabb ) :
         id_( id ), process_( process ), state_( Set<SUID>::emptySet() ), aabb_( aabb ) {}
      NeighborhoodReconstructionBlock( const BlockID& id, const uint_t process, const AABBReconstruction& aabb ) :
         id_( id ), process_( process ), state_( Set<SUID>::emptySet() ) { aabb( aabb_, id ); }
      NeighborhoodReconstructionBlock( const BlockID& id, const uint_t process, const Set<SUID>& state, const AABBReconstruction& aabb ) :
         id_( id ), process_( process ), state_( state ) { aabb( aabb_, id ); }
      bool operator==( const NeighborhoodReconstructionBlock& rhs ) const { return id_ == rhs.id_; }
      const BlockID&   getId()      const { return id_; }
            uint_t     getProcess() const { return process_; }
      const Set<SUID>& getState()   const { return state_; }
      const AABB&      getAABB()    const { return aabb_; }
   private:
      BlockID   id_;
      uint_t    process_;
      Set<SUID> state_;
      AABB      aabb_;
   };



   template< typename BLOCK >
   class NeighborhoodReconstruction {
   public:
      NeighborhoodReconstruction( const AABB& domain, const bool xPeriodic, const bool yPeriodic, const bool zPeriodic ) :
         domain_( domain ), xPeriodic_( xPeriodic ), yPeriodic_( yPeriodic ), zPeriodic_( zPeriodic ) {}
      void operator()( BLOCK* block, const std::vector< NeighborhoodReconstructionBlock >& neighbors ) const {
         BlockReconstruction::reconstructNeighborhood( block, neighbors, domain_, xPeriodic_, yPeriodic_, zPeriodic_ ); }
   private:
      const AABB domain_;
      const bool xPeriodic_, yPeriodic_, zPeriodic_;
   };



   template< typename BLOCK >
   class NeighborhoodSectionReconstruction {
   public:
      NeighborhoodSectionReconstruction( const AABB& domain, const uint_t xSize, const uint_t ySize, const uint_t zSize, const bool xPeriodic,
                                         const bool yPeriodic, const bool zPeriodic, const uint_t treeIdDigits ) :
         domain_( domain ), xSize_( xSize ), ySize_( ySize ), zSize_( zSize ), xPeriodic_( xPeriodic ),
         yPeriodic_( yPeriodic ), zPeriodic_( zPeriodic ), treeIdDigits_( treeIdDigits ) {}
      void operator()( BLOCK* block ) const { BlockReconstruction::reconstructNeighborhoodSections( block, domain_, xSize_, ySize_, zSize_,
                                                                                          xPeriodic_, yPeriodic_, zPeriodic_, treeIdDigits_ ); }
   private:
      const AABB   domain_;
      const uint_t xSize_, ySize_, zSize_;
      const bool   xPeriodic_, yPeriodic_, zPeriodic_;
      const uint_t treeIdDigits_;
   };



   static uint_t reconstructAABB( AABB& aabb, const BlockID& blockId, const AABB& domain, const uint_t xSize, const uint_t ySize,
                                  const uint_t zSize, const uint_t treeIdDigits );

   template< typename BLOCK >
   static void reconstructNeighborhood( BLOCK* block, const std::vector< NeighborhoodReconstructionBlock >& neighbors, const AABB& domain,
                                        const bool xPeriodic, const bool yPeriodic, const bool zPeriodic );

   template< typename BLOCK >
   static void reconstructNeighborhoodSections( BLOCK* block, const AABB& domain, const uint_t xSize, const uint_t ySize, const uint_t zSize,
                                                const bool xPeriodic, const bool yPeriodic, const bool zPeriodic, const uint_t treeIdDigits );
};



template< typename BLOCK >
void BlockReconstruction::reconstructNeighborhood( BLOCK* block, const std::vector< NeighborhoodReconstructionBlock >& neighbors,
                                                   const AABB& domain, const bool xPeriodic, const bool yPeriodic, const bool zPeriodic ) {

   std::vector< real_t >                                      neighborhoodSectionBlockCenters;
   std::set< const NeighborhoodReconstructionBlock* >         neighborhood;
   std::map< const NeighborhoodReconstructionBlock*, uint_t > neighborhoodIndex;
   std::vector< uint_t >                                      neighborhoodSectionBlocks[26];

   block->clearNeighborhood();

   for( uint_t n = 0; n != 26; ++n ) {

      constructNeighborhoodSectionBlockCenters( n, block->getAABB(), neighborhoodSectionBlockCenters );

      WALBERLA_ASSERT_EQUAL( neighborhoodSectionBlockCenters.size() % 3, 0 );

      for( uint_t p = 0; p != neighborhoodSectionBlockCenters.size(); p += 3 ) {

         real_t x = neighborhoodSectionBlockCenters[p];
         real_t y = neighborhoodSectionBlockCenters[p+1];
         real_t z = neighborhoodSectionBlockCenters[p+2];

         // treat periodicity
         if( x <  domain.xMin() && xPeriodic ) x = domain.xMax() - domain.xMin() + x;
         if( x >= domain.xMax() && xPeriodic ) x = domain.xMin() - domain.xMax() + x;
         if( y <  domain.yMin() && yPeriodic ) y = domain.yMax() - domain.yMin() + y;
         if( y >= domain.yMax() && yPeriodic ) y = domain.yMin() - domain.yMax() + y;
         if( z <  domain.zMin() && zPeriodic ) z = domain.zMax() - domain.zMin() + z;
         if( z >= domain.zMax() && zPeriodic ) z = domain.zMin() - domain.zMax() + z;

         for( uint_t i = 0; i != neighbors.size(); ++i ) {
            if( neighbors[i].getAABB().contains( x, y, z ) ) {

               const NeighborhoodReconstructionBlock* neighbor = &(neighbors[i]);
               uint_t index = 0;

               if( neighborhood.insert( neighbor ).second ) {

                  index = block->getNeighborhoodSize();
                  neighborhoodIndex[ neighbor ] = index;

                  block->addNeighbor( neighbor->getId(), neighbor->getProcess(), neighbor->getState() );
               }
               else index = neighborhoodIndex[ neighbor ];

               if( neighborhoodSectionBlocks[n].empty() || neighborhoodSectionBlocks[n].back() != index )
                  neighborhoodSectionBlocks[n].push_back( index );

               break;
            }
         }
      }

#ifndef NDEBUG
      for( uint_t v = 0; v != neighborhoodSectionBlocks[n].size(); ++v )
         for( uint_t w = v+1; w != neighborhoodSectionBlocks[n].size(); ++w )
            WALBERLA_ASSERT_UNEQUAL( neighborhoodSectionBlocks[n][v], neighborhoodSectionBlocks[n][w] );
      if( !neighborhoodSectionBlocks[n].empty() &&
            block->getNeighbor( neighborhoodSectionBlocks[n].back() ).getId().getUsedBits() > block->getId().getUsedBits() )
      {
         WALBERLA_ASSERT_EQUAL( neighborhoodSectionBlocks[n].size(), getBlockMaxNeighborhoodSectionSize(n) );
      }
      else
         WALBERLA_ASSERT( neighborhoodSectionBlocks[n].empty() || neighborhoodSectionBlocks[n].size() == 1 );
#endif

      neighborhoodSectionBlockCenters.clear();
   }

   for( uint_t n = 0; n != 26; ++n ) {

      block->clearNeighborhoodSection(n);

      for( uint_t i = 0; i != neighborhoodSectionBlocks[n].size(); ++i )
         block->addNeighbor( n, neighborhoodSectionBlocks[n][i] );
   }

#ifndef NDEBUG
   std::set< const typename BLOCK::NeighborBlock* > blockNeighborhood;
   for( uint_t i = 0; i != block->getNeighborhoodSize(); ++i )
      blockNeighborhood.insert( &(block->getNeighbor(i)) );
   for( uint_t i = 0; i != 26; ++i )
      for( uint_t j = 0; j != block->getNeighborhoodSectionSize(i); ++j )
         WALBERLA_ASSERT( blockNeighborhood.find( &(block->getNeighbor(i,j)) ) != blockNeighborhood.end() );
#endif
}



template< typename BLOCK >
void BlockReconstruction::reconstructNeighborhoodSections( BLOCK* block, const AABB& domain, const uint_t xSize, const uint_t ySize,
                                                           const uint_t zSize, const bool xPeriodic, const bool yPeriodic,
                                                           const bool zPeriodic, const uint_t treeIdDigits ) {

   std::vector< AABB > aabb( block->getNeighborhoodSize() );

   for( uint_t i = 0; i != aabb.size(); ++i )
      reconstructAABB( aabb[i], block->getNeighborId(i), domain, xSize, ySize, zSize, treeIdDigits );

   std::vector< real_t > neighborhoodSectionBlockCenters;
   std::vector< uint_t > neighborhoodSectionBlocks;

   for( uint_t n = 0; n != 26; ++n ) {

      constructNeighborhoodSectionBlockCenters( n, block->getAABB(), neighborhoodSectionBlockCenters );

      WALBERLA_ASSERT_EQUAL( neighborhoodSectionBlockCenters.size() % 3, 0 );

      for( uint_t p = 0; p != neighborhoodSectionBlockCenters.size(); p += 3 ) {

         real_t x = neighborhoodSectionBlockCenters[p];
         real_t y = neighborhoodSectionBlockCenters[p+1];
         real_t z = neighborhoodSectionBlockCenters[p+2];

         // treat periodicity
         if( x <  domain.xMin() && xPeriodic ) x = domain.xMax() - domain.xMin() + x;
         if( x >= domain.xMax() && xPeriodic ) x = domain.xMin() - domain.xMax() + x;
         if( y <  domain.yMin() && yPeriodic ) y = domain.yMax() - domain.yMin() + y;
         if( y >= domain.yMax() && yPeriodic ) y = domain.yMin() - domain.yMax() + y;
         if( z <  domain.zMin() && zPeriodic ) z = domain.zMax() - domain.zMin() + z;
         if( z >= domain.zMax() && zPeriodic ) z = domain.zMin() - domain.zMax() + z;

         for( uint_t i = 0; i != aabb.size(); ++i ) {
            if( aabb[i].contains( x, y, z ) && ( neighborhoodSectionBlocks.empty() || neighborhoodSectionBlocks.back() != i ) )
            {
               neighborhoodSectionBlocks.push_back(i);
               break;
            }
         }
      }

#ifndef NDEBUG
      for( uint_t v = 0; v != neighborhoodSectionBlocks.size(); ++v )
         for( uint_t w = v+1; w != neighborhoodSectionBlocks.size(); ++w )
            WALBERLA_ASSERT_UNEQUAL( neighborhoodSectionBlocks[v], neighborhoodSectionBlocks[w] );
      if( !neighborhoodSectionBlocks.empty() &&
          block->getNeighborId( neighborhoodSectionBlocks.back() ).getUsedBits() > block->getId().getUsedBits() )
      {
         WALBERLA_ASSERT_EQUAL( neighborhoodSectionBlocks.size(), getBlockMaxNeighborhoodSectionSize(n) );
      }
      else
         WALBERLA_ASSERT( neighborhoodSectionBlocks.empty() || neighborhoodSectionBlocks.size() == 1 );
#endif

      block->clearNeighborhoodSection(n);

      for( uint_t i = 0; i != neighborhoodSectionBlocks.size(); ++i )
         block->addNeighbor( n, neighborhoodSectionBlocks[i] );

      neighborhoodSectionBlocks.clear();
      neighborhoodSectionBlockCenters.clear();
   }
}



} // namespace blockforest
} // namespace walberla


