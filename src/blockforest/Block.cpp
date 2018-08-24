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
//! \file Block.cpp
//! \ingroup blockforest
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#include "Block.h"
#include "BlockForest.h"
#include "PhantomBlock.h"
#include "SetupBlock.h"

#include "core/mpi/MPIManager.h"

#include <map>


namespace walberla {
namespace blockforest {



Block::NeighborBlock::NeighborBlock( const BlockForest & forest, const BlockID & id, const uint_t process, const Set<SUID> & state ) :
         id_( id ), process_( process ), state_( state ), aabb_( forest.getAABBFromBlockId(id) )
{}



Block::Block( BlockForest & forest, const SetupBlock * const block ) :

   IBlock( forest, block->getAABB(), block->getId().getID() ),

   forest_( forest ), id_( block->getId() ), level_( block->getLevel() ), targetLevel_( block->getLevel() )
{
   setState( block->getState() );

   std::map< BlockID, uint_t > neighborhoodMapping;

   for( uint_t i = 0; i != block->getNeighborhoodSize(); ++i )
   {
      neighborhood_.emplace_back( forest, block->getNeighborId(i), block->getNeighborProcess(i), block->getNeighborState(i) );
      neighborhoodMapping[ block->getNeighborId(i) ] = i;
   }

   for( uint_t i = 0; i != 26; ++i ) {
      for( uint_t j = 0; j != block->getNeighborhoodSectionSize(i); ++j )
      {
         WALBERLA_ASSERT( neighborhoodMapping.find( block->getNeighborId(i,j) ) != neighborhoodMapping.end() );

         neighborhoodSection_[i].push_back( &(neighborhood_[ neighborhoodMapping[block->getNeighborId(i,j)] ]) );
      }
   }
}



Block::Block( BlockForest & forest, const BlockID & id, const AABB & aabb, const Set<SUID> & state, const uint_t level,
              const BlockReconstruction::NeighborhoodReconstruction< Block > & neighborhoodReconstruction,
              const std::vector< BlockReconstruction::NeighborhoodReconstructionBlock > & neighbors ) :

   IBlock( forest, aabb, id.getID() ),

   forest_( forest ), id_( id ), level_( level ), targetLevel_( level )
{
   setState( state );

   neighborhoodReconstruction( this, neighbors );
}



Block::Block( BlockForest & forest, const PhantomBlock & phantom ) :

   IBlock( forest, phantom.getAABB(), phantom.getId().getID() ),

   forest_( forest ), id_( phantom.getId() ), level_( phantom.getLevel() ), targetLevel_( phantom.getLevel() )
{
   setState( phantom.getState() );

   resetNeighborhood( phantom );
}



Block::Block( BlockForest & forest, const BlockID & id, const AABB & aabb, const uint_t level, mpi::RecvBuffer & buffer,
              const std::function< uint_t ( const uint_t ) > & processMapping ) :

   IBlock( forest, aabb, id.getID() ),

   forest_( forest ), id_( id ), level_( level ), targetLevel_( level )
{
   Set<SUID> state;
   buffer >> state;
   setState( state );

   uint_t size(0);
   buffer >> size;

   if( processMapping )
   {
      for( uint_t i = uint_t(0); i != size; ++i )
      {
         BlockID   nId;
         Set<SUID> nState;
         uint_t    nProcess;
         buffer >> nId >> nState >> nProcess;
         addNeighbor( nId, processMapping(nProcess), nState );
         WALBERLA_ASSERT_LESS( processMapping(nProcess), uint_c( mpi::MPIManager::instance()->numProcesses() ) );
      }
   }
   else
   {
      for( uint_t i = uint_t(0); i != size; ++i )
      {
         BlockID   nId;
         Set<SUID> nState;
         uint_t    nProcess;
         buffer >> nId >> nState >> nProcess;
         addNeighbor( nId, nProcess, nState );
         WALBERLA_ASSERT_LESS( nProcess, uint_c( mpi::MPIManager::instance()->numProcesses() ) );
      }
   }

   for( uint_t i = uint_t(0); i != uint_t(26); ++i )
   {
      buffer >> size;
      for( uint_t j = uint_t(0); j != size; ++j )
      {
         uint8_t index(0);
         buffer >> index;
         WALBERLA_ASSERT_LESS( index, neighborhood_.size() );
         neighborhoodSection_[i].push_back( &(neighborhood_[index]) );
      }
   }
}



void Block::toBuffer( mpi::SendBuffer & buffer ) const
{
   buffer << this->getState()
          << uint_c( neighborhood_.size() );

   for( uint_t i = uint_t(0); i != neighborhood_.size(); ++i )
   {
      buffer << neighborhood_[i].getId()
             << neighborhood_[i].getState()
             << neighborhood_[i].getProcess();
   }

   for( uint_t i = uint_t(0); i != uint_t(26); ++i )
   {
      buffer << uint_c( neighborhoodSection_[i].size() );
      for( uint_t j = uint_t(0); j != neighborhoodSection_[i].size(); ++j )
      {
         uint8_t index(0);
         for( uint_t k = uint_t(0); k != neighborhood_.size() && neighborhoodSection_[i][j] != &(neighborhood_[k]); ++k, ++index ) {}
         WALBERLA_ASSERT_LESS( index, neighborhood_.size() );
         buffer << index;
      }
   }
}



uint_t Block::getProcess() const {

   return forest_.getProcess();
}



/// The following members are not used for checking if two Block objects are equal: forest_
bool Block::equal( const IBlock* rhs ) const {

   const Block* block = dynamic_cast< const Block* >( rhs );

   if( block != rhs )
      return false;

   if( id_ != block->id_ || level_ != block->level_ )
      return false;

   for( uint_t i = 0; i != 26; ++i ) {
      if( neighborhoodSection_[i].size() != block->neighborhoodSection_[i].size() )
         return false;
      for( uint_t j = 0; j != neighborhoodSection_[i].size(); ++j )
         if( *(neighborhoodSection_[i][j]) != *(block->neighborhoodSection_[i][j]) )
            return false;
   }

   if( neighborhood_.size() != block->neighborhood_.size() )
      return false;

   for( uint_t i = 0; i != neighborhood_.size(); ++i )
      if( neighborhood_[i] != block->neighborhood_[i] )
         return false;

   return true;
}



void Block::resetNeighborhood( const PhantomBlock & phantom )
{
   std::map< BlockID, uint_t > neighborhoodMapping;

   neighborhood_.clear();
   for( uint_t i = 0; i != phantom.getNeighborhoodSize(); ++i )
   {
      neighborhood_.emplace_back( forest_, phantom.getNeighborId(i), phantom.getNeighborProcess(i), phantom.getNeighborState(i) );
      neighborhoodMapping[ phantom.getNeighborId(i) ] = i;
   }

   for( uint_t i = 0; i != 26; ++i )
   {
      neighborhoodSection_[i].clear();
      for( uint_t j = 0; j != phantom.getNeighborhoodSectionSize(i); ++j )
      {
         WALBERLA_ASSERT( neighborhoodMapping.find( phantom.getNeighborId(i,j) ) != neighborhoodMapping.end() );

         neighborhoodSection_[i].push_back( &(neighborhood_[ neighborhoodMapping[phantom.getNeighborId(i,j)] ]) );
      }
   }
}



} // namespace blockforest
} // namespace walberla
