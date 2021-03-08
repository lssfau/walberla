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
//! \file PhantomBlockForest.h
//! \ingroup blockforest
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "BlockID.h"
#include "PhantomBlock.h"

#include <map>
#include <vector>



namespace walberla {
namespace blockforest {



class BlockForest;



class PhantomBlockForest
{
public:

   using BlockStateDeterminationFunction = std::function<Set<SUID> (const std::vector<std::pair<BlockID, Set<SUID>>> &, const BlockID &)>;

   using PhantomBlockDataAssignmentFunction = std::function<void (std::vector<std::pair<const PhantomBlock *, walberla::any>> &, const PhantomBlockForest &)>;

   /// \param iteration execution counter of this callback
   /// \return should the callback rerun after phantom block migration?
   using MigrationPreparationFunction = std::function<bool (std::vector<std::pair<const PhantomBlock *, uint_t>> &, std::set<uint_t> &, const PhantomBlockForest &, const uint_t)>; // = load balancing

   using PhantomBlockDataPackFunction = std::function<void (mpi::SendBuffer &, const PhantomBlock &)>;
   using PhantomBlockDataUnpackFunction = std::function<void (mpi::RecvBuffer &, const PhantomBlock &, walberla::any &)>;



   PhantomBlockForest( BlockForest & blockforest );

   const BlockForest & getBlockForest() const { return blockforest_; }

   uint_t getDepth()          const { return depth_; }
   uint_t getNumberOfLevels() const { return depth_ + 1; }

   uint_t getNumberOfBlocks() const { return blocks_.size(); }

   bool blockExists( const BlockID & id ) const { return blocks_.find(id) != blocks_.end(); }

   inline shared_ptr< const PhantomBlock > getBlock( const BlockID & id ) const;
   inline shared_ptr<       PhantomBlock > getBlock( const BlockID & id );

   const std::map< BlockID, shared_ptr< PhantomBlock > > & getBlockMap() const { return blocks_; }

   const std::vector< uint_t > & getNeighborhood() const { return neighborhood_; }
   const std::vector< uint_t > & getNeighboringProcesses() const { return getNeighborhood(); }

   void initialize( const BlockStateDeterminationFunction & function, const bool recalculateDepth );
   void assignBlockData( const PhantomBlockDataAssignmentFunction & function );
   bool calculateMigrationInformation( const MigrationPreparationFunction & function, const uint_t iteration );
   void migrate( const PhantomBlockDataPackFunction & packBlockData, const PhantomBlockDataUnpackFunction & unpackBlockData );
   
private:

   void updateNeighborhood();
   void prepareMigration();
   
   
   
   BlockForest & blockforest_;

   uint_t depth_; // might be different to 'blockforest_' !

   std::map< BlockID, shared_ptr< PhantomBlock > > blocks_;
   
   std::vector< uint_t > neighborhood_; // neighbor processes (not entirely reconstructable from 'blocks_' -> empty buffer processes!)
   std::set< uint_t > processesToRecvFrom_;

}; // class PhantomBlockForest



inline shared_ptr< const PhantomBlock > PhantomBlockForest::getBlock( const BlockID & id ) const
{ 
   auto it = blocks_.find( id );
   if( it != blocks_.end() )
      return it->second;
   return shared_ptr< PhantomBlock >();
}



inline shared_ptr< PhantomBlock > PhantomBlockForest::getBlock( const BlockID & id )
{
   auto it = blocks_.find( id );
   if( it != blocks_.end() )
      return it->second;
   return shared_ptr< PhantomBlock >();
}



} // namespace blockforest

using PhantomBlockForest = blockforest::PhantomBlockForest;

} // namespace walberla
