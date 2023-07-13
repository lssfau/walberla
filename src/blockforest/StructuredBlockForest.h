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
//! \file StructuredBlockForest.h
//! \ingroup blockforest
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "BlockForest.h"
#include "domain_decomposition/StructuredBlockStorage.h"

#include <functional>

namespace walberla {
namespace blockforest {



class StructuredBlockForest : public StructuredBlockStorage {

public:

   inline StructuredBlockForest( const shared_ptr< BlockForest >& blockForest,
                                 const uint_t blockXCells, const uint_t blockYCells, const uint_t blockZCells );

   const BlockForest& getBlockForest() const { return *blockForest_; }
         BlockForest& getBlockForest()       { return *blockForest_; }

   const shared_ptr<       BlockForest > & getBlockForestPointer()       { return blockForest_; }
         shared_ptr< const BlockForest >   getBlockForestPointer() const { return blockForest_; }

   // implementation of pure virtual functions of base class StructuredBlockStorage

   using StructuredBlockStorage::blockExists;
   using StructuredBlockStorage::blockExistsLocally;
   using StructuredBlockStorage::blockExistsRemotely;

          bool blockExists        ( const Cell& cell, const uint_t level = 0 ) const override;
   inline bool blockExistsLocally ( const Cell& cell, const uint_t level = 0 ) const override;
   inline bool blockExistsRemotely( const Cell& cell, const uint_t level = 0 ) const override;

   void getBlockID( IBlockID& id, const Cell& cell, const uint_t level = 0 ) const override;

   inline uint_t getLevel( const IBlock& block ) const override;

   using StructuredBlockStorage::getNumberOfXCells;
   using StructuredBlockStorage::getNumberOfYCells;
   using StructuredBlockStorage::getNumberOfZCells;

#ifdef NDEBUG
   uint_t getNumberOfXCells( const IBlock& /*block*/ ) const override { return blockCells_[0]; }
   uint_t getNumberOfYCells( const IBlock& /*block*/ ) const override { return blockCells_[1]; }
   uint_t getNumberOfZCells( const IBlock& /*block*/ ) const override { return blockCells_[2]; }
#else
   uint_t getNumberOfXCells( const IBlock& block ) const override { WALBERLA_ASSERT_EQUAL( &(getBlockStorage()), &(block.getBlockStorage()) ) return blockCells_[0]; }
   uint_t getNumberOfYCells( const IBlock& block ) const override { WALBERLA_ASSERT_EQUAL( &(getBlockStorage()), &(block.getBlockStorage()) ) return blockCells_[1]; }
   uint_t getNumberOfZCells( const IBlock& block ) const override { WALBERLA_ASSERT_EQUAL( &(getBlockStorage()), &(block.getBlockStorage()) ) return blockCells_[2]; }
#endif

   using StructuredBlockStorage::getNumberOfCells;

   inline uint_t getNumberOfCells( const IBlock& block, const uint_t index ) const override;

   // direct access to all member functions which are special to of BlockForest (-> for documentation of these functions see class BlockForest)

   uint_t getProcess()        const { return blockForest_->getProcess(); }
   uint_t getProcessIdBytes() const { return blockForest_->getProcessIdBytes(); }

   uint_t getXSize() const { return blockForest_->getXSize(); }
   uint_t getYSize() const { return blockForest_->getYSize(); }
   uint_t getZSize() const { return blockForest_->getZSize(); }
   uint_t getSize( const uint_t index ) const { return blockForest_->getSize( index ); }

   real_t getRootBlockXSize() const { return blockForest_->getRootBlockXSize(); }
   real_t getRootBlockYSize() const { return blockForest_->getRootBlockYSize(); }
   real_t getRootBlockZSize() const { return blockForest_->getRootBlockZSize(); }

   bool storesUniformBlockGrid() const { return blockForest_->storesUniformBlockGrid(); }

   uint_t getDepth()          const { return blockForest_->getDepth(); }
   uint_t getNumberOfLevels() const { return blockForest_->getNumberOfLevels(); }
   uint_t getTreeIdDigits()   const { return blockForest_->getTreeIdDigits(); }
   uint_t getBlockIdBytes()   const { return blockForest_->getBlockIdBytes(); }

   uint_t getNumberOfBlocks() const { return blockForest_->getNumberOfBlocks(); }
   uint_t getNumberOfBlocks( const uint_t level ) const { return blockForest_->getNumberOfBlocks( level ); }

   using StructuredBlockStorage::getBlocks;

   void getBlocks( std::vector< const Block* >& blocks, const uint_t level ) const { blockForest_->getBlocks( blocks, level ); }
   void getBlocks( std::vector<       Block* >& blocks, const uint_t level )       { blockForest_->getBlocks( blocks, level ); }

   const Block* getRootBlock( const uint_t x, const uint_t y, const uint_t z ) const { return blockForest_->getRootBlock(x,y,z); }
         Block* getRootBlock( const uint_t x, const uint_t y, const uint_t z )       { return blockForest_->getRootBlock(x,y,z); }

   bool rootBlockExists        ( const uint_t x, const uint_t y, const uint_t z ) const { return blockForest_->rootBlockExists        (x,y,z); }
   bool rootBlockExistsLocally ( const uint_t x, const uint_t y, const uint_t z ) const { return blockForest_->rootBlockExistsLocally (x,y,z); }
   bool rootBlockExistsRemotely( const uint_t x, const uint_t y, const uint_t z ) const { return blockForest_->rootBlockExistsRemotely(x,y,z); }

   void getRootBlockAABB       ( AABB&      aabb,  const uint_t x, const uint_t y, const uint_t z ) const { blockForest_->getRootBlockAABB       ( aabb,  x, y, z ); }
   void getRootBlockState      ( Set<SUID>& state, const uint_t x, const uint_t y, const uint_t z ) const { blockForest_->getRootBlockState      ( state, x, y, z ); }
   void getRootBlockProcessRank( uint_t&    rank,  const uint_t x, const uint_t y, const uint_t z ) const { blockForest_->getRootBlockProcessRank( rank,  x, y, z ); }

   const BlockForest::BlockInformation& getBlockInformation() const { return blockForest_->getBlockInformation(); }

   uint_t getLevelFromBlockId( const BlockID& id ) const { return blockForest_->getLevelFromBlockId( id ); }
   uint_t getAABBFromBlockId( AABB& aabb, const BlockID& id ) const { return blockForest_->getAABBFromBlockId( aabb, id ); }

   void getForestCoordinates   ( uint_t& x, uint_t& y, uint_t& z, const BlockID& id ) const { blockForest_->getForestCoordinates   (x,y,z,id); }
   void getRootBlockCoordinates( uint_t& x, uint_t& y, uint_t& z, const BlockID& id ) const { blockForest_->getRootBlockCoordinates(x,y,z,id); }
   void getRootBlockID( BlockID& id, const uint_t x, const uint_t y, const uint_t z ) const { blockForest_->getRootBlockID( id, x, y, z ); }

   bool insertBuffersIntoProcessNetwork() const { return blockForest_->insertBuffersIntoProcessNetwork(); }

   const std::vector< uint_t >& getNeighborhood() const { return blockForest_->getNeighborhood(); }

   internal::BlockDataHandlingAdder addBlockData( const std::string & identifier = std::string() ) { return blockForest_->addBlockData( identifier ); }
   
   template< typename T >
   inline BlockDataID addBlockData( const shared_ptr< T > & dataHandling,
                                    const std::string & identifier          = std::string(),
                                    const Set<SUID> & requiredSelectors     = Set<SUID>::emptySet(),
                                    const Set<SUID> & incompatibleSelectors = Set<SUID>::emptySet() )
   { return blockForest_->addBlockData( dataHandling, identifier, requiredSelectors, incompatibleSelectors ); }
   
   template< typename T >
   inline BlockDataID addBlockData( std::function< T* ( IBlock* const block ) > function,
                                    const std::string& identifier = std::string(),
                                    const Set<SUID>& requiredSelectors     = Set<SUID>::emptySet(),
                                    const Set<SUID>& incompatibleSelectors = Set<SUID>::emptySet() )
   { return blockForest_->addBlockData( function, identifier, requiredSelectors, incompatibleSelectors ); }

   template< typename T >
   inline BlockDataID loadBlockData( const std::string & file, const shared_ptr< T > & dataHandling,
                                     const std::string & identifier          = std::string(),
                                     const Set<SUID> & requiredSelectors     = Set<SUID>::emptySet(),
                                     const Set<SUID> & incompatibleSelectors = Set<SUID>::emptySet() )
   { return blockForest_->loadBlockData( file, dataHandling, identifier, requiredSelectors, incompatibleSelectors ); }
   
   inline void refresh() { blockForest_->refresh(); }

   // new, additional, cell-related functions special to StructuredBlockForest

   uint_t getNumberOfXCellsPerBlock() const { return blockCells_[0]; }
   uint_t getNumberOfYCellsPerBlock() const { return blockCells_[1]; }
   uint_t getNumberOfZCellsPerBlock() const { return blockCells_[2]; }

   uint_t getNumberOfCellsPerBlock( const uint_t index ) const { WALBERLA_ASSERT_LESS( index, 3 ); return blockCells_[ index ]; }

protected:

   inline bool equal( const StructuredBlockStorage* rhs ) const override;

   // helper class for 'StructuredBlockForest::addCellBoundingBoxesAsBlockData'
   class CellBoundingBoxHandling : public AlwaysInitializeBlockDataHandling< CellInterval >
   {
   public:
      CellBoundingBoxHandling( const StructuredBlockForest & forest ) : forest_( forest ) {}
      CellInterval * initialize( IBlock * const block ) override { return forest_.initializeCellBoundingBox( block ); }
   private:
      const StructuredBlockForest & forest_;
   };
   friend class CellBoundingBoxHandling;

   //using StructuredBlockStorage::initializeCellBoundingBox;
   inline BlockDataID addCellBoundingBoxesAsBlockData( const std::string & identifier ) override;

private:

   static inline std::vector< uint_t > constructNumberOfCells( uint_t cells, const uint_t levels );
   static inline void resetCellDecompositionInStorage( StructuredBlockForest & storage );



   shared_ptr< BlockForest > blockForest_;

   uint_t blockCells_[3];

}; // class StructuredBlockForest



inline StructuredBlockForest::StructuredBlockForest( const shared_ptr< BlockForest >& blockForest,
                                                     const uint_t blockXCells, const uint_t blockYCells, const uint_t blockZCells ) :

   StructuredBlockStorage( blockForest,
                           constructNumberOfCells( blockXCells * blockForest->getXSize(), blockForest->getNumberOfLevels() ),
                           constructNumberOfCells( blockYCells * blockForest->getYSize(), blockForest->getNumberOfLevels() ),
                           constructNumberOfCells( blockZCells * blockForest->getZSize(), blockForest->getNumberOfLevels() ) ),

   blockForest_( blockForest ) {

   blockForest_->addRefreshCallbackFunctionBeforeBlockDataIsUnpacked(
            BlockForest::RefreshCallbackWrappper( std::bind( resetCellDecompositionInStorage, std::ref(*this) ) ) );

   blockCells_[0] = blockXCells;
   blockCells_[1] = blockYCells;
   blockCells_[2] = blockZCells;
}



inline bool StructuredBlockForest::blockExistsLocally( const Cell& cell, const uint_t level ) const {

   return getBlock( cell, level ) != nullptr;
}



inline bool StructuredBlockForest::blockExistsRemotely( const Cell& cell, const uint_t level ) const {

   return blockExists( cell, level ) && !blockExistsLocally( cell, level );
}



inline uint_t StructuredBlockForest::getLevel( const IBlock& block ) const {

   WALBERLA_ASSERT_EQUAL( dynamic_cast< const Block* >( &block ), &block )

   return static_cast< const Block* >( &block )->getLevel();
}



#ifdef NDEBUG
inline uint_t StructuredBlockForest::getNumberOfCells( const IBlock& /*block*/, const uint_t index ) const {
#else
inline uint_t StructuredBlockForest::getNumberOfCells( const IBlock& block, const uint_t index ) const {
#endif

   WALBERLA_ASSERT_EQUAL( &(getBlockStorage()), &(block.getBlockStorage()) )
   WALBERLA_ASSERT_LESS( index, 3 )

   return blockCells_[ index ];
}



inline bool StructuredBlockForest::equal( const StructuredBlockStorage* rhs ) const {

   const StructuredBlockForest* forest = dynamic_cast< const StructuredBlockForest* >( rhs );

   if( forest != rhs )
      return false;

   // blockForest_: not checked since there is already an equality check in the base class StructuredBlockStorage

   return blockCells_[0] == forest->blockCells_[0] && blockCells_[1] == forest->blockCells_[1] && blockCells_[2] == forest->blockCells_[2];
}



inline BlockDataID StructuredBlockForest::addCellBoundingBoxesAsBlockData( const std::string & identifier )
{
   return addBlockData( walberla::make_shared< CellBoundingBoxHandling >( *this ), identifier );
}



inline std::vector< uint_t > StructuredBlockForest::constructNumberOfCells( uint_t cells, const uint_t levels ) {

   std::vector< uint_t > cellsVector( 1, cells );

   for( uint_t i = 1; i < levels; ++i ) {
      cells *= 2;
      cellsVector.push_back( cells );
   }

   return cellsVector;
}



inline void StructuredBlockForest::resetCellDecompositionInStorage( StructuredBlockForest & storage )
{
   storage.resetCellDecomposition( constructNumberOfCells( storage.blockCells_[0] * storage.getXSize(), storage.getNumberOfLevels() ),
                                   constructNumberOfCells( storage.blockCells_[1] * storage.getYSize(), storage.getNumberOfLevels() ),
                                   constructNumberOfCells( storage.blockCells_[2] * storage.getZSize(), storage.getNumberOfLevels() ) );
}



} // namespace blockforest

using blockforest::StructuredBlockForest;

} // namespace walberla

