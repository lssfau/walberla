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
//! \file BlockForest.h
//! \ingroup blockforest
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "Block.h"
#include "BlockDataHandling.h"
#include "BlockID.h"
#include "PhantomBlockForest.h"
#include "Types.h"

#include "core/debug/Debug.h"
#include "core/math/Uint.h"
#include "core/timing/TimingPool.h"

#include "domain_decomposition/BlockStorage.h"

#include <map>
#include <vector>


namespace walberla {
namespace blockforest {



class SetupBlockForest;



class BlockForest : public BlockStorage
{
public:

   using RefreshMinTargetLevelDeterminationFunction = std::function<void (std::vector<std::pair<const Block *, uint_t>> &, std::vector<const Block *> &, const BlockForest &)>;

   using RefreshCallbackFunction = std::function<void (BlockForest &, const PhantomBlockForest &)>;

   using SnapshotCreationFunction = std::function<void (std::vector<uint_t> &, std::vector<uint_t> &)>;
   using SnapshotRestoreFunction = std::function<uint_t (const uint_t)>;
   using SnapshotRestoreCallbackFunction = std::function<void ()>;

   enum FileIOMode { MPI_PARALLEL, MASTER_SLAVE, SERIALIZED_DISTRIBUTED };



   class RefreshFunctor {
   public:
      RefreshFunctor( BlockForest & forest, const uint_t checkFrequency = uint_t(1) ) :
         forest_( forest ), executionCounter_( uint_t(0) ), checkFrequency_( checkFrequency ) {}
      void operator()()
      {
         ++executionCounter_;
         if( checkFrequency_ == uint_t(0) || ( executionCounter_ - uint_c(1) ) % checkFrequency_ != 0 )
            return;
         forest_.refresh();
      }
   private:
      BlockForest & forest_;
      uint_t executionCounter_;
      uint_t checkFrequency_;
   };



   // So that existing void-void functors can be registered as refresh callback functions
   class RefreshCallbackWrappper
   {
   public:
      using Functor_T = std::function<void ()>;
      RefreshCallbackWrappper( const Functor_T & functor ) : functor_( functor ) {}
      void operator()( BlockForest &, const PhantomBlockForest & ) { functor_(); }
   private:
      Functor_T functor_;
   };



   class SnapshotCreationFunctor {
   public:
      SnapshotCreationFunctor( BlockForest & forest, const SnapshotCreationFunction & function,
                               const uint_t checkFrequency = uint_t(1) ) :
         forest_( forest ), executionCounter_( uint_t(0) ), checkFrequency_( checkFrequency ), function_( function ) {}
      void operator()()
      {
         ++executionCounter_;
         if( checkFrequency_ == uint_t(0) || ( executionCounter_ - uint_c(1) ) % checkFrequency_ != 0 )
            return;
         std::vector<uint_t> sendTo;
         std::vector<uint_t> recvFrom;
         function_( sendTo, recvFrom );
         forest_.createSnapshot( sendTo, recvFrom );
      }
   private:
      BlockForest & forest_;
      uint_t executionCounter_;
      uint_t checkFrequency_;
      SnapshotCreationFunction function_;
   };



   class BlockInformation;
   friend class BlockInformation;
   class BlockInformation {

      friend class BlockForest;

   private:

      struct Node {

         Node() : process_( uint_c(0) ), state_( Set<SUID>::emptySet() ) {}
         Node( const uint_t process, const Set<SUID>& state ) : process_( process ), state_( state ) {}

         bool operator==( const Node& rhs ) const {
            if( process_ != rhs.process_ || state_ != rhs.state_ || children_.size() != rhs.children_.size() )
               return false;
            for( uint_t i = 0; i != children_.size(); ++i )
               if( *(children_[i]) != *(rhs.children_[i]) )
                  return false;
            return true;
         }
         bool operator!=( const Node& rhs ) const { return !operator==( rhs ); }

         void setChild( const uint_t index, shared_ptr< Node > const child )
         {
            if( children_.size() <= index )
               children_.resize( index+1 );
            children_[index] = child;
         }

         uint_t    process_;
         Set<SUID> state_;
         std::vector< shared_ptr< Node > >  children_;
      };

      void clear() { nodes_.clear(); }

      bool operator==( const BlockInformation& rhs ) const;
      bool operator!=( const BlockInformation& rhs ) const { return !operator==( rhs ); }

   public:

      BlockInformation( const BlockForest & forest ) : forest_( forest ) {}

      bool active() const { return !nodes_.empty(); }

      void getAllBlocks( std::vector< shared_ptr< IBlockID > >& blocks ) const;

      // blocks

      bool getProcess( uint_t& process, const real_t x, const real_t y, const real_t z ) const
         { const Node* node = getNode( x, y, z ); if( node ) { process = node->process_; return true; } return false; }
      bool getProcess( uint_t& process, const BlockID& id ) const
         { const Node* node = getNode( id ); if( node ) { process = node->process_; return true; } return false; }

      bool getState( Set<SUID>& state, const real_t x, const real_t y, const real_t z ) const
         { const Node* node = getNode( x, y, z ); if( node ) { state = node->state_; return true; } return false; }
      bool getState( Set<SUID>& state, const BlockID& id ) const
         { const Node* node = getNode( id ); if( node ) { state = node->state_; return true; } return false; }

      bool exists( const real_t x, const real_t y, const real_t z ) const { return getNode(x,y,z) != nullptr; }
      bool exists( const BlockID& id ) const { return getNode( id ) != nullptr; }

      bool existsRemotely( const real_t x, const real_t y, const real_t z ) const
         { const Node* node = getNode( x, y, z ); return ( node != nullptr && node->process_ != forest_.getProcess() ); }
      bool existsRemotely( const BlockID& id ) const
         { const Node* node = getNode( id ); return ( node != nullptr && node->process_ != forest_.getProcess() ); }

      bool getId( BlockID& id, const real_t x, const real_t y, const real_t z ) const;

      // root blocks

      bool getRootBlockProcess( uint_t& process, const uint_t x, const uint_t y, const uint_t z ) const
         { const Node* node = getRootNode(x,y,z); if( node ) { process = node->process_; return true; } return false; }

      bool getRootBlockState( Set<SUID>& state, const uint_t x, const uint_t y, const uint_t z ) const
         { const Node* node = getRootNode(x,y,z); if( node ) { state = node->state_; return true; } return false; }

      bool rootBlockExists( const uint_t x, const uint_t y, const uint_t z ) const { return getRootNode(x,y,z) != nullptr; }

      bool rootBlockExistsRemotely( const uint_t x, const uint_t y, const uint_t z ) const
         { const Node* node = getRootNode(x,y,z); return ( node != nullptr && node->process_ != forest_.getProcess() ); }

   private:

      const Node * getNode( const real_t x, const real_t y, const real_t z ) const;
      const Node * getNode( const BlockID & id ) const;

      const Node * getRootNode( const uint_t x, const uint_t y, const uint_t z ) const {
         const uint_t index =  z * forest_.getYSize() * forest_.getXSize() + y * forest_.getXSize() + x;
         if( index >= nodes_.size() )
            return nullptr;
         return nodes_[ index ].get();
      }

      const BlockForest &                forest_;
      std::vector< shared_ptr< Node > >  nodes_;
   };



   BlockForest( const uint_t process, const SetupBlockForest& forest, const bool keepGlobalBlockInformation = false );
   BlockForest( const uint_t process, const char* const filename, const bool broadcastFile = true, const bool keepGlobalBlockInformation = false );

   ~BlockForest() override = default;

   uint_t getProcess()        const { return process_; }
   uint_t getProcessIdBytes() const { return processIdBytes_; }

   uint_t getXSize() const { return size_[0]; } ///< number of coarse blocks on the initial grid (= number of octree root blocks) in x direction
   uint_t getYSize() const { return size_[1]; } ///< number of coarse blocks on the initial grid (= number of octree root blocks) in y direction
   uint_t getZSize() const { return size_[2]; } ///< number of coarse blocks on the initial grid (= number of octree root blocks) in z direction
   uint_t getSize( const uint_t index ) const { WALBERLA_ASSERT_LESS( index, 3 ); return size_[index]; }

   real_t getRootBlockXSize() const { return ( domain_.xMax() - domain_.xMin() ) / real_c( size_[0] ); }
   real_t getRootBlockYSize() const { return ( domain_.yMax() - domain_.yMin() ) / real_c( size_[1] ); }
   real_t getRootBlockZSize() const { return ( domain_.zMax() - domain_.zMin() ) / real_c( size_[2] ); }

   bool storesUniformBlockGrid() const { return depth_ == 0; }

   uint_t getTreeIdDigits()   const { return treeIdDigits_; }
   uint_t getBlockIdBytes()   const { uint_t const bits = treeIdDigits_ + 3 * depth_; return (bits >> 3) + (( bits & 7 ) ? uint_c(1) : uint_c(0)); }
   
   uint_t getDepth()          const { return depth_; }
   uint_t getNumberOfLevels() const { return depth_ + uint_t(1); }
   
   bool limitedDepth() const
   {
#ifdef WALBERLA_BLOCKFOREST_PRIMITIVE_BLOCKID
      return true;
#else
      return false;
#endif
   }
   bool limitedLevels() const { return limitedDepth(); }
   
   uint_t getMaxDepth() const
   {
#ifdef WALBERLA_BLOCKFOREST_PRIMITIVE_BLOCKID
      return ( math::UINT_BITS - treeIdDigits_ ) / uint_t(3);
#else
      return std::numeric_limits< uint_t >::max();
#endif
   }
   uint_t getMaxLevels() const { return limitedLevels() ? (getMaxDepth() + uint_t(1)) : std::numeric_limits< uint_t >::max(); }



   inline uint_t getNumberOfBlocks() const { WALBERLA_ASSERT_EQUAL( BlockStorage::getNumberOfBlocks(), blocks_.size() ); return blocks_.size(); }
   inline uint_t getNumberOfBlocks( const uint_t level ) const;

   inline const std::map< BlockID, shared_ptr< Block > > & getBlockMap() const { return blocks_; }
   
   using BlockStorage::getBlocks;

   inline void getBlocks( std::vector< const Block* >& blocks, const uint_t level ) const;
   inline void getBlocks( std::vector<       Block* >& blocks, const uint_t level );

   inline void getBlocksContainedWithinAABB( std::vector< const IBlock* >& blocks, const AABB& aabb ) const override;
   inline void getBlocksContainedWithinAABB( std::vector<       IBlock* >& blocks, const AABB& aabb ) override;

   inline void getBlocksOverlappedByAABB( std::vector< const IBlock* >& blocks, const AABB& aabb ) const override;
   inline void getBlocksOverlappedByAABB( std::vector<       IBlock* >& blocks, const AABB& aabb ) override;

   using BlockStorage::getBlock;

   inline const Block* getBlock( const IBlockID& id ) const override;
   inline       Block* getBlock( const IBlockID& id ) override;

   inline const Block* getBlock( const real_t x, const real_t y, const real_t z ) const override;
   inline       Block* getBlock( const real_t x, const real_t y, const real_t z ) override;

   inline const Block* getRootBlock( const uint_t x, const uint_t y, const uint_t z ) const;
   inline       Block* getRootBlock( const uint_t x, const uint_t y, const uint_t z );



   bool containsGlobalBlockInformation() const override { return blockInformation_->active(); }

   inline void getAllBlocks( std::vector< shared_ptr< IBlockID > >& blocks ) const override;

   inline bool blockExists        ( const real_t x, const real_t y, const real_t z ) const override;
   inline bool blockExistsLocally ( const real_t x, const real_t y, const real_t z ) const override;
   inline bool blockExistsRemotely( const real_t x, const real_t y, const real_t z ) const override;

   inline bool blockExists        ( const IBlockID& id ) const override;
   inline bool blockExistsLocally ( const IBlockID& id ) const override;
   inline bool blockExistsRemotely( const IBlockID& id ) const override;

   inline bool rootBlockExists        ( const uint_t x, const uint_t y, const uint_t z ) const;
   inline bool rootBlockExistsLocally ( const uint_t x, const uint_t y, const uint_t z ) const;
   inline bool rootBlockExistsRemotely( const uint_t x, const uint_t y, const uint_t z ) const;

   void getBlockID( IBlockID& id, const real_t x, const real_t y, const real_t z ) const override;
   void getAABB       ( AABB&      aabb,  const IBlockID& id ) const override;
   void getState      ( Set<SUID>& state, const IBlockID& id ) const override;
   void getProcessRank( uint_t&    rank,  const IBlockID& id ) const override;

   void getRootBlockAABB       ( AABB&      aabb,  const uint_t x, const uint_t y, const uint_t z ) const;
   void getRootBlockState      ( Set<SUID>& state, const uint_t x, const uint_t y, const uint_t z ) const;
   void getRootBlockProcessRank( uint_t&    rank,  const uint_t x, const uint_t y, const uint_t z ) const;

   const BlockInformation & getBlockInformation() const { return *blockInformation_; }


   inline uint_t getLevel( const IBlock& block ) const override;
   inline uint_t getLevelFromBlockId( const BlockID& id ) const;
   inline uint_t getAABBFromBlockId( AABB& aabb, const BlockID& id ) const;
   inline AABB   getAABBFromBlockId( const BlockID& id ) const;

          void getForestCoordinates   ( uint_t& x, uint_t& y, uint_t& z, const BlockID& id ) const;
          void getRootBlockCoordinates( uint_t& x, uint_t& y, uint_t& z, const BlockID& id ) const { return getForestCoordinates( x, y, z, id ); }
   inline void getRootBlockID( BlockID& id, const uint_t x, const uint_t y, const uint_t z ) const;



   bool insertBuffersIntoProcessNetwork() const { return insertBuffersIntoProcessNetwork_; }

   const std::vector< uint_t > & getNeighborhood() const { return neighborhood_; }
   const std::vector< uint_t > & getNeighboringProcesses() const override { return getNeighborhood(); }

   std::map< uint_t, std::vector< Vector3<real_t> > > getNeighboringProcessOffsets() const override;


   
   internal::BlockDataHandlingAdder addBlockData( const std::string & identifier = std::string() ) { return internal::BlockDataHandlingAdder( *this, identifier ); }

   template< typename T >
   inline BlockDataID addBlockData( const shared_ptr< T > & dataHandling,
                                    const std::string & identifier          = std::string(),
                                    const Set<SUID> & requiredSelectors     = Set<SUID>::emptySet(),
                                    const Set<SUID> & incompatibleSelectors = Set<SUID>::emptySet() );
                                    
   template< typename T >
   inline BlockDataID addBlockData( std::function< T* ( IBlock* const block ) > function,
                                    const std::string & identifier          = std::string(),
                                    const Set<SUID> & requiredSelectors     = Set<SUID>::emptySet(),
                                    const Set<SUID> & incompatibleSelectors = Set<SUID>::emptySet() )
   { return BlockStorage::addBlockData( function, identifier, requiredSelectors, incompatibleSelectors ); }
                                    
   BlockDataID addBlockData( const domain_decomposition::internal::SelectableBlockDataHandlingWrapper & dataHandling, const std::string& identifier = std::string() )
   { return BlockStorage::addBlockData( dataHandling, identifier ); }
   
   template< typename T >
   inline BlockDataID loadBlockData( const std::string & file, const shared_ptr< T > & dataHandling,
                                     const std::string & identifier          = std::string(),
                                     const Set<SUID> & requiredSelectors     = Set<SUID>::emptySet(),
                                     const Set<SUID> & incompatibleSelectors = Set<SUID>::emptySet() );
                                    
   BlockDataID loadBlockData( const std::string & file,
                              const domain_decomposition::internal::SelectableBlockDataHandlingWrapper & dataHandling, const std::string & identifier = std::string() )
   { return BlockStorage::loadBlockData( file, dataHandling, identifier ); }


   // AMR Pipeline
   // 1) distributed block level adjustment (callbacks are called that determine if a block level changes)
   //       - how to register this functor? setRefreshMinTargetLevelDeterminationFunction()
   //       - can be switched on/off via function recalculateBlockLevelsInRefresh()  (is on by default)
   // 2) create proxy/phantom blocks which are moved around instead of the real blocks because cheaper
   // 3) dynamic load balancing (using phantom blocks)
   //       - load balancing strategy, can be configured via setRefreshPhantomBlockMigrationPreparationFunction()
   //         e.g. diffusive or space filling curve strategy
   // 4) data migration / refinement / coarsening

   /// Triggers AMR Pipeline
   void refresh();

   /// Functor that calls refresh with given frequency
   RefreshFunctor getRefreshFunctor( const uint_t checkFrequency = uint_t(1) ) { return RefreshFunctor( *this, checkFrequency ); }
   
   /// Modification stamp is changed when refresh moves a block or refines/coarsens at least one block
   /// however, stamp may (in rare cases) also change if block structure was not altered
   /// i.e. if modification stamp stays the same the block structure was not changed (no block movements, no level changes)
   uint_t getModificationStamp() const { return modificationStamp_; }

   /// Returns timing pool that contains measurement for sub-steps of AMR pipeline
   const WcTimingPool & getRefreshTiming() const { return refreshTiming_; }
   void clearRefreshTiming() { refreshTiming_.clear(); }

   /// switches AMR pipeline step 1 on/off
   bool recalculateBlockLevelsInRefresh() const { return recalculateBlockLevelsInRefresh_; }
   void recalculateBlockLevelsInRefresh( const bool r ) { recalculateBlockLevelsInRefresh_ = r; }

   /// by default rebalancing can be skipped if no block level has changed
   /// with this switch one can force the execution of the rebalancing step
   void alwaysRebalanceInRefresh( const bool b ) { alwaysRebalanceInRefresh_ = b; }
   bool alwaysRebalanceInRefresh() const { return alwaysRebalanceInRefresh_; }

   /// callback for AMR pipeline step 1
   /// guaranteed minimal block level
   void setRefreshMinTargetLevelDeterminationFunction( const RefreshMinTargetLevelDeterminationFunction & f ) { refreshMinTargetLevelDeterminationFunction_ = f; }

   /// allow multiple AMR passes (possibility to refine more than once within one time step)
   /// if true, one all-to-all reduction with a boolean value is performed during refresh
   bool allowMultipleRefreshCycles() const { return allowMultipleRefreshCycles_; }
   void allowMultipleRefreshCycles( const bool m ) { allowMultipleRefreshCycles_ = m; }

   /// call setRefreshMinTargetLevelDeterminationFunction callback again after AMR internal refinement to keep 2:1 balance
   /// if true, one all-to-all reduction with a boolean value is performed during refresh
   bool reevaluateMinTargetLevelsAfterForcedRefinement() const { return reevaluateMinTargetLevelsAfterForcedRefinement_; }
   void reevaluateMinTargetLevelsAfterForcedRefinement( const bool c ) { reevaluateMinTargetLevelsAfterForcedRefinement_ = c; }

   /// exit pipeline after first AMR step if no block level has changed after setRefreshMinTargetLevelDeterminationFunction callback
   /// can be overwritten by alwaysRebalanceInRefresh
   /// if true, one all-to-all reduction with a boolean value is performed during refresh
   bool checkForEarlyOutInRefresh() const { return checkForEarlyOutInRefresh_; }
   void checkForEarlyOutInRefresh( const bool c ) { checkForEarlyOutInRefresh_ = c; }

   /// exit pipeline after first AMR step if no block level has changed
   /// can be overwritten by alwaysRebalanceInRefresh
   /// if true, one all-to-all reduction with a boolean value is performed during refresh
   bool checkForLateOutInRefresh() const { return checkForLateOutInRefresh_; }
   void checkForLateOutInRefresh( const bool c ) { checkForLateOutInRefresh_ = c; }
   
   /// maximal depth of block structure is allowed to change?
   /// if true, one all-to-all reduction with an unsigned int value is performed after refresh
   bool allowRefreshChangingDepth() const { return allowChangingDepth_; }
   void allowRefreshChangingDepth( const bool c ) { allowChangingDepth_ = c; }
   
   /// if true, one all-to-all reduction with a boolean value is performed during refresh
   bool checkForEarlyOutAfterLoadBalancing() const { return checkForEarlyOutAfterLoadBalancing_; }
   void checkForEarlyOutAfterLoadBalancing( const bool c ) { checkForEarlyOutAfterLoadBalancing_ = c; }

   /// callback which determines the state id (SUID) of the new block during refinement/coarsening
   void setRefreshBlockStateDeterminationFunction( const PhantomBlockForest::BlockStateDeterminationFunction & f ) { refreshBlockStateDeterminationFunction_ = f; }
   /// callback to assign arbitrary data to a phantom block (used only during load balancing), e.g. weights
   void setRefreshPhantomBlockDataAssignmentFunction( const PhantomBlockForest::PhantomBlockDataAssignmentFunction & f ) { refreshPhantomBlockDataAssignmentFunction_ = f; }
   /// load balancing algorithm
   void setRefreshPhantomBlockMigrationPreparationFunction( const PhantomBlockForest::MigrationPreparationFunction & f ) { refreshPhantomBlockMigrationPreparationFunction_ = f; }
   void setRefreshPhantomBlockDataPackFunction( const PhantomBlockForest::PhantomBlockDataPackFunction & f ) { refreshPhantomBlockDataPackFunction_ = f; }
   void setRefreshPhantomBlockDataUnpackFunction( const PhantomBlockForest::PhantomBlockDataUnpackFunction & f ) { refreshPhantomBlockDataUnpackFunction_ = f; }

   inline bool loadBalancingFunctionRegistered() const { return static_cast<bool>(refreshPhantomBlockMigrationPreparationFunction_); }
   /// get number of "setRefreshPhantomBlockMigrationPreparationFunction" calls
   inline uint_t phantomBlockMigrationIterations() const { return phantomBlockMigrationIterations_; }

   inline uint_t    addRefreshCallbackFunctionBeforeBlockDataIsPacked( const RefreshCallbackFunction & f );
   inline void   removeRefreshCallbackFunctionBeforeBlockDataIsPacked( const uint_t handle );

   inline uint_t    addRefreshCallbackFunctionBeforeBlockDataIsUnpacked( const RefreshCallbackFunction & f );
   inline void   removeRefreshCallbackFunctionBeforeBlockDataIsUnpacked( const uint_t handle );

   inline uint_t    addRefreshCallbackFunctionAfterBlockDataIsUnpacked( const RefreshCallbackFunction & f );
   inline void   removeRefreshCallbackFunctionAfterBlockDataIsUnpacked( const uint_t handle );



   void  createSnapshot( const std::vector<uint_t> & sendTo, const std::vector<uint_t> & recvFrom );
   void restoreSnapshot( const SnapshotRestoreFunction & processMapping, const bool rebelance = true );

   SnapshotCreationFunctor getSnapshotCreationFunctor( const SnapshotCreationFunction & function,
                                                       const uint_t checkFrequency = uint_t(1) ) { return SnapshotCreationFunctor( *this, function, checkFrequency ); }

   inline uint_t    addCallbackFunctionAfterBlockDataIsRestored( const SnapshotRestoreCallbackFunction & f );
   inline void   removeCallbackFunctionAfterBlockDataIsRestored( const uint_t handle );



   /// Block states are reduced using MPI
   void saveToFile( const std::string & filename, FileIOMode fileIOMode = MPI_PARALLEL ) const;
   /// Block states are passed by argument (ATTENTION: 'blockStates' must be identical for every process!)
   void saveToFile( const std::string & filename, const Set<SUID> & blockStates, FileIOMode fileIOMode = MPI_PARALLEL ) const;

protected:

   bool equal( const BlockStorage* rhs ) const override;
   
   void addBlockData( IBlock * const block, const BlockDataID & index, domain_decomposition::internal::BlockData * const data )
   { BlockStorage::addBlockData( block, index, data ); }

private:

   void constructBlockInformation( const SetupBlockForest & forest );
   void constructBlockInformation( const std::vector< BlockID > & ids, const std::vector< shared_ptr< BlockInformation::Node > > & nodes );
   void constructBlockInformation();

   void registerRefreshTimer();
   bool determineBlockTargetLevels( bool & additionalRefreshCycleRequired, bool & rerun );
   void update( PhantomBlockForest & phantomForest );

   void saveToFile( const std::string & filename, FileIOMode fileIOMode,
                    const std::map< SUID, std::vector< bool > > & suidMap, const uint_t suidBytes ) const;
   void storeFileHeader( std::vector< uint8_t > & data, uint_t & offset ) const;



   uint_t process_;
   uint_t processIdBytes_;

   uint_t size_[3];     // number of coarse blocks on the initial grid (= number of octree root blocks) in each direction
   uint_t depth_;       // depth := number of levels - 1
   uint_t treeIdDigits_;

   std::map< BlockID, shared_ptr< Block > > blocks_;

   bool insertBuffersIntoProcessNetwork_;
   std::vector< uint_t > neighborhood_; // neighbor processes (not entirely reconstructable from 'blocks_' -> empty buffer processes!)

   shared_ptr< BlockInformation > blockInformation_;



   uint_t modificationStamp_;
   WcTimingPool refreshTiming_;
   
   bool recalculateBlockLevelsInRefresh_;
   bool alwaysRebalanceInRefresh_;

   RefreshMinTargetLevelDeterminationFunction refreshMinTargetLevelDeterminationFunction_;

   bool allowMultipleRefreshCycles_;                     // if true, one all-to-all reduction with a boolean value is performed during refresh
   bool reevaluateMinTargetLevelsAfterForcedRefinement_; // if true, one all-to-all reduction with a boolean value is performed during refresh
   bool checkForEarlyOutInRefresh_;                      // if true, one all-to-all reduction with a boolean value is performed during refresh
   bool checkForLateOutInRefresh_;                       // if true, one all-to-all reduction with a boolean value is performed during refresh
   bool allowChangingDepth_;                             // if true, one all-to-all reduction with an unsigned int value is performed after refresh
   bool checkForEarlyOutAfterLoadBalancing_;             // if true, one all-to-all reduction with a boolean value is performed during refresh
   
   PhantomBlockForest::BlockStateDeterminationFunction refreshBlockStateDeterminationFunction_;
   PhantomBlockForest::PhantomBlockDataAssignmentFunction refreshPhantomBlockDataAssignmentFunction_;
   PhantomBlockForest::MigrationPreparationFunction refreshPhantomBlockMigrationPreparationFunction_;
   PhantomBlockForest::PhantomBlockDataPackFunction refreshPhantomBlockDataPackFunction_;
   PhantomBlockForest::PhantomBlockDataUnpackFunction refreshPhantomBlockDataUnpackFunction_;
   
   uint_t phantomBlockMigrationIterations_;

   std::map< uint_t, RefreshCallbackFunction > callbackBeforeBlockDataIsPacked_;
   std::map< uint_t, RefreshCallbackFunction > callbackBeforeBlockDataIsUnpacked_;
   std::map< uint_t, RefreshCallbackFunction > callbackAfterBlockDataIsUnpacked_;

   uint_t nextCallbackBeforeBlockDataIsPackedHandle_;
   uint_t nextCallbackBeforeBlockDataIsUnpackedHandle_;
   uint_t nextCallbackAfterBlockDataIsUnpackedHandle_;



   bool snapshotExists_;

   uint_t snapshotDepth_;
   uint_t snapshotBlockDataItems_;

   std::map< uint_t, mpi::RecvBuffer > snapshot_;

   std::map< uint_t, SnapshotRestoreCallbackFunction > callbackAfterBlockDataIsRestored_;
   uint_t nextCallbackAfterBlockDataIsRestoredHandle_;



#ifndef NDEBUG
public:
   void checkBlockInformationConsistency( const SetupBlockForest& forest ) const;
#endif
};



inline uint_t BlockForest::getNumberOfBlocks( const uint_t level ) const {

   uint_t numberOfBlocks = 0;

   for( auto it = blocks_.begin(); it != blocks_.end(); ++it )
      if( it->second->getLevel() == level ) ++numberOfBlocks;

   return numberOfBlocks;
}



inline void BlockForest::getBlocks( std::vector< const Block* >& blocks, const uint_t level ) const {

   for( auto it = blocks_.begin(); it != blocks_.end(); ++it )
      if( it->second->getLevel() == level ) blocks.push_back( it->second.get() );
}



inline void BlockForest::getBlocks( std::vector< Block* >& blocks, const uint_t level ) {

   for( auto it = blocks_.begin(); it != blocks_.end(); ++it )
      if( it->second->getLevel() == level ) blocks.push_back( it->second.get() );
}



inline void BlockForest::getBlocksContainedWithinAABB( std::vector< const IBlock* >& blocks, const AABB& aabb ) const {

   for( auto it = blocks_.begin(); it != blocks_.end(); ++it )
      if( aabb.contains( it->second->getAABB() ) ) blocks.push_back( it->second.get() );
}



inline void BlockForest::getBlocksContainedWithinAABB( std::vector< IBlock* >& blocks, const AABB& aabb ) {

   for( auto it = blocks_.begin(); it != blocks_.end(); ++it )
      if( aabb.contains( it->second->getAABB() ) ) blocks.push_back( it->second.get() );
}



inline void BlockForest::getBlocksOverlappedByAABB( std::vector< const IBlock* >& blocks, const AABB& aabb ) const {

   for( auto it = blocks_.begin(); it != blocks_.end(); ++it )
      if( it->second->getAABB().intersects( aabb ) ) blocks.push_back( it->second.get() );
}



inline void BlockForest::getBlocksOverlappedByAABB( std::vector< IBlock* >& blocks, const AABB& aabb ) {

   for( auto it = blocks_.begin(); it != blocks_.end(); ++it )
      if( it->second->getAABB().intersects( aabb ) ) blocks.push_back( it->second.get() );
}



inline const Block* BlockForest::getBlock( const IBlockID& id ) const {

   WALBERLA_ASSERT_EQUAL( dynamic_cast< const BlockID* >( &id ), &id )

   auto it = blocks_.find( *static_cast< const BlockID* >( &id ) );

   if( it != blocks_.end() )
      return it->second.get();

   return nullptr;
}



inline Block* BlockForest::getBlock( const IBlockID& id ) {

   WALBERLA_ASSERT_EQUAL( dynamic_cast< const BlockID* >( &id ), &id )

   auto it = blocks_.find( *static_cast< const BlockID* >( &id ) );

   if( it != blocks_.end() )
      return it->second.get();

   return nullptr;
}



inline const Block* BlockForest::getBlock( const real_t x, const real_t y, const real_t z ) const {

   for( auto it = blocks_.begin(); it != blocks_.end(); ++it )
      if( it->second->getAABB().contains(x,y,z) ) return it->second.get();

   return nullptr;
}



inline Block* BlockForest::getBlock( const real_t x, const real_t y, const real_t z ) {

   for( auto it = blocks_.begin(); it != blocks_.end(); ++it )
      if( it->second->getAABB().contains(x,y,z) ) return it->second.get();

   return nullptr;
}



inline const Block* BlockForest::getRootBlock( const uint_t x, const uint_t y, const uint_t z ) const {

   BlockID id;
   getRootBlockID( id, x, y, z );

   return getBlock( id );
}



inline Block* BlockForest::getRootBlock( const uint_t x, const uint_t y, const uint_t z ) {

   BlockID id;
   getRootBlockID( id, x, y, z );

   return getBlock( id );
}



inline void BlockForest::getAllBlocks( std::vector< shared_ptr< IBlockID > >& blocks ) const {

   if( blockInformation_->active() )
      return blockInformation_->getAllBlocks( blocks );

   for( auto block = begin(); block != end(); ++block ) {
      const IBlockID& id = block->getId();
      WALBERLA_ASSERT_EQUAL( dynamic_cast< const BlockID* >( &id ), &id )
      blocks.push_back( walberla::make_shared< BlockID >( *static_cast< const BlockID* >( &id ) ) );
   }
}



inline bool BlockForest::blockExists( const real_t x, const real_t y, const real_t z ) const {

   if( blockInformation_->active() )
      return blockInformation_->exists(x,y,z);

   return getBlock(x,y,z) != nullptr;
}



inline bool BlockForest::blockExistsLocally( const real_t x, const real_t y, const real_t z ) const {

   return getBlock(x,y,z) != nullptr;
}



inline bool BlockForest::blockExistsRemotely( const real_t x, const real_t y, const real_t z ) const {

   if( blockInformation_->active() )
      return blockInformation_->existsRemotely(x,y,z);

   return getBlock(x,y,z) == nullptr;
}



inline bool BlockForest::blockExists( const IBlockID& id ) const {

   WALBERLA_ASSERT_EQUAL( dynamic_cast< const BlockID* >( &id ), &id )

   if( blockInformation_->active() )
      return blockInformation_->exists( *static_cast< const BlockID* >( &id ) );

   return getBlock( id ) != nullptr;
}



inline bool BlockForest::blockExistsLocally( const IBlockID& id ) const {

   return getBlock( id ) != nullptr;
}



inline bool BlockForest::blockExistsRemotely( const IBlockID& id ) const {

   WALBERLA_ASSERT_EQUAL( dynamic_cast< const BlockID* >( &id ), &id )

   if( blockInformation_->active() )
      return blockInformation_->existsRemotely( *static_cast< const BlockID* >( &id ) );

   return getBlock( id ) == nullptr;
}



inline bool BlockForest::rootBlockExists( const uint_t x, const uint_t y, const uint_t z ) const {

   if( blockInformation_->active() )
      return blockInformation_->rootBlockExists(x,y,z);

   return getRootBlock(x,y,z) != nullptr;
}



inline bool BlockForest::rootBlockExistsLocally( const uint_t x, const uint_t y, const uint_t z ) const {

   return getRootBlock(x,y,z) != nullptr;
}



inline bool BlockForest::rootBlockExistsRemotely( const uint_t x, const uint_t y, const uint_t z ) const {

   if( blockInformation_->active() )
      return blockInformation_->rootBlockExistsRemotely(x,y,z);

   return getRootBlock(x,y,z) == nullptr;
}


inline uint_t BlockForest::getLevel( const IBlock& block ) const {

   WALBERLA_ASSERT_EQUAL( dynamic_cast< const Block* >( &block ), &block )

   return static_cast< const Block* >( &block )->getLevel();
}



inline uint_t BlockForest::getLevelFromBlockId( const BlockID& id ) const {

   WALBERLA_ASSERT_GREATER_EQUAL( id.getUsedBits(), treeIdDigits_ )
   WALBERLA_ASSERT_EQUAL( ( id.getUsedBits() - treeIdDigits_ ) % 3, 0 )

   return ( id.getUsedBits() - treeIdDigits_ ) / 3;
}



inline uint_t BlockForest::getAABBFromBlockId( AABB& aabb, const BlockID& id ) const {

   BlockReconstruction::AABBReconstruction const aabbReconstruction( domain_, size_[0], size_[1], size_[2], treeIdDigits_ );

   return aabbReconstruction( aabb, id );
}



inline AABB BlockForest::getAABBFromBlockId( const BlockID& id ) const {

   AABB aabb;
   getAABBFromBlockId( aabb, id );
   return aabb;
}



inline void BlockForest::getRootBlockID( BlockID& id, const uint_t x, const uint_t y, const uint_t z ) const {

   WALBERLA_ASSERT_LESS( x, size_[0] )
   WALBERLA_ASSERT_LESS( y, size_[1] )
   WALBERLA_ASSERT_LESS( z, size_[2] )

   id.clear();
   id.init( z * size_[1] * size_[0] + y * size_[0] + x, uint_c(1) << ( treeIdDigits_ - 1 ) );
}



template< typename T >
inline BlockDataID BlockForest::addBlockData( const shared_ptr< T > & dataHandling, const std::string & identifier,
                                              const Set<SUID> & requiredSelectors, const Set<SUID> & incompatibleSelectors )
{
   //static_assert( std::is_base_of< BlockDataHandling<typename T::value_type>, T >::value );

   auto downcast = dynamic_pointer_cast< blockforest::BlockDataHandling<typename T::value_type> >( dataHandling );

   if( downcast )
   {
      domain_decomposition::internal::SelectableBlockDataHandlingWrapper sbdhw(
               walberla::make_shared< internal::BlockDataHandlingHelper<typename T::value_type> >( downcast ),
               requiredSelectors, incompatibleSelectors, identifier );

      return addBlockData( sbdhw, identifier );
   }

   domain_decomposition::internal::SelectableBlockDataHandlingWrapper sbdhw(
            walberla::make_shared< domain_decomposition::internal::BlockDataHandlingHelper<typename T::value_type> >( dataHandling ),
            requiredSelectors, incompatibleSelectors, identifier );

   return addBlockData( sbdhw, identifier );
}



template< typename T >
inline BlockDataID BlockForest::loadBlockData( const std::string & file, const shared_ptr< T > & dataHandling, const std::string & identifier,
                                               const Set<SUID> & requiredSelectors, const Set<SUID> & incompatibleSelectors )
{
   auto downcast = dynamic_pointer_cast< blockforest::BlockDataHandling<typename T::value_type> >( dataHandling );

   if( downcast )
   {
      domain_decomposition::internal::SelectableBlockDataHandlingWrapper sbdhw(
               walberla::make_shared< internal::BlockDataHandlingHelper<typename T::value_type> >( downcast ),
               requiredSelectors, incompatibleSelectors, identifier );

      return loadBlockData( file, sbdhw, identifier );
   }

   domain_decomposition::internal::SelectableBlockDataHandlingWrapper sbdhw(
            walberla::make_shared< domain_decomposition::internal::BlockDataHandlingHelper<typename T::value_type> >( dataHandling ),
            requiredSelectors, incompatibleSelectors, identifier );

   return loadBlockData( file, sbdhw, identifier );
}



inline uint_t BlockForest::addRefreshCallbackFunctionBeforeBlockDataIsPacked( const RefreshCallbackFunction & f )
{
   callbackBeforeBlockDataIsPacked_.insert( callbackBeforeBlockDataIsPacked_.end(), std::make_pair( nextCallbackBeforeBlockDataIsPackedHandle_, f ) );
   ++nextCallbackBeforeBlockDataIsPackedHandle_;
   return nextCallbackBeforeBlockDataIsPackedHandle_ - uint_t(1);
}

inline void BlockForest::removeRefreshCallbackFunctionBeforeBlockDataIsPacked( const uint_t handle )
{
   callbackBeforeBlockDataIsPacked_.erase( handle );
}



inline uint_t BlockForest::addRefreshCallbackFunctionBeforeBlockDataIsUnpacked( const RefreshCallbackFunction & f )
{
   callbackBeforeBlockDataIsUnpacked_.insert( callbackBeforeBlockDataIsUnpacked_.end(), std::make_pair( nextCallbackBeforeBlockDataIsUnpackedHandle_, f ) );
   ++nextCallbackBeforeBlockDataIsUnpackedHandle_;
   return nextCallbackBeforeBlockDataIsUnpackedHandle_ - uint_t(1);
}

inline void BlockForest::removeRefreshCallbackFunctionBeforeBlockDataIsUnpacked( const uint_t handle )
{
   callbackBeforeBlockDataIsUnpacked_.erase( handle );
}



inline uint_t BlockForest::addRefreshCallbackFunctionAfterBlockDataIsUnpacked( const RefreshCallbackFunction & f )
{
   callbackAfterBlockDataIsUnpacked_.insert( callbackAfterBlockDataIsUnpacked_.end(), std::make_pair( nextCallbackAfterBlockDataIsUnpackedHandle_, f ) );
   ++nextCallbackAfterBlockDataIsUnpackedHandle_;
   return nextCallbackAfterBlockDataIsUnpackedHandle_ - uint_t(1);
}

inline void BlockForest::removeRefreshCallbackFunctionAfterBlockDataIsUnpacked( const uint_t handle )
{
   callbackAfterBlockDataIsUnpacked_.erase( handle );
}



inline uint_t BlockForest::addCallbackFunctionAfterBlockDataIsRestored( const SnapshotRestoreCallbackFunction & f )
{
   callbackAfterBlockDataIsRestored_.insert( callbackAfterBlockDataIsRestored_.end(), std::make_pair( nextCallbackAfterBlockDataIsRestoredHandle_, f ) );
   ++nextCallbackAfterBlockDataIsRestoredHandle_;
   return nextCallbackAfterBlockDataIsRestoredHandle_ - uint_t(1);
}

inline void BlockForest::removeCallbackFunctionAfterBlockDataIsRestored( const uint_t handle )
{
   callbackAfterBlockDataIsRestored_.erase( handle );
}



////////////////////
// Useful helpers //
////////////////////



class MinTargetLevelDeterminationFunctions
{
public:

   using MinTargetLevelDeterminationFunction = blockforest::BlockForest::RefreshMinTargetLevelDeterminationFunction;

   void add( const MinTargetLevelDeterminationFunction & function )
   {
      functions_.push_back( function );
   }

   void operator()( std::vector< std::pair< const Block *, uint_t > > & minTargetLevels,
                    std::vector< const Block * > & blocksAlreadyMarkedForRefinement,
                    const blockforest::BlockForest & forest )
   {
      for( auto function = functions_.begin(); function != functions_.end(); ++function )
         (*function)( minTargetLevels, blocksAlreadyMarkedForRefinement, forest );
   }

private:

   std::vector< MinTargetLevelDeterminationFunction > functions_;

}; // class MinTargetLevelDeterminationFunctions



class CombinedMinTargetLevelDeterminationFunctions
{
public:

   using MinTargetLevelDeterminationFunction = blockforest::BlockForest::RefreshMinTargetLevelDeterminationFunction;

   CombinedMinTargetLevelDeterminationFunctions(const std::function<uint_t(const std::vector<uint_t> &)> & targetLevelReductionFct = [](const std::vector<uint_t> & t){ return *std::max_element(t.begin(), t.end());})
   : targetLevelReductionFct_( targetLevelReductionFct )
   {

   }

   void add( const MinTargetLevelDeterminationFunction & function )
   {
      functions_.push_back( function );
   }

   void operator()( std::vector< std::pair< const Block *, uint_t > > & minTargetLevels,
                    std::vector< const Block * > & blocksAlreadyMarkedForRefinement,
                    const blockforest::BlockForest & forest )
   {
      const uint_t numberOfBlocks = minTargetLevels.size();

      std::vector< std::vector< std::pair< const Block *, uint_t > > > minTargetLevelsPerFunction(functions_.size(), minTargetLevels);

      // evaluate the different determination functions
      auto iter = minTargetLevelsPerFunction.begin();
      for( auto function = functions_.begin(); function != functions_.end(); ++function, ++iter )
      {
         (*function)( *iter, blocksAlreadyMarkedForRefinement, forest );
         WALBERLA_ASSERT_EQUAL(iter->size(), numberOfBlocks, "Number of blocks has changed during min target level determination!")
      }

      // combine the outcome of the different functions into a single target level
      std::vector<uint_t> targetLevels(functions_.size());
      for( uint_t block = 0; block < numberOfBlocks; ++block )
      {
         for( uint_t fct = 0; fct < functions_.size(); ++fct)
         {
            WALBERLA_ASSERT_EQUAL(minTargetLevelsPerFunction[fct][block].first->getId(), minTargetLevels[block].first->getId())
            targetLevels[fct] = minTargetLevelsPerFunction[fct][block].second;
         }
         minTargetLevels[block].second = targetLevelReductionFct_(targetLevels);
      }

   }

private:

   std::vector< MinTargetLevelDeterminationFunction > functions_;
   std::function<uint_t(const std::vector<uint_t> &)> targetLevelReductionFct_;

}; // class CombinedMinTargetLevelDeterminationFunctions



} // namespace blockforest

using blockforest::BlockForest;

} // namespace walberla
