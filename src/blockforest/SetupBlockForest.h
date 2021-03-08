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
//! \file SetupBlockForest.h
//! \ingroup blockforest
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "GlobalLoadBalancing.h"
#include "SetupBlock.h"
#include "Types.h"

#include "core/NonCopyable.h"
#include "core/debug/Debug.h"
#include "core/selectable/SetSelectableObject.h"
#include "core/math/AABB.h"
#include "core/uid/SUID.h"

#include <functional>
#include <set>
#include <string>
#include <vector>


namespace walberla {
namespace blockforest {


class SetupBlockForest : private NonCopyable {

public:

   using TargetProcessAssignmentFunction = std::function<uint_t (SetupBlockForest &, const uint_t, const memory_t)>; // returns number of processes (may be lower than numberOfProcesses)



   class RootBlockAABB {
      friend class SetupBlockForest;
   public:

      AABB operator()( const uint_t index ) const; // index = treeIndex

   private:

      RootBlockAABB( const AABB& domain, const real_t rootBlockXSize, const real_t rootBlockYSize, const real_t rootBlockZSize,
                     const uint_t xSize, const uint_t ySize, const uint_t zSize ) : domain_( domain ),
         rootBlockXSize_( rootBlockXSize ), rootBlockYSize_( rootBlockYSize ), rootBlockZSize_( rootBlockZSize ),
         xSize_( xSize ), ySize_( ySize ), zSize_( zSize ) {}

      const AABB   domain_;
      const real_t rootBlockXSize_, rootBlockYSize_, rootBlockZSize_;
      const uint_t xSize_, ySize_, zSize_;
   };



   // Do not use a vector of bool's! Due to the implementation of this vector in the standard library, parallel access to a
   // vector of bool's - even on different elements - is not thread-safe!
   using RootBlockExclusionFunction = std::function<void (std::vector<uint8_t> &, const RootBlockAABB &)>;

   using RefinementSelectionFunction = std::function<void (SetupBlockForest &)>;
   using WorkloadMemorySUIDAssignmentFunction = std::function<void (SetupBlockForest &)>;

   using CommunicationPairs = std::vector<std::pair<const SetupBlock *, const SetupBlock *>>;
   using CommunicationWeights = std::vector<real_t>;
   using CommunicationWeightFunction = std::function<void (const CommunicationPairs &, CommunicationWeights &)>;

   inline static void NullCommunicationWeightFunction( const CommunicationPairs &, CommunicationWeights & )
   {
   }


   class const_iterator;
   friend class const_iterator;

   class iterator;
   friend class iterator;
   class iterator {

      friend class const_iterator;
      friend class SetupBlockForest;

   public:

      iterator( const iterator& it )  = default;

      iterator& operator++() { WALBERLA_ASSERT_NOT_NULLPTR( block_ ); block_ = forest_->getNextBlock( block_ ); return *this; } // prefix ++X
      iterator  operator++(int) { iterator it( *this ); operator++(); return it; };                                             // postfix X++

      bool operator==( const iterator& rhs ) const { return block_ == rhs.block_; }
      bool operator!=( const iterator& rhs ) const { return block_ != rhs.block_; }

      SetupBlock* get() { return block_; }

      SetupBlock& operator*()  { WALBERLA_ASSERT_NOT_NULLPTR( block_ ); return *block_; }
      SetupBlock* operator->() { WALBERLA_ASSERT_NOT_NULLPTR( block_ ); return  block_; }

   private:

      iterator( SetupBlockForest* const forest, SetupBlock* const block ) : forest_( forest ), block_( block ) {}

      SetupBlockForest* const forest_;
      SetupBlock* block_;
   };


   class const_iterator {

      friend class SetupBlockForest;

   public:

      const_iterator( const       iterator& it ) : forest_( it.forest_ ), block_( it.block_ ) {}
      const_iterator( const const_iterator& it )  = default;

      const_iterator& operator++() { WALBERLA_ASSERT_NOT_NULLPTR( block_ ); block_ = forest_->getNextBlock( block_ ); return *this; } // prefix ++X
      const_iterator  operator++(int) { const_iterator it( *this ); operator++(); return it; };                                       // postfix X++

      bool operator==( const const_iterator& rhs ) const { return block_ == rhs.block_; }
      bool operator!=( const const_iterator& rhs ) const { return block_ != rhs.block_; }

      const SetupBlock* get() const { return block_; }

      const SetupBlock& operator*()  const { WALBERLA_ASSERT_NOT_NULLPTR( block_ ); return *block_; }
      const SetupBlock* operator->() const { WALBERLA_ASSERT_NOT_NULLPTR( block_ ); return  block_; }

   private:

      const_iterator( const SetupBlockForest* const forest, const SetupBlock* const block ) : forest_( forest ), block_( block ) {}

      const SetupBlockForest* const forest_;
      const SetupBlock* block_;
   };



   inline  SetupBlockForest();
   inline ~SetupBlockForest();

   const AABB& getDomain() const { return domain_; }

   real_t getRootBlockXSize() const { return rootBlockSize_[0]; }
   real_t getRootBlockYSize() const { return rootBlockSize_[1]; }
   real_t getRootBlockZSize() const { return rootBlockSize_[2]; }
   real_t getRootBlockSize( const uint_t index ) const { WALBERLA_ASSERT_LESS( index, 3 ); return rootBlockSize_[index]; }

   uint_t getXSize() const { return size_[0]; }
   uint_t getYSize() const { return size_[1]; }
   uint_t getZSize() const { return size_[2]; }
   uint_t getSize( const uint_t index ) const { WALBERLA_ASSERT_LESS( index, 3 ); return size_[index]; }

   bool isXPeriodic() const { return periodic_[0]; }
   bool isYPeriodic() const { return periodic_[1]; }
   bool isZPeriodic() const { return periodic_[2]; }
   bool isPeriodic( const uint_t index ) const { WALBERLA_ASSERT_LESS( index, 3 ); return periodic_[index]; }

   uint_t getDepth()              const { return depth_; }
   uint_t getNumberOfLevels()     const { return depth_ + uint_t(1); }
   uint_t getMinLevel()           const;
   uint_t getMaxLevel()           const;
   uint_t getTreeIdDigits()       const { return treeIdDigits_; }
   uint_t getBlockIdBytes()       const { uint_t bits = treeIdDigits_ + 3 * depth_; return (bits >> 3) + (( bits & 7 ) ? uint_c(1) : uint_c(0)); }

   uint_t getNumberOfTrees()      const { return forest_.size(); }
   uint_t getNumberOfRootBlocks() const { return numberOfRootBlocks_; }
   uint_t getNumberOfBlocks()     const { return numberOfBlocks_; }
   uint_t getNumberOfBlocks( const uint_t level ) const;

   inline const_iterator begin() const;
   inline const_iterator end()   const { return const_iterator( this, nullptr ); }

   inline iterator begin();
   inline iterator end() { return iterator( this, nullptr ); }

   const SetupBlock* getFirstBlock() const;
         SetupBlock* getFirstBlock();

   const SetupBlock* getNextBlock( const SetupBlock* block ) const;
         SetupBlock* getNextBlock( const SetupBlock* block );

   const SetupBlock* getTree( const uint_t treeIndex ) const { return getRootBlock( treeIndex ); }
         SetupBlock* getTree( const uint_t treeIndex )       { return getRootBlock( treeIndex ); }

   const SetupBlock* getRootBlock( const uint_t treeIndex ) const { WALBERLA_ASSERT_LESS( treeIndex, forest_.size() ); return forest_[treeIndex]; }
         SetupBlock* getRootBlock( const uint_t treeIndex )       { WALBERLA_ASSERT_LESS( treeIndex, forest_.size() ); return forest_[treeIndex]; }

   inline const SetupBlock* getRootBlock( const uint_t x, const uint_t y, const uint_t z ) const;
   inline       SetupBlock* getRootBlock( const uint_t x, const uint_t y, const uint_t z );

   const SetupBlock* getBlock( const BlockID& id ) const;

   inline const SetupBlock* getBlock( const real_t px, const real_t py, const real_t pz ) const;
   inline       SetupBlock* getBlock( const real_t px, const real_t py, const real_t pz );

   void getBlocks( std::vector< const SetupBlock* >& blocks ) const;
   void getBlocks( std::vector<       SetupBlock* >& blocks );

   void getBlocks( std::vector< const SetupBlock* >& blocks, const uint_t level ) const;
   void getBlocks( std::vector<       SetupBlock* >& blocks, const uint_t level );

   void getMortonOrder ( std::vector< SetupBlock* >& blocks ) { getBlocks( blocks ); }
   void getHilbertOrder( std::vector< SetupBlock* >& blocks );

   void getProcessSpecificBlocks( std::vector< const SetupBlock* >& blocks, const uint_t process ) const;

   void getBlocksOverlappedByAABB( std::vector< SetupBlock* >& blocks, const AABB& aabb );

   void getBlocks( std::vector< SetupBlock* >& blocks, const uint_t xmin, const uint_t ymin, const uint_t zmin,   // min incl.
                                                       const uint_t xmax, const uint_t ymax, const uint_t zmax ); // max excl.



   inline uint_t mapForestCoordinatesToTreeIndex( const uint_t x, const uint_t y, const uint_t z ) const;

   inline static void mapTreeIndexToForestCoordinates( const uint_t treeIndex, const uint_t xSize, const uint_t ySize, uint_t& x, uint_t& y, uint_t& z );
   inline        void mapTreeIndexToForestCoordinates( const uint_t treeIndex, uint_t& x, uint_t& y, uint_t& z ) const;

   void mapPointToPeriodicDomain( real_t& px, real_t& py, real_t& pz ) const;

   uint_t mapPointToTreeIndex( const real_t px, const real_t py, const real_t pz ) const;

   void mapAABBToBoundingForestCoordinates( const AABB& aabb, uint_t (&min)[3], uint_t (&max)[3] ) const;

   static void getRootBlockAABB( AABB& aabb, const AABB& domain,
                                 const real_t rootBlockXSize, const real_t rootBlockYSize, const real_t rootBlockZSize,
                                 const uint_t xSize, const uint_t ySize, const uint_t zSize,
                                 const uint_t x, const uint_t y, const uint_t z );

   inline void getRootBlockAABB( AABB& aabb, const uint_t x, const uint_t y, const uint_t z ) const;
   inline void getRootBlockAABB( AABB& aabb, const uint_t treeIndex ) const;



   /// Returns true if the block 'block' is located at the lower x-axis border of the domain
   inline bool atDomainXMinBorder( const SetupBlock & block ) const;
   /// Returns true if the block 'block' is located at the upper x-axis border of the domain
   inline bool atDomainXMaxBorder( const SetupBlock & block ) const;
   /// Returns true if the block 'block' is located at the lower y-axis border of the domain
   inline bool atDomainYMinBorder( const SetupBlock & block ) const;
   /// Returns true if the block 'block' is located at the upper y-axis border of the domain
   inline bool atDomainYMaxBorder( const SetupBlock & block ) const;
   /// Returns true if the block 'block' is located at the lower z-axis border of the domain
   inline bool atDomainZMinBorder( const SetupBlock & block ) const;
   /// Returns true if the block 'block' is located at the upper z-axis border of the domain
   inline bool atDomainZMaxBorder( const SetupBlock & block ) const;

   inline bool atDomainMinBorder( const uint_t index, const SetupBlock & block ) const;
   inline bool atDomainMaxBorder( const uint_t index, const SetupBlock & block ) const;



   void init( const AABB& domain, const uint_t xSize, const uint_t ySize, const uint_t zSize,
              const bool xPeriodic, const bool yPeriodic, const bool zPeriodic, const Set<SUID>& selector = Set<SUID>::emptySet() );



   // PROCESS DISTRIBUTION FUNCTIONS

   void assignAllBlocksToRootProcess();

   void balanceLoad( const TargetProcessAssignmentFunction & function,
                     const uint_t numberOfProcesses, const real_t minBufferProcessesFraction = real_t(0),
                     const memory_t perProcessMemoryLimit = memory_t(0),
                     const bool reorderProcessesByBFS = false, const bool insertBufferProcesses = false );
                     
   void balanceLoad( const TargetProcessAssignmentFunction & function,
                     const uint_t numberOfProcesses, const uint_t numberOfBufferProcesses,
                     const memory_t perProcessMemoryLimit = memory_t(0),
                     const bool reorderProcessesByBFS = false, const bool insertBufferProcesses = false );

   // calculate process distribution, target = specific number of processes
   void calculateProcessDistribution_Default( const uint_t        numberOfProcesses,
                                              const memory_t      memoryLimit,
                                              const std::string&  sfcMethod               = std::string( "hilbert" ),
                                              const uint_t        sfcIterations           = 10,
                                              const bool          sortByLevel             = false,
                                              const GlobalLoadBalancing::MetisConfiguration< SetupBlock >& metisConfig =
                                                 GlobalLoadBalancing::MetisConfiguration< SetupBlock >(),
                                              const bool          reorderProcessesByBFS   = false,
                                              const bool          insertBufferProcesses   = false,
                                              const real_t        bufferProcessesFraction = real_c(0) );

   // calculate process distribution, target = specific number of processes
   void calculateProcessDistribution_LevelwiseMetis( const uint_t   numberOfProcesses,
                                                     const bool     reorderProcessesByBFS   = false,
                                                     const CommunicationWeightFunction & communicationWeightFunction = NullCommunicationWeightFunction );

   void calculateProcessDistribution_Greedy( const uint_t   numberOfProcesses,
                                             const memory_t memoryLimit,
                                             const bool     reorderProcessesByBFS   = false,
                                             const bool     insertBufferProcesses   = false,
                                             const real_t   bufferProcessesFraction = real_c(0) );



   uint_t getNumberOfProcesses      () const { return numberOfProcesses_; }
   uint_t getNumberOfWorkerProcesses() const { return numberOfProcesses_ - numberOfBufferProcesses_; }
   uint_t getNumberOfBufferProcesses() const { return numberOfBufferProcesses_; }
   uint_t getProcessIdBytes         () const { uint_t bits = uintMSBPosition( numberOfProcesses_ - 1 );
                                               return ( bits >> 3 ) + ( ( bits & 7 ) ? uint_c(1) : uint_c(0) ); }

   inline bool insertBuffersIntoProcessNetwork() const { return insertBuffersIntoProcessNetwork_; }
   inline bool isWorkerProcess( const uint_t process ) const;
   inline bool isBufferProcess( const uint_t process ) const;



   inline void addRootBlockExclusionFunction( RootBlockExclusionFunction function,
                                              const Set<SUID>& requiredSelectors     = Set<SUID>::emptySet(),
                                              const Set<SUID>& incompatibleSelectors = Set<SUID>::emptySet(),
                                              const std::string& identifier = std::string() );

   inline void addRefinementSelectionFunction( RefinementSelectionFunction function,
                                               const Set<SUID>& requiredSelectors     = Set<SUID>::emptySet(),
                                               const Set<SUID>& incompatibleSelectors = Set<SUID>::emptySet(),
                                               const std::string& identifier = std::string() );

   inline void addWorkloadMemorySUIDAssignmentFunction( WorkloadMemorySUIDAssignmentFunction function,
                                                        const Set<SUID>& requiredSelectors     = Set<SUID>::emptySet(),
                                                        const Set<SUID>& incompatibleSelectors = Set<SUID>::emptySet(),
                                                        const std::string& identifier = std::string() );



   void saveToFile( const char* const filename ) const;

   void writeVTKOutput( const std::string & filestem ) const;
   void writeCSV( const std::string & filestem ) const;



          void        toStream( std::ostream& os ) const;
   inline std::string toString() const;



private:

   void createForest( const Set<SUID>& selector );

   inline void initWorkloadMemorySUID( const Set<SUID>& selector );

   inline void updateNeighborhood( std::set< SetupBlock* >& blocksToUpdate );
          void updateNeighborhood( std::vector< SetupBlock* >& blocks );

   void createNeighborhood();

   static SetupBlock* mapPointToBlock( SetupBlock* const block, const real_t px, const real_t py, const real_t pz );

   void balanceLoadHelper( const TargetProcessAssignmentFunction & function, const uint_t numberOfProcesses, const uint_t numberOfBufferProcesses,
                           const memory_t perProcessMemoryLimit, const bool reorderProcessesByBFS, const bool insertBufferProcesses );
   void calculateProcessDistributionFinalization( const bool reorderProcessesByBFS = false, const bool insertBufferProcesses = false );



   std::vector< SetupBlock* > forest_; // == coarse grid

   uint_t numberOfRootBlocks_;
   uint_t numberOfBlocks_;

   AABB   domain_;           // the simulation space/region
   real_t rootBlockSize_[3]; // [(domain x width) / size_[0]], etc. [equivalent to the size of each root block's bounding box]
   uint_t size_[3];          // number of coarse blocks on the initial grid (= number of octree root blocks) in each direction
   bool   periodic_[3];

   uint_t  depth_; // depth := number of levels - 1
   uint_t  treeIdDigits_;

   uint_t  numberOfProcesses_;
   uint_t  numberOfBufferProcesses_;
   bool    insertBuffersIntoProcessNetwork_;

   std::vector< std::vector< SetupBlock* > > blockDistribution_;



   selectable::SetSelectableObject< RootBlockExclusionFunction, SUID >                     rootBlockExclusionFunctions_;
   selectable::SetSelectableObject< RefinementSelectionFunction, SUID >                   refinementSelectionFunctions_;
   selectable::SetSelectableObject< WorkloadMemorySUIDAssignmentFunction, SUID > workloadMemorySUIDAssignmentFunctions_;



#ifndef NDEBUG
   void checkNeighborhoodConsistency() const;
#endif
};



inline SetupBlockForest::SetupBlockForest() :

   numberOfRootBlocks_( 0 ), numberOfBlocks_( 0 ),
   depth_( 0 ), treeIdDigits_( 0 ), numberOfProcesses_( 0 ), numberOfBufferProcesses_( 0 ),
   insertBuffersIntoProcessNetwork_( false ) {

   rootBlockSize_[0] = rootBlockSize_[1] = rootBlockSize_[2] = real_c(0);
   size_[0]          = size_[1]          = size_[2]          = 0;
   periodic_[0]      = periodic_[1]      = periodic_[2]      = false;
}



inline SetupBlockForest::~SetupBlockForest() {

   for( uint_t i = 0; i != forest_.size(); ++i )
   {
      if( forest_[i] != nullptr ) delete forest_[i];
   }
}



inline SetupBlockForest::const_iterator SetupBlockForest::begin() const {

   const SetupBlock* block = getFirstBlock();

   if( block == nullptr )
      return end();

   return SetupBlockForest::const_iterator( this, block );
}



inline SetupBlockForest::iterator SetupBlockForest::begin() {

   SetupBlock* block = getFirstBlock();

   if( block == nullptr )
      return end();

   return SetupBlockForest::iterator( this, block );
}



inline const SetupBlock* SetupBlockForest::getRootBlock( const uint_t x, const uint_t y, const uint_t z ) const {

   return getRootBlock( mapForestCoordinatesToTreeIndex(x,y,z) );
}



inline SetupBlock* SetupBlockForest::getRootBlock( const uint_t x, const uint_t y, const uint_t z ) {

   return getRootBlock( mapForestCoordinatesToTreeIndex(x,y,z) );
}



inline const SetupBlock* SetupBlockForest::getBlock( const real_t px, const real_t py, const real_t pz ) const {

   if( !domain_.contains( px, py, pz ) )
      return nullptr;

   SetupBlock* block = forest_[ mapPointToTreeIndex( px, py, pz ) ];
   if( block == nullptr ) return nullptr;

   return mapPointToBlock( block, px, py, pz );
}



inline SetupBlock* SetupBlockForest::getBlock( const real_t px, const real_t py, const real_t pz ) {

   if( !domain_.contains( px, py, pz ) )
      return nullptr;

   SetupBlock* block = forest_[ mapPointToTreeIndex( px, py, pz ) ];
   if( block == nullptr ) return nullptr;

   return mapPointToBlock( block, px, py, pz );
}



inline uint_t SetupBlockForest::mapForestCoordinatesToTreeIndex( const uint_t x, const uint_t y, const uint_t z ) const {

   WALBERLA_ASSERT_LESS( x, size_[0] );
   WALBERLA_ASSERT_LESS( y, size_[1] );
   WALBERLA_ASSERT_LESS( z, size_[2] );

   return z * size_[1] * size_[0] + y * size_[0] + x;
}



inline void SetupBlockForest::mapTreeIndexToForestCoordinates( const uint_t treeIndex, const uint_t xSize, const uint_t ySize,
                                                               uint_t& x, uint_t& y, uint_t& z )
{
              z = treeIndex / ( xSize * ySize );
   uint_t index = treeIndex % ( xSize * ySize );
              y = index / xSize;
              x = index % xSize;
}



inline void SetupBlockForest::mapTreeIndexToForestCoordinates( const uint_t treeIndex, uint_t& x, uint_t& y, uint_t& z ) const {

   WALBERLA_ASSERT_LESS( treeIndex, size_[0] * size_[1] * size_[2] );

   mapTreeIndexToForestCoordinates( treeIndex, size_[0], size_[1], x, y, z );
}



inline void SetupBlockForest::getRootBlockAABB( AABB& aabb, const uint_t x, const uint_t y, const uint_t z ) const
{
   getRootBlockAABB( aabb, domain_, rootBlockSize_[0], rootBlockSize_[1], rootBlockSize_[2], size_[0], size_[1], size_[2], x, y, z );
}



inline void SetupBlockForest::getRootBlockAABB( AABB& aabb, const uint_t treeIndex ) const {

   uint_t x,y,z;
   mapTreeIndexToForestCoordinates( treeIndex, x, y, z );
   getRootBlockAABB( aabb, x, y, z );
}



inline bool SetupBlockForest::atDomainXMinBorder( const SetupBlock & block ) const
{
   const AABB & blockAABB = block.getAABB();
   return realIsEqual( blockAABB.xMin(), domain_.xMin(), real_c( 1.0E-6 ) * ( blockAABB.xMax() - blockAABB.xMin() ) );
}



inline bool SetupBlockForest::atDomainXMaxBorder( const SetupBlock & block ) const
{
   const AABB & blockAABB = block.getAABB();
   return realIsEqual( blockAABB.xMax(), domain_.xMax(), real_c( 1.0E-6 ) * ( blockAABB.xMax() - blockAABB.xMin() ) );
}



inline bool SetupBlockForest::atDomainYMinBorder( const SetupBlock & block ) const
{
   const AABB & blockAABB = block.getAABB();
   return realIsEqual( blockAABB.yMin(), domain_.yMin(), real_c( 1.0E-6 ) * ( blockAABB.yMax() - blockAABB.yMin() ) );
}



inline bool SetupBlockForest::atDomainYMaxBorder( const SetupBlock & block ) const
{
   const AABB & blockAABB = block.getAABB();
   return realIsEqual( blockAABB.yMax(), domain_.yMax(), real_c( 1.0E-6 ) * ( blockAABB.yMax() - blockAABB.yMin() ) );
}



inline bool SetupBlockForest::atDomainZMinBorder( const SetupBlock & block ) const
{
   const AABB & blockAABB = block.getAABB();
   return realIsEqual( blockAABB.zMin(), domain_.zMin(), real_c( 1.0E-6 ) * ( blockAABB.zMax() - blockAABB.zMin() ) );
}



inline bool SetupBlockForest::atDomainZMaxBorder( const SetupBlock & block ) const
{
   const AABB & blockAABB = block.getAABB();
   return realIsEqual( blockAABB.zMax(), domain_.zMax(), real_c( 1.0E-6 ) * ( blockAABB.zMax() - blockAABB.zMin() ) );
}



inline bool SetupBlockForest::atDomainMinBorder( const uint_t index, const SetupBlock & block ) const
{
   WALBERLA_ASSERT_LESS( index, uint_t(3) );
   const AABB & blockAABB = block.getAABB();
   return realIsEqual( blockAABB.min( index ), domain_.min( index ), real_c( 1.0E-6 ) * ( blockAABB.max( index ) - blockAABB.min( index) ) );
}



inline bool SetupBlockForest::atDomainMaxBorder( const uint_t index, const SetupBlock & block ) const
{
   WALBERLA_ASSERT_LESS( index, uint_t(3) );
   const AABB & blockAABB = block.getAABB();
   return realIsEqual( blockAABB.max( index ), domain_.max( index ), real_c( 1.0E-6 ) * ( blockAABB.max( index ) - blockAABB.min( index ) ) );
}



inline void SetupBlockForest::initWorkloadMemorySUID( const Set<SUID>& selector ) {

   std::vector< WorkloadMemorySUIDAssignmentFunction > workloadMemorySUIDAssignmentFunctions;

   workloadMemorySUIDAssignmentFunctions_.get( workloadMemorySUIDAssignmentFunctions, selector );

   if( workloadMemorySUIDAssignmentFunctions.empty() ) {
      WALBERLA_LOG_PROGRESS( "Initializing SetupBlockForest: No callback functions for assigning workload, memory requirements, and SUIDs to blocks found.\n"
                             "                               Default values (zero/nothing) will be assigned ..." );
   }
   else {
      WALBERLA_LOG_PROGRESS( "Initializing SetupBlockForest: Assigning workload, memory requirements, and SUIDs to blocks ..." );
   }

   for( uint_t i = 0; i != workloadMemorySUIDAssignmentFunctions.size(); ++i )
      workloadMemorySUIDAssignmentFunctions[i]( *this );
}



inline void SetupBlockForest::updateNeighborhood( std::set< SetupBlock* >& blocksToUpdate ) {

   std::vector< SetupBlock* > blocks;

   for( std::set< SetupBlock* >::iterator it = blocksToUpdate.begin(); it != blocksToUpdate.end(); ++it )
      blocks.push_back( *it );

   updateNeighborhood( blocks );
}



inline bool SetupBlockForest::isWorkerProcess( const uint_t process ) const
{
   WALBERLA_ASSERT_LESS( process, numberOfProcesses_ );
   WALBERLA_ASSERT_LESS( process, blockDistribution_.size() );

   return !(blockDistribution_[ process ].empty());
}



inline bool SetupBlockForest::isBufferProcess( const uint_t process ) const {

   WALBERLA_ASSERT_LESS( process, numberOfProcesses_ );
   WALBERLA_ASSERT_LESS( process, blockDistribution_.size() );

   return blockDistribution_[ process ].empty();
}



inline void SetupBlockForest::addRootBlockExclusionFunction( RootBlockExclusionFunction function,
                                                             const Set<SUID>& requiredSelectors, const Set<SUID>& incompatibleSelectors,
                                                             const std::string& identifier ) {

   rootBlockExclusionFunctions_.add( function, requiredSelectors, incompatibleSelectors, identifier );
}



inline void SetupBlockForest::addRefinementSelectionFunction( RefinementSelectionFunction function,
                                                              const Set<SUID>& requiredSelectors, const Set<SUID>& incompatibleSelectors,
                                                              const std::string& identifier ) {

   refinementSelectionFunctions_.add( function, requiredSelectors, incompatibleSelectors, identifier );
}



inline void SetupBlockForest::addWorkloadMemorySUIDAssignmentFunction( WorkloadMemorySUIDAssignmentFunction function,
                                                                       const Set<SUID>& requiredSelectors, const Set<SUID>& incompatibleSelectors,
                                                                       const std::string& identifier ) {

   workloadMemorySUIDAssignmentFunctions_.add( function, requiredSelectors, incompatibleSelectors, identifier );
}



inline std::string SetupBlockForest::toString() const
{
   std::ostringstream oss;
   toStream( oss );
   return oss.str();
}



//////////////////////
// Global Functions //
//////////////////////



inline std::ostream & operator<<( std::ostream & os, const SetupBlockForest & forest )
{
   forest.toStream( os );
   return os;
}



////////////////////
// Useful helpers //
////////////////////



class RefinementSelectionFunctions
{
public:

   using RefinementSelectionFunction = blockforest::SetupBlockForest::RefinementSelectionFunction;

   void add( const RefinementSelectionFunction & function )
   {
      function_.push_back( function );
   }

   void operator()( SetupBlockForest & forest )
   {
      for( auto function = function_.begin(); function != function_.end(); ++function )
         (*function)( forest );
   }

private:

   std::vector< RefinementSelectionFunction > function_;

}; // class RefinementSelectionFunctions



} // namespace blockforest

using blockforest::SetupBlockForest;

} // namespace walberla
