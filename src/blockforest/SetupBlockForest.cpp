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
//! \file SetupBlockForest.cpp
//! \ingroup blockforest
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#include "BlockForestFile.h"
#include "BlockNeighborhoodConstruction.h"
#include "BlockNeighborhoodSection.h"
#include "HilbertCurveConstruction.h"
#include "OutputColor.h"
#include "SetupBlockForest.h"
#include "Utility.h"
#include "loadbalancing/Cartesian.h"
#include "loadbalancing/StaticCurve.h"

#include "core/Abort.h"
#include "core/OpenMP.h"
#include "core/debug/CheckFunctions.h"
#include "core/logging/Logging.h"
#include "core/math/Sample.h"

#include "domain_decomposition/MapPointToPeriodicDomain.h"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stack>


namespace walberla {
namespace blockforest {



AABB SetupBlockForest::RootBlockAABB::operator()( const uint_t index ) const // index = treeIndex
{
   AABB aabb;

   uint_t x,y,z;
   SetupBlockForest::mapTreeIndexToForestCoordinates( index, xSize_, ySize_, x, y, z );
   SetupBlockForest::getRootBlockAABB( aabb, domain_, rootBlockXSize_, rootBlockYSize_, rootBlockZSize_, xSize_, ySize_, zSize_, x, y, z );

   return aabb;
}



uint_t SetupBlockForest::getNumberOfBlocks( const uint_t level ) const
{
   if( depth_ == uint_t(0) )
      return getNumberOfBlocks();

   uint_t count( uint_t(0) );

   for( uint_t i = 0; i != forest_.size(); ++i ) {

      if( forest_[i] == nullptr )
         continue;

      std::stack< SetupBlock* > stack;

      stack.push( forest_[i] );

      while( !stack.empty() ) {

         SetupBlock* const block = stack.top();
         stack.pop();

         WALBERLA_ASSERT_NOT_NULLPTR( block );

         if( block->hasChildren() ) {
            for( uint_t c = 8; c-- != 0; )
               stack.push( block->getChild(c) );
         }
         else if( block->getLevel() == level )
            ++count;
      }
   }

   return count;
}



const SetupBlock* SetupBlockForest::getFirstBlock() const {

   SetupBlock* block = nullptr;
   for( uint_t i = 0; i != forest_.size() && block == nullptr; ++i )
      block = forest_[i];

   if( block == nullptr )
      return nullptr;

   while( block->hasChildren() )
      block = block->getChild(0);

   return block;
}



SetupBlock* SetupBlockForest::getFirstBlock() {

   SetupBlock* block = nullptr;
   for( uint_t i = 0; i != forest_.size() && block == nullptr; ++i )
      block = forest_[i];

   if( block == nullptr )
      return nullptr;

   while( block->hasChildren() )
      block = block->getChild(0);

   return block;
}



const SetupBlock* SetupBlockForest::getNextBlock( const SetupBlock* block ) const {

   if( block == nullptr )
      return nullptr;

   // ASCEND

   uint_t child = 0;

   bool loop = true;
   while( block->hasFather() && loop ) {
      if( block->getId().getBranchId() != 7 ) { // if( block_ is NOT the rightmost child of its father )
         loop  = false;
         child = block->getId().getBranchId() + 1;
      }
      block = block->getFather();
   }

   // SWITCH TO THE NEXT TREE

   if( child == 0 ) {

      uint_t treeIndex = block->getId().getTreeIndex() + uint_c(1);

      WALBERLA_ASSERT_LESS( treeIndex-1 ,forest_.size() );
      WALBERLA_ASSERT_EQUAL( block, forest_[ treeIndex-1 ] );

      while( treeIndex < forest_.size() && forest_[ treeIndex ] == nullptr ) ++treeIndex;

      if( treeIndex == forest_.size() )
         return nullptr;

      block = forest_[ treeIndex ];
   }

   // DESCEND

   while( block->hasChildren() ) {
      block = block->getChild( child );
      child  = 0;
   }

   return block;
}



SetupBlock* SetupBlockForest::getNextBlock( const SetupBlock* block ) {

   if( block == nullptr )
      return nullptr;

   // ASCEND

   uint_t child = 0;

   bool loop = true;
   while( block->hasFather() && loop ) {
      if( block->getId().getBranchId() != 7 ) { // if( block_ is NOT the rightmost child of its father )
         loop  = false;
         child = block->getId().getBranchId() + 1;
      }
      block = block->getFather();
   }

   // SWITCH TO THE NEXT TREE

   if( child == 0 ) {

      uint_t treeIndex = block->getId().getTreeIndex() + uint_c(1);

      WALBERLA_ASSERT_LESS( treeIndex-1, forest_.size() );
      WALBERLA_ASSERT_EQUAL( block, forest_[ treeIndex-1 ] );

      while( treeIndex < forest_.size() && forest_[ treeIndex ] == nullptr ) ++treeIndex;

      if( treeIndex == forest_.size() )
         return nullptr;

      block = forest_[ treeIndex ];
   }

   // DESCEND

   while( block->hasChildren() ) {
      block = block->getChild( child );
      child  = 0;
   }

   return const_cast< SetupBlock* >( block );
}



const SetupBlock* SetupBlockForest::getBlock( const BlockID& id ) const {

   WALBERLA_ASSERT_GREATER_EQUAL( id.getUsedBits(), treeIdDigits_ );
   WALBERLA_ASSERT_EQUAL( ( id.getUsedBits() - treeIdDigits_ ) % 3, 0 );

   BlockID blockId( id );

   const uint_t levels = ( blockId.getUsedBits() - treeIdDigits_ ) / 3;

   std::vector< uint_t > branchId( levels );

   for( uint_t i = levels; i-- != 0; ) {
      branchId[i] = blockId.getBranchId();
      blockId.removeBranchId();
   }

   const uint_t index = blockId.getTreeIndex();

   WALBERLA_ASSERT_LESS( index, forest_.size() );
   WALBERLA_ASSERT_NOT_NULLPTR( forest_[index] );

   SetupBlock* block = forest_[index];

   for( uint_t i = 0; i != levels; ++i ) {
      WALBERLA_ASSERT( block->hasChildren() );
      WALBERLA_ASSERT_NOT_NULLPTR( block->getChild( branchId[i] ) );
      block = block->getChild( branchId[i] );
   }

   return block;
}



void SetupBlockForest::getBlocks( std::vector< const SetupBlock* >& blocks ) const {

   // ATTENTION: the vector 'blocks' is not emptied

   for( uint_t i = 0; i != forest_.size(); ++i ) {

      if( forest_[i] == nullptr )
         continue;

      // depth-first search

      std::stack< SetupBlock* > stack;

      stack.push( forest_[i] );

      while( !stack.empty() ) {

         SetupBlock* const block = stack.top();
         stack.pop();

         WALBERLA_ASSERT_NOT_NULLPTR( block );

         if( block->hasChildren() ) {
            for( uint_t c = 8; c-- != 0; )
               stack.push( block->getChild(c) );
         }
         else blocks.push_back( block );
      }
   }
}



void SetupBlockForest::getBlocks( std::vector< SetupBlock* >& blocks ) {

   // ATTENTION: the vector 'blocks' is not emptied

   for( uint_t i = 0; i != forest_.size(); ++i ) {

      if( forest_[i] == nullptr )
         continue;

      // depth-first search

      std::stack< SetupBlock* > stack;

      stack.push( forest_[i] );

      while( !stack.empty() ) {

         SetupBlock* const block = stack.top();
         stack.pop();

         WALBERLA_ASSERT_NOT_NULLPTR( block );

         if( block->hasChildren() ) {
            for( uint_t c = 8; c-- != 0; )
               stack.push( block->getChild(c) );
         }
         else blocks.push_back( block );
      }
   }
}



void SetupBlockForest::getBlocks( std::vector< const SetupBlock* >& blocks, const uint_t level ) const {

   // ATTENTION: the vector 'blocks' is not emptied

   for( uint_t i = 0; i != forest_.size(); ++i ) {

      if( forest_[i] == nullptr )
         continue;

      std::stack< SetupBlock* > stack;

      stack.push( forest_[i] );

      while( !stack.empty() ) {

         SetupBlock* const block = stack.top();
         stack.pop();

         WALBERLA_ASSERT_NOT_NULLPTR( block );

         if( block->hasChildren() ) {
            for( uint_t c = 8; c-- != 0; )
               stack.push( block->getChild(c) );
         }
         else if( block->getLevel() == level )
            blocks.push_back( block );
      }
   }
}



void SetupBlockForest::getBlocks( std::vector< SetupBlock* >& blocks, const uint_t level ) {

   // ATTENTION: the vector 'blocks' is not emptied

   for( uint_t i = 0; i != forest_.size(); ++i ) {

      if( forest_[i] == nullptr )
         continue;

      std::stack< SetupBlock* > stack;

      stack.push( forest_[i] );

      while( !stack.empty() ) {

         SetupBlock* const block = stack.top();
         stack.pop();

         WALBERLA_ASSERT_NOT_NULLPTR( block );

         if( block->hasChildren() ) {
            for( uint_t c = 8; c-- != 0; )
               stack.push( block->getChild(c) );
         }
         else if( block->getLevel() == level )
            blocks.push_back( block );
      }
   }
}



void SetupBlockForest::getHilbertOrder( std::vector< SetupBlock* >& blocks ) {

   // ATTENTION: the vector 'blocks' is not emptied

   uint_t y = 0;
   uint_t x = 0;

   uint_t yLoopEnd = size_[1];
   uint_t xLoopEnd = size_[0];

   for( uint_t z = 0; z != size_[2]; ++z ) {
      while( y != yLoopEnd ) {
         if( yLoopEnd == 0 ) --y;
         while( x != xLoopEnd ) {
            if( xLoopEnd == 0 ) --x;

            WALBERLA_ASSERT_LESS( z*size_[0]*size_[1] + y*size_[0] + x, forest_.size() );

            SetupBlock* root = forest_[ z*size_[0]*size_[1] + y*size_[0] + x ];

            if( root != nullptr ) {

               std::stack< SetupBlock* > stack;
               std::stack< uint_t > orientation;

               stack.push( root );
               orientation.push( 0 );

               while( !stack.empty() ) {

                  SetupBlock* const block = stack.top();
                  uint_t index = orientation.top();

                  stack.pop();
                  orientation.pop();

                  WALBERLA_ASSERT_NOT_NULLPTR( block );

                  if( block->hasChildren() ) {
                     for( uint_t c = 8; c-- != 0; ) {
                        stack.push( block->getChild( hilbertOrder[index][c] ) );
                        orientation.push( hilbertOrientation[index][c] );
                     }
                  }
                  else blocks.push_back( block );
               }
            }
            if( xLoopEnd != 0 ) ++x;
         }
         WALBERLA_ASSERT_EQUAL( x, xLoopEnd );
         xLoopEnd = ( xLoopEnd == 0 ) ? size_[0] : 0;
         if( yLoopEnd != 0 ) ++y;
      }
      WALBERLA_ASSERT_EQUAL( y, yLoopEnd );
      yLoopEnd = ( yLoopEnd == 0 ) ? size_[1] : 0;
   }
}



void SetupBlockForest::getProcessSpecificBlocks( std::vector< const SetupBlock* >& blocks, const uint_t process ) const {

   // ATTENTION: the vector 'blocks' is not emptied

   WALBERLA_ASSERT_LESS( process, blockDistribution_.size() );

   const std::vector< SetupBlock* >& processBlocks = blockDistribution_[ process ];

   for( uint_t i = 0; i < processBlocks.size(); ++i )
      blocks.push_back( processBlocks[i] );
}



void SetupBlockForest::getBlocksOverlappedByAABB( std::vector< SetupBlock* >& blocks, const AABB& aabb ) {

   // ATTENTION: the vector 'blocks' is not emptied

   if( aabb.xMin() >= domain_.xMax() || aabb.xMax() <= domain_.xMin() ||
       aabb.yMin() >= domain_.yMax() || aabb.yMax() <= domain_.yMin() ||
       aabb.zMin() >= domain_.zMax() || aabb.zMax() <= domain_.zMin() )
      return;

   uint_t min[3], max[3];

   mapAABBToBoundingForestCoordinates( aabb, min, max );

   std::stack< SetupBlock* > stack;

   for( uint_t z = max[2]; z-- != min[2]; ) {
      for( uint_t y = max[1]; y-- != min[1]; ) {
         for( uint_t x = max[0]; x-- != min[0]; ) {
            SetupBlock* const block = forest_[ mapForestCoordinatesToTreeIndex(x,y,z) ];
            if( block != nullptr && block->getAABB().intersects( aabb ) )
               stack.push( block );
         }
      }
   }

   while( !stack.empty() ) {

      SetupBlock* const block = stack.top();
      stack.pop();

      WALBERLA_ASSERT( block->getAABB().intersects( aabb ) );

      if( block->hasChildren() ) {

         for( uint_t i = 8; i-- != 0; )
            if( block->getChild(i)->getAABB().intersects( aabb ) )
               stack.push( block->getChild(i) );
      }
      else blocks.push_back( block );
   }
}



void SetupBlockForest::getBlocks( std::vector< SetupBlock* >& blocks, const uint_t xmin, const uint_t ymin, const uint_t zmin,    // min incl.
                                                                      const uint_t xmax, const uint_t ymax, const uint_t zmax ) { // max excl.

   WALBERLA_ASSERT_LESS_EQUAL( xmin, xmax ); WALBERLA_ASSERT_LESS_EQUAL( xmax, size_[0] );
   WALBERLA_ASSERT_LESS_EQUAL( ymin, ymax ); WALBERLA_ASSERT_LESS_EQUAL( ymax, size_[1] );
   WALBERLA_ASSERT_LESS_EQUAL( zmin, zmax ); WALBERLA_ASSERT_LESS_EQUAL( zmax, size_[2] );

   std::stack< SetupBlock* > stack;

   for( uint_t z = zmax; z-- != zmin; ) {
      for( uint_t y = ymax; y-- != ymin; ) {
         for( uint_t x = xmax; x-- != xmin; ) {
            SetupBlock* const block = forest_[ mapForestCoordinatesToTreeIndex(x,y,z) ];
            if( block != nullptr ) stack.push( block );
         }
      }
   }

   while( !stack.empty() ) {

      SetupBlock* const block = stack.top();
      stack.pop();

      if( block->hasChildren() ) {
         for( uint_t i = 8; i-- != 0; )
            stack.push( block->getChild(i) );
      }
      else blocks.push_back( block );
   }
}



void SetupBlockForest::mapPointToPeriodicDomain( real_t & px, real_t & py, real_t & pz ) const
{
   std::array< bool, 3 > periodic;
   periodic[0] = periodic_[0];
   periodic[1] = periodic_[1];
   periodic[2] = periodic_[2];

   domain_decomposition::mapPointToPeriodicDomain( periodic, domain_, px, py, pz );
}



uint_t SetupBlockForest::mapPointToTreeIndex( const real_t px, const real_t py, const real_t pz ) const {

   WALBERLA_ASSERT( domain_.contains( px, py, pz ) );

   uint_t x = static_cast< uint_t >( ( px - domain_.xMin() ) / rootBlockSize_[0] );
   uint_t y = static_cast< uint_t >( ( py - domain_.yMin() ) / rootBlockSize_[1] );
   uint_t z = static_cast< uint_t >( ( pz - domain_.zMin() ) / rootBlockSize_[2] );

   if( x >= size_[0] ) x = size_[0] - 1; // shouldn't happen, ...
   if( y >= size_[1] ) y = size_[1] - 1; // ... but might happen due to ...
   if( z >= size_[2] ) z = size_[2] - 1; // ... floating point inaccuracy?

   return mapForestCoordinatesToTreeIndex( x, y, z );
}



void SetupBlockForest::mapAABBToBoundingForestCoordinates( const AABB& aabb, uint_t (&min)[3], uint_t (&max)[3] ) const {

   // ATTENTION: min[3] incl., max[3] excl.

   WALBERLA_ASSERT_LESS( aabb.xMin(), domain_.xMax() ); WALBERLA_ASSERT_GREATER( aabb.xMax(), domain_.xMin() );
   WALBERLA_ASSERT_LESS( aabb.yMin(), domain_.yMax() ); WALBERLA_ASSERT_GREATER( aabb.yMax(), domain_.yMin() );
   WALBERLA_ASSERT_LESS( aabb.zMin(), domain_.zMax() ); WALBERLA_ASSERT_GREATER( aabb.zMax(), domain_.zMin() );

   min[0] = min[1] = min[2] = 0;
   max[0] = size_[0] - 1;
   max[1] = size_[1] - 1;
   max[2] = size_[2] - 1;

   for( uint_t i = 0; i != 3; ++i ) {

      if( aabb.min(i) > domain_.min(i) ) {

         min[i] = static_cast< uint_t >( ( aabb.min(i) - domain_.min(i) ) / rootBlockSize_[i] );

         if( min[i] >= size_[i] ) min[i] = size_[i] - 1; // shouldn't happen, but might happen due to floating point inaccuracy?

         // shouldn't happen, but might happen due to floating point inaccuracy?
         SetupBlock* block = forest_[ mapForestCoordinatesToTreeIndex( min[0], min[1], min[2] ) ];
         if( block != nullptr ) {
            const AABB& tree = block->getAABB();
            if( aabb.min(i) < tree.min(i) ) --min[i];
            else if( aabb.min(i) >= tree.max(i) ) ++min[i];
         }
      }

      if( aabb.max(i) < domain_.max(i) ) {

         max[i] = static_cast< uint_t >( ( aabb.max(i) - domain_.min(i) ) / rootBlockSize_[i] );

         if( max[i] >= size_[i] ) max[i] = size_[i] - 1; // shouldn't happen, but might happen due to floating point inaccuracy?

         // shouldn't happen, but might happen due to floating point inaccuracy?
         SetupBlock* block = forest_[ mapForestCoordinatesToTreeIndex( max[0], max[1], max[2] ) ];
         if( block != nullptr ) {
            const AABB& tree = block->getAABB();
            if( aabb.max(i) <= tree.min(i) ) --max[i];
            else if( aabb.max(i) > tree.max(i) ) ++max[i];
         }
      }
   }

   ++max[0]; ++max[1]; ++max[2];
}



void SetupBlockForest::getRootBlockAABB( AABB & aabb, const AABB & domain,
                                         const real_t rootBlockXSize, const real_t rootBlockYSize, const real_t rootBlockZSize,
                                         const uint_t xSize, const uint_t ySize, const uint_t zSize,
                                         const uint_t x, const uint_t y, const uint_t z )
{
   WALBERLA_ASSERT_LESS( x, xSize );
   WALBERLA_ASSERT_LESS( y, ySize );
   WALBERLA_ASSERT_LESS( z, zSize );

   aabb.initMinMaxCorner( domain.xMin() + static_cast< real_t >(  x  ) * rootBlockXSize,
                          domain.yMin() + static_cast< real_t >(  y  ) * rootBlockYSize,
                          domain.zMin() + static_cast< real_t >(  z  ) * rootBlockZSize,
                          ( x+1 == xSize ) ? domain.xMax() : domain.xMin() + static_cast< real_t >( x+1 ) * rootBlockXSize,
                          ( y+1 == ySize ) ? domain.yMax() : domain.yMin() + static_cast< real_t >( y+1 ) * rootBlockYSize,
                          ( z+1 == zSize ) ? domain.zMax() : domain.zMin() + static_cast< real_t >( z+1 ) * rootBlockZSize );
}



void SetupBlockForest::init( const AABB& domain, const uint_t xSize, const uint_t ySize, const uint_t zSize,
                             const bool xPeriodic, const bool yPeriodic, const bool zPeriodic, const Set<SUID>& selector )
{
   WALBERLA_LOG_PROGRESS( "Initializing SetupBlockForest:" <<
                          "\n - AABB: " << domain <<
                          "\n - forest size (root blocks / blocks on the initial grid): " << xSize << " x " << ySize << " x " << zSize <<
                          "\n - periodicity: " << std::boolalpha << xPeriodic << " x " << yPeriodic << " x " << zPeriodic );

   if( xSize * ySize * zSize == uint_c(0) )
      WALBERLA_ABORT( "Initializing SetupBlockForest failed: xSize (= " << xSize << ") * "
                                                            "ySize (= " << ySize << ") * zSize (= " << zSize << ") == 0!" );

   if( !( ( xSize <= std::numeric_limits< std::vector< SetupBlock* >::size_type >::max() ) &&
          ( ySize <= std::numeric_limits< std::vector< SetupBlock* >::size_type >::max() / xSize ) &&
          ( zSize <= std::numeric_limits< std::vector< SetupBlock* >::size_type >::max() / ( xSize * ySize ) ) ) )
      WALBERLA_ABORT( "Initializing SetupBlockForest failed: You requested too many blocks "
                      "(xSize (= " << xSize << ") * ySize (= " << ySize << ") * zSize (= " << zSize << "))!" );

   if( !( ( xSize <= std::numeric_limits< uint_t >::max() ) &&
          ( ySize <= std::numeric_limits< uint_t >::max() / xSize ) &&
          ( zSize <= std::numeric_limits< uint_t >::max() / ( xSize * ySize ) ) ) )
      WALBERLA_ABORT( "Initializing SetupBlockForest failed: You requested too many blocks "
                      "(xSize (= " << xSize << ") * ySize (= " << ySize << ") * zSize (= " << zSize << "))!" );

   if( static_cast< uint_t >( std::numeric_limits< uint_t >::digits ) < 1 + uintMSBPosition( xSize * ySize * zSize - 1 ) )
      WALBERLA_ABORT( "Initializing SetupBlockForest failed: You requested too many blocks "
                      "(xSize (= " << xSize << ") * ySize (= " << ySize << ") * zSize (= " << zSize << "))!" );

   if( !forest_.empty() ) {
      for( uint_t i = 0; i != forest_.size(); ++i ) {
         if( forest_[i] != nullptr ) delete forest_[i];
      }
      forest_.clear();
   }

   // initialization of data members

   domain_ = domain;

   rootBlockSize_[0] = ( domain.xMax() - domain.xMin() ) / static_cast< real_t >( xSize );
   rootBlockSize_[1] = ( domain.yMax() - domain.yMin() ) / static_cast< real_t >( ySize );
   rootBlockSize_[2] = ( domain.zMax() - domain.zMin() ) / static_cast< real_t >( zSize );

   size_[0] = xSize;
   size_[1] = ySize;
   size_[2] = zSize;

   periodic_[0] = xPeriodic;
   periodic_[1] = yPeriodic;
   periodic_[2] = zPeriodic;

   depth_ = 0;

   const uint_t maxIndex = xSize * ySize * zSize - 1;
   uint_t treeIdMarker = 1;
   while( treeIdMarker <= maxIndex )
      treeIdMarker <<= 1;

   treeIdDigits_ = uintMSBPosition( treeIdMarker );

   // determine which root blocks are not needed and can be excluded

   const uint_t size = xSize * ySize * zSize;

   std::vector<uint8_t> excludeBlock( size, 0 );

   std::vector< RootBlockExclusionFunction > rootBlockExclusionFunctions;

   rootBlockExclusionFunctions_.get( rootBlockExclusionFunctions, selector );

   RootBlockAABB rootBlockAABB( domain, rootBlockSize_[0], rootBlockSize_[1], rootBlockSize_[2], xSize, ySize, zSize );

   if( rootBlockExclusionFunctions.size() > 0 )
      WALBERLA_LOG_PROGRESS( "Initializing SetupBlockForest: Calling root block exclusion callback functions ..." );

   for( uint_t i = 0; i != rootBlockExclusionFunctions.size(); ++i )
      rootBlockExclusionFunctions[i]( excludeBlock, rootBlockAABB );

   // creation of all root blocks

   WALBERLA_LOG_PROGRESS( "Initializing SetupBlockForest: Allocating root blocks ..." );

   forest_.resize( size, nullptr );
   numberOfRootBlocks_ = uint_c(0);

   AABB aabb;

   for( uint_t z = 0; z != zSize; ++z ) {
      for( uint_t y = 0; y != ySize; ++y ) {
         for( uint_t x = 0; x != xSize; ++x ) {

            const uint_t treeIndex = z * ySize * xSize + y * xSize + x;

            if( excludeBlock[treeIndex] == 0 )
            {
               getRootBlockAABB( aabb, x, y, z );

               forest_[ treeIndex ] = new SetupBlock( nullptr, BlockID( treeIndex, treeIdMarker ),
                                                      aabb.xMin(), aabb.yMin(), aabb.zMin(), aabb.xMax(), aabb.yMax(), aabb.zMax(), 0 );
               ++numberOfRootBlocks_;
            }
         }
      }
   }

   numberOfBlocks_ = numberOfRootBlocks_;

   // neighborhood setup (with respect to periodicity)

   WALBERLA_LOG_PROGRESS( "Initializing SetupBlockForest: Setting up neighborhood information for each block (with respect to periodicity) ..." );

#ifdef _OPENMP
   const int izSize = int_c( zSize );
   #pragma omp parallel for schedule(static)
   for( int iz = 0; iz < izSize; ++iz ) {
      uint_t z = uint_c( iz );
#else
   for( uint_t z = 0; z < zSize; ++z ) {
#endif
      for( uint_t y = 0; y != ySize; ++y ) {
         for( uint_t x = 0; x != xSize; ++x ) {

            const uint_t treeIndex = z * ySize * xSize + y * xSize + x;

            WALBERLA_ASSERT_LESS( treeIndex, forest_.size() );

            if( forest_[ treeIndex ] != nullptr ) {

               for( uint_t w = 0; w != 3; ++w ) {
                  for( uint_t v = 0; v != 3; ++v ) {
                     for( uint_t u = 0; u != 3; ++u ) {

                        if( ( u == 1 && v == 1 && w == 1 ) ||
                            ( !xPeriodic && ( (u == 0 && x == 0) || (u == 2 && x == xSize-1) ) ) ||
                            ( !yPeriodic && ( (v == 0 && y == 0) || (v == 2 && y == ySize-1) ) ) ||
                            ( !zPeriodic && ( (w == 0 && z == 0) || (w == 2 && z == zSize-1) ) ) ) continue;

                        uint_t n = w * 9 + v * 3 + u;
                        if( n > 13 ) --n;

                        uint_t nIndex = 0;

                        if( u == 0 && x == 0 ) nIndex += xSize-1;
                        else if( !(u == 2 && x == xSize-1) ) nIndex += ( x + u ) - 1;

                        if( v == 0 && y == 0 ) nIndex += ( ySize-1 ) * xSize;
                        else if( !(v == 2 && y == ySize-1) ) nIndex += (( y + v ) - 1) * xSize;

                        if( w == 0 && z == 0 ) nIndex += ( zSize-1 ) * xSize * ySize;
                        else if( !(w == 2 && z == zSize-1) ) nIndex += (( z + w ) - 1) * xSize * ySize;

                        WALBERLA_ASSERT_LESS( n, 26 );
                        WALBERLA_ASSERT_LESS( nIndex, forest_.size() );

                        if( forest_[ nIndex ] != nullptr )
                           forest_[ treeIndex ]->addNeighbor( n, forest_[ nIndex ] );
                     }
                  }
               }

               forest_[ treeIndex ]->assembleNeighborhood();
            }
         }
      }
   }

   createForest( selector );

   initWorkloadMemorySUID( selector );

   WALBERLA_LOG_PROGRESS( "Initializing SetupBlockForest: finished!\nThe following block structure has been created:\n" << *this );
}



void SetupBlockForest::createForest( const Set<SUID>& selector ) {

#ifndef NDEBUG
   checkNeighborhoodConsistency();
#endif

   std::vector< RefinementSelectionFunction > refinementSelectionFunctions;

   refinementSelectionFunctions_.get( refinementSelectionFunctions, selector );

   if( refinementSelectionFunctions.empty() )
      return;

   WALBERLA_LOG_PROGRESS( "Initializing SetupBlockForest: There is at least one refinement selection function.\n"
                          "                               Creating the forest (iteratively) ..." );

   bool loop = true;
   while( loop ) {

      // MARK ALL BLOCKS THAT NEED TO BE SPLIT

      for( uint_t i = 0; i < refinementSelectionFunctions.size(); ++i )
         refinementSelectionFunctions[i]( *this );

      // GET ALL BLOCKS (= ALL LEAVES)

      std::vector< SetupBlock* > blocks;
      getBlocks( blocks );

      // MARK ADDITIONAL BLOCKS IN ORDER TO ENSURE 2:1 SIZE RATIO OF NEIGHBORING BLOCKS

      bool markerLoop = true;
      while( markerLoop ) {
         markerLoop = false;

         const int blockssize = int_c( blocks.size() );
#ifdef _OPENMP
         #pragma omp parallel for schedule(static)
#endif
         for( int i = 0; i < blockssize; ++i ) { // mark every block that is not yet marked, but has a smaller
            SetupBlock* const block = blocks[ uint_c(i) ]; // neighboring block that is marked to be split
            for( uint_t n = 0; n < 26 && !block->isMarked(); ++n ) {
               if( block->neighborhoodSectionHasSmallerBlocks(n) ) {
                  for( uint_t c = 0; c < block->getNeighborhoodSectionSize(n) && !block->isMarked(); ++c ) {
                     if( block->getNeighbor(n,c)->isMarked() ) {
                        block->setMarker( true );
                        markerLoop = true;
                     }
                  }
               }
            }
         }
      }

      // IDENTIFY ALL BLOCKS THAT ARE MARKED TO BE SPLIT

      std::vector< SetupBlock* > blocksToSplit;
      for( uint_t i = 0; i < blocks.size(); ++i )
         if( blocks[i]->isMarked() ) blocksToSplit.push_back( blocks[i] );

      numberOfBlocks_ += blocksToSplit.size() * 7;

      if( blocksToSplit.empty() ) {
         WALBERLA_LOG_PROGRESS( "Initializing SetupBlockForest: No blocks marked for refinement, aborting refinement process ..." );
      }
      else if( blocksToSplit.size() == 1 ) {
         WALBERLA_LOG_PROGRESS( "Initializing SetupBlockForest: One block marked for refinement, splitting block ..." );
      }
      else {
         WALBERLA_LOG_PROGRESS( "Initializing SetupBlockForest: " << blocksToSplit.size() << " blocks marked for refinement, splitting blocks ..." );
      }

      // SPLIT BLOCKS

      int blocksToSplitsize = int_c( blocksToSplit.size() );
#ifdef _OPENMP
      #pragma omp parallel for schedule(static)
#endif
      for( int i = 0; i < blocksToSplitsize; ++i )
      {
         WALBERLA_LOG_DETAIL( "Initializing SetupBlockForest: Splitting block with ID " << blocksToSplit[ uint_c(i) ]->getId() <<
                              "\n                               - AABB: " << blocksToSplit[ uint_c(i) ]->getAABB() <<
                              "\n                               - level: " << blocksToSplit[ uint_c(i) ]->getLevel() );
         blocksToSplit[ uint_c(i) ]->split();
      }

      // IDENTIFY ALL BLOCKS THAT NEED TO UPDATE THEIR NEIGHBORHOOD & ADAPT 'depth_' DATA MEMBER

      std::set< SetupBlock* > blocksToUpdate;

      for( uint_t i = 0; i < blocksToSplit.size(); ++i ) {

         SetupBlock* const block = blocksToSplit[i];

         for( uint_t c = 0; c != 8; ++c ) {
            blocksToUpdate.insert( block->getChild(c) );
            depth_ = ( block->getChild(c)->getLevel() > depth_ ) ? block->getChild(c)->getLevel() : depth_;
         }

         for( uint_t n = 0; n != 26; ++n )
            for( uint_t c = 0; c != block->getNeighborhoodSectionSize(n); ++c )
               if( !block->getNeighbor(n,c)->hasChildren() )
                  blocksToUpdate.insert( block->getNeighbor(n,c) );
      }

      // UPDATE BLOCK NEIGHBORHOODS

      if( !blocksToSplit.empty() )
         WALBERLA_LOG_PROGRESS( "Initializing SetupBlockForest: Updating block neighborhood information ..." );

      updateNeighborhood( blocksToUpdate );

      // RESET BLOCK MARKERS

      blocksToSplitsize = int_c( blocksToSplit.size() );
#ifdef _OPENMP
      #pragma omp parallel for schedule(static)
#endif
      for( int i = 0; i < blocksToSplitsize; ++i )
         blocksToSplit[ uint_c(i) ]->setMarker( false );

      // IF NO BLOCKS WERE MARKED TO BE SPLIT, THE FOREST CREATION ALGORITHM IS FINISHED

      loop = !blocksToSplit.empty();

#ifndef NDEBUG
      checkNeighborhoodConsistency();
#endif
   }

   WALBERLA_LOG_PROGRESS( "Initializing SetupBlockForest: Creating block forest finished." );
}



void SetupBlockForest::updateNeighborhood( std::vector< SetupBlock* >& blocks ) {

   const int blockssize = int_c( blocks.size() );
#ifdef _OPENMP
   #pragma omp parallel for schedule(static)
#endif
   for( int i = 0; i < blockssize; ++i ) {
      SetupBlock* const block = blocks[ uint_c(i)  ];

      const std::vector< SetupBlock* >& neighborhood = ( block->hasFather() ) ? block->getFather()->getNeighborhood() :
                                                                                block->getNeighborhood();
      std::vector< real_t >      neighborhoodSectionBlockCenters;
      std::vector< SetupBlock* > neighborhoodSectionBlocks;

      for( uint_t n = 0; n != 26; ++n ) {

         constructNeighborhoodSectionBlockCenters( n, block->getAABB(), neighborhoodSectionBlockCenters );

         WALBERLA_ASSERT_EQUAL( neighborhoodSectionBlockCenters.size() % 3, 0 );

         for( uint_t p = 0; p != neighborhoodSectionBlockCenters.size(); p += 3 ) {

            real_t x = neighborhoodSectionBlockCenters[p];
            real_t y = neighborhoodSectionBlockCenters[p+1];
            real_t z = neighborhoodSectionBlockCenters[p+2];

            // treat periodicity
            if( x <  domain_.xMin() && periodic_[0] ) x = domain_.xMax() - domain_.xMin() + x;
            if( x >= domain_.xMax() && periodic_[0] ) x = domain_.xMin() - domain_.xMax() + x;
            if( y <  domain_.yMin() && periodic_[1] ) y = domain_.yMax() - domain_.yMin() + y;
            if( y >= domain_.yMax() && periodic_[1] ) y = domain_.yMin() - domain_.yMax() + y;
            if( z <  domain_.zMin() && periodic_[2] ) z = domain_.zMax() - domain_.zMin() + z;
            if( z >= domain_.zMax() && periodic_[2] ) z = domain_.zMin() - domain_.zMax() + z;

            SetupBlock* neighbor = nullptr;

            for( uint_t j = 0; j != neighborhood.size() && neighbor == nullptr; ++j ) {
               if( neighborhood[j]->getAABB().contains( x, y, z ) )
                  neighbor = mapPointToBlock( neighborhood[j], x, y, z );
            }
            if( neighbor == nullptr && block->hasFather() && block->getFather()->getAABB().contains( x, y, z ) )
               neighbor = mapPointToBlock( block->getFather(), x, y, z );

            if( neighborhoodSectionBlocks.empty() || neighborhoodSectionBlocks.back() != neighbor )
               neighborhoodSectionBlocks.push_back( neighbor );
         }

#ifndef NDEBUG
         for( uint_t v = 0; v != neighborhoodSectionBlocks.size(); ++v )
            for( uint_t w = v+1; w != neighborhoodSectionBlocks.size(); ++w )
               WALBERLA_ASSERT_UNEQUAL( neighborhoodSectionBlocks[v], neighborhoodSectionBlocks[w] );
         WALBERLA_ASSERT( !neighborhoodSectionBlocks.empty() );
#endif

         block->clearNeighborhoodSection(n);
         if( neighborhoodSectionBlocks.back() != nullptr ) {

#ifndef NDEBUG
            if( neighborhoodSectionBlocks.back()->getLevel() > block->getLevel() )
            {
               WALBERLA_ASSERT_EQUAL( neighborhoodSectionBlocks.size(), getBlockMaxNeighborhoodSectionSize(n) );
            }
            else
            {
               WALBERLA_ASSERT_EQUAL( neighborhoodSectionBlocks.size(), 1 );
            }
#endif
            for( uint_t j = 0; j != neighborhoodSectionBlocks.size(); ++j )
               block->addNeighbor( n, neighborhoodSectionBlocks[j] );
         }

         neighborhoodSectionBlocks.clear();
         neighborhoodSectionBlockCenters.clear();
      }
      block->assembleNeighborhood();
   }
}



void SetupBlockForest::createNeighborhood() {

   std::vector< SetupBlock* > blocks;
   getBlocks( blocks );

   const int blockssize = int_c( blocks.size() );
#ifdef _OPENMP
   #pragma omp parallel for schedule(static)
#endif
   for( int i = 0; i < blockssize; ++i ) {
      SetupBlock* const block = blocks[ uint_c(i)  ];

      std::vector< real_t >      neighborhoodSectionBlockCenters;
      std::vector< SetupBlock* > neighborhoodSectionBlocks;

      for( uint_t n = 0; n != 26; ++n ) {

         constructNeighborhoodSectionBlockCenters( n, block->getAABB(), neighborhoodSectionBlockCenters );

         WALBERLA_ASSERT_EQUAL( neighborhoodSectionBlockCenters.size() % 3, 0 );

         for( uint_t p = 0; p != neighborhoodSectionBlockCenters.size(); p += 3 ) {

            real_t x = neighborhoodSectionBlockCenters[p];
            real_t y = neighborhoodSectionBlockCenters[p+1];
            real_t z = neighborhoodSectionBlockCenters[p+2];

            // treat periodicity
            if( x <  domain_.xMin() && periodic_[0] ) x = domain_.xMax() - domain_.xMin() + x;
            if( x >= domain_.xMax() && periodic_[0] ) x = domain_.xMin() - domain_.xMax() + x;
            if( y <  domain_.yMin() && periodic_[1] ) y = domain_.yMax() - domain_.yMin() + y;
            if( y >= domain_.yMax() && periodic_[1] ) y = domain_.yMin() - domain_.yMax() + y;
            if( z <  domain_.zMin() && periodic_[2] ) z = domain_.zMax() - domain_.zMin() + z;
            if( z >= domain_.zMax() && periodic_[2] ) z = domain_.zMin() - domain_.zMax() + z;

            SetupBlock* neighbor = getBlock( x, y, z );

            if( neighborhoodSectionBlocks.empty() || neighborhoodSectionBlocks.back() != neighbor )
               neighborhoodSectionBlocks.push_back( neighbor );
         }

#ifndef NDEBUG
         for( uint_t v = 0; v != neighborhoodSectionBlocks.size(); ++v )
            for( uint_t w = v+1; w != neighborhoodSectionBlocks.size(); ++w )
               WALBERLA_ASSERT_UNEQUAL( neighborhoodSectionBlocks[v], neighborhoodSectionBlocks[w] );
         WALBERLA_ASSERT( !neighborhoodSectionBlocks.empty() );
#endif

         block->clearNeighborhoodSection(n);
         if( neighborhoodSectionBlocks.back() != nullptr ) {

#ifndef NDEBUG
            if( neighborhoodSectionBlocks.back()->getLevel() > block->getLevel() )
            {
               WALBERLA_ASSERT_EQUAL( neighborhoodSectionBlocks.size(), getBlockMaxNeighborhoodSectionSize(n) );
            }
            else
            {
               WALBERLA_ASSERT_EQUAL( neighborhoodSectionBlocks.size(), 1 );
            }
#endif
            for( uint_t j = 0; j != neighborhoodSectionBlocks.size(); ++j )
               block->addNeighbor( n, neighborhoodSectionBlocks[j] );
         }

         neighborhoodSectionBlocks.clear();
         neighborhoodSectionBlockCenters.clear();
      }
      block->assembleNeighborhood();
   }
}



SetupBlock* SetupBlockForest::mapPointToBlock( SetupBlock* const block, const real_t px, const real_t py, const real_t pz ) {

   WALBERLA_ASSERT( block->getAABB().contains( px, py, pz ) );

   if( !block->hasChildren() )
      return block;

   uint_t branchId = 0;

   const AABB& aabb = block->getChild(0)->getAABB();

   if( px >= aabb.xMax() ) ++branchId;
   if( py >= aabb.yMax() ) branchId += 2;
   if( pz >= aabb.zMax() ) branchId += 4;

   return mapPointToBlock( block->getChild( branchId ), px, py, pz );
}



void SetupBlockForest::assignAllBlocksToRootProcess()
{
   WALBERLA_LOG_PROGRESS( "Balancing SetupBlockForest: Assigning all blocks to the root process ..." );

   numberOfProcesses_       = 1;
   numberOfBufferProcesses_ = 0;

   std::vector< SetupBlock* > blocks;
   getBlocks( blocks );

   for( uint_t i = 0; i != blocks.size(); ++i )
      blocks[i]->assignTargetProcess( uint_c(0) );

   calculateProcessDistributionFinalization();
}



void SetupBlockForest::balanceLoad( const TargetProcessAssignmentFunction & function,
                                    const uint_t numberOfProcesses, const real_t minBufferProcessesFraction,
                                    const memory_t perProcessMemoryLimit,
                                    const bool reorderProcessesByBFS, const bool insertBufferProcesses )
{
   WALBERLA_LOG_PROGRESS( "Balancing SetupBlockForest: Creating a process distribution for " << numberOfProcesses_ << " process(es) ..." );

   if( minBufferProcessesFraction < real_t(0) || !(minBufferProcessesFraction < real_t(1)) )
      WALBERLA_ABORT( "Load balancing failed: \'buffer processes fraction\' must be in [0,1). "
                      "The value you provided was \'" << minBufferProcessesFraction << "\'." );
   
   // numberOfProcesses = numberOfWorkerProcesses + numberOfBufferProcesses
   //
   //  numberOfBufferProcesses                 numberOfBufferProcesses
   // ------------------------- = --------------------------------------------------- = bufferProcessesFraction
   //    numberOfProcesses         numberOfWorkerProcesses + numberOfBufferProcesses
   //
   // ==> numberOfWorkerProcesses = ( 1 - bufferProcessesFraction ) * numberOfProcesses
   //
   // integer = cast< integer >( 0.5 + floating point )  [for correct rounding]

   const uint_t numberOfWorkerProcesses = uint_c( real_c(0.5) + ( real_t(1) - minBufferProcessesFraction ) * real_c( numberOfProcesses ) );
   WALBERLA_CHECK_LESS_EQUAL( numberOfWorkerProcesses, numberOfProcesses );
   const uint_t numberOfBufferProcesses = numberOfProcesses - numberOfWorkerProcesses;
   
   balanceLoadHelper( function, numberOfProcesses, numberOfBufferProcesses, perProcessMemoryLimit, reorderProcessesByBFS, insertBufferProcesses );
}



void SetupBlockForest::balanceLoad( const TargetProcessAssignmentFunction & function,
                                    const uint_t numberOfProcesses, const uint_t numberOfBufferProcesses,
                                    const memory_t perProcessMemoryLimit,
                                    const bool reorderProcessesByBFS, const bool insertBufferProcesses )
{
   WALBERLA_LOG_PROGRESS( "Balancing SetupBlockForest: Creating a process distribution for " << numberOfProcesses_ << " process(es) ..." );
   
   balanceLoadHelper( function, numberOfProcesses, numberOfBufferProcesses, perProcessMemoryLimit, reorderProcessesByBFS, insertBufferProcesses );
}



uint_t SetupBlockForest::getMinLevel() const
{
   std::vector< const SetupBlock* > blocks;
   getBlocks( blocks );

   WALBERLA_ASSERT( !blocks.empty() );

   uint_t minLevel = blocks.front()->getLevel();

   for( auto blockIt = blocks.begin() + 1; blockIt != blocks.end(); ++blockIt )
   {
      uint_t level = ( *blockIt )->getLevel();
      if( level < minLevel )
         minLevel = level;
   }

   return minLevel;
}

uint_t SetupBlockForest::getMaxLevel() const
{
   std::vector< const SetupBlock* > blocks;
   getBlocks( blocks );

   WALBERLA_ASSERT( !blocks.empty() );

   uint_t maxLevel = blocks.front()->getLevel();

   for( auto blockIt = blocks.begin() + 1; blockIt != blocks.end(); ++blockIt )
   {
      uint_t level = ( *blockIt )->getLevel();
      if( level > maxLevel )
         maxLevel = level;
   }

   return maxLevel;
}



void SetupBlockForest::calculateProcessDistribution_Default( const uint_t        numberOfProcesses,
                                                             const memory_t      memoryLimit,
                                                             const std::string&  sfcMethod               /* = std::string( "hilbert" ) */,
                                                             const uint_t        sfcIterations           /* = 10 */,
                                                             const bool          sortByLevel             /* = false */,
                                                             const GlobalLoadBalancing::MetisConfiguration< SetupBlock >& metisConfig,
                                                             const bool          reorderProcessesByBFS   /* = false */,
                                                             const bool          insertBufferProcesses   /* = false */,
                                                             const real_t        bufferProcessesFraction /* = real_c(0) */ )
{
   WALBERLA_LOG_PROGRESS( "Balancing SetupBlockForest: Creating a process distribution for " << numberOfProcesses << " process(es) ..." );

   // error checks

   if( numberOfProcesses == 0 )
      WALBERLA_ABORT( "Load balancing failed: \'numberOfProcesses\' must be greater than 0!" );

   if( memoryLimit <= numeric_cast< memory_t >(0) )
      WALBERLA_ABORT( "Load balancing failed: You must provide a per process memory limit greater than 0!\n"
                      "                       (The memory limit you provided was \'" << memoryLimit << "\')"      );

   if( sfcMethod != "hilbert" && sfcMethod != "morton" )
      WALBERLA_ABORT( "Load balancing failed: SFC method \"" << sfcMethod << "\" unavailable "
                      "(the only available methods are \"hilbert\" and \"morton\")" )

   if( bufferProcessesFraction < real_c(0) || bufferProcessesFraction >= real_c(1) )
      WALBERLA_ABORT( "Load balancing failed: \'bufferProcessesFraction\' must be in [0,1). "
                      "The value you provided was \'" << bufferProcessesFraction << "\'." );

   // get all blocks (either in morton or in hilbert order)

   numberOfProcesses_ = numberOfProcesses;

   std::vector< SetupBlock* > blocks;

   if( sfcMethod == "hilbert" )
      getHilbertOrder( blocks );
   else
      getMortonOrder( blocks );

   // calculate process distribution

   // SOME IMPORTANT FORMULAS USED IN THE FOLLOWING LINES OF CODE
   //
   // numberOfProcesses = numberOfWorkerProcesses + numberOfBufferProcesses
   //
   //  numberOfBufferProcesses                 numberOfBufferProcesses
   // ------------------------- = --------------------------------------------------- = bufferProcessesFraction
   //    numberOfProcesses         numberOfWorkerProcesses + numberOfBufferProcesses
   //
   // ==> numberOfWorkerProcesses = ( 1 - bufferProcessesFraction ) * numberOfProcesses
   //
   // integer = cast< integer >( 0.5 + floating point )  [for correct rounding]

   uint_t numberOfWorkerProcesses = uint_c( real_c(0.5) + ( real_c(1) - bufferProcessesFraction ) * real_c( numberOfProcesses ) );

   WALBERLA_ASSERT_LESS_EQUAL( numberOfWorkerProcesses, numberOfProcesses );

   if( sortByLevel )
      numberOfWorkerProcesses = GlobalLoadBalancing::balanceSorted( blocks, sfcIterations, memoryLimit, metisConfig, numberOfWorkerProcesses );
   else
      numberOfWorkerProcesses = GlobalLoadBalancing::balance( blocks, sfcIterations, memoryLimit, metisConfig, numberOfWorkerProcesses );

   if( numberOfWorkerProcesses == 0 )
      WALBERLA_ABORT( "Load balancing failed: A distribution to " << numberOfProcesses << " processes given a memory limit of \"" << memoryLimit <<
                      "\" is impossible.\n                       (Are the memory coefficients correctly assigned to all blocks via "
                      "a callback function that was registered with \"addWorkloadMemorySUIDAssignmentFunction()\"?)" );

   numberOfBufferProcesses_ = numberOfProcesses - numberOfWorkerProcesses;

   calculateProcessDistributionFinalization( reorderProcessesByBFS, insertBufferProcesses );
}



void SetupBlockForest::calculateProcessDistribution_LevelwiseMetis( const uint_t numberOfProcesses,
   const bool reorderProcessesByBFS /*= false*/,
   const CommunicationWeightFunction & communicationWeightFunction /*= NullCommunicationWeightFunction*/ )
{
   
#ifndef WALBERLA_BUILD_WITH_METIS
   WALBERLA_ABORT( "You are trying to balance your SetupBlockForest using METIS, but you did not compile with METIS "
                   "support. Make sure that METIS is found during the CMake configuration process!");

   WALBERLA_UNUSED( numberOfProcesses );
   WALBERLA_UNUSED( reorderProcessesByBFS );
   WALBERLA_UNUSED( communicationWeightFunction );
#else
   WALBERLA_LOG_PROGRESS( "Balancing SetupBlockForest: Creating a process distribution for " << numberOfProcesses << " process(es) ..." );

   // error checks

   if( numberOfProcesses == 0 )
      WALBERLA_ABORT( "Load balancing failed: \'numberOfProcesses\' must be greater than 0!" );

   // get all blocks

   numberOfProcesses_ = numberOfProcesses;

   const uint_t minLevel = getMinLevel();
   const uint_t maxLevel = getMaxLevel();

   for( uint_t level = minLevel; level <= maxLevel; ++level )
   {
      std::vector< SetupBlock* > blocks;
      getBlocks( blocks, level );
      if( blocks.size() > numberOfProcesses_ )
      {
         blockforest::GlobalLoadBalancing::metis2( blocks, numberOfProcesses_, communicationWeightFunction );
      }
      else
      {
         uint_t targetProcess = 0;
         for( auto blockIt = blocks.begin(); blockIt != blocks.end(); ++blockIt )
         {
            (*blockIt)->assignTargetProcess( targetProcess++ );
         }

      }
      
   }

   numberOfBufferProcesses_ = 0;

   calculateProcessDistributionFinalization( reorderProcessesByBFS, false );

#endif
}



// sort functors (used in "SetupBlockForest::calculateProcessDistribution_Greedy")
class SetSorter {
public:
   SetSorter( const std::vector< workload_t >& workload ) : workload_( workload ) {}
   bool operator()( const uint_t& lhs, const uint_t& rhs ) const
   { WALBERLA_ASSERT_LESS( lhs, workload_.size() ); WALBERLA_ASSERT_LESS( rhs, workload_.size() ); return workload_[lhs] < workload_[rhs]; }
private:
   const std::vector< workload_t >& workload_;
};

struct BlockSorter {
   bool operator() ( const SetupBlock * const i, const SetupBlock * const j ) { return ( i->getWorkload() > j->getWorkload() ); }
};

void SetupBlockForest::calculateProcessDistribution_Greedy( const uint_t   numberOfProcesses,
                                                            const memory_t memoryLimit,
                                                            const bool     reorderProcessesByBFS   /* = false */,
                                                            const bool     insertBufferProcesses   /* = false */,
                                                            const real_t   bufferProcessesFraction /* = real_c(0) */ )
{
   WALBERLA_LOG_PROGRESS( "Balancing SetupBlockForest: Creating a process distribution for " << numberOfProcesses << " process(es) ..." );

   // error checks

   if( numberOfProcesses == 0 )
      WALBERLA_ABORT( "Load balancing failed: \'numberOfProcesses\' must be greater than 0!" );

   if( memoryLimit <= numeric_cast< memory_t >(0) )
      WALBERLA_ABORT( "Load balancing failed: You must provide a per process memory limit greater than 0!\n"
                      "                       (The memory limit you provided was \'" << memoryLimit << "\')"      );

   if( bufferProcessesFraction < real_c(0) || bufferProcessesFraction >= real_c(1) )
      WALBERLA_ABORT( "Load balancing failed: \'bufferProcessesFraction\' must be in [0,1). "
                      "The value you provided was \'" << bufferProcessesFraction << "\'." );

   // get all blocks

   numberOfProcesses_ = numberOfProcesses;

   std::vector< SetupBlock* > blocks;
   getBlocks( blocks );

   // calculate process distribution

   // SOME IMPORTANT FORMULAS USED IN THE FOLLOWING LINES OF CODE
   //
   // numberOfProcesses = numberOfWorkerProcesses + numberOfBufferProcesses
   //
   //  numberOfBufferProcesses                 numberOfBufferProcesses
   // ------------------------- = --------------------------------------------------- = bufferProcessesFraction
   //    numberOfProcesses         numberOfWorkerProcesses + numberOfBufferProcesses
   //
   // ==> numberOfWorkerProcesses = ( 1 - bufferProcessesFraction ) * numberOfProcesses
   //
   // integer = cast< integer >( 0.5 + floating point )  [for correct rounding]

   uint_t numberOfWorkerProcesses = uint_c( real_c(0.5) + ( real_c(1) - bufferProcessesFraction ) * real_c( numberOfProcesses ) );

   WALBERLA_ASSERT_LESS_EQUAL( numberOfWorkerProcesses, numberOfProcesses );

   std::sort( blocks.begin(), blocks.end(), BlockSorter() );

   std::vector< workload_t > workload( numberOfWorkerProcesses, workload_c(0) );
   std::vector< memory_t   > memory  ( numberOfWorkerProcesses, memory_c(0) );

   SetSorter sorter( workload );
   std::multiset< uint_t, SetSorter > distributition( sorter );
   for( uint_t i = 0; i < numberOfWorkerProcesses; ++i )
      distributition.insert( distributition.end(), i );

   numberOfWorkerProcesses = 0;
   for( auto block = blocks.begin(); block != blocks.end(); ++block )
   {
      auto process = distributition.begin();
      while( process != distributition.end() )
      {
         if( memory[ *process ] + (*block)->getMemory() <= memoryLimit )
         {
            workload[ *process ] += (*block)->getWorkload();
            memory[ *process ]   += (*block)->getMemory();
            (*block)->assignTargetProcess( *process );
            WALBERLA_ASSERT_LESS_EQUAL( *process, numberOfWorkerProcesses );
            numberOfWorkerProcesses = std::max( numberOfWorkerProcesses, *process + uint_c(1) );
            break;
         }
         ++process;
      }
      if( process == distributition.end() )
         WALBERLA_ABORT( "Load balancing failed: A distribution to " << numberOfProcesses << " processes given a memory limit of \"" << memoryLimit <<
                         "\" is impossible.\n                       (Are the memory coefficients correctly assigned to all blocks via "
                         "a callback function that was registered with \"addWorkloadMemorySUIDAssignmentFunction()\"?)" );
      uint_t p = *process;
      distributition.erase( process );
      distributition.insert( p );
   }

   WALBERLA_ASSERT_LESS_EQUAL( numberOfWorkerProcesses, numberOfProcesses );
   numberOfBufferProcesses_ = numberOfProcesses - numberOfWorkerProcesses;

   calculateProcessDistributionFinalization( reorderProcessesByBFS, insertBufferProcesses );
}



void SetupBlockForest::balanceLoadHelper( const TargetProcessAssignmentFunction & function,
                                          const uint_t numberOfProcesses, const uint_t numberOfBufferProcesses,
                                          const memory_t perProcessMemoryLimit,
                                          const bool reorderProcessesByBFS, const bool insertBufferProcesses )
{
   if( !function )
      WALBERLA_ABORT( "Load balancing failed: the load balancing callback function is empty!" );
   
   if( numberOfProcesses == 0 )
      WALBERLA_ABORT( "Load balancing failed: \'number of processes\' must be greater than 0!" );

   if( numberOfBufferProcesses >= numberOfProcesses )
      WALBERLA_ABORT( "Load balancing failed: The number of \'buffer\' processes must be smaller than the total number of processes.\n"
                      "                       The values you provided: " << numberOfProcesses << " processes and " << numberOfBufferProcesses << " \'buffer\' processes." );

   numberOfProcesses_ = numberOfProcesses;

   uint_t numberOfWorkerProcesses = numberOfProcesses - numberOfBufferProcesses;
   WALBERLA_CHECK_LESS_EQUAL( numberOfWorkerProcesses, numberOfProcesses );

   const uint_t returnedNumberOfWorkerProcesses = function( *this, numberOfWorkerProcesses, perProcessMemoryLimit );
   WALBERLA_CHECK_LESS_EQUAL( returnedNumberOfWorkerProcesses, numberOfWorkerProcesses );
   numberOfWorkerProcesses = returnedNumberOfWorkerProcesses;

   numberOfBufferProcesses_ = numberOfProcesses - numberOfWorkerProcesses;

   std::vector< SetupBlock* > blocks;
   getBlocks( blocks );

   // make sure that every process from '0' to 'numberOfWorkerProcesses-1' holds at least one block
   // and every block is assigned to a process in [0,numberOfWorkerProcesses-1]

   std::vector< bool > processHasBlocks( numberOfWorkerProcesses, false );

   for( auto block = blocks.begin(); block != blocks.end(); ++block )
   {
      WALBERLA_CHECK_LESS( (*block)->getTargetProcess(), numberOfWorkerProcesses );
      processHasBlocks[ (*block)->getTargetProcess() ] = true;
   }

   for( auto it = processHasBlocks.begin(); it != processHasBlocks.end(); ++it )
      WALBERLA_CHECK( *it );

   // make sure that the per process memory limit is satisfied

   if( perProcessMemoryLimit > memory_t(0) )
   {
      std::vector< memory_t > memory( numberOfWorkerProcesses, memory_c(0) );
      for( auto block = blocks.begin(); block != blocks.end(); ++block )
      {
         const uint_t targetProcess = (*block)->getTargetProcess();
         if( memory[ targetProcess ] + (*block)->getMemory() > perProcessMemoryLimit )
            WALBERLA_ABORT( "Load balancing failed: A distribution to " << numberOfProcesses << " processes given a memory limit of \"" << perProcessMemoryLimit <<
                            "\" is impossible.\n                       (Are the memory coefficients correctly assigned to all blocks via "
                            "a callback function that was registered with \"addWorkloadMemorySUIDAssignmentFunction()\"?)" );
         memory[ targetProcess ] += (*block)->getMemory();
      }
   }

   calculateProcessDistributionFinalization( reorderProcessesByBFS, insertBufferProcesses );
}



void SetupBlockForest::calculateProcessDistributionFinalization( const bool reorderProcessesByBFS /* = false */,
                                                                 const bool insertBufferProcesses /* = false */ )
{
   std::vector< SetupBlock* > blocks;
   getBlocks( blocks );

   // breadth-first search process reordering

   if( reorderProcessesByBFS )
   {
      WALBERLA_LOG_PROGRESS( "Balancing SetupBlockForest: Executing breadth-first search process reordering ..." );

      std::vector< std::vector< uint_t > > processNeighbors( getNumberOfWorkerProcesses() );

      GlobalLoadBalancing::prepareProcessReordering( blocks, processNeighbors );

      GlobalLoadBalancing::reorderProcessesByBFS( blocks, processNeighbors );
   }

   // adapt process distribution to buffer processes

   insertBuffersIntoProcessNetwork_ = insertBufferProcesses;

   if( insertBuffersIntoProcessNetwork_ )
   {
      WALBERLA_LOG_PROGRESS( "Balancing SetupBlockForest: Inserting buffer processes into process network ..." );

      std::vector< uint_t > bufferProcesses;

      const uint_t div = getNumberOfWorkerProcesses() / ( numberOfBufferProcesses_ + 1 );
      const uint_t mod = getNumberOfWorkerProcesses() % ( numberOfBufferProcesses_ + 1 );

      bufferProcesses.resize( numberOfBufferProcesses_, div );

      for( uint_t i = 0; i < mod; ++i )
         ++bufferProcesses[i];

      for( uint_t i = 1; i < numberOfBufferProcesses_; ++i )
         bufferProcesses[i] += bufferProcesses[i-1];

      const int blockssize = int_c( blocks.size() );
#ifdef _OPENMP
      #pragma omp parallel for schedule(static)
#endif
      for( int i = 0; i < blockssize; ++i ) {

         uint_t offset = 0;

         const uint_t process = blocks[ uint_c(i) ]->getTargetProcess();
         for( uint_t j = 0; j < numberOfBufferProcesses_ && process >= bufferProcesses[j]; ++j )
               ++offset;

         blocks[ uint_c(i)  ]->assignTargetProcess( blocks[ uint_c(i)  ]->getTargetProcess() + offset );
      }

#ifndef NDEBUG
      for( uint_t i = 0; i < numberOfBufferProcesses_; ++i )
         WALBERLA_ASSERT_LESS( bufferProcesses[i] + i, numberOfProcesses_ ); // bufferProcesses[i] + i == process IDs of buffer processes
#endif
   }

   // setup process to block mapping

   blockDistribution_.clear();
   blockDistribution_.resize( numberOfProcesses_ );
   for( uint_t i = 0; i != blocks.size(); ++i )
      blockDistribution_[ blocks[i]->getProcess() ].push_back( blocks[i] );

#ifndef NDEBUG
   std::vector< bool > processHasBlocks( numberOfProcesses_, false );

   for( uint_t i = 0; i != blocks.size(); ++i ) {
      WALBERLA_ASSERT_LESS( blocks[i]->getProcess(), numberOfProcesses_ );
      processHasBlocks[ blocks[i]->getProcess() ] = true;
   }

   for( uint_t i = 0; i != numberOfProcesses_; ++i )
      WALBERLA_ASSERT( processHasBlocks[i] == isWorkerProcess(i) );
#endif

   WALBERLA_LOG_PROGRESS( "Balancing SetupBlockForest: process distribution to " << numberOfProcesses_ << " process(es) finished!\n"
                          "- number of worker processes:       " << getNumberOfWorkerProcesses() << "\n" <<
                          "- number of empty buffer processes: " << numberOfBufferProcesses_ << "\n" <<
                          "- buffer processes are inserted into the process network: " << ( insertBuffersIntoProcessNetwork_ ? "yes\n" : "no\n" ) <<
                          "The resulting block structure looks like as follows:\n" << *this );
}



//**********************************************************************************************************************
/// \brief
///
/// For a description of the file format see BlockForestFile.h \see BlockForestFile.h
//**********************************************************************************************************************

void SetupBlockForest::saveToFile( const char* const filename ) const {

   std::ofstream file( filename, std::ofstream::binary );

   // HEADER

   uint_t offset = 0;
   std::vector< uint8_t > buffer( internal::FILE_HEADER_SIZE );

   // domain AABB

   for( uint_t i = 0; i != 3; ++i )
      offset += realToByteArray( domain_.min(i), buffer, offset );

   for( uint_t i = 0; i != 3; ++i )
      offset += realToByteArray( domain_.max(i), buffer, offset );

   // number of coarse/root blocks in each direction

   for( uint_t i = 0; i != 3; ++i ) {
      uintToByteArray( size_[i], buffer, offset, 4 );
      offset += 4;
   }

   // domain periodicity

   for( uint_t i = 0; i != 3; ++i ) {
      uintToByteArray( periodic_[i] ? uint_c(1) : uint_c(0), buffer, offset, 1 );
      ++offset;
   }

   // block forest depth (= number of levels - 1)

   uintToByteArray( depth_, buffer, offset, 1 );
   ++offset;

   // treeIdDigits (= number of bits used for storing the tree ID [tree ID marker + tree index])

   uintToByteArray( treeIdDigits_, buffer, offset, 1 );
   ++offset;

   // processIdBytes (= number of bytes required for storing process IDs)

   uintToByteArray( getProcessIdBytes(), buffer, offset, 1 );
   ++offset;

   // insertBuffersIntoProcessNetwork?

   uintToByteArray( insertBuffersIntoProcessNetwork_ ? uint_c(1) : uint_c(0), buffer, offset, 1 );
   ++offset;

   // number of processes

   uintToByteArray( numberOfProcesses_, buffer, offset, 4 );
   offset += 4;

   file.write( reinterpret_cast< const char* >( &(buffer[0]) ), numeric_cast< std::streamsize >( buffer.size() ) );
   buffer.clear();
   offset = 0;

   // SUID MAPPING

   Set<SUID> suids;

   for( uint_t i = 0; i != numberOfProcesses_; ++i ) {
      for( uint_t j = 0; j != blockDistribution_[i].size(); ++j )
         suids += blockDistribution_[i][j]->getState();
   }

   // number of SUIDs

   WALBERLA_CHECK_LESS( suids.size(), uint_c(256), "When saving the block structure to file, only 255 different SUIDs (block states) are allowed!" );

   buffer.resize(1);
   uintToByteArray( suids.size(), buffer, 0, 1 );

   file.write( reinterpret_cast< const char* >( &(buffer[0]) ), numeric_cast< std::streamsize >( buffer.size() ) );
   buffer.clear();

   std::map< SUID, std::vector< bool > > suidMap;

   const uint_t suidBytes = ( ( suids.size() % 8 == 0 ) ? ( suids.size() / 8 ) : ( suids.size() / 8 + 1 ) );

   // for every SUID ...

   uint_t i = 0;
   for( Set<SUID>::const_iterator it = suids.begin(); it != suids.end(); ++it ) {

      std::vector< bool > suidBoolVec( 8 * suidBytes );
      suidBoolVec[i] =  true;
      suidMap[ *it ] = suidBoolVec;

      // length of its identifier string

      const uint_t length = it->getIdentifier().length();
      WALBERLA_CHECK_LESS( length, 256, "SUID identifiers are allowed to consist of 255 characters at most when saving the block structure to file!" );

      buffer.resize( 1 + length );
      uintToByteArray( length, buffer, 0, 1 );

      // the identifier string

      const char* str = it->getIdentifier().c_str();
      for( uint_t j = 0; j != length; ++j )
         buffer[1+j] = *reinterpret_cast< const uint8_t* >( str + j );

      file.write( reinterpret_cast< const char* >( &(buffer[0]) ), numeric_cast< std::streamsize >( buffer.size() ) );
      buffer.clear();

      ++i;
   }

   // BLOCK DATA

   const uint_t blockIdBytes   = getBlockIdBytes();
   const uint_t processIdBytes = getProcessIdBytes();

   // for each process ...

   for( i = 0; i != numberOfProcesses_; ++i )
   {
      // number of blocks (can be '0' -> buffer process!)

      buffer.resize(2);
      uintToByteArray( uint_c( blockDistribution_[i].size() ), buffer, 0, 2 );

      file.write( reinterpret_cast< const char* >( &(buffer[0]) ), numeric_cast< std::streamsize >( buffer.size() ) );
      buffer.clear();
      offset = 0;

      std::set< uint_t > neighbors;

      // if there are blocks on this process ...

      if( !blockDistribution_[i].empty() ) {

         WALBERLA_ASSERT( isWorkerProcess(i) );

         buffer.resize( blockDistribution_[i].size() * ( blockIdBytes + suidBytes ) );

         // for each block ...

         for( uint_t j = 0; j != blockDistribution_[i].size(); ++j ) {
            const SetupBlock* const block = blockDistribution_[i][j];

            // block ID

            block->getId().toByteArray( buffer, offset, blockIdBytes );
            offset += blockIdBytes;

            // block state (SUID set)

            if( suidBytes > 0 ) {
               std::vector< bool > suidBoolVec( 8 * suidBytes );

               const Set<SUID>& state = block->getState();
               for( Set<SUID>::const_iterator suid = state.begin(); suid != state.end(); ++suid ) {
                  WALBERLA_ASSERT( suidMap.find( *suid ) != suidMap.end() );
                  //Elementwise OR of all elements
                  for (uint_t k = 0;k  < suidBoolVec.size(); ++k) {
                     suidBoolVec[k] = suidBoolVec[k] | suidMap.find( *suid )->second[k];
                  }
               }

               boolVectorToByteArray( suidBoolVec, buffer, offset );
               offset += suidBytes;
            }

            for( uint_t k = 0; k != block->getNeighborhoodSize(); ++k )
               if( block->getNeighbor(k)->getProcess() != i )
                  neighbors.insert( block->getNeighbor(k)->getProcess() );
         }

         file.write( reinterpret_cast< const char* >( &(buffer[0]) ), numeric_cast< std::streamsize >( buffer.size() ) );
         buffer.clear();
         offset = 0;
      }
      else
      {
         WALBERLA_ASSERT( isBufferProcess(i) );
         WALBERLA_ASSERT_GREATER( i, 0 );
         if( ( i + 1 ) < numberOfProcesses_ && insertBuffersIntoProcessNetwork_ )
            WALBERLA_ASSERT( !( isBufferProcess( i - 1 ) && isWorkerProcess( i + 1 ) ) );

         if( insertBuffersIntoProcessNetwork_ )
            neighbors.insert( i - 1 );
      }

      // process neighborhood (= all neighboring processes)

      if( insertBuffersIntoProcessNetwork_ && ( i + 1 ) < numberOfProcesses_ && isBufferProcess( i + 1 ) )
         neighbors.insert( i + 1 );

      buffer.resize( 2 + neighbors.size() * processIdBytes );

      uintToByteArray( uint_c( neighbors.size() ), buffer, offset, 2 );
      offset += 2;

      for( std::set< uint_t >::iterator it = neighbors.begin(); it != neighbors.end(); ++it ) {
         uintToByteArray( *it, buffer, offset, processIdBytes );
         offset += processIdBytes;
      }

      file.write( reinterpret_cast< const char* >( &(buffer[0]) ), numeric_cast< std::streamsize >( buffer.size() ) );
      buffer.clear();
      offset = 0;
   }

   file.close();
}


void SetupBlockForest::writeVTKOutput( const std::string & filestem ) const
{
   std::ostringstream oss;
   oss << filestem << ".vtk";

   std::ofstream outfile( oss.str().c_str() );

   outfile << "# vtk DataFile Version 3.0\n"
           << "SetupBlockForest\n"
           << "ASCII\n\n"
           << "DATASET UNSTRUCTURED_GRID\n\n";

   uint_t allocatedBlocks = 0;

   std::vector< const SetupBlock* > blocks;
   getBlocks( blocks );

   for( uint_t i = 0; i != blocks.size(); ++i )
      if( blocks[i]->getMemory() > 0 ) ++allocatedBlocks;

   outfile << "POINTS " << ( 8 * allocatedBlocks ) << " " << typeToString<real_t>() << "\n\n";

   for( uint_t i = 0; i != blocks.size(); ++i ) {
      const SetupBlock* const block = blocks[i];
      if( block->getMemory() > 0 ) {
         for( uint_t z = 0; z != 2; ++z ) {
            for( uint_t y = 0; y != 2; ++y ) {
               for( uint_t x = 0; x != 2; ++x ) {
                  outfile << ( ( x == 0 ) ? block->getAABB().xMin() : block->getAABB().xMax() ) << " "
                          << ( ( y == 0 ) ? block->getAABB().yMin() : block->getAABB().yMax() ) << " "
                          << ( ( z == 0 ) ? block->getAABB().zMin() : block->getAABB().zMax() ) << "\n";
               }
            }
         }
         outfile << std::endl;
      }
   }

   outfile << "\n\nCELLS " << allocatedBlocks << " " << ( 9 * allocatedBlocks ) << "\n\n";

   for( uint_t i = 0, c = 0; i != allocatedBlocks; ++i ) {

      outfile << "8";

      for( uint_t j = 0; j != 8; ++j, ++c )
         outfile << " " << c;

      outfile << std::endl;
   }

   outfile << "\n\nCELL_TYPES " << allocatedBlocks << "\n\n";

   for( uint_t i = 0; i != allocatedBlocks; ++i )
      outfile << "11\n";

   outfile << "\n\nCELL_DATA " << allocatedBlocks;

   outfile << "\n\nSCALARS workload double 1"
           <<   "\nLOOKUP_TABLE default\n";

   for( uint_t i = 0; i != blocks.size(); ++i ) {
      if( blocks[i]->getMemory() > 0 )
         outfile << blocks[i]->getWorkload() << "\n";
   }

   outfile << "\n\nSCALARS process int 1"
           <<   "\nLOOKUP_TABLE default\n";

   for( uint_t i = 0; i != blocks.size(); ++i ) {
      if( blocks[i]->getMemory() > 0 )
         outfile << blocks[i]->getProcess() << "\n";
   }

#ifdef WALBERLA_BLOCKFOREST_PRIMITIVE_BLOCKID

   outfile << "\n\nSCALARS blockId unsigned_long 1"
           <<   "\nLOOKUP_TABLE default\n";

   for( uint_t i = 0; i != blocks.size(); ++i ) {
      if( blocks[i]->getMemory() > 0 )
         outfile << blocks[i]->getId().getTreeId() << "\n";
   }

#endif // WALBERLA_BLOCKFOREST_PRIMITIVE_BLOCKID

   outfile << "\n\nSCALARS level double 1"
           <<   "\nLOOKUP_TABLE colors\n";

   for( uint_t i = 0; i != blocks.size(); ++i ) {
      if( blocks[i]->getMemory() > 0 )
         outfile << (( depth_ == 0 ) ? 0 : ( static_cast< double >( blocks[i]->getLevel() ) / static_cast< double >( depth_ ) )) << "\n";
   }

   outfile << "\n\nLOOKUP_TABLE colors " << ( depth_ + 1 );

   for( uint_t i = 0; i != depth_ + 1; ++i ) {

      outfile << "\n";
      writeColor( outfile, i, PARAVIEW );
   }

   outfile << std::endl;
   outfile.close();
}



void SetupBlockForest::writeCSV( const std::string & filestem ) const
{
   std::ostringstream oss;
   oss << filestem << ".csv";

   // each entry in the csv file: process-ID, number of blocks, number of blocks on level 0, number of blocks on level 1, ...

   std::ofstream outfile( oss.str().c_str() );

   uint_t c = uint_t(0);

   for( auto blocks = blockDistribution_.begin(); blocks != blockDistribution_.end(); ++blocks, ++c )
   {
      outfile << c << "," << blocks->size();
      for( uint_t i = 0; i <= depth_; ++i )
      {
         uint_t n = uint_t(0);
         for( auto block = blocks->begin(); block != blocks->end(); ++block )
            if( (*block)->getLevel() == i ) ++n;
         outfile << "," << n;
      }
      outfile << "\n";
   }

   outfile << std::endl;
   outfile.close();
}



void SetupBlockForest::toStream( std::ostream & os ) const
{
   uint_t discardedRootBlocks = uint_t(0);
   for( auto block = forest_.begin(); block != forest_.end(); ++block )
      if( *block == NULL ) ++discardedRootBlocks;

   os << "- AABB: " << domain_ << "\n"
      << "- initial decomposition: " << size_[0] << " x " << size_[1] << " x " << size_[2] << " (= forest size)\n"
      << std::boolalpha << "- periodicity: " << periodic_[0] << " x " << periodic_[1] << " x " << periodic_[2] << "\n"
      << "- number of blocks discarded from the initial grid: " << discardedRootBlocks << " (= "
                                                                << ( real_t(100) * real_c(discardedRootBlocks) / real_c(size_[0]*size_[1]*size_[2]) ) << " %)\n"
      << "- number of levels: " << getNumberOfLevels() << "\n"
      << "- tree ID digits: " << treeIdDigits_ << " (-> block ID bytes = " << getBlockIdBytes() << ")\n"
      << "- total number of blocks: " << numberOfBlocks_ << "\n";

   if( numberOfProcesses_ == uint_t(0) )
   {
      os << "- blocks have not yet been distributed to processes";

      if( depth_ > uint_t(0) )
      {
         os << "\n- distribution of space/memory/work to different grid levels:\n";

         std::vector< real_t > space         ( depth_ + 2, real_t(0) );
         std::vector< real_t > memory        ( depth_ + 2, real_t(0) );
         std::vector< real_t > workload      ( depth_ + 2, real_t(0) );
         std::vector< uint_t > numberOfBlocks( depth_ + 2, uint_t(0) );

         for( auto block = begin(); block != end(); ++block )
         {
            const auto level = block->getLevel();

            space[ level ]          += block->getAABB().volume();
            memory[ level ]         += real_c( block->getMemory() );
            workload[ level ]       += real_c( block->getWorkload() );
            numberOfBlocks[ level ] += uint_t(1);

            space.back()          += block->getAABB().volume();
            memory.back()         += real_c( block->getMemory() );
            workload.back()       += real_c( block->getWorkload() );
            numberOfBlocks.back() += uint_t(1);
         }

         WALBERLA_ASSERT_EQUAL( numberOfBlocks_, numberOfBlocks.back() );

         for( uint_t l = uint_t(0); l <= depth_; ++l )
         {
            os << "   + level " << l << ":\n"
               << "      - " << numberOfBlocks[l] << " blocks ...\n"
               << "      - ... cover " << ( real_t(100) * space[l] / space.back() ) << " % of the total simulation space\n"
               << "      - ... account for " << ( real_t(100) * memory[l] / memory.back() ) << " % of the total memory foot print\n"
               << "      - ... generate " << ( real_t(100) * workload[l] / workload.back() ) << " % of the total workload";
            if( l != depth_ ) os << "\n";
         }
      }
   }
   else
   {
      std::vector< math::Sample > blocks  ( depth_ + 2 );
      std::vector< math::Sample > memory  ( depth_ + 2 );
      std::vector< math::Sample > workload( depth_ + 2 );

      WALBERLA_ASSERT_EQUAL( numberOfProcesses_, blockDistribution_.size() );

      for( auto process = blockDistribution_.begin(); process != blockDistribution_.end(); ++process )
      {
         std::vector< uint_t >     processBlocks  ( depth_ + 2, uint_t(0) );
         std::vector< memory_t >   processMemory  ( depth_ + 2, memory_t(0) );
         std::vector< workload_t > processWorkload( depth_ + 2, workload_t(0) );

         for( auto block = process->begin(); block != process->end(); ++block )
         {
            const auto level = (*block)->getLevel();

            processBlocks[level]   += uint_t(1);
            processMemory[level]   += (*block)->getMemory();
            processWorkload[level] += (*block)->getWorkload();

            processBlocks.back()   += uint_t(1);
            processMemory.back()   += (*block)->getMemory();
            processWorkload.back() += (*block)->getWorkload();
         }

         for( uint_t i = 0; i < depth_ + 2; ++i )
         {
            blocks[i].insert( real_c( processBlocks[i] ) );
            memory[i].insert( real_c( processMemory[i] ) );
            workload[i].insert( real_c( processWorkload[i] ) );
         }
      }

      os << "- number of processes: " << numberOfProcesses_ << " (" << getNumberOfWorkerProcesses() << " worker process(es) / "
                                                                    << getNumberOfBufferProcesses() << " empty buffer process(es))\n"
         << "- buffer processes are inserted into the process network: " << ( insertBuffersIntoProcessNetwork_ ? "yes\n" : "no\n" )
         << "- process ID bytes: " << getProcessIdBytes() << "\n"
         << "- blocks/memory/workload per process:\n"
         << "   + blocks:\n"
         << "      - min       = " << blocks.back().min() << "\n"
         << "      - max       = " << blocks.back().max() << "\n"
         << "      - avg       = " << blocks.back().mean() << "\n"
         << "      - stdDev    = " << blocks.back().stdDeviation() << "\n"
         << "      - relStdDev = " << blocks.back().relativeStdDeviation() << "\n"
         << "   + memory:\n"
         << "      - min       = " << memory.back().min() << "\n"
         << "      - max       = " << memory.back().max() << "\n"
         << "      - avg       = " << memory.back().mean() << "\n"
         << "      - stdDev    = " << memory.back().stdDeviation() << "\n"
         << "      - relStdDev = " << memory.back().relativeStdDeviation() << "\n"
         << "   + workload:\n"
         << "      - min       = " << workload.back().min() << "\n"
         << "      - max       = " << workload.back().max() << "\n"
         << "      - avg       = " << workload.back().mean() << "\n"
         << "      - stdDev    = " << workload.back().stdDeviation() << "\n"
         << "      - relStdDev = " << workload.back().relativeStdDeviation();

      if( depth_ > uint_t(0) )
      {
         os << "\n- distribution of space/memory/work to different grid levels:\n";

         std::vector< real_t > space( depth_ + 2, real_t(0) );
         std::vector< uint_t > numberOfBlocks( depth_ + 2, uint_t(0) );

         for( auto block = begin(); block != end(); ++block )
         {
            const auto level = block->getLevel();

            space[ level ] += block->getAABB().volume();
            numberOfBlocks[ level ] += uint_t(1);

            space.back() += block->getAABB().volume();
            numberOfBlocks.back() += uint_t(1);
         }

         WALBERLA_ASSERT_EQUAL( numberOfBlocks_, numberOfBlocks.back() );

         for( uint_t l = uint_t(0); l <= depth_; ++l )
         {
            os << "   + level " << l << ":\n"
               << "      - " << numberOfBlocks[l] << " blocks ...\n"
               << "      - ... cover " << ( real_t(100) * space[l] / space.back() ) << " % of the total simulation space\n"
               << "      - ... account for " << ( real_t(100) * memory[l].sum() / memory.back().sum() ) << " % of the total memory foot print\n"
               << "      - ... generate " << ( real_t(100) * workload[l].sum() / workload.back().sum() ) << " % of the total work load\n"
               << "      - blocks per process:\n"
               << "         + min       = " << blocks[l].min() << "\n"
               << "         + max       = " << blocks[l].max() << "\n"
               << "         + avg       = " << blocks[l].mean() << "\n"
               << "         + stdDev    = " << blocks[l].stdDeviation() << "\n"
               << "         + relStdDev = " << blocks[l].relativeStdDeviation() << "\n"
               << "      - memory per process:\n"
               << "         + min       = " << memory[l].min() << "\n"
               << "         + max       = " << memory[l].max() << "\n"
               << "         + avg       = " << memory[l].mean() << "\n"
               << "         + stdDev    = " << memory[l].stdDeviation() << "\n"
               << "         + relStdDev = " << memory[l].relativeStdDeviation() << "\n"
               << "      - workload per process:\n"
               << "         + min       = " << workload[l].min() << "\n"
               << "         + max       = " << workload[l].max() << "\n"
               << "         + avg       = " << workload[l].mean() << "\n"
               << "         + stdDev    = " << workload[l].stdDeviation() << "\n"
               << "         + relStdDev = " << workload[l].relativeStdDeviation();
            if( l != depth_ ) os << "\n";
         }

         std::vector< const SetupBlock* > fineBlocks;
         getBlocks( fineBlocks, depth_ );

         real_t minSpace    = fineBlocks[0]->getAABB().volume();
         real_t maxSpace    = fineBlocks[0]->getAABB().volume();
         real_t minMemory   = real_c( fineBlocks[0]->getMemory() );
         real_t maxMemory   = real_c( fineBlocks[0]->getMemory() );
         real_t minWorkload = real_c( fineBlocks[0]->getWorkload() );
         real_t maxWorkload = real_c( fineBlocks[0]->getWorkload() );

         for( uint_t i = 1; i < fineBlocks.size(); ++i )
         {
            minSpace    = std::min( minSpace,    fineBlocks[i]->getAABB().volume() );
            maxSpace    = std::max( maxSpace,    fineBlocks[i]->getAABB().volume() );
            minMemory   = std::min( minMemory,   real_c( fineBlocks[i]->getMemory() ) );
            maxMemory   = std::max( maxMemory,   real_c( fineBlocks[i]->getMemory() ) );
            minWorkload = std::min( minWorkload, real_c( fineBlocks[i]->getWorkload() ) );
            maxWorkload = std::max( maxWorkload, real_c( fineBlocks[i]->getWorkload() ) );
         }

         WALBERLA_ASSERT( realIsEqual( minSpace, maxSpace ) );

         if( realIsEqual( minMemory, maxMemory ) || realIsEqual( minWorkload, maxWorkload ) )
         {
            os << "\n- using a uniform decomposition with a resolution equal to the finest level, one would ...\n";

            if( realIsEqual( minMemory, maxMemory ) )
            {
               const real_t memoryPerFineSpace = minMemory / minSpace;
               const real_t everythingFineMemory = space.back() * memoryPerFineSpace;
               os << "   + ... need " << ( everythingFineMemory / memory.back().sum() ) << " times the memory";
            }

            if( realIsEqual( minWorkload, maxWorkload ) )
            {
               const real_t workloadPerFineSpace = minWorkload / minSpace;
               const real_t everythingFineWorkload = space.back() * workloadPerFineSpace;
               os << "\n   + ... generate " << ( everythingFineWorkload / workload.back().sum() ) << " times the workload";
            }
         }
      }
   }
}



#ifndef NDEBUG

void SetupBlockForest::checkNeighborhoodConsistency() const {

   std::vector< const SetupBlock* > blocks;
   getBlocks( blocks );

   const int blockssize = int_c( blocks.size() );
#ifdef _OPENMP
   #pragma omp parallel for schedule(static)
#endif
   for( int i = 0; i < blockssize; ++i ) {

      const SetupBlock* const block = blocks[uint_c(i)];

      std::vector< real_t > neighborhoodSectionBlockCenters;

      for( uint_t n = 0; n != 26; ++n ) {

         std::vector< bool > hit( block->getNeighborhoodSectionSize(n), false );

         constructNeighborhoodSectionBlockCenters( n, block->getAABB(), neighborhoodSectionBlockCenters );

         WALBERLA_ASSERT_EQUAL( neighborhoodSectionBlockCenters.size() % 3, 0 );

         for( uint_t p = 0; p != neighborhoodSectionBlockCenters.size(); p += 3 ) {

            real_t x = neighborhoodSectionBlockCenters[p];
            real_t y = neighborhoodSectionBlockCenters[p+1];
            real_t z = neighborhoodSectionBlockCenters[p+2];

            // treat periodicity
            if( x <  domain_.xMin() && periodic_[0] ) x = domain_.xMax() - domain_.xMin() + x;
            if( x >= domain_.xMax() && periodic_[0] ) x = domain_.xMin() - domain_.xMax() + x;
            if( y <  domain_.yMin() && periodic_[1] ) y = domain_.yMax() - domain_.yMin() + y;
            if( y >= domain_.yMax() && periodic_[1] ) y = domain_.yMin() - domain_.yMax() + y;
            if( z <  domain_.zMin() && periodic_[2] ) z = domain_.zMax() - domain_.zMin() + z;
            if( z >= domain_.zMax() && periodic_[2] ) z = domain_.zMin() - domain_.zMax() + z;

            bool noHit = true;
            for( uint_t c = 0; c != block->getNeighborhoodSectionSize(n) && noHit; ++c ) {
               if( block->getNeighbor(n,c)->getAABB().contains(x,y,z) ) {
                  hit[c] = true;
                  noHit = false;
               }
            }

            // either one neighbor must be hit OR the block is located at the border of the (non-periodic) simulation domain
            if( noHit )
               WALBERLA_ASSERT_NULLPTR( getBlock(x,y,z) );
         }

         // every neighbor must be hit by at least one point
         for( uint_t c = 0; c != block->getNeighborhoodSectionSize(n); ++c )
            WALBERLA_ASSERT( hit[c] );

         neighborhoodSectionBlockCenters.clear();
      }
   }
}

#endif



} // namespace blockforest
} // namespace walberla
