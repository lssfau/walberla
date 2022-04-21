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
//! \file BlockForestInitialization.cpp
//! \ingroup mesh
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//
//======================================================================================================================

#include "BlockForestInitialization.h"

#include "blockforest/Initialization.h"
#include "blockforest/loadbalancing/StaticCurve.h"
#include "blockforest/loadbalancing/StaticParMetis.h"

#include "waLBerlaDefinitions.h"

#include "core/math/IntegerFactorization.h"
#include "core/math/Primes.h"

#include <map>

namespace walberla {
namespace mesh {

static inline uint_t uintAbsDiff( const uint_t x, const uint_t y )
{
   return x > y ? x - y : y - x;
}

static inline void compareAABB( const AABB & oldAABB, const AABB & newAABB )
{
   if(oldAABB != newAABB)
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Domain AABB has been adapted to requested block size:\n" \
                                 "Old AABB: " << oldAABB   << "\n" \
                                 "New AABB: " << newAABB << "\n" )
   }
}

ComplexGeometryBlockforestCreator::ComplexGeometryBlockforestCreator( const AABB & aabb, const blockforest::SetupBlockForest::RootBlockExclusionFunction & rootBlockExclusionFunction )
   : aabb_(aabb), maxIterations_(25), acceptableRelativeError_( real_t(0.1) ), maxBlockSkewness_(2.0),
     processMemoryLimit_( real_t( 0.0 ) ), periodicity_( false, false, false ), rootBlockExclusionFunction_ ( rootBlockExclusionFunction ),
     workloadMemorySUIDAssignmentFunction_( blockforest::uniformWorkloadAndMemoryAssignment ),
#ifdef WALBERLA_BUILD_WITH_PARMETIS
     targetProcessAssignmentFunction_( blockforest::StaticLevelwiseParMetis() )
#else
     targetProcessAssignmentFunction_( blockforest::StaticLevelwiseCurveBalanceWeighted() )
#endif

{
}

shared_ptr<SetupBlockForest> ComplexGeometryBlockforestCreator::createSetupBlockForest( const uint_t targetNumRootBlocks, const uint_t numProcesses ) const
{
   std::set< uint_t > blockGridSizeTested;

   uint_t sizeBlockGrid = targetNumRootBlocks * uint_t(2);

   uint_t          bestSizeBlockGrid = sizeBlockGrid;
   uint_t          bestNumRootBlocks = std::numeric_limits<uint_t>::max();
   Vector3<uint_t> bestSizeBlockGrid3D( std::numeric_limits<uint_t>::max(), std::numeric_limits<uint_t>::max(), std::numeric_limits<uint_t>::max() );

   enum Dir { UP, DOWN };

   Dir dir = UP;

   uint_t i = 0;
   while( i < maxIterations_ )
   {
      while( blockGridSizeTested.find( sizeBlockGrid ) != blockGridSizeTested.end() )
      {
         switch(dir)
         {
         case UP:   ++sizeBlockGrid; break;
         case DOWN: --sizeBlockGrid; break;
         }

         if( sizeBlockGrid == uint_t(0) )
         {
            sizeBlockGrid = uint_t(1);
            dir = UP;
         }
      }

      const Vector3<uint_t> sizeBlockGrid3D = math::getFactors3D( sizeBlockGrid, aabb_.sizes() );

      if( maxBlockSkewness_ * real_c(sizeBlockGrid3D.min() ) < real_c( sizeBlockGrid3D.max() ) )
      {
         blockGridSizeTested.insert( sizeBlockGrid );
         continue;
      }

      const uint_t numRootBlocks = findNumBlocks( sizeBlockGrid3D );

      blockGridSizeTested.insert( sizeBlockGrid );

      if( uintAbsDiff( numRootBlocks, targetNumRootBlocks ) < uintAbsDiff( bestNumRootBlocks, targetNumRootBlocks ) )
      {
         bestNumRootBlocks   = numRootBlocks;
         bestSizeBlockGrid   = sizeBlockGrid;
         bestSizeBlockGrid3D = sizeBlockGrid3D;
      }

      const real_t factor = real_c( numRootBlocks ) / real_c( targetNumRootBlocks );

      if( std::fabs( factor - real_t(1) ) < acceptableRelativeError_ )
      {
         break;
      }
      if( numRootBlocks < targetNumRootBlocks )
      {
         sizeBlockGrid = uint_c( real_c( sizeBlockGrid ) / factor + real_t(0.5) );
         dir = UP;
      }
      else if( numRootBlocks > targetNumRootBlocks )
      {
         sizeBlockGrid = uint_c( real_c( sizeBlockGrid ) / factor + real_t(0.5) );
         dir = DOWN;
      }
      else
      {
         break;
      }

      ++i;
   }

   WALBERLA_LOG_INFO_ON_ROOT( "Using a block grid of size " << bestSizeBlockGrid3D << " (" << bestSizeBlockGrid << ") resulting in " << bestNumRootBlocks << " root blocks." )

   auto setupBlockForest = make_shared<SetupBlockForest>();
   setupBlockForest->addRootBlockExclusionFunction( rootBlockExclusionFunction_ );
   setupBlockForest->addWorkloadMemorySUIDAssignmentFunction( workloadMemorySUIDAssignmentFunction_ );

   setupBlockForest->init( aabb_, bestSizeBlockGrid3D[0], bestSizeBlockGrid3D[1], bestSizeBlockGrid3D[2], periodicity_[0], periodicity_[1], periodicity_[2] );

   setupBlockForest->balanceLoad( targetProcessAssignmentFunction_, numProcesses, uint_t(0), processMemoryLimit_, true, false );

   return setupBlockForest;
}



 shared_ptr<SetupBlockForest> ComplexGeometryBlockforestCreator::createSetupBlockForest( const Vector3<real_t> & blockSize, const uint_t numProcesses ) const
{
   Vector3<uint_t> numBlocks;

   for( uint_t i = uint_t(0); i < uint_t(3); ++i )
      numBlocks[i] = uint_c( std::ceil( aabb_.size( i ) / blockSize[i] ) );

   AABB newAABB( real_t(0), real_t(0), real_t(0),
                 real_c( numBlocks[0] ) * blockSize[0], real_c( numBlocks[1] ) * blockSize[1], real_c( numBlocks[2] ) * blockSize[2] );

   newAABB.translate( aabb_.center() - newAABB.center() );
   compareAABB(aabb_, newAABB);

   auto setupBlockForest = make_shared<SetupBlockForest>();
   setupBlockForest->addRootBlockExclusionFunction( rootBlockExclusionFunction_ );
   setupBlockForest->addWorkloadMemorySUIDAssignmentFunction( workloadMemorySUIDAssignmentFunction_ );

   setupBlockForest->init( newAABB, numBlocks[0], numBlocks[1], numBlocks[2], periodicity_[0], periodicity_[1], periodicity_[2] );

   setupBlockForest->balanceLoad( targetProcessAssignmentFunction_, numProcesses, uint_t(0), processMemoryLimit_, true, false );

   return setupBlockForest;
}

shared_ptr<BlockForest> ComplexGeometryBlockforestCreator::createBlockForest( const uint_t targetNumRootBlocks ) const
{
   shared_ptr< blockforest::SetupBlockForest> setupBlockForest = createSetupBlockForest( targetNumRootBlocks );

   auto blockForest = make_shared< blockforest::BlockForest >( MPIManager::instance()->rank(), *setupBlockForest );

   return blockForest;
}

shared_ptr<BlockForest> ComplexGeometryBlockforestCreator::createBlockForest( const Vector3<real_t> & blockSize ) const
{
   shared_ptr< blockforest::SetupBlockForest> setupBlockForest = createSetupBlockForest( blockSize );

   auto blockForest = make_shared< blockforest::BlockForest >( MPIManager::instance()->rank(), *setupBlockForest );

   return blockForest;
}

uint_t ComplexGeometryBlockforestCreator::findNumBlocks( const Vector3<uint_t> & numRootBlocks3D ) const
{
   WALBERLA_LOG_DEVEL_ON_ROOT( "Testing block grid " << numRootBlocks3D )

   SetupBlockForest setupBlockForest;
   setupBlockForest.addRootBlockExclusionFunction( rootBlockExclusionFunction_ );

   setupBlockForest.init( aabb_, numRootBlocks3D[0], numRootBlocks3D[1], numRootBlocks3D[2], periodicity_[0], periodicity_[1], periodicity_[2] );

   WALBERLA_LOG_DEVEL_ON_ROOT( "Testing block grid " << numRootBlocks3D << " resulted in " << setupBlockForest.getNumberOfBlocks() )

   return uint_c( setupBlockForest.getNumberOfBlocks() );

}















ComplexGeometryStructuredBlockforestCreator::ComplexGeometryStructuredBlockforestCreator( const AABB & aabb, const Vector3<real_t> & cellSize, const blockforest::SetupBlockForest::RootBlockExclusionFunction & rootBlockExclusionFunction )
   : aabb_(aabb), cellSize_( cellSize ), maxIterations_(25), acceptableRelativeError_( real_t(0.1) ),
     processMemoryLimit_( real_t( 0.0 ) ), periodicity_( false, false, false ),
     rootBlockExclusionFunction_ ( rootBlockExclusionFunction ),
     workloadMemorySUIDAssignmentFunction_( blockforest::uniformWorkloadAndMemoryAssignment ),
#ifdef WALBERLA_BUILD_WITH_PARMETIS
     targetProcessAssignmentFunction_( blockforest::StaticLevelwiseParMetis() )
#else
     targetProcessAssignmentFunction_( blockforest::StaticLevelwiseCurveBalanceWeighted() )
#endif
{
}


shared_ptr<SetupBlockForest> ComplexGeometryStructuredBlockforestCreator::createSetupBlockForest( const uint_t targetNumRootBlocks, const uint_t numProcesses ) const
{

   Vector3<uint_t> numCells = domainSizeCells();
   real_t domainVolume = real_c( numCells[0] ) * real_c( numCells[1] ) * real_c( numCells[2] );

   std::set< uint_t > blockSizeTested;

   uint_t blockSize = uint_c( std::pow( domainVolume / ( real_t(2) * real_t( targetNumRootBlocks ) ), real_t(1) / real_t(3) ) + 0.5 );

   uint_t          bestBlockSize = blockSize;
   uint_t          bestNumRootBlocks = std::numeric_limits<uint_t>::max();

   enum Dir { UP, DOWN };

   Dir dir = UP;

   uint_t i = 0;
   while( i < maxIterations_ )
   {
      while( blockSizeTested.find( blockSize ) != blockSizeTested.end() )
      {
         switch(dir)
         {
         case UP:   ++blockSize; break;
         case DOWN: --blockSize; break;
         }

         if( blockSize == uint_t(0) )
         {
            blockSize = uint_t(1);
            dir = UP;
         }
      }

      const uint_t numRootBlocks = findNumBlocks( Vector3<uint_t>( blockSize, blockSize, blockSize ) );

      blockSizeTested.insert( blockSize );

      if( uintAbsDiff( numRootBlocks, targetNumRootBlocks ) < uintAbsDiff( bestNumRootBlocks, targetNumRootBlocks ) )
      {
         bestNumRootBlocks   = numRootBlocks;
         bestBlockSize       = blockSize;
      }

      const real_t factor = real_c( numRootBlocks ) / real_c( targetNumRootBlocks );

      if( std::fabs( factor - real_t(1) ) < acceptableRelativeError_ )
      {
         break;
      }
      if( numRootBlocks < targetNumRootBlocks )
      {
         blockSize = uint_c( real_c( blockSize ) * std::pow( factor, real_t(1) / real_t(3) ) + real_t(0.5) );
         dir = UP;
      }
      else if( numRootBlocks > targetNumRootBlocks )
      {
         blockSize = uint_c( real_c( blockSize ) * std::pow( factor, real_t(1) / real_t(3) ) + real_t(0.5) );
         dir = DOWN;
      }
      else
      {
         break;
      }

      ++i;
   }

   WALBERLA_LOG_INFO_ON_ROOT( "Using a block of size " << bestBlockSize << " resulting in " << bestNumRootBlocks << " root blocks." )

   Vector3<uint_t> numBlocks( uint_c( std::ceil( real_c( numCells[0] ) / real_c( bestBlockSize ) ) ),
                              uint_c( std::ceil( real_c( numCells[1] ) / real_c( bestBlockSize ) ) ),
                              uint_c( std::ceil( real_c( numCells[2] ) / real_c( bestBlockSize ) ) ) );

   AABB newAABB( real_t(0), real_t(0), real_t(0),
                 real_c( numBlocks[0] * bestBlockSize ) * cellSize_[0],
                 real_c( numBlocks[1] * bestBlockSize ) * cellSize_[1],
                 real_c( numBlocks[2] * bestBlockSize ) * cellSize_[2] );

   newAABB.translate( aabb_.center() - newAABB.center() );
   compareAABB(aabb_, newAABB);

   auto setupBlockForest = make_shared<SetupBlockForest>();
   setupBlockForest->addRootBlockExclusionFunction( rootBlockExclusionFunction_ );
   setupBlockForest->addWorkloadMemorySUIDAssignmentFunction( workloadMemorySUIDAssignmentFunction_ );

   if( refinementSelectionFunction_ )
      setupBlockForest->addRefinementSelectionFunction( refinementSelectionFunction_ );

   setupBlockForest->init( newAABB, numBlocks[0], numBlocks[1], numBlocks[2], periodicity_[0], periodicity_[1], periodicity_[2] );

   setupBlockForest->balanceLoad( targetProcessAssignmentFunction_, numProcesses, uint_t(0), processMemoryLimit_, true, false );

   return setupBlockForest;
}



shared_ptr<SetupBlockForest> ComplexGeometryStructuredBlockforestCreator::createSetupBlockForest( const Vector3<uint_t> & blockSize, const uint_t numProcesses ) const
{
   Vector3<uint_t> numCells = domainSizeCells();

   Vector3<uint_t> numBlocks( uint_c( std::ceil( real_c( numCells[0] ) / real_c( blockSize[0] ) ) ),
                              uint_c( std::ceil( real_c( numCells[1] ) / real_c( blockSize[1] ) ) ),
                              uint_c( std::ceil( real_c( numCells[2] ) / real_c( blockSize[2] ) ) ) );

   AABB newAABB( real_t(0), real_t(0), real_t(0),
                 real_c( numBlocks[0] * blockSize[0] ) * cellSize_[0],
                 real_c( numBlocks[1] * blockSize[1] ) * cellSize_[1],
                 real_c( numBlocks[2] * blockSize[2] ) * cellSize_[2] );

   newAABB.translate( aabb_.center() - newAABB.center() );
   compareAABB(aabb_, newAABB);

   auto setupBlockForest = make_shared<SetupBlockForest>();
   setupBlockForest->addRootBlockExclusionFunction( rootBlockExclusionFunction_ );
   setupBlockForest->addWorkloadMemorySUIDAssignmentFunction( workloadMemorySUIDAssignmentFunction_ );

   if( refinementSelectionFunction_ )
      setupBlockForest->addRefinementSelectionFunction( refinementSelectionFunction_ );

   setupBlockForest->init( newAABB, numBlocks[0], numBlocks[1], numBlocks[2], periodicity_[0], periodicity_[1], periodicity_[2] );

   setupBlockForest->balanceLoad( targetProcessAssignmentFunction_, numProcesses, uint_t(0), processMemoryLimit_, true, false );

   return setupBlockForest;
}

shared_ptr<SetupBlockForest> ComplexGeometryStructuredBlockforestCreator::createSetupBlockForest( const Vector3<uint_t> & cellsPerBlock, const Vector3<uint_t> & numBlocks ) const
{
   uint_t numProcesses = numBlocks[0] * numBlocks[1] * numBlocks[2];


   AABB newAABB(real_t(0), real_t(0), real_t(0), real_c(numBlocks[0] * cellsPerBlock[0]) * cellSize_[0],
                real_c(numBlocks[1] * cellsPerBlock[1]) * cellSize_[1], real_c(numBlocks[2] * cellsPerBlock[2]) * cellSize_[2]);

   newAABB.translate( aabb_.center() - newAABB.center() );
   compareAABB(aabb_, newAABB);

   auto setupBlockForest = make_shared<SetupBlockForest>();
   setupBlockForest->addRootBlockExclusionFunction( rootBlockExclusionFunction_ );
   setupBlockForest->addWorkloadMemorySUIDAssignmentFunction( workloadMemorySUIDAssignmentFunction_ );

   if( refinementSelectionFunction_ )
      setupBlockForest->addRefinementSelectionFunction( refinementSelectionFunction_ );

   setupBlockForest->init( newAABB, numBlocks[0], numBlocks[1], numBlocks[2], periodicity_[0], periodicity_[1], periodicity_[2] );

   setupBlockForest->balanceLoad( targetProcessAssignmentFunction_, numProcesses, uint_t(0), processMemoryLimit_, true, false );

   return setupBlockForest;
}

shared_ptr<StructuredBlockForest> ComplexGeometryStructuredBlockforestCreator::createStructuredBlockForest( const uint_t targetNumRootBlocks ) const
{
   shared_ptr< blockforest::SetupBlockForest> setupBlockForest = createSetupBlockForest( targetNumRootBlocks );

   Vector3< real_t > blockSize( setupBlockForest->getRootBlockXSize(),
                                setupBlockForest->getRootBlockYSize(),
                                setupBlockForest->getRootBlockZSize() );

   Vector3<uint_t> blockSizeCells( uint_c( blockSize[0] / cellSize_[0] + real_t(0.5) ), uint_c( blockSize[1] / cellSize_[1] + real_t(0.5) ), uint_c( blockSize[2] / cellSize_[2] + real_t(0.5) ) );

   WALBERLA_ASSERT_FLOAT_EQUAL( blockSize[0] / cellSize_[0], real_c( blockSizeCells[0] ) )
   WALBERLA_ASSERT_FLOAT_EQUAL( blockSize[1] / cellSize_[1], real_c( blockSizeCells[1] ) )
   WALBERLA_ASSERT_FLOAT_EQUAL( blockSize[2] / cellSize_[2], real_c( blockSizeCells[2] ) )

   auto blockForest = make_shared< blockforest::BlockForest >( MPIManager::instance()->rank(), *setupBlockForest );

   auto structuredBlockForest = make_shared< blockforest::StructuredBlockForest >( blockForest, blockSizeCells[0], blockSizeCells[1], blockSizeCells[2] );

   structuredBlockForest->createCellBoundingBoxes();

   return structuredBlockForest;
}

shared_ptr<StructuredBlockForest> ComplexGeometryStructuredBlockforestCreator::createStructuredBlockForest( const Vector3<uint_t> & blockSize ) const
{
   shared_ptr< blockforest::SetupBlockForest> setupBlockForest = createSetupBlockForest( blockSize );

   auto blockForest = make_shared< blockforest::BlockForest >( MPIManager::instance()->rank(), *setupBlockForest );

   auto structuredBlockForest = make_shared< blockforest::StructuredBlockForest >( blockForest, blockSize[0], blockSize[1], blockSize[2] );

   structuredBlockForest->createCellBoundingBoxes();

   return structuredBlockForest;
}

shared_ptr<StructuredBlockForest> ComplexGeometryStructuredBlockforestCreator::createStructuredBlockForest( const Vector3<uint_t> & cellsPerBlock, const Vector3<uint_t> & numBlocks ) const
{
   shared_ptr< blockforest::SetupBlockForest> setupBlockForest = createSetupBlockForest( cellsPerBlock, numBlocks );

   auto blockForest = make_shared< blockforest::BlockForest >( MPIManager::instance()->rank(), *setupBlockForest );

   auto structuredBlockForest = make_shared< blockforest::StructuredBlockForest >( blockForest, cellsPerBlock[0], cellsPerBlock[1], cellsPerBlock[2] );

   structuredBlockForest->createCellBoundingBoxes();

   return structuredBlockForest;
}

uint_t ComplexGeometryStructuredBlockforestCreator::findNumBlocks( const Vector3<uint_t> & blockSize ) const
{
   WALBERLA_LOG_DEVEL_ON_ROOT( "Testing block size " << blockSize )

   Vector3<uint_t> numCells( uint_c( std::ceil( aabb_.xSize() / cellSize_[0] ) ),
                             uint_c( std::ceil( aabb_.ySize() / cellSize_[1] ) ),
                             uint_c( std::ceil( aabb_.zSize() / cellSize_[2] ) ) );

   Vector3<uint_t> numBlocks( uint_c( std::ceil( real_c( numCells[0] ) / real_c( blockSize[0] ) ) ),
                              uint_c( std::ceil( real_c( numCells[1] ) / real_c( blockSize[1] ) ) ),
                              uint_c( std::ceil( real_c( numCells[2] ) / real_c( blockSize[2] ) ) ) );

   AABB newAABB( real_t(0), real_t(0), real_t(0),
                 real_c( numBlocks[0] * blockSize[0] ) * cellSize_[0],
                 real_c( numBlocks[1] * blockSize[1] ) * cellSize_[1],
                 real_c( numBlocks[2] * blockSize[2] ) * cellSize_[2] );

   newAABB.translate( aabb_.center() - newAABB.center() );


   SetupBlockForest setupBlockForest;
   setupBlockForest.addRootBlockExclusionFunction( rootBlockExclusionFunction_ );

   if( refinementSelectionFunction_ )
      setupBlockForest.addRefinementSelectionFunction( refinementSelectionFunction_ );

   setupBlockForest.init( newAABB, numBlocks[0], numBlocks[1], numBlocks[2], periodicity_[0], periodicity_[1], periodicity_[2] );

   WALBERLA_LOG_DEVEL_ON_ROOT( "Testing block size " << blockSize << " resulted in " << setupBlockForest.getNumberOfBlocks() )

   return setupBlockForest.getNumberOfBlocks();

}

} // namespace mesh
} // namespace walberla
