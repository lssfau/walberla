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
//! \file BlockWorkloadMemory.h
//! \ingroup mesh
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//
//======================================================================================================================

#pragma once

#include "blockforest/SetupBlockForest.h"
#include "blockforest/Types.h"

#include "core/DataTypes.h"
#include "core/cell/CellInterval.h"
#include "core/mpi/MPIManager.h"
#include "core/mpi/Reduce.h"
#include "core/uid/SUID.h"

#include <random>
#include <vector>

namespace walberla {
namespace mesh {

template< typename DistanceObject >
class MeshWorkloadMemory
{
public:
   typedef blockforest::memory_t           memory_t;
   typedef blockforest::workload_t         workload_t;
   typedef typename DistanceObject::Scalar Scalar;

   MeshWorkloadMemory( const shared_ptr< DistanceObject > & distanceObject, const Vector3< Scalar > & cellSize ) : distanceObject_( distanceObject ), cellSize_( cellSize ) { defaultInit(); }
   MeshWorkloadMemory( const shared_ptr< DistanceObject > & distanceObject, const Scalar cellSize ) : distanceObject_( distanceObject ), cellSize_( cellSize, cellSize, cellSize ) { defaultInit(); }

   void setInsideCellWorkload ( const workload_t workload ) { insideCellWorkload_ = workload;  }
   void setOutsideCellWorkload( const workload_t workload ) { outsideCellWorkload_ = workload; }

   void setInsideCellMemoryConsumption ( const memory_t memoryConsumption ) { insideCellMemoryConsumption_ = memoryConsumption;  }
   void setOutsideCellMemoryConsumption( const memory_t memoryConsumption ) { outsideCellMemoryConsumption_ = memoryConsumption; }

   void setForceZeroMemoryOnZeroWorkload( const bool value ) { forceZeroMemoryOnZeroWorkload_ = value; }
   void setZeroWorkloadSUID( const SUID & suid ) { zeroWorkloadSUID_ = make_shared<SUID>( suid ); }

   void operator()( blockforest::SetupBlockForest & forest ) const;

private:

   void defaultInit();
   void countCells( const math::GenericAABB< Scalar > & aabb, const uint_t level, uint_t & numCellsInside, uint_t & numCellsOutside ) const;

   shared_ptr< DistanceObject > distanceObject_;
   Vector3< typename DistanceObject::Scalar > cellSize_;

   memory_t outsideCellMemoryConsumption_;
   memory_t insideCellMemoryConsumption_;

   workload_t outsideCellWorkload_;
   workload_t insideCellWorkload_;

   shared_ptr< SUID > zeroWorkloadSUID_;
   bool               forceZeroMemoryOnZeroWorkload_;
};

template< typename DistanceObject >
MeshWorkloadMemory< DistanceObject > makeMeshWorkloadMemory( const shared_ptr< DistanceObject > & distanceObject, const Vector3< typename DistanceObject::Scalar > & cellSize ) { return MeshWorkloadMemory< DistanceObject >( distanceObject, cellSize ); }


template< typename DistanceObject >
MeshWorkloadMemory< DistanceObject > makeMeshWorkloadMemory( const shared_ptr< DistanceObject > & distanceObject, const typename DistanceObject::Scalar cellSize ) { return MeshWorkloadMemory< DistanceObject >( distanceObject, cellSize ); }

template< typename DistanceObject >
inline void walberla::mesh::MeshWorkloadMemory<DistanceObject>::defaultInit()
{
   setInsideCellWorkload ( blockforest::memory_c(1) );
   setOutsideCellWorkload( blockforest::memory_c(0) );

   setInsideCellMemoryConsumption ( blockforest::memory_c(1) );
   setOutsideCellMemoryConsumption( blockforest::memory_c(1) );
}


template< typename DistanceObject >
void walberla::mesh::MeshWorkloadMemory<DistanceObject>::countCells( const math::GenericAABB< Scalar > & aabb, const uint_t level, uint_t & numCellsInside, uint_t & numCellsOutside ) const
{
   typedef math::GenericAABB< Scalar > Box;

   numCellsInside = uint_t(0);
   numCellsOutside = uint_t(0);

    if( aabb.empty() )
      return;

   std::queue< Box > boxQueue;
   boxQueue.push( aabb );

   Vector3<Scalar> levelCellSize = cellSize_ / numeric_cast<Scalar>( uint_t(1) << level );

   while( !boxQueue.empty() )
   {
      const Box & curAabb = boxQueue.front();

      if(curAabb.empty())
      {
         boxQueue.pop();
         continue;
      }

      const CellInterval ci( cell_idx_c( ( curAabb.xMin() - aabb.xMin() ) / levelCellSize[0] + Scalar(0.5) ),
                             cell_idx_c( ( curAabb.yMin() - aabb.yMin() ) / levelCellSize[1] + Scalar(0.5) ),
                             cell_idx_c( ( curAabb.zMin() - aabb.zMin() ) / levelCellSize[2] + Scalar(0.5) ),
                             cell_idx_c( ( curAabb.xMax() - aabb.xMin() ) / levelCellSize[0] - Scalar(0.5) ),
                             cell_idx_c( ( curAabb.yMax() - aabb.yMin() ) / levelCellSize[1] - Scalar(0.5) ),
                             cell_idx_c( ( curAabb.zMax() - aabb.zMin() ) / levelCellSize[2] - Scalar(0.5) ) );

      WALBERLA_ASSERT( !ci.empty(), ci );

      const Scalar sqSignedDistance = distanceObject_->sqSignedDistance( toOpenMesh( curAabb.center() ) );

      if( ci.numCells() == uint_t(1) )
      {
         boxQueue.pop();
         sqSignedDistance < Scalar(0) ? ++numCellsInside : ++numCellsOutside;
         continue;
      }

      const Scalar circumRadius = curAabb.sizes().length() * Scalar(0.5);
      const Scalar sqCircumRadius = circumRadius * circumRadius;

      if( sqSignedDistance < -sqCircumRadius )
      {
         boxQueue.pop();
         numCellsInside += ci.numCells();
         continue; // clearly the box is fully covered by the mesh
      }

      if( sqSignedDistance > sqCircumRadius )
      {
         boxQueue.pop();
         numCellsOutside += ci.numCells();
         continue; // clearly the box is fully outside of the mesh
      }

      const auto &    min = curAabb.minCorner();
      const auto &    max = curAabb.maxCorner();

      const Vector3< Scalar > center( curAabb.xMin() + numeric_cast<Scalar>( ci.xSize() / uint_t(2) ) * levelCellSize[0],
                                      curAabb.yMin() + numeric_cast<Scalar>( ci.ySize() / uint_t(2) ) * levelCellSize[1],
                                      curAabb.zMin() + numeric_cast<Scalar>( ci.zSize() / uint_t(2) ) * levelCellSize[2] );

      boxQueue.push( Box::createFromMinMaxCorner(    min[0],    min[1],    min[2], center[0], center[1], center[2] ) );
      boxQueue.push( Box::createFromMinMaxCorner(    min[0],    min[1], center[2], center[0], center[1],    max[2] ) );
      boxQueue.push( Box::createFromMinMaxCorner(    min[0], center[1],    min[2], center[0],    max[1], center[2] ) );
      boxQueue.push( Box::createFromMinMaxCorner(    min[0], center[1], center[2], center[0],    max[1],    max[2] ) );
      boxQueue.push( Box::createFromMinMaxCorner( center[0],    min[1],    min[2],    max[0], center[1], center[2] ) );
      boxQueue.push( Box::createFromMinMaxCorner( center[0],    min[1], center[2],    max[0], center[1],    max[2] ) );
      boxQueue.push( Box::createFromMinMaxCorner( center[0], center[1],    min[2],    max[0],    max[1], center[2] ) );
      boxQueue.push( Box::createFromMinMaxCorner( center[0], center[1], center[2],    max[0],    max[1],    max[2] ) );

      boxQueue.pop();
   }
}


template< typename DistanceObject >
inline void walberla::mesh::MeshWorkloadMemory<DistanceObject>::operator()( blockforest::SetupBlockForest & forest ) const
{
   std::vector< blockforest::SetupBlock* > blocks;
   forest.getBlocks( blocks );

   const uint_t numBlocks    = uint_c( blocks.size() );
   const uint_t numProcesses = uint_c( MPIManager::instance()->numProcesses() );
   const uint_t chunkSize    = uint_c( std::ceil( real_c( numBlocks ) / real_c( numProcesses ) ) );

   const uint_t rank    = uint_c( MPIManager::instance()->rank() );
   const int chunkBegin = int_c( rank * chunkSize );
   const int chunkEnd   = std::min( int_c( ( rank + 1 ) * chunkSize ), int_c( numBlocks ) );

   std::vector<size_t> shuffle( numBlocks );
   for( size_t i = 0; i < shuffle.size(); ++i )
   {
      shuffle[i] = i;
   }

   //old compilers might need this
   //std::srand( 42 );
   //std::random_shuffle( shuffle.begin(), shuffle.end() );

   std::mt19937 g( 42 );
   std::shuffle( shuffle.begin(), shuffle.end(), g );

   std::vector<uint_t> insideCells( numBlocks, uint_t( 0 ) );
   std::vector<uint_t> outsideCells( numBlocks, uint_t( 0 ) );

   #ifdef _OPENMP
   #pragma omp parallel for schedule( dynamic )
   #endif
   for( int i = chunkBegin; i < chunkEnd; ++i )
   {
      const size_t is = numeric_cast<size_t>( i );
      const math::AABB & aabb = blocks[ shuffle[ is ] ]->getAABB();
      const uint_t level = blocks[ shuffle[ is ] ]->getLevel();
      countCells( aabb, level, insideCells[ shuffle[ is ] ], outsideCells[ shuffle[ is ] ] );
   }

   allReduceInplace( insideCells, mpi::MAX );
   allReduceInplace( outsideCells, mpi::MAX );

   for( uint_t i = 0; i != blocks.size(); ++i )
   {
      blockforest::workload_t workload = blockforest::workload_c( insideCells[i] ) * insideCellWorkload_ + blockforest::workload_c( outsideCells[i] ) * outsideCellWorkload_;
      blocks[i]->setWorkload( workload );

      blockforest::memory_t memory = blockforest::memory_t( 0 );
      
      if( !forceZeroMemoryOnZeroWorkload_ || !floatIsEqual( workload, blockforest::workload_t(0) ) )
         memory = blockforest::memory_c( insideCells[i] ) * insideCellMemoryConsumption_ + blockforest::memory_c( outsideCells[i] ) * outsideCellMemoryConsumption_;

      blocks[i]->setMemory( memory );

      if( zeroWorkloadSUID_ && floatIsEqual( workload, blockforest::workload_t(0) ) )
         blocks[i]->setState( *zeroWorkloadSUID_ );
   }
}

} // namespace mesh
} // namespace walberla

