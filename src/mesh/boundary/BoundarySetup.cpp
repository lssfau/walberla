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
//! \file BoundarySetup.cpp
//! \ingroup mesh
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//
//======================================================================================================================

#include "BoundarySetup.h"

#include "field/AddToStorage.h"
#include "field/vtk/VTKWriter.h"

namespace walberla {
namespace mesh {

void BoundarySetup::divideAndPushCellInterval( const CellInterval & ci, std::queue< CellInterval > & outputQueue )
{
   WALBERLA_ASSERT( !ci.empty() );

   Cell newMax( ci.xMin() + std::max( cell_idx_c( ci.xSize() ) / cell_idx_t(2) - cell_idx_t(1), cell_idx_t(0) ),
                ci.yMin() + std::max( cell_idx_c( ci.ySize() ) / cell_idx_t(2) - cell_idx_t(1), cell_idx_t(0) ),
                ci.zMin() + std::max( cell_idx_c( ci.zSize() ) / cell_idx_t(2) - cell_idx_t(1), cell_idx_t(0) ) );

   WALBERLA_ASSERT( ci.contains( newMax ) );

   Cell newMin( newMax[0] + cell_idx_c( 1 ), newMax[1] + cell_idx_c( 1 ), newMax[2] + cell_idx_c( 1 ) );

   outputQueue.push( CellInterval( ci.xMin(), ci.yMin(), ci.zMin(), newMax[0], newMax[1], newMax[2] ) );
   if( newMin[2] <= ci.zMax() )
      outputQueue.push( CellInterval( ci.xMin(), ci.yMin(), newMin[2], newMax[0], newMax[1], ci.zMax() ) );
   if( newMin[1] <= ci.yMax() )
   {
      outputQueue.push( CellInterval( ci.xMin(), newMin[1], ci.zMin(), newMax[0], ci.yMax(), newMax[2]) );
      if( newMin[2] <= ci.zMax() )
         outputQueue.push( CellInterval( ci.xMin(), newMin[1], newMin[2], newMax[0], ci.yMax(), ci.zMax() ) );
   }
   if( newMin[0] <= ci.xMax() )
   {
      outputQueue.push( CellInterval( newMin[0], ci.yMin(), ci.zMin(), ci.xMax(), newMax[1], newMax[2] ) );
      if( newMin[2] <= ci.zMax() )
         outputQueue.push( CellInterval( newMin[0], ci.yMin(), newMin[2], ci.xMax(), newMax[1], ci.zMax() ) );
      if( newMin[1] <= ci.yMax() )
      {
         outputQueue.push( CellInterval( newMin[0], newMin[1], ci.zMin(), ci.xMax(), ci.yMax(), newMax[2] ) );
         if( newMin[2] <= ci.zMax() )
            outputQueue.push( CellInterval( newMin[0], newMin[1], newMin[2], ci.xMax(), ci.yMax(), ci.zMax() ) );
      }
   }
}

void BoundarySetup::allocateOrResetVoxelizationField()
{
   if( voxelizationFieldId_ )
   {
      for( auto blockIt = structuredBlockStorage_->begin(); blockIt != structuredBlockStorage_->end(); ++blockIt )
      {
         VoxelizationField * voxelizationField = blockIt->getData< VoxelizationField >( *voxelizationFieldId_ );
         voxelizationField->setWithGhostLayer( uint8_t(0) );
      }
   }
   else
   {
      voxelizationFieldId_ = make_shared< BlockDataID >( field::addToStorage< VoxelizationField >( structuredBlockStorage_, "flag field", uint8_t(0), field::zyxf, numGhostLayers_ ) );
   }

   WALBERLA_ASSERT_NOT_NULLPTR( voxelizationFieldId_ );
}

void BoundarySetup::deallocateVoxelizationField()
{
   if( voxelizationFieldId_ )
   {
      structuredBlockStorage_->clearBlockData( *voxelizationFieldId_ );
      voxelizationFieldId_.reset();
   }
}

void BoundarySetup::voxelize()
{
   allocateOrResetVoxelizationField();

   for( auto blockIt = structuredBlockStorage_->begin(); blockIt != structuredBlockStorage_->end(); ++blockIt )
   {
      VoxelizationField * voxelizationField = blockIt->getData< VoxelizationField >( *voxelizationFieldId_ );

      WALBERLA_ASSERT_NOT_NULLPTR( voxelizationField );
      WALBERLA_ASSERT_EQUAL( numGhostLayers_, voxelizationField->nrOfGhostLayers() );

      CellInterval blockCi = voxelizationField->xyzSizeWithGhostLayer();
      structuredBlockStorage_->transformBlockLocalToGlobalCellInterval( blockCi, *blockIt );

      std::queue< CellInterval > ciQueue;
      ciQueue.push( blockCi );

      while( !ciQueue.empty() )
      {
         const CellInterval & curCi = ciQueue.front();

         WALBERLA_ASSERT( !curCi.empty(), "Cell Interval: " << curCi );

         AABB curAABB = structuredBlockStorage_->getAABBFromCellBB( curCi, structuredBlockStorage_->getLevel( *blockIt ) );

         WALBERLA_ASSERT( !curAABB.empty(), "AABB: " << curAABB );

         Vector3<real_t> cellCenter = curAABB.center();
         structuredBlockStorage_->mapToPeriodicDomain( cellCenter );
         const real_t sqSignedDistance = distanceFunction_( cellCenter );

         if( curCi.numCells() == uint_t(1) )
         {
            if( ( sqSignedDistance < real_t(0) ) )
            {
               Cell localCell;
               structuredBlockStorage_->transformGlobalToBlockLocalCell( localCell, *blockIt, curCi.min() );
               voxelizationField->get( localCell ) = uint8_t(1);
            }

            ciQueue.pop();
            continue;
         }

         const real_t circumRadius = curAABB.sizes().length() * real_t(0.5);
         const real_t sqCircumRadius = circumRadius * circumRadius;

         if( sqSignedDistance < -sqCircumRadius )
         {
            // clearly the cell interval is fully covered by the mesh
            CellInterval localCi;
            structuredBlockStorage_->transformGlobalToBlockLocalCellInterval( localCi, *blockIt, curCi );
            std::fill( voxelizationField->beginSliceXYZ( localCi ), voxelizationField->end(), uint8_t(1) );

            ciQueue.pop();
            continue;
         }

         if( sqSignedDistance > sqCircumRadius )
         {
            // clearly the cell interval is fully outside of the mesh

            ciQueue.pop();
            continue;
         }

         WALBERLA_ASSERT_GREATER( curCi.numCells(), uint_t(1) );
         divideAndPushCellInterval( curCi, ciQueue );

         ciQueue.pop();
      }
   }
}


void BoundarySetup::writeVTKVoxelfile( const std::string & identifier, bool writeGhostLayers, const std::string & baseFolder, const std::string & executionFolder )
{
   WALBERLA_ASSERT_NOT_NULLPTR( voxelizationFieldId_ );
   WALBERLA_ASSERT_NOT_NULLPTR( structuredBlockStorage_ );

   field::createVTKOutput< VoxelizationField >( *voxelizationFieldId_, *structuredBlockStorage_, identifier, uint_t(1), writeGhostLayers ? numGhostLayers_ : uint_t(0), false, baseFolder, executionFolder )();
}


} // namespace mesh
} // namespace walberla