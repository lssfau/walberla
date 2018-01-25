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
//! \file BoundarySetup.h
//! \ingroup mesh
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//
//======================================================================================================================

#pragma once

#include "BoundaryInfo.h"

#include "core/DataTypes.h"
#include "core/cell/CellInterval.h"

#include "domain_decomposition/StructuredBlockStorage.h"
#include "domain_decomposition/BlockDataID.h"

#include "field/GhostLayerField.h"
#include "field/FlagField.h"

#include "stencil/D3Q27.h"

#include <functional>

#include <queue>

namespace walberla {
namespace mesh {

class BoundarySetup
{
public:
   typedef field::GhostLayerField< uint8_t, uint_t(1) > VoxelizationField;

   typedef std::function< real_t ( const Vector3< real_t > &) > DistanceFunction;

   enum Location { INSIDE, OUTSIDE };

   BoundarySetup( const shared_ptr< StructuredBlockStorage > & structuredBlockStorage, const DistanceFunction & distanceFunction, const uint_t numGhostLayers )
      : structuredBlockStorage_( structuredBlockStorage ), distanceFunction_( distanceFunction ), numGhostLayers_( numGhostLayers ), cellVectorChunkSize_( size_t(1000) )
   {
      voxelize();
   }

   ~BoundarySetup() { deallocateVoxelizationField(); }

   template< typename BoundaryHandlingType >
   void setDomainCells( const BlockDataID boundaryHandlingID, const Location domainLocation );

   template< typename BoundaryHandlingType, typename BoundaryFunction, typename Stencil = stencil::D3Q27 >
   void setBoundaries( const BlockDataID boundaryHandlingID, const BoundaryFunction & boundaryFunction, Location boundaryLocation );


   template<typename FlagField_T>
   void setFlag( const BlockDataID flagFieldID, field::FlagUID flagUID, Location boundaryLocation );

   void writeVTKVoxelfile( const std::string & identifier = "voxelization", bool writeGhostLayers = false, const std::string & baseFolder = std::string("vtk_out"),
                           const std::string & executionFolder = std::string("voxelization") );

private:

   void divideAndPushCellInterval( const CellInterval & ci, std::queue< CellInterval > & outputQueue );

   void allocateOrResetVoxelizationField();
   void deallocateVoxelizationField();

   void voxelize();

   shared_ptr< StructuredBlockStorage >       structuredBlockStorage_;
   shared_ptr< BlockDataID >                  voxelizationFieldId_;
   DistanceFunction                           distanceFunction_;  // function providing the squared signed distance to an object
   uint_t                                     numGhostLayers_;
   size_t                                     cellVectorChunkSize_;
};



template< typename BoundaryHandlingType >
void BoundarySetup::setDomainCells( const BlockDataID boundaryHandlingId, const Location domainLocation )
{
   for( auto blockIt = structuredBlockStorage_->begin(); blockIt != structuredBlockStorage_->end(); ++blockIt )
   {
      BoundaryHandlingType * boundaryHandling  = blockIt->getData< BoundaryHandlingType >( boundaryHandlingId  );
      WALBERLA_CHECK_NOT_NULLPTR( boundaryHandling, "boundaryHandlingId invalid!" );

      const VoxelizationField * voxelizationField = blockIt->getData< VoxelizationField >( *voxelizationFieldId_ );
      WALBERLA_ASSERT_NOT_NULLPTR( voxelizationField );
      WALBERLA_CHECK_LESS_EQUAL( numGhostLayers_, boundaryHandling->getFlagField()->nrOfGhostLayers(), "You want to use mesh boundary setup with " \
                                 << numGhostLayers_ << " but your flag field has only " << boundaryHandling->getFlagField()->nrOfGhostLayers() << " ghost layers!" );

      const uint8_t domainValue = domainLocation == INSIDE ? uint8_t(1) : uint8_t(0);

      std::vector<Cell> domainCells;

      for( cell_idx_t z = -cell_idx_c( numGhostLayers_ ); z < cell_idx_c( voxelizationField->zSize() + numGhostLayers_ ); ++z )
         for( cell_idx_t y = -cell_idx_c( numGhostLayers_ ); y < cell_idx_c( voxelizationField->ySize() + numGhostLayers_ ); ++y )
            for( cell_idx_t x = -cell_idx_c( numGhostLayers_ ); x < cell_idx_c( voxelizationField->xSize() + numGhostLayers_ ); ++x )
            {
               if( voxelizationField->get( x, y, z ) == domainValue )
                  domainCells.push_back( Cell(x,y,z) );

               if( domainCells.size() > cellVectorChunkSize_ )
               {
                  boundaryHandling->setDomain( domainCells.begin(), domainCells.end() );
                  domainCells.clear();
               }
            }

      boundaryHandling->setDomain( domainCells.begin(), domainCells.end() );
   }
}


template<typename FlagField_T>
void BoundarySetup::setFlag( const BlockDataID flagFieldID, field::FlagUID flagUID, Location boundaryLocation )
{
   for( auto blockIt = structuredBlockStorage_->begin(); blockIt != structuredBlockStorage_->end(); ++blockIt )
   {
      FlagField_T * flagField  = blockIt->getData< FlagField_T >( flagFieldID );
      WALBERLA_CHECK_NOT_NULLPTR( flagField, "flagFieldID invalid!" );
      auto flag = flagField->getFlag(flagUID);

      const VoxelizationField * voxelizationField = blockIt->getData< VoxelizationField >( *voxelizationFieldId_ );
      WALBERLA_ASSERT_NOT_NULLPTR( voxelizationField );
      WALBERLA_CHECK_LESS_EQUAL( numGhostLayers_, flagField->nrOfGhostLayers(), "You want to use mesh boundary setup with " \
                                 << numGhostLayers_ << " but your flag field has only " << flagField->nrOfGhostLayers() << " ghost layers!" );

      const uint8_t domainValue = boundaryLocation == INSIDE ? uint8_t(0) : uint8_t(1);

      for( cell_idx_t z = -cell_idx_c( numGhostLayers_ ); z < cell_idx_c( voxelizationField->zSize() + numGhostLayers_ ); ++z )
         for( cell_idx_t y = -cell_idx_c( numGhostLayers_ ); y < cell_idx_c( voxelizationField->ySize() + numGhostLayers_ ); ++y )
            for( cell_idx_t x = -cell_idx_c( numGhostLayers_ ); x < cell_idx_c( voxelizationField->xSize() + numGhostLayers_ ); ++x )
            {
               if( voxelizationField->get( x, y, z ) != domainValue )
                  flagField->addFlag(x, y, z, flag);
            }
   }
}


template< typename BoundaryHandlingType, typename BoundaryFunction, typename Stencil >
void BoundarySetup::setBoundaries( const BlockDataID boundaryHandlingID, const BoundaryFunction & boundaryFunction, Location boundaryLocation )
{
   for( auto blockIt = structuredBlockStorage_->begin(); blockIt != structuredBlockStorage_->end(); ++blockIt )
   {
      BoundaryHandlingType * boundaryHandling  = blockIt->getData< BoundaryHandlingType >( boundaryHandlingID  );
      WALBERLA_CHECK_NOT_NULLPTR( boundaryHandling, "boundaryHandlingId invalid!" );

      const VoxelizationField * voxelizationField = blockIt->getData< VoxelizationField    >( *voxelizationFieldId_ );
      WALBERLA_ASSERT_NOT_NULLPTR( voxelizationField );
      WALBERLA_CHECK_LESS_EQUAL( numGhostLayers_, boundaryHandling->getFlagField()->nrOfGhostLayers(), "You want to use mesh boundary setup with " \
                                 << numGhostLayers_ << " but your flag field has only " << boundaryHandling->getFlagField()->nrOfGhostLayers() << " ghost layers!" );

      const uint8_t domainValue   = boundaryLocation == INSIDE ? uint8_t(0) : uint8_t(1);

      const CellInterval blockCi = voxelizationField->xyzSizeWithGhostLayer();

      for( cell_idx_t z = -cell_idx_c( numGhostLayers_ ); z < cell_idx_c( voxelizationField->zSize() + numGhostLayers_ ); ++z )
         for( cell_idx_t y = -cell_idx_c( numGhostLayers_ ); y < cell_idx_c( voxelizationField->ySize() + numGhostLayers_ ); ++y )
            for( cell_idx_t x = -cell_idx_c( numGhostLayers_ ); x < cell_idx_c( voxelizationField->xSize() + numGhostLayers_ ); ++x )
            {
               if( voxelizationField->get( x, y, z ) == domainValue )
                  continue;

               for(auto dirIt = Stencil::beginNoCenter(); dirIt != Stencil::end(); ++dirIt)
               {
                  const Cell cell( x, y, z );
                  const Cell neighborCell = cell + *dirIt;
                  if( blockCi.contains( neighborCell ) && voxelizationField->get( neighborCell ) == domainValue )
                  {
                     const Vector3< real_t > cellCenter = structuredBlockStorage_->getBlockLocalCellCenter( *blockIt, cell );
                     const BoundaryInfo & bi = boundaryFunction( cellCenter );
                     const auto boundaryMask = boundaryHandling->getBoundaryMask( bi.getUid() );

                     boundaryHandling->setBoundary( boundaryMask, cell.x(), cell.y(), cell.z(), *( bi.getConfig() ) );

                     break;
                  }
               }
            }
   }
}



} // namespace mesh
} // namespace walberla