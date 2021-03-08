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
//! \file BodyAndVolumeFractionMapping.cpp
//! \ingroup pe_coupling
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//
//======================================================================================================================

#include "BodyAndVolumeFractionMapping.h"

#include "geometry/bodies/BodyOverlapFunctions.h"

#include "pe/rigidbody/BodyIterators.h"

#include "pe_coupling/geometry/PeOverlapFraction.h"
#include "pe_coupling/mapping/BodyBBMapping.h"

namespace walberla {
namespace pe_coupling {


void mapPSMBodyAndVolumeFraction( const pe::BodyID body, IBlock & block, StructuredBlockStorage & blockStorage,
                                  const BlockDataID bodyAndVolumeFractionFieldID )
{
   using BodyAndVolumeFraction_T = std::pair<pe::BodyID, real_t>;
   using BodyAndVolumeFractionField_T = GhostLayerField<std::vector<BodyAndVolumeFraction_T>, 1>;

   BodyAndVolumeFractionField_T * bodyAndVolumeFractionField = block.getData< BodyAndVolumeFractionField_T >( bodyAndVolumeFractionFieldID );
   WALBERLA_ASSERT_NOT_NULLPTR( bodyAndVolumeFractionField );

   // get bounding box of body
   CellInterval cellBB = getCellBB( body, block, blockStorage, bodyAndVolumeFractionField->nrOfGhostLayers() );

   uint_t level = blockStorage.getLevel( block );
   Vector3<real_t> dxVec(blockStorage.dx(level), blockStorage.dy(level), blockStorage.dz(level));

   for( auto cellIt = cellBB.begin(); cellIt != cellBB.end(); ++cellIt )
   {
      Cell cell( *cellIt );

      // get the cell's center
      Vector3<real_t> cellCenter;
      cellCenter = blockStorage.getBlockLocalCellCenter( block, cell );

      const real_t fraction = overlapFractionPe( *body, cellCenter, dxVec );

      // if the cell intersected with the body, store a pointer to that body and the corresponding volume fraction in the field
      if( fraction > real_t(0) )
      {
         bodyAndVolumeFractionField->get(cell).emplace_back( body, fraction );
      }
   }
}


void BodyAndVolumeFractionMapping::initialize()
{
   for( auto blockIt = blockStorage_->begin(); blockIt != blockStorage_->end(); ++blockIt )
   {
      BodyAndVolumeFractionField_T * bodyAndVolumeFractionField = blockIt->getData< BodyAndVolumeFractionField_T >( bodyAndVolumeFractionFieldID_ );

      if( updatedBodyAndVolumeFractionField_ == nullptr )
      {
         // hold internally an identical field for swapping
         updatedBodyAndVolumeFractionField_ = shared_ptr<BodyAndVolumeFractionField_T>( bodyAndVolumeFractionField->cloneUninitialized() );

         auto xyzFieldSize = updatedBodyAndVolumeFractionField_->xyzSize();
         for( auto fieldIt = xyzFieldSize.begin(); fieldIt != xyzFieldSize.end(); ++fieldIt )
         {
            (updatedBodyAndVolumeFractionField_->get(*fieldIt)).clear();
         }
      }

      // clear the field
      auto xyzFieldSize = bodyAndVolumeFractionField->xyzSize();
      for( auto fieldIt = xyzFieldSize.begin(); fieldIt != xyzFieldSize.end(); ++fieldIt )
      {
         (bodyAndVolumeFractionField->get(*fieldIt)).clear();
      }

      for( auto bodyIt = pe::BodyIterator::begin( *blockIt, bodyStorageID_); bodyIt != pe::BodyIterator::end(); ++bodyIt )
      {
         if( mappingBodySelectorFct_(bodyIt.getBodyID()) )
         {
            mapPSMBodyAndVolumeFraction(bodyIt.getBodyID(), *blockIt, *blockStorage_, bodyAndVolumeFractionFieldID_ );
            lastUpdatedPositionMap_.insert( std::pair< walberla::id_t, Vector3< real_t > >( bodyIt->getSystemID(), bodyIt->getPosition() ) );
         }
      }

      for( auto bodyIt = globalBodyStorage_->begin(); bodyIt != globalBodyStorage_->end(); ++bodyIt )
      {
         if( mappingBodySelectorFct_(bodyIt.getBodyID()) )
         {
            mapPSMBodyAndVolumeFraction(bodyIt.getBodyID(), *blockIt, *blockStorage_, bodyAndVolumeFractionFieldID_);
            lastUpdatedPositionMap_.insert( std::pair< walberla::id_t, Vector3< real_t > >( bodyIt->getSystemID(), bodyIt->getPosition() ) );
         }
      }

   }
}

void BodyAndVolumeFractionMapping::update()
{
   std::map< walberla::id_t, Vector3< real_t > > tempLastUpdatedPositionMap;

   for( auto blockIt = blockStorage_->begin(); blockIt != blockStorage_->end(); ++blockIt )
   {
      BodyAndVolumeFractionField_T * bodyAndVolumeFractionField = blockIt->getData< BodyAndVolumeFractionField_T >( bodyAndVolumeFractionFieldID_ );

      for( auto bodyIt = pe::BodyIterator::begin( *blockIt, bodyStorageID_); bodyIt != pe::BodyIterator::end(); ++bodyIt )
      {
         if( mappingBodySelectorFct_(bodyIt.getBodyID()) )
         {
            updatePSMBodyAndVolumeFraction(bodyIt.getBodyID(), *blockIt, bodyAndVolumeFractionField, tempLastUpdatedPositionMap);
         }
      }

      for( auto bodyIt = globalBodyStorage_->begin(); bodyIt != globalBodyStorage_->end(); ++bodyIt )
      {
         if( mappingBodySelectorFct_(bodyIt.getBodyID()) )
         {
            updatePSMBodyAndVolumeFraction(bodyIt.getBodyID(), *blockIt, bodyAndVolumeFractionField, tempLastUpdatedPositionMap);
         }
      }

      bodyAndVolumeFractionField->swapDataPointers( *updatedBodyAndVolumeFractionField_ );

      auto xyzFieldSize = updatedBodyAndVolumeFractionField_->xyzSize();
      for( auto fieldIt = xyzFieldSize.begin(); fieldIt != xyzFieldSize.end(); ++fieldIt )
      {
         (updatedBodyAndVolumeFractionField_->get(*fieldIt)).clear();
      }
   }

   lastUpdatedPositionMap_ = tempLastUpdatedPositionMap;
}

void BodyAndVolumeFractionMapping::updatePSMBodyAndVolumeFraction( pe::BodyID body, IBlock & block,
                                                                   BodyAndVolumeFractionField_T * oldBodyAndVolumeFractionField,
                                                                   std::map< walberla::id_t, Vector3< real_t > > & tempLastUpdatedPositionMap )
{

   WALBERLA_ASSERT_NOT_NULLPTR( oldBodyAndVolumeFractionField );

   // estimate traveled distance since last volume fraction update
   real_t traveledSquaredDistance( real_t(0) );
   auto mapBodyIt = lastUpdatedPositionMap_.find( body->getSystemID() );
   if( mapBodyIt != lastUpdatedPositionMap_.end() )
   {
      // body found and traveled distance can be estimated
      Vector3<real_t> distanceVec = body->getPosition() - mapBodyIt->second;
      traveledSquaredDistance = distanceVec.sqrLength();
   } else
   {
      // body was not found in map -> is a new body, so no information is known
      traveledSquaredDistance = std::numeric_limits<real_t>::max();
   }

   // get bounding box of body
   CellInterval cellBB = getCellBB( body, block, *blockStorage_, oldBodyAndVolumeFractionField->nrOfGhostLayers() );

   uint_t level = blockStorage_->getLevel( block );
   Vector3<real_t> dxVec(blockStorage_->dx(level), blockStorage_->dy(level), blockStorage_->dz(level));

   // if body has not moved (specified by some epsilon), just reuse old fraction values
   if( body->getLinearVel().sqrLength()  < velocityUpdatingEpsilonSquared_ &&
       body->getAngularVel().sqrLength() < velocityUpdatingEpsilonSquared_ &&
       traveledSquaredDistance < positionUpdatingEpsilonSquared_ )
   {      
      for( cell_idx_t z = cellBB.zMin(); z <= cellBB.zMax(); ++z )
      {
         for( cell_idx_t y = cellBB.yMin(); y <= cellBB.yMax(); ++y )
         {
            for( cell_idx_t x = cellBB.xMin(); x <= cellBB.xMax(); ++x )
            {

               auto oldVec = oldBodyAndVolumeFractionField->get(x,y,z);
               for( auto pairIt = oldVec.begin(); pairIt != oldVec.end(); ++pairIt ) 
               {
                  if( pairIt->first == body )
                  {
                     updatedBodyAndVolumeFractionField_->get(x,y,z).push_back(*pairIt);
                     break;
                  }
               }
            }
         }
      }
      tempLastUpdatedPositionMap.insert( std::pair< walberla::id_t, Vector3< real_t > >( body->getSystemID(), mapBodyIt->second ) );

   } else
   {
      // else body has moved significantly or is new to this block, thus the volume fraction has to be calculated
      for( cell_idx_t z = cellBB.zMin(); z <= cellBB.zMax(); ++z )
      {
         for( cell_idx_t y = cellBB.yMin(); y <= cellBB.yMax(); ++y )
         {
            for( cell_idx_t x = cellBB.xMin(); x <= cellBB.xMax(); ++x )
            {
               // get the cell's center
               Vector3<real_t> cellCenter;
               cellCenter = blockStorage_->getBlockLocalCellCenter( block, Cell(x,y,z) );

               const real_t fraction = overlapFractionPe( *body, cellCenter, dxVec, superSamplingDepth_ );

               // if the cell intersected with the body, store a pointer to that body and the corresponding volume fraction in the field
               if( fraction > real_t(0) )
               {
                  updatedBodyAndVolumeFractionField_->get(x,y,z).emplace_back( body, fraction );
               }
            }
         }
      }

      tempLastUpdatedPositionMap.insert( std::pair< walberla::id_t, Vector3< real_t > >( body->getSystemID(), body->getPosition() ) );
   }
}


} // namespace pe_coupling
} // namespace walberla

