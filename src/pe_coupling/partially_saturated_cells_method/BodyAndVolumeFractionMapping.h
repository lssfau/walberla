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
//! \file BodyAndVolumeFractionMapping.h
//! \ingroup pe_coupling
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//
//======================================================================================================================

#pragma once

#include "domain_decomposition/StructuredBlockStorage.h"

#include "field/GhostLayerField.h"

#include "geometry/bodies/BodyOverlapFunctions.h"

#include "pe/rigidbody/BodyIterators.h"
#include "pe/Types.h"

#include "pe_coupling/geometry/PeOverlapFraction.h"
#include "pe_coupling/mapping/BodyBBMapping.h"
#include "pe_coupling/utility/BodySelectorFunctions.h"

#include <functional>

namespace walberla {
namespace pe_coupling {

/*!\brief Computes the solid volume fraction for each cell intersecting with the specific pe body
 *
 * Depending on the implementation of the "overlapFraction" for the respective pe body type (sphere, ellipsoid, etc. ),
 * super sampling or an analytical formula is used.
 *
 * The information is stored as a pair of a pointer to a PE body and the corresponding solid volume fraction.
 * As several bodies could intersect with one cell, the pairs are stored in a vector with the size of the amount of intersecting bodies.
 *
 */
void mapPSMBodyAndVolumeFraction( const pe::BodyID body, IBlock & block, StructuredBlockStorage & blockStorage,
                                  const BlockDataID bodyAndVolumeFractionFieldID )
{
   typedef std::pair< pe::BodyID, real_t >                              BodyAndVolumeFraction_T;
   typedef GhostLayerField< std::vector< BodyAndVolumeFraction_T >, 1 > BodyAndVolumeFractionField_T;

   BodyAndVolumeFractionField_T * bodyAndVolumeFractionField = block.getData< BodyAndVolumeFractionField_T >( bodyAndVolumeFractionFieldID );
   WALBERLA_ASSERT_NOT_NULLPTR( bodyAndVolumeFractionField );

   // get bounding box of body
   CellInterval cellBB = getCellBB( body, block, blockStorage, bodyAndVolumeFractionField->nrOfGhostLayers() );

   for( auto cellIt = cellBB.begin(); cellIt != cellBB.end(); ++cellIt )
   {
      Cell cell( *cellIt );

      // get the cell's center
      Vector3<real_t> cellCenter;
      cellCenter = blockStorage.getBlockLocalCellCenter( block, cell );

      const real_t fraction = overlapFractionPe( *body, cellCenter, blockStorage.dx( blockStorage.getLevel( block ) ) );

      // if the cell intersected with the body, store a pointer to that body and the corresponding volume fraction in the field
      if( fraction > real_t(0) )
      {
         bodyAndVolumeFractionField->get(cell).push_back( BodyAndVolumeFraction_T( body, fraction ) );
      }
   }
}

/*!\brief Mapping class that can be used inside the timeloop to update the volume fractions and body mappings
 *
 * Upon construction, it uses the free mapping function to initialize the fraction and body mapping field.
 *
 * Whether or not a body is considered by this mapping depends on the return value of 'mappingBodySelectorFct'.
 *
 * Successive calls try to update this field instead of completely recalculating all volume fractions which is expensive.
 * To do so, the update functionality determines whether bodies have very small velocities, both translational and rotational,
 * (limited with velocityUpdatingEpsilon) or have not traveled very far (limited by positionUpdatingEpsilon) since the last volume fraction recalculation.
 * If both criteria are fulfilled, the volume fraction information is simply copied from the old data field.
 * If not, the recalculation and remapping is carried out.
 * This functionality is a pure performance optimization and might affect the accuracy of the simulation if the limits are too large.
 *
 * Similarly, the argument superSamplingDepth can be used to limit the default depth (=4) of the super sampling approach
 * which is applied to approximate the volume fraction in each cell when no analytical formula exists. If the depth is smaller,
 * the approximation will be less accurate but faster.
 *
 * WARNING: Use these functionalities with care! Extensive use might result in wrong results or crashing simulations.
 */
class BodyAndVolumeFractionMapping
{
public:

   typedef std::pair< pe::BodyID, real_t >                              BodyAndVolumeFraction_T;
   typedef GhostLayerField< std::vector< BodyAndVolumeFraction_T >, 1 > BodyAndVolumeFractionField_T;

   BodyAndVolumeFractionMapping( const shared_ptr<StructuredBlockStorage> & blockStorage,
                                 const shared_ptr<pe::BodyStorage> & globalBodyStorage,
                                 const BlockDataID & bodyStorageID,
                                 const BlockDataID & bodyAndVolumeFractionFieldID,
                                 const std::function<bool(pe::BodyID)> & mappingBodySelectorFct = selectRegularBodies,
                                 const real_t velocityUpdatingEpsilon = real_t(0),
                                 const real_t positionUpdatingEpsilon = real_t(0),
                                 const uint_t superSamplingDepth = uint_t(4) )
   : blockStorage_( blockStorage), globalBodyStorage_( globalBodyStorage ), bodyStorageID_( bodyStorageID ),
     bodyAndVolumeFractionFieldID_( bodyAndVolumeFractionFieldID ), mappingBodySelectorFct_( mappingBodySelectorFct ),
     velocityUpdatingEpsilonSquared_( velocityUpdatingEpsilon * velocityUpdatingEpsilon ),
     positionUpdatingEpsilonSquared_( positionUpdatingEpsilon * positionUpdatingEpsilon ),
     superSamplingDepth_( superSamplingDepth )
   {
      initialize();
   }

   void operator()()
   {
      update();
   }

   void initialize();
   void update();

private:

   void updatePSMBodyAndVolumeFraction( pe::BodyID body, IBlock & block,
                                        BodyAndVolumeFractionField_T * oldBodyAndVolumeFractionField,
                                        std::map< walberla::id_t, Vector3< real_t > > & tempLastUpdatedPositionMap );

   shared_ptr<StructuredBlockStorage> blockStorage_;
   shared_ptr<pe::BodyStorage> globalBodyStorage_;

   const BlockDataID bodyStorageID_;
   const BlockDataID bodyAndVolumeFractionFieldID_;

   std::function<bool(pe::BodyID)> mappingBodySelectorFct_;

   shared_ptr<BodyAndVolumeFractionField_T> updatedBodyAndVolumeFractionField_;
   std::map< walberla::id_t, Vector3< real_t > > lastUpdatedPositionMap_;

   const real_t velocityUpdatingEpsilonSquared_;
   const real_t positionUpdatingEpsilonSquared_;
   const uint_t superSamplingDepth_;
};

void BodyAndVolumeFractionMapping::initialize()
{
   for( auto blockIt = blockStorage_->begin(); blockIt != blockStorage_->end(); ++blockIt )
   {
      BodyAndVolumeFractionField_T * bodyAndVolumeFractionField = blockIt->getData< BodyAndVolumeFractionField_T >( bodyAndVolumeFractionFieldID_ );

      if( updatedBodyAndVolumeFractionField_ == NULL )
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
            mapPSMBodyAndVolumeFraction( bodyIt.getBodyID(), *blockIt, *blockStorage_, bodyAndVolumeFractionFieldID_ );
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

               const real_t fraction = overlapFractionPe( *body, cellCenter, blockStorage_->dx( blockStorage_->getLevel( block ) ), superSamplingDepth_ );

               // if the cell intersected with the body, store a pointer to that body and the corresponding volume fraction in the field
               if( fraction > real_t(0) )
               {
                  updatedBodyAndVolumeFractionField_->get(x,y,z).push_back( BodyAndVolumeFraction_T( body, fraction ) );
               }
            }
         }
      }

      tempLastUpdatedPositionMap.insert( std::pair< walberla::id_t, Vector3< real_t > >( body->getSystemID(), body->getPosition() ) );
   }
}


} // namespace pe_coupling
} // namespace walberla

