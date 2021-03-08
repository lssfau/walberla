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
#include "pe/Types.h"
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
                                  const BlockDataID bodyAndVolumeFractionFieldID );

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

   using BodyAndVolumeFraction_T = std::pair<pe::BodyID, real_t>;
   using BodyAndVolumeFractionField_T = GhostLayerField<std::vector<BodyAndVolumeFraction_T>, 1>;

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



} // namespace pe_coupling
} // namespace walberla

