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
//! \file SolidVolumeFractionFieldEvaluator.h
//! \ingroup pe_coupling
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/DataTypes.h"

#include "domain_decomposition/StructuredBlockStorage.h"

#include "field/distributors/all.h"
#include "field/GhostLayerField.h"

#include "pe/rigidbody/BodyIterators.h"

#include "pe_coupling/utility/BodySelectorFunctions.h"

#include <functional>

namespace walberla {
namespace pe_coupling {
namespace discrete_particle_methods {

/*!\brief Evaluator of the solid volume fraction field.
 *
 * Updates the solid volume fraction field. Includes firstly removing all old entries from the field, then remapping
 * the local bodies' volumes to the cells.
 * Potentially writes to the ghost layer of the solid volume fraction field which thus requires a special
 * PullReducion communication afterwards.
 *
 * Depending on the used Distributor_T, the resulting solid volume fraction field will vary:
 * - NearestNeighborDistributor:
 *       Corresponds to simply assigning the total body volume to the cell containing the body center.
 * - KernelDistributor:
 *       The body's volume is not directly assigned to one cell only but spread over the neighboring cells as well.
 *       See Finn, Li, Apte - "Particle based modelling and simulation of natural sand dynamics in the wave bottom boundary layer" (2016)
 *       for the application, even though different kernel was used there.
 *
 * Whether or not a body gets treated by the evaluator depends on the return value of 'dpmBodySelectorFct'.
 *
 * For more infos on distributors, see src/field/distributors.
 */
template< typename FlagField_T, template <typename,typename> class Distributor_T >
class SolidVolumeFractionFieldEvaluator
{
public:

   using ScalarField_T = GhostLayerField<real_t, 1>;
   using ScalarDistributor_T = Distributor_T<ScalarField_T, FlagField_T>;

   SolidVolumeFractionFieldEvaluator( const shared_ptr<StructuredBlockStorage> & blockStorage,
                                      const BlockDataID & solidVolumeFractionFieldID, const BlockDataID & bodyStorageID,
                                      const BlockDataID & flagFieldID, const Set< FlagUID > & domainFlags,
                                      const std::function<bool(pe::BodyID)> & dpmBodySelectorFct = selectRegularBodies )
   : blockStorage_( blockStorage ), solidVolumeFractionFieldID_( solidVolumeFractionFieldID ),
     bodyStorageID_( bodyStorageID ), dpmBodySelectorFct_( dpmBodySelectorFct)
   {
      scalarDistributorID_ = field::addDistributor< ScalarDistributor_T, FlagField_T >( blockStorage, solidVolumeFractionFieldID, flagFieldID, domainFlags );
   }

   void operator()( IBlock * const block )
   {
      ScalarField_T* svfField = block->getData< ScalarField_T >( solidVolumeFractionFieldID_ );

      // reset field
      WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ( svfField,
         svfField->get(x,y,z) = real_c(0);
      )

      ScalarDistributor_T * scalarDistributor = block->getData< ScalarDistributor_T >( scalarDistributorID_ );

      // assign the local bodies' volume to the cell, depending on the chosen Distributor_T
      for( auto bodyIt = pe::LocalBodyIterator::begin(*block, bodyStorageID_); bodyIt != pe::LocalBodyIterator::end(); ++bodyIt )
      {
         if(!dpmBodySelectorFct_(bodyIt.getBodyID())) continue;

         real_t bodyVolume = bodyIt->getVolume();
         const Vector3<real_t> bodyPosition = bodyIt->getPosition();

         scalarDistributor->distribute( bodyPosition, &bodyVolume );
      }
   }

private:

   shared_ptr<StructuredBlockStorage> blockStorage_;
   BlockDataID solidVolumeFractionFieldID_;
   BlockDataID bodyStorageID_;

   std::function<bool(pe::BodyID)> dpmBodySelectorFct_;

   BlockDataID scalarDistributorID_;
};


} // namespace discrete_particle_methods
} // namespace pe_coupling
} // namespace walberla
