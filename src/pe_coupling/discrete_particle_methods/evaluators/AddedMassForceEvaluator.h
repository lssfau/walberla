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
//! \file AddedMassForceEvaluator.h
//! \ingroup pe_coupling
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/DataTypes.h"
#include "core/debug/Debug.h"

#include "domain_decomposition/StructuredBlockStorage.h"

#include "field/interpolators/all.h"
#include "field/distributors/all.h"
#include "field/GhostLayerField.h"

#include "lbm/field/MacroscopicValueCalculation.h"
#include "lbm/field/PdfField.h"
#include "lbm/lattice_model/ForceModel.h"

#include "pe/rigidbody/BodyIterators.h"
#include "pe/rigidbody/RigidBody.h"
#include "pe/Types.h"
#include "pe/Materials.h"

#include "pe_coupling/utility/BodySelectorFunctions.h"

#include "BodyVelocityTimeDerivativeEvaluator.h"

#include "stencil/Directions.h"

#include <functional>

namespace walberla {
namespace pe_coupling {
namespace discrete_particle_methods {

/*!\brief Evaluator of the added (virtual) mass force based on the given added mass correlation.
 *
 * The added mass force is evaluated for each body (requires interpolation of several quantities from the fluid fields)
 * and then applied onto them.
 * The corresponding reaction force is added to the fluid via the given distributor.
 *
 * Note that the fluid velocity (and its derivatives) has to be the fluid-phase velocity (and not the volume-averaged velocity).
 *
 * Whether or not a body gets treated by the evaluator depends on the return value of 'dpmBodySelectorFct'.
 *
 * For more infos on interpolators, see field interpolators in src/field/interpolators.
 * For more infos on distributors, see src/field/distributors.
 */
template< typename FlagField_T, template<typename,typename> class FieldInterpolator_T, template<typename,typename> class Distributor_T >
class AddedMassForceEvaluator
{  

public:

   using Vec3Field_T = GhostLayerField<Vector3<real_t>, 1>;
   using Vec3FieldInterpolator_T = FieldInterpolator_T<Vec3Field_T, FlagField_T>;
   using ForceDistributor_T = Distributor_T<Vec3Field_T, FlagField_T>;

   AddedMassForceEvaluator( const shared_ptr<StructuredBlockStorage> & blockStorage,
                            const BlockDataID & forceFieldID, const BlockDataID & bodyStorageID,
                            const BlockDataID & flagFieldID, const Set< FlagUID > & domainFlags,
                            const BlockDataID & velocityTimeDerivativeFieldID,
                            const std::function<Vector3<real_t> ( const Vector3<real_t> &, const Vector3<real_t> &, real_t, real_t )> & addedMassForceCorrelationFunction,
                            const shared_ptr< BodyVelocityTimeDerivativeEvaluator > & bodyVelocityTimeDerivativeEvaluator,
                            const std::function<bool(pe::BodyID)> & dpmBodySelectorFct = selectRegularBodies )
   : blockStorage_( blockStorage ), bodyStorageID_( bodyStorageID ),
     addedMassForceCorrelationFunction_( addedMassForceCorrelationFunction ),
     bodyVelocityTimeDerivativeEvaluator_( bodyVelocityTimeDerivativeEvaluator ),
     dpmBodySelectorFct_( dpmBodySelectorFct)
   {
      velocityTimeDerivativeFieldInterpolatorID_ = field::addFieldInterpolator< Vec3FieldInterpolator_T, FlagField_T >( blockStorage, velocityTimeDerivativeFieldID, flagFieldID, domainFlags );
      forceDistributorID_ = field::addDistributor< ForceDistributor_T, FlagField_T >( blockStorage, forceFieldID, flagFieldID, domainFlags );
   }

   void operator()();

private:

   shared_ptr<StructuredBlockStorage> blockStorage_;
   const BlockDataID bodyStorageID_;

   std::function<Vector3<real_t> ( const Vector3<real_t> &, const Vector3<real_t> &, real_t, real_t )> addedMassForceCorrelationFunction_;

   shared_ptr< BodyVelocityTimeDerivativeEvaluator > bodyVelocityTimeDerivativeEvaluator_;

   std::function<bool(pe::BodyID)> dpmBodySelectorFct_;

   BlockDataID velocityTimeDerivativeFieldInterpolatorID_;
   BlockDataID forceDistributorID_;

};

template< typename FlagField_T, template<typename,typename> class FieldInterpolator_T, template<typename,typename> class Distributor_T >
void AddedMassForceEvaluator< FlagField_T, FieldInterpolator_T, Distributor_T >
::operator()()
{

   for( auto blockIt = blockStorage_->begin(); blockIt != blockStorage_->end(); ++blockIt )
   {

      Vec3FieldInterpolator_T * velocityTimeDerivativeInterpolator = blockIt->getData< Vec3FieldInterpolator_T >( velocityTimeDerivativeFieldInterpolatorID_ );
      ForceDistributor_T * forceDistributor                        = blockIt->getData< ForceDistributor_T >( forceDistributorID_ );

      for( auto bodyIt = pe::LocalBodyIterator::begin(*blockIt, bodyStorageID_); bodyIt != pe::LocalBodyIterator::end(); ++bodyIt )
      {
         if(!dpmBodySelectorFct_(bodyIt.getBodyID())) continue;

         Vector3<real_t> forceOnFluid( real_t(0) );

         Vector3<real_t> bodyPosition = bodyIt->getPosition();
         Vector3<real_t> bodyVelocity = bodyIt->getLinearVel();
         real_t bodyVolume = bodyIt->getVolume();
         real_t fluidDensity( real_t(1) );

         // evaluate added (virtual) mass force
         Vector3<real_t> timeDerivativeFluidVelocity( real_t(0) );
         Vector3<real_t> timeDerivativeBodyVelocity( real_t(0) );
         bodyVelocityTimeDerivativeEvaluator_->get(timeDerivativeBodyVelocity, bodyVelocity, bodyIt->getSystemID() );
         velocityTimeDerivativeInterpolator->get( bodyPosition, &timeDerivativeFluidVelocity );
         Vector3<real_t> addedMassForce = addedMassForceCorrelationFunction_( timeDerivativeFluidVelocity, timeDerivativeBodyVelocity, bodyVolume, fluidDensity );
         WALBERLA_ASSERT( !math::isnan(addedMassForce[0]) && !math::isnan(addedMassForce[1]) && !math::isnan(addedMassForce[2]),
                          "NaN found in added mass force for body at position " << bodyPosition );

         bodyIt->addForce( addedMassForce );
         forceOnFluid += ( -addedMassForce );

         // set/distribute force on fluid
         forceDistributor->distribute( bodyPosition, &forceOnFluid );
      }
   }

}

} // namespace discrete_particle_methods
} // namespace pe_coupling
} // namespace walberla
