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
//! \file AveragedInteractionForceFieldToForceFieldAdder.h
//! \ingroup pe_coupling
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/math/Vector3.h"

#include "domain_decomposition/StructuredBlockStorage.h"

#include "field/GhostLayerField.h"

namespace walberla {
namespace pe_coupling {
namespace discrete_particle_methods {

using math::Vector3;

/*!\brief Adds an accumulated force field value to another force field, applying continuous averaging.
 *
 * The force field contains all external force, etc acting on the fluid, except the fluid-solid interaction forces
 * from the DPS coupling.
 * The interaction force field only contains the interaction forces, simply added up over several subcycling timesteps,
 * with the maximum number given by maximumAveragingSteps.
 * This class will average the interaction force over the currently carried out sub cycles (counted via an internal counter)
 * and add this averaged quantity to the force field, to get rid of possible oscillations in the interaction force
 * and have a valid force field in each sub cycling iteration.
 *
 * Can also be used to simply add the current interaction force field to the force field if maximumAveragingSteps = 1.
 */
class AveragedInteractionForceFieldToForceFieldAdder
{
public:
   using ForceField_T = GhostLayerField<Vector3<real_t>, 1>;

   AveragedInteractionForceFieldToForceFieldAdder( const shared_ptr<StructuredBlockStorage> & blockStorage, const BlockDataID & forceFieldID,
                                                   const ConstBlockDataID & interactionForceFieldID, uint_t maximumAveragingSteps )
      : blockStorage_( blockStorage), forceFieldID_( forceFieldID ), interactionForceFieldID_( interactionForceFieldID ),
        maximumAveragingSteps_( maximumAveragingSteps ), averagingCounter_( uint_t(0) )
   {
   }

   void operator()()
   {
      ++averagingCounter_;

      const real_t averagingFactor = real_t(1)/real_c(averagingCounter_);


      for( auto blockIt = blockStorage_->begin(); blockIt != blockStorage_->end(); ++blockIt )
      {
         ForceField_T* forceField = blockIt->getData< ForceField_T >( forceFieldID_ );
         const ForceField_T* interactionForceField = blockIt->getData< ForceField_T >( interactionForceFieldID_ );

         WALBERLA_FOR_ALL_CELLS_XYZ( forceField,
            forceField->get(x,y,z) += averagingFactor * interactionForceField->get(x,y,z);
         );
      }

      if( averagingCounter_ >= maximumAveragingSteps_ )
      {
         resetAveragingCounter();
      }

   }

   void resetAveragingCounter( )
   {
      averagingCounter_ = uint_t(0);
   }

private:
   shared_ptr<StructuredBlockStorage> blockStorage_;
   const BlockDataID forceFieldID_;
   const ConstBlockDataID interactionForceFieldID_;
   uint_t maximumAveragingSteps_;
   uint_t averagingCounter_;
};


} // namespace discrete_particle_methods
} // namespace pe_coupling
} // namespace walberla
