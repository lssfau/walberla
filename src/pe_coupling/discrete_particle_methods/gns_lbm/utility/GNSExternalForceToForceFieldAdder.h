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
//! \file GNSExternalForceToForceFieldAdder.h
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

/*!\brief Applies a given constant external force to a force field.
 *
 * Utility function when using external forces with methods that couple over the force field.
 * GNS speciality: the external force is additionally multiplied by the fluid volume fraction before it is added.
 */
class GNSExternalForceToForceFieldAdder
{
public:
   using ForceField_T = GhostLayerField<Vector3<real_t>, 1>;
   using ScalarField_T = GhostLayerField<real_t, 1>;

   GNSExternalForceToForceFieldAdder( const BlockDataID & forceFieldID, const ConstBlockDataID & solidVolumeFractionFieldID, const Vector3<real_t> & externalForce )
      : forceFieldID_( forceFieldID ), solidVolumeFractionFieldID_( solidVolumeFractionFieldID ), externalForce_( externalForce )
   {
   }

   void operator()( IBlock * const block )
   {
      ForceField_T* forceField = block->getData< ForceField_T >( forceFieldID_ );
      const ScalarField_T* solidVolumeFractionField = block->getData< ScalarField_T >( solidVolumeFractionFieldID_ );

      //TODO include dx force scaling
      WALBERLA_FOR_ALL_CELLS_XYZ( forceField,
         forceField->get(x,y,z) += ( real_t(1) - solidVolumeFractionField->get(x,y,z) ) * externalForce_;
      );
   }

   void reset( const Vector3<real_t> & newExternalForce )
   {
      externalForce_ = newExternalForce;
   }

private:
   const BlockDataID forceFieldID_;
   const ConstBlockDataID solidVolumeFractionFieldID_;
   Vector3<real_t> externalForce_;
};


} // namespace discrete_particle_methods
} // namespace pe2_coupling
} // namespace walberla
