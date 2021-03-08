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
//! \file ExternalForceToForceFieldAdder.h
//! \ingroup pe_coupling
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/math/Vector3.h"

#include "field/GhostLayerField.h"

namespace walberla {
namespace pe_coupling {
namespace discrete_particle_methods {

using math::Vector3;

/*!\brief Applies a given constant external force to a force field.
 *
 * Utility function when using external forces with methods that couple over the force field
 * (all discrete particle methods, e.g.) since the force field is typically reset in each iteration.
 */
class ExternalForceToForceFieldAdder
{
public:
   using ForceField_T = GhostLayerField<Vector3<real_t>, 1>;

   ExternalForceToForceFieldAdder( const BlockDataID & forceFieldID, const Vector3<real_t> & externalForce )
      : forceFieldID_( forceFieldID ), externalForce_( externalForce )
   {
   }

   void operator()( IBlock * const block )
   {
      ForceField_T* forceField = block->getData< ForceField_T >( forceFieldID_ );

      WALBERLA_FOR_ALL_CELLS_XYZ( forceField,
         //TODO include level dependent dx force scaling
         forceField->get(x,y,z) += externalForce_;
      );
   }

   void reset( const Vector3<real_t> & newExternalForce )
   {
      externalForce_ = newExternalForce;
   }

private:
   const BlockDataID forceFieldID_;
   Vector3<real_t> externalForce_;
};


} // namespace discrete_particle_methods
} // namespace pe_coupling
} // namespace walberla
