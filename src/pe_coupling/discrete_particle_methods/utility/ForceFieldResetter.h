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
//! \file ForceFieldResetter.h
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

/*!\brief Resets the values currently stored in the force field to (0,0,0).
 */
class ForceFieldResetter
{
public:
   using ForceField_T = GhostLayerField<Vector3<real_t>, 1>;

   ForceFieldResetter( const BlockDataID & forceFieldID )
      : forceFieldID_( forceFieldID )
   {
   }

   void operator()( IBlock * const block )
   {
      ForceField_T* forceField = block->getData<ForceField_T>(forceFieldID_);
      WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ(forceField,
                                                       forceField->get(x,y,z).reset();
      )
   }

private:
   const BlockDataID forceFieldID_;
};


} // namespace discrete_particle_methods
} // namespace pe_coupling
} // namespace walberla
