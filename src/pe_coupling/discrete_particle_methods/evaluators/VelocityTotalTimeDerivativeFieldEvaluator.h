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
//! \file VelocityTotalTimeDerivativeFieldEvaluator.h
//! \ingroup pe_coupling
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/math/Vector3.h"
#include "core/math/Matrix3.h"

#include "domain_decomposition/StructuredBlockStorage.h"

#include "field/GhostLayerField.h"

namespace walberla {
namespace pe_coupling {
namespace discrete_particle_methods {

using math::Vector3;


/*!\brief Evaluator of the total derivative of the fluid velocity.
 *
 * This evaluates Du/Dt = du/dt + u * grad(u).
 * The gradient of the velocity, grad(u), has to be provided by a field, e.g. evaluated with VelocityGradientFieldEvaluator.h.
 * The time derivative, du/dt, is approximated by a simple backward difference: du/dt = (u_new - u_old ) / deltaT.
 * This requires two velocity field.
 *
 */
class VelocityTotalTimeDerivativeFieldEvaluator
{
public:
   using VelocityField_T = GhostLayerField<Vector3<real_t>, 1>;
   using TensorField_T = GhostLayerField<Matrix3<real_t>, 1>;

   VelocityTotalTimeDerivativeFieldEvaluator( const BlockDataID & totalTimeDerivativeVelocityFieldID,
                                              const ConstBlockDataID & currentVelocityFieldID,
                                              const ConstBlockDataID & formerVelocityFieldID,
                                              const ConstBlockDataID & velocityGradientFieldID,
                                              const real_t & deltaT = real_t(1) )
      : totalTimeDerivativeVelocityFieldID_( totalTimeDerivativeVelocityFieldID ), currentVelocityFieldID_( currentVelocityFieldID ),
        formerVelocityFieldID_( formerVelocityFieldID ), velocityGradientFieldID_( velocityGradientFieldID ), deltaTinv_( real_t(1) / deltaT )
   { }

   void operator()(IBlock * const block)
   {
      VelocityField_T* totalTimeDerivativeVelocityField = block->getData< VelocityField_T >( totalTimeDerivativeVelocityFieldID_ );
      const VelocityField_T* currentVelocityField       = block->getData< VelocityField_T >( currentVelocityFieldID_ );
      const VelocityField_T* formerVelocityField        = block->getData< VelocityField_T >( formerVelocityFieldID_ );
      const TensorField_T*   velocityGradientField      = block->getData< TensorField_T >( velocityGradientFieldID_ );

      // simple backward difference approximation of Du/Dt
      WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ( totalTimeDerivativeVelocityField,
         Vector3<real_t> vel = currentVelocityField->get(x,y,z);

         // grad(u) =
         // | du1/dx1 du2/dx1 du3/dx1 |   | 0 1 2 |   | 0,0  0,1  0,2 |
         // | du1/dx2 du2/dx2 du3/dx2 | = | 3 4 5 | = | 1,0  1,1  1,2 |
         // | du1/dx3 du2/dx3 du3/dx3 |   | 6 7 8 |   | 2,0  2,1  2,2 |
         Matrix3<real_t> gradVel = velocityGradientField->get(x,y,z);

         // evaluate u * grad(u) = u_i (d u_j / d x_i)
         Vector3<real_t> velGradVel( gradVel(0,0) * vel[0] + gradVel(1,0) * vel[1] + gradVel(2,0) * vel[2],
                                     gradVel(0,1) * vel[0] + gradVel(1,1) * vel[1] + gradVel(2,1) * vel[2],
                                     gradVel(0,2) * vel[0] + gradVel(1,2) * vel[1] + gradVel(2,2) * vel[2] );

         totalTimeDerivativeVelocityField->get(x,y,z) = ( vel - formerVelocityField->get(x,y,z) ) * deltaTinv_ + velGradVel;
      );
   }

   void resetDeltaT( const real_t & deltaT )
   {
      deltaTinv_ = real_t(1) / deltaT;
   }

private:
   const BlockDataID totalTimeDerivativeVelocityFieldID_;
   const ConstBlockDataID currentVelocityFieldID_;
   const ConstBlockDataID formerVelocityFieldID_;
   const ConstBlockDataID velocityGradientFieldID_;
   real_t deltaTinv_;
};


} // namespace discrete_particle_methods
} // namespace pe_coupling
} // namespace walberla
