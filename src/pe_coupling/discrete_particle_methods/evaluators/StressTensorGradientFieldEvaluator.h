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
//! \file StressTensorGradientFieldEvaluator.h
//! \ingroup pe_coupling
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/math/Vector3.h"
#include "core/math/Matrix3.h"

#include "domain_decomposition/StructuredBlockStorage.h"

#include "field/GhostLayerField.h"

#include "stencil/Directions.h"

namespace walberla {
namespace pe_coupling {
namespace discrete_particle_methods {

using math::Vector3;

/*!\brief Evaluator of the stress tensor gradient field, given a velocity field.
 *
 */
template <typename LatticeModel_T, typename BoundaryHandling_T>
class StressTensorGradientFieldEvaluator
{
public:
   using TensorField_T = GhostLayerField<Matrix3<real_t>, 1>;
   using VectorField_T = GhostLayerField<Vector3<real_t>, 1>;
   using Stencil_T = typename LatticeModel_T::Stencil;

   StressTensorGradientFieldEvaluator( const BlockDataID & stressTensorGradientFieldID,
                                       const ConstBlockDataID & velocityGradientFieldID,
                                       const ConstBlockDataID & boundaryHandlingID,
                                       const real_t & dynamicFluidViscosity )
      : stressTensorGradientFieldID_( stressTensorGradientFieldID ),
        velocityGradientFieldID_( velocityGradientFieldID ),
        boundaryHandlingID_( boundaryHandlingID ),
        dynamicFluidViscosity_( dynamicFluidViscosity )
   { }

   void operator()(IBlock * const block)
   {
      VectorField_T* stressTensorGradientField = block->getData< VectorField_T >( stressTensorGradientFieldID_ );
      const TensorField_T* velocityGradientField     = block->getData< TensorField_T >( velocityGradientFieldID_ );
      const BoundaryHandling_T * boundaryHandling    = block->getData< BoundaryHandling_T >( boundaryHandlingID_ );

      WALBERLA_FOR_ALL_CELLS_XYZ( stressTensorGradientField,
         if( boundaryHandling->isDomain(x,y,z) )
         {
            stressTensorGradientField->get(x,y,z) = getStressTensorGradient( Cell(x,y,z), velocityGradientField, boundaryHandling );
         }
      );
   }

   void resetViscosity( real_t newDynamicFluidViscosity )
   {
      dynamicFluidViscosity_ = newDynamicFluidViscosity;
   }

private:

   // calculates grad( nu * ( ( grad(u) ) + (grad(u))**T ) ), i.e. the gradient of the stress tensor
   Vector3<real_t> getStressTensorGradient( const Cell & cell, const TensorField_T* velocityGradientField, const BoundaryHandling_T * boundaryHandling )
   {

      std::vector< Matrix3< real_t > > stressTensorValues( Stencil_T::Size, Matrix3< real_t >(real_t(0)) );


      Matrix3< real_t > velGradientInCenterCell = velocityGradientField->get( cell );
      Matrix3< real_t > stressTensorInCenterCell = ( velGradientInCenterCell + velGradientInCenterCell.getTranspose() ) * dynamicFluidViscosity_;

      for( auto dir = Stencil_T::beginNoCenter(); dir != Stencil_T::end(); ++dir)
      {
         // check if boundary treatment is necessary
         if( !boundaryHandling->isDomain( cell + *dir ) )
         {
            // copy from center cell
            stressTensorValues[ *dir ] = stressTensorInCenterCell;
         } else {
            Matrix3< real_t > velGradient = velocityGradientField->get( cell + *dir );
            stressTensorValues[ *dir ] = ( velGradient + velGradient.getTranspose() ) * dynamicFluidViscosity_;
         }
      }

      // obtain the gradient of the tensor using the gradient formula
      // See: Ramadugu et al - "Lattice differential operators for computational physics" (2013)
      // with T = c_s**2
      const real_t inv_c_s_sqr = real_t(3);
      Vector3<real_t> gradStressTensor( real_t(0) );
      for( auto dir = Stencil_T::beginNoCenter(); dir != Stencil_T::end(); ++dir)
      {
         Vector3<real_t> latticeVel( real_c(dir.cx()), real_c(dir.cy()), real_c(dir.cz()) );
         Matrix3<real_t> & tempTau = stressTensorValues[ *dir ];
         Vector3<real_t> tau1( tempTau[0], tempTau[3], tempTau[6] );
         Vector3<real_t> tau2( tempTau[1], tempTau[4], tempTau[7] );
         Vector3<real_t> tau3( tempTau[2], tempTau[5], tempTau[8] );

         gradStressTensor[0] += LatticeModel_T::w[ dir.toIdx() ] * latticeVel * tau1;
         gradStressTensor[1] += LatticeModel_T::w[ dir.toIdx() ] * latticeVel * tau2;
         gradStressTensor[2] += LatticeModel_T::w[ dir.toIdx() ] * latticeVel * tau3;

      }
      gradStressTensor *= inv_c_s_sqr;

      return gradStressTensor;

   }

   const BlockDataID stressTensorGradientFieldID_;
   const ConstBlockDataID velocityGradientFieldID_;
   const ConstBlockDataID boundaryHandlingID_;
   real_t dynamicFluidViscosity_;
};


} // namespace discrete_particle_methods
} // namespace pe_coupling
} // namespace walberla
