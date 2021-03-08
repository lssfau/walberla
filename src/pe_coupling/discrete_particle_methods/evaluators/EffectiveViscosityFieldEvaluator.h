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
//! \file EffectiveViscosityFieldEvaluator.h
//! \ingroup pe_coupling
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//
//======================================================================================================================

#pragma once

#include "lbm/lattice_model/CollisionModel.h"

#include "field/GhostLayerField.h"

namespace walberla {
namespace pe_coupling {
namespace discrete_particle_methods {

// correlations for the viscosity

// no change in viscosity
real_t calculateUnchangedEffectiveViscosity( real_t fluidViscosity, real_t /*porosity*/ )
{
   return fluidViscosity;
}

// see: Fattahi, Waluga, Wohlmuth - "Large scale lattice Boltzmann simulation for the coupling of free and porous media flow"
real_t calculateRescaledEffectiveViscosity( real_t fluidViscosity, real_t porosity )
{
   return fluidViscosity / porosity;
}

// see: J. R. Finn, M. Li, S. V. Apte - "Particle based modelling and simulation of natural sand dynamics in the wave
// bottom boundary layer", Journal of Fluid Mechanics 796 (2016) 340–385. doi:10.1017/jfm.2016.246.
real_t calculateEilersEffectiveViscosity( real_t fluidViscosity, real_t porosity )
{
   const real_t closePackingFraction = real_t(0.64);
   const real_t intrinsicViscosity = real_t(2.5); //for monosized spheres
   const real_t temp = real_t(1) + real_t(0.5) * intrinsicViscosity * ( real_t(1) - porosity ) / ( porosity  / closePackingFraction );
   return fluidViscosity * temp * temp;
}

/*!\brief Evaluates the (local) effective viscosity.
 *
 * The effective viscosity (using the provided correlation) is evaluated for each cell and stored,
 * as the respective omega value, in a field, that can be given to a suitable LBM sweep.
 *
 */
class EffectiveViscosityFieldEvaluator
{
public:
   using ScalarField_T = GhostLayerField<real_t, 1>;

   EffectiveViscosityFieldEvaluator( const BlockDataID & omegaFieldID, const ConstBlockDataID & solidVolumeFractionFieldID,
                                     const real_t & fluidViscosity,
                                     const std::function<real_t (real_t, real_t)> & effectiveViscosityFunc )
      : omegaFieldID_( omegaFieldID ), solidVolumeFractionFieldID_( solidVolumeFractionFieldID ), fluidViscosity_( fluidViscosity ),
        effectiveViscosityFunc_( effectiveViscosityFunc )
   {
   }

   void operator()( IBlock * const block )
   {
      ScalarField_T* omegaField = block->getData<ScalarField_T>(omegaFieldID_);
      const ScalarField_T* svfField   = block->getData<ScalarField_T>(solidVolumeFractionFieldID_);

      WALBERLA_FOR_ALL_CELLS_XYZ(omegaField,
          const real_t porosity = real_t(1) - svfField->get(x,y,z);
          WALBERLA_ASSERT_FLOAT_UNEQUAL(porosity, real_t(0));

          real_t effectiveViscosity = effectiveViscosityFunc_(fluidViscosity_, porosity);

          //TODO: add level dependence for omega calculation
          omegaField->get(x,y,z) = lbm::collision_model::omegaFromViscosity(effectiveViscosity);
      )
   }

   void resetFluidViscosity( real_t newFluidViscosity )
   {
      fluidViscosity_ = newFluidViscosity;
   }

private:
   const BlockDataID omegaFieldID_;
   const ConstBlockDataID solidVolumeFractionFieldID_;
   real_t fluidViscosity_;
   std::function<real_t (real_t, real_t)> effectiveViscosityFunc_;
};


} // namespace discrete_particle_methods
} // namespace pe_coupling
} // namespace walberla
