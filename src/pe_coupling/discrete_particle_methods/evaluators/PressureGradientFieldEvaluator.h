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
//! \file PressureGradientFieldEvaluator.h
//! \ingroup pe_coupling
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/math/Vector3.h"

#include "domain_decomposition/StructuredBlockStorage.h"

#include "field/GhostLayerField.h"

#include "stencil/Directions.h"

namespace walberla {
namespace pe_coupling {
namespace discrete_particle_methods {

using math::Vector3;

/*!\brief Evaluator of the pressure gradient field, given a pressure field.
 *
 */
template <typename LatticeModel_T, typename BoundaryHandling_T>
class PressureGradientFieldEvaluator
{
public:
   using VectorField_T = GhostLayerField<Vector3<real_t>, 1>;
   using ScalarField_T = GhostLayerField<real_t, 1>;
   using Stencil_T = typename LatticeModel_T::Stencil;

   PressureGradientFieldEvaluator( const BlockDataID & pressureGradientFieldID, const ConstBlockDataID & pressureFieldID,
                                   const ConstBlockDataID & boundaryHandlingID )
      : pressureGradientFieldID_( pressureGradientFieldID ), pressureFieldID_( pressureFieldID ),
        boundaryHandlingID_( boundaryHandlingID )
   { }

   void operator()(IBlock * const block)
   {
      VectorField_T* pressureGradientField    = block->getData< VectorField_T >( pressureGradientFieldID_ );
      const ScalarField_T* pressureField            = block->getData< ScalarField_T >( pressureFieldID_ );
      const BoundaryHandling_T * boundaryHandling   = block->getData< BoundaryHandling_T >( boundaryHandlingID_ );

      WALBERLA_FOR_ALL_CELLS_XYZ( pressureGradientField,
         if( boundaryHandling->isDomain(x,y,z) )
         {
            pressureGradientField->get(x,y,z) = getPressureGradient( Cell(x,y,z), pressureField, boundaryHandling );
         }
      );
   }

private:

   Vector3<real_t> getPressureGradient( const Cell & cell, const ScalarField_T * pressureField, const BoundaryHandling_T * boundaryHandling )
   {

      // temporarily store pressure values of surrounding cells
      std::vector<real_t> pressureValues( Stencil_T::Size, real_t(0) );

      // store pressure value in center cell to potentially apply Neumann like boundary conditions
      real_t pressureInCenterCell = pressureField->get( cell );

      for( auto dir = Stencil_T::beginNoCenter(); dir != Stencil_T::end(); ++dir)
      {
         // check if boundary treatment is necessary
         if( !boundaryHandling->isDomain( cell + *dir ) )
         {
            // copy from center cell
            pressureValues[ *dir ] = pressureInCenterCell;
         } else 
         {
            pressureValues[ *dir ] = pressureField->get( cell + *dir );

         }
      }

      // apply gradient formula
      // uses LBM weighting to determine gradient
      // See: Ramadugu et al - "Lattice differential operators for computational physics" (2013)
      // with T = c_s**2
      const real_t inv_c_s_sqr = real_c(3);
      Vector3<real_t> gradient( real_c(0) );
      for( auto dir = Stencil_T::beginNoCenter(); dir != Stencil_T::end(); ++dir)
      {
         real_t pressure = pressureValues[ *dir ];
         gradient[0] += LatticeModel_T::w[ dir.toIdx() ] * real_c(dir.cx()) * pressure;
         gradient[1] += LatticeModel_T::w[ dir.toIdx() ] * real_c(dir.cy()) * pressure;
         gradient[2] += LatticeModel_T::w[ dir.toIdx() ] * real_c(dir.cz()) * pressure;
      }
      gradient *= inv_c_s_sqr;

      return gradient;
   }

   const BlockDataID pressureGradientFieldID_;
   const ConstBlockDataID pressureFieldID_;
   const ConstBlockDataID boundaryHandlingID_;
};


} // namespace discrete_particle_methods
} // namespace pe_coupling
} // namespace walberla
