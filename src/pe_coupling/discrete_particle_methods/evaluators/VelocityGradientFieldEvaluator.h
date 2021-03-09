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
//! \file VelocityGradientFieldEvaluator.h
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


/*!\brief Evaluator of the velocity gradient field, given a velocity field.
 *
 */
template <typename LatticeModel_T, typename BoundaryHandling_T>
class VelocityGradientFieldEvaluator
{
public:
   using TensorField_T = GhostLayerField<Matrix3<real_t>, 1>;
   using VectorField_T = GhostLayerField<Vector3<real_t>, 1>;
   using Stencil_T = typename LatticeModel_T::Stencil;

   VelocityGradientFieldEvaluator( const BlockDataID & velocityGradientFieldID,
                                   const ConstBlockDataID & velocityFieldID,
                                   const ConstBlockDataID & boundaryHandlingID )
      : velocityGradientFieldID_( velocityGradientFieldID ),
        velocityFieldID_( velocityFieldID ),
        boundaryHandlingID_( boundaryHandlingID )
   { }

   void operator()(IBlock * const block)
   {
      const VectorField_T* velocityField = block->getData< VectorField_T >( velocityFieldID_ );
      TensorField_T* velocityGradientField = block->getData< TensorField_T >( velocityGradientFieldID_ );
      const BoundaryHandling_T * boundaryHandling = block->getData< BoundaryHandling_T >( boundaryHandlingID_ );

      Matrix3<real_t> velocityGradient( real_t(0) );

      WALBERLA_FOR_ALL_CELLS_XYZ( velocityGradientField,
         if( boundaryHandling->isDomain(x,y,z) )
         {
            getVelocityGradient( Cell(x,y,z), velocityField, boundaryHandling, velocityGradient );
            velocityGradientField->get(x,y,z) = velocityGradient;
         }
      );
   }

private:

   // calculates grad(u)
   // grad(u) =
   // | du1/dx1 du2/dx1 du3/dx1 |   | 0 1 2 |   | 0,0  0,1  0,2 |
   // | du1/dx2 du2/dx2 du3/dx2 | = | 3 4 5 | = | 1,0  1,1  1,2 |
   // | du1/dx3 du2/dx3 du3/dx3 |   | 6 7 8 |   | 2,0  2,1  2,2 |
   void getVelocityGradient( const Cell & cell, const VectorField_T * velocityField, const BoundaryHandling_T * boundaryHandling, Matrix3<real_t> & velocityGradient )
   {

      std::vector< Vector3<real_t> > velocityValues( Stencil_T::Size, Vector3<real_t>(real_t(0)) );

      Vector3<real_t> velocityInCenterCell = velocityField->get( cell );

      for( auto dir = Stencil_T::beginNoCenter(); dir != Stencil_T::end(); ++dir)
      {
         // check if boundary treatment is necessary
         if( !boundaryHandling->isDomain( cell + *dir ) )
         {
            // copy from center cell
            velocityValues[ *dir ] = velocityInCenterCell;
         } else {
            velocityValues[ *dir ] = velocityField->get( cell + *dir );
         }
      }

      // obtain the matrix grad(u) with the help of the gradient formula from
      // See: Ramadugu et al - Lattice differential operators for computational physics (2013)
      // with T = c_s**2
      const real_t inv_c_s_sqr = real_t(3);
      velocityGradient = real_t(0);
      for( auto dir = Stencil_T::beginNoCenter(); dir != Stencil_T::end(); ++dir)
      {
         real_t cx = real_c(dir.cx());
         real_t cy = real_c(dir.cy());
         real_t cz = real_c(dir.cz());

         // grad(ux)
         real_t ux = velocityValues[ *dir ][0];
         velocityGradient[ 0 ] += LatticeModel_T::w[ dir.toIdx() ] * cx * ux;
         velocityGradient[ 3 ] += LatticeModel_T::w[ dir.toIdx() ] * cy * ux;
         velocityGradient[ 6 ] += LatticeModel_T::w[ dir.toIdx() ] * cz * ux;

         // grad(uy)
         real_t uy = velocityValues[ *dir ][1];
         velocityGradient[ 1 ] += LatticeModel_T::w[ dir.toIdx() ] * cx * uy;
         velocityGradient[ 4 ] += LatticeModel_T::w[ dir.toIdx() ] * cy * uy;
         velocityGradient[ 7 ] += LatticeModel_T::w[ dir.toIdx() ] * cz * uy;

         // grad(uz)
         real_t uz = velocityValues[ *dir ][2];
         velocityGradient[ 2 ] += LatticeModel_T::w[ dir.toIdx() ] * cx * uz;
         velocityGradient[ 5 ] += LatticeModel_T::w[ dir.toIdx() ] * cy * uz;
         velocityGradient[ 8 ] += LatticeModel_T::w[ dir.toIdx() ] * cz * uz;

      }
      velocityGradient *= inv_c_s_sqr;

   }

   const BlockDataID velocityGradientFieldID_;
   const ConstBlockDataID velocityFieldID_;
   const ConstBlockDataID boundaryHandlingID_;
};


} // namespace discrete_particle_methods
} // namespace pe_coupling
} // namespace walberla
