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
//! \file VelocityCurlFieldEvaluator.h
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

/*!\brief Evaluator of the curl of a given velocity field.
 *
 */
template <typename LatticeModel_T, typename BoundaryHandling_T>
class VelocityCurlFieldEvaluator
{
public:
   using VectorField_T = GhostLayerField<Vector3<real_t>, 1>;
   using Stencil_T = typename LatticeModel_T::Stencil;

   VelocityCurlFieldEvaluator( const BlockDataID & velocityCurlFieldID, const ConstBlockDataID & velocityFieldID,
                               const ConstBlockDataID & boundaryHandlingID )
      : velocityCurlFieldID_( velocityCurlFieldID ), velocityFieldID_( velocityFieldID ),
        boundaryHandlingID_( boundaryHandlingID )
   { }

   void operator()(IBlock * const block)
   {
      const VectorField_T* velocityField            = block->getData< VectorField_T >( velocityFieldID_ );
      VectorField_T* velocityCurlField        = block->getData< VectorField_T >( velocityCurlFieldID_ );
      const BoundaryHandling_T * boundaryHandling   = block->getData< BoundaryHandling_T >( boundaryHandlingID_ );

      WALBERLA_FOR_ALL_CELLS_XYZ( velocityCurlField,
         if( boundaryHandling->isDomain(x,y,z) )
         {
            velocityCurlField->get(x,y,z) = getVelocityCurl( Cell(x,y,z), velocityField, boundaryHandling );
         }
      );
   }

private:

   Vector3<real_t> getVelocityCurl( const Cell & cell, const VectorField_T * velocityField, const BoundaryHandling_T * boundaryHandling )
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

      // apply curl formula
      // See: Ramadugu et al - Lattice differential operators for computational physics (2013)
      // with T = c_s**2
      const real_t inv_c_s_sqr = real_t(3);
      Vector3<real_t> curl( real_t(0) );
      for( auto dir = Stencil_T::beginNoCenter(); dir != Stencil_T::end(); ++dir)
      {
         Vector3<real_t> latticeVel( real_c(dir.cx()), real_c(dir.cy()), real_c(dir.cz()) );
         curl += LatticeModel_T::w[ dir.toIdx() ] * ( latticeVel % velocityValues[ *dir ] );
      }
      curl *= inv_c_s_sqr;

      return curl;
   }

   const BlockDataID velocityCurlFieldID_;
   const ConstBlockDataID velocityFieldID_;
   const ConstBlockDataID boundaryHandlingID_;
};


} // namespace discrete_particle_methods
} // namespace pe_coupling
} // namespace walberla
