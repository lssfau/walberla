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
//! \file ExtrapolationDirectionFinder.h
//! \ingroup lbm_mesapd_coupling
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//
//======================================================================================================================


#pragma once

#include "domain_decomposition/StructuredBlockStorage.h"

#include "field/Field.h"

#include "mesa_pd/data/IAccessor.h"

#include "stencil/Directions.h"
#include "stencil/D3Q6.h"
#include "stencil/D3Q27.h"



namespace walberla {
namespace lbm_mesapd_coupling {


//**************************************************************************************************************************************
/*!
 *   \brief Classes to be used with some Reconstructor variants to find an extrapolation direction
 *
 *   The goal is to find a suitable extrapolation direction which is a three dimensional vector with cell index offsets to the current cell.
 *   Each ExtrapolationDirectionFinder class must exactly implement the member function
 *   template< typename ParticleAccessor_T >
 *   Vector3<cell_idx_t> operator()( IBlock * const block, const cell_idx_t & x, const cell_idx_t & y, const cell_idx_t & z,
 *                                   const size_t particleIdx, const ParticleAccessor_T & ac) const
 *   that returns an appropriate extrapolation direction (i.e. a certain lattice direction).
 *
 *   Different variants are available:
 *
 *   - FlagFieldNormalExtrapolationDirectionFinder:
 *     Uses information from the direct (N,S,E,W,T,B) neighbors' flags to find the extrapolation direction
 *     by checking for the domain (here: fluid) flag.
 *     This can lead to incorrect directions in close packing scenarios or in the vicinity of other boundaries.
 *
 *   - SphereNormalExtrapolationDirectionFinder:
 *     Calculates the local normal of the particle which is in its current form only correct for spherical shapes.
 *     The corresponding direction with cell offsets is then determined based on this normal.
 *
 */
//**************************************************************************************************************************************

// finds lattice direction in given stencil that corresponds best to the provided direction
template < typename Stencil_T >
Vector3<cell_idx_t> findCorrespondingLatticeDirection( const Vector3<real_t> & direction )
{
   stencil::Direction correspondingDirection = stencil::C;
   real_t innerProduct = real_t(0);
   for( auto d = Stencil_T::beginNoCenter(); d != Stencil_T::end(); ++d )
   {
      // compute inner product <dir,c_i>
      real_t temporaryInnerProduct = direction[0] * stencil::cNorm[0][*d] + direction[1] * stencil::cNorm[1][*d] + direction[2] * stencil::cNorm[2][*d];
      if( temporaryInnerProduct > innerProduct )
      {
         // largest inner product is the most fitting lattice direction
         innerProduct = temporaryInnerProduct;
         correspondingDirection = *d;
      }
   }
   return Vector3<cell_idx_t>(stencil::cx[correspondingDirection],
                              stencil::cy[correspondingDirection],
                              stencil::cz[correspondingDirection]);
}


template< typename BoundaryHandling_T >
class FlagFieldNormalExtrapolationDirectionFinder
{
public:

   FlagFieldNormalExtrapolationDirectionFinder( const shared_ptr<StructuredBlockStorage> & blockStorage, const BlockDataID & boundaryHandlingID )
   : blockStorage_( blockStorage ), boundaryHandlingID_( boundaryHandlingID )
   {}

   template< typename ParticleAccessor_T >
   Vector3<cell_idx_t> operator()( IBlock * const block, const cell_idx_t & x, const cell_idx_t & y, const cell_idx_t & z,
                                   const size_t /*particleIdx*/, const ParticleAccessor_T & /*ac*/) const
   {
      BoundaryHandling_T * boundaryHandling = block->getData< BoundaryHandling_T >( boundaryHandlingID_ );

      WALBERLA_ASSERT_NOT_NULLPTR( boundaryHandling );

      Vector3<cell_idx_t> extrapolationDirection;

      for( auto directNeighborDir = stencil::D3Q6::begin(); directNeighborDir != stencil::D3Q6::end(); ++directNeighborDir)
      {
         if( boundaryHandling->isDomain( Cell( x + directNeighborDir.cx(), y + directNeighborDir.cy(), z + directNeighborDir.cz() ) ) )
         {
            extrapolationDirection[0] += directNeighborDir.cx();
            extrapolationDirection[1] += directNeighborDir.cy();
            extrapolationDirection[2] += directNeighborDir.cz();
         }
      }

      return extrapolationDirection;
   }

private:
   shared_ptr<StructuredBlockStorage> blockStorage_;
   const BlockDataID boundaryHandlingID_;
};


class SphereNormalExtrapolationDirectionFinder
{
public:

   explicit SphereNormalExtrapolationDirectionFinder( const shared_ptr<StructuredBlockStorage> & blockStorage )
   : blockStorage_( blockStorage )
   {}

   template< typename ParticleAccessor_T >
   Vector3<cell_idx_t> operator()( IBlock * const block, const cell_idx_t & x, const cell_idx_t & y, const cell_idx_t & z,
                                   const size_t particleIdx, const ParticleAccessor_T & ac) const
   {
      static_assert(std::is_base_of<mesa_pd::data::IAccessor, ParticleAccessor_T>::value, "Provide a valid accessor as template");

      real_t cx, cy, cz;
      blockStorage_->getBlockLocalCellCenter( *block, Cell(x,y,z), cx, cy, cz );

      Vector3<real_t> particleCenterPosition = ac.getPosition(particleIdx);
      WALBERLA_ASSERT( !math::isnan(particleCenterPosition) );

      Vector3<real_t> particleNormal( cx - particleCenterPosition[0], cy - particleCenterPosition[1], cz - particleCenterPosition[2] );

      return findCorrespondingLatticeDirection< stencil::D3Q27 >( particleNormal );
   }

private:
   shared_ptr<StructuredBlockStorage> blockStorage_;

};

} // namespace lbm_mesapd_coupling
} // namespace walberla
