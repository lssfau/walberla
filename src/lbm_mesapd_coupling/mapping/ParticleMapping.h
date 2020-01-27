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
//! \file ParticleMapping.h
//! \ingroup lbm_mesapd_coupling
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//
//======================================================================================================================

#pragma once

#include "ParticleBoundingBox.h"

#include "core/debug/Debug.h"
#include "core/math/Vector3.h"
#include "domain_decomposition/StructuredBlockStorage.h"
#include "field/FlagField.h"

#include "mesa_pd/common/Contains.h"
#include "mesa_pd/data/IAccessor.h"
#include "mesa_pd/kernel/SingleCast.h"

#include <functional>

namespace walberla {
namespace lbm_mesapd_coupling {

/*
 * Kernel that can be used to map MESA_PD particles into the domain with a given LBM boundary flag.
 * In coupled simulations, usually carried out for walls or fixed particles.
 *
 */
template< typename BoundaryHandling_T>
class ParticleMappingKernel
{
public:
   ParticleMappingKernel(const shared_ptr<StructuredBlockStorage> & blockStorage, const BlockDataID & boundaryHandlingID)
   : blockStorage_(blockStorage), boundaryHandlingID_(boundaryHandlingID){}

   template<typename ParticleAccessor_T>
   void operator()(const size_t particleIdx, const ParticleAccessor_T& ac, const FlagUID & obstacle )
   {
      static_assert(std::is_base_of<mesa_pd::data::IAccessor, ParticleAccessor_T>::value, "Provide a valid accessor as template");

      for( auto blockIt = blockStorage_->begin(); blockIt != blockStorage_->end(); ++blockIt )
      {
         mapParticleOnBlock( particleIdx, ac, *blockIt, obstacle );
      }
   }

private:

   template< typename ParticleAccessor_T >
   void mapParticleOnBlock( const size_t particleIdx, const ParticleAccessor_T& ac,
                           IBlock & block, const FlagUID & obstacle )
   {
      static_assert(std::is_base_of<mesa_pd::data::IAccessor, ParticleAccessor_T>::value, "Provide a valid accessor as template");

      WALBERLA_ASSERT_EQUAL( &block.getBlockStorage(), &(blockStorage_->getBlockStorage()) );

      BoundaryHandling_T * boundaryHandling = block.getData< BoundaryHandling_T >( boundaryHandlingID_);
      WALBERLA_ASSERT_NOT_NULLPTR( boundaryHandling );

      auto * flagField = boundaryHandling->getFlagField();
      WALBERLA_ASSERT_NOT_NULLPTR( flagField );
      WALBERLA_ASSERT( flagField->flagExists( obstacle ) );

      const auto obstacleFlag = flagField->getFlag( obstacle );

      auto cellBB = getParticleCellBB(particleIdx, ac, block, *blockStorage_, flagField->nrOfGhostLayers() );

      if( cellBB.empty() ) return;

      mesa_pd::kernel::SingleCast singleCast;
      mesa_pd::ContainsPointFunctor containsPointFctr;

      Vector3<real_t> startCellCenter = blockStorage_->getBlockLocalCellCenter( block, cellBB.min() );
      auto blockLevel = blockStorage_->getLevel(block);
      const real_t dx = blockStorage_->dx( blockLevel );
      const real_t dy = blockStorage_->dy( blockLevel );
      const real_t dz = blockStorage_->dz( blockLevel );

      real_t cz = startCellCenter[2];
      for( cell_idx_t z = cellBB.zMin(); z <= cellBB.zMax(); ++z )
      {
         real_t cy = startCellCenter[1];
         for( cell_idx_t y = cellBB.yMin(); y <= cellBB.yMax(); ++y )
         {
            real_t cx = startCellCenter[0];
            for( cell_idx_t x = cellBB.xMin(); x <= cellBB.xMax(); ++x )
            {
               if( singleCast(particleIdx, ac, containsPointFctr, ac, Vector3<real_t>(cx,cy,cz)) )
               {
                  boundaryHandling->forceBoundary(obstacleFlag, x, y, z);
               }
               cx += dx;
            }
            cy += dy;
         }
         cz += dz;
      }
   }

   shared_ptr<StructuredBlockStorage> blockStorage_;
   BlockDataID boundaryHandlingID_;
};


} // namespace lbm_mesapd_coupling
} // namespace walberla
