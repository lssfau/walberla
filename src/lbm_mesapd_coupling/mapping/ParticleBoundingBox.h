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
//! \file ParticleBoundingBox.h
//! \ingroup lbm_mesapd_coupling
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/debug/Debug.h"

#include "domain_decomposition/StructuredBlockStorage.h"

#include "mesa_pd/common/AABBConversion.h"
#include "mesa_pd/data/Flags.h"
#include "mesa_pd/data/IAccessor.h"

namespace walberla {
namespace lbm_mesapd_coupling {

/*
 * Obtain a block-local cell bounding box from a given AABB (e.g. the particle's AABB)
 * If the given AABB is (partly) infinite, AABBIsInfinite should be set to true (e.g. for infinite particles)
 */
CellInterval getCellBBFromAABB( const math::AABB & aabb, bool AABBIsInfinite,
                                const IBlock & block, StructuredBlockStorage & blockStorage,
                                const uint_t numberOfGhostLayersToInclude)
{

   CellInterval cellBB;

   if( AABBIsInfinite )
   {
      auto level = blockStorage.getLevel(block);
      const real_t dx = blockStorage.dx(level);
      const real_t dy = blockStorage.dy(level);
      const real_t dz = blockStorage.dz(level);
      Vector3<real_t> aabbExtensionByGhostLayers(real_c(numberOfGhostLayersToInclude) * dx,
                                                 real_c(numberOfGhostLayersToInclude) * dy,
                                                 real_c(numberOfGhostLayersToInclude) * dz);
      auto extendedBlockAABB = blockStorage.getAABB(block.getId()).getExtended( aabbExtensionByGhostLayers );

      // intersect the (partly) infinite aabb with the block AABB, extended by its ghost layers
      // then determine the cell bounding box of the intersection
      blockStorage.getCellBBFromAABB( cellBB, aabb.getIntersection( extendedBlockAABB ), level );

      // if (partly) infinite aabb does not intersect with the extended block AABB, return an empty interval
      if( cellBB.empty() ) return CellInterval();
   }
   else
   {
      blockStorage.getCellBBFromAABB( cellBB, aabb, blockStorage.getLevel(block) );
   }

   WALBERLA_ASSERT( !cellBB.empty() );

   cellBB.xMin() -= cell_idx_t(1); cellBB.yMin() -= cell_idx_t(1); cellBB.zMin() -= cell_idx_t(1);
   cellBB.xMax() += cell_idx_t(1); cellBB.yMax() += cell_idx_t(1); cellBB.zMax() += cell_idx_t(1);

   CellInterval blockBB = blockStorage.getBlockCellBB( block );

   cell_idx_t layers = cell_idx_c( numberOfGhostLayersToInclude );

   blockBB.xMin() -= layers; blockBB.yMin() -= layers; blockBB.zMin() -= layers;
   blockBB.xMax() += layers; blockBB.yMax() += layers; blockBB.zMax() += layers;

   cellBB.intersect( blockBB );

   blockStorage.transformGlobalToBlockLocalCellInterval( cellBB, block );

   return cellBB;

}

template< typename ParticleAccessor_T >
CellInterval getParticleCellBB(const size_t particleIdx, const ParticleAccessor_T& ac,
                               const IBlock & block, StructuredBlockStorage & blockStorage,
                               const uint_t numberOfGhostLayersToInclude)
{
   return getCellBBFromAABB(mesa_pd::getParticleAABB(particleIdx, ac),
                            mesa_pd::data::particle_flags::isSet( ac.getFlags(particleIdx), mesa_pd::data::particle_flags::INFINITE),
                            block, blockStorage, numberOfGhostLayersToInclude );
}


} // namespace lbm_mesapd_coupling
} // namespace walberla
