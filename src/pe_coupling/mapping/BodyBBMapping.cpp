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
//! \file BodyBBMapping.cpp
//! \ingroup pe_coupling
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include "BodyBBMapping.h"

#include "pe/rigidbody/RigidBody.h"

namespace walberla {
namespace pe_coupling {

CellInterval getCellBB( const pe::ConstBodyID body, const IBlock & block, StructuredBlockStorage & blockStorage,
                        const uint_t numberOfGhostLayersToInclude )
{

   WALBERLA_ASSERT_NOT_NULLPTR( body );

   CellInterval cellBB;

   if( body->isFinite() )
   {
      blockStorage.getCellBBFromAABB( cellBB, body->getAABB(), blockStorage.getLevel(block) );
   }
   else
   {
      // if body is infinite (global), its AABB is also infinite which then requires special treatment

      auto level = blockStorage.getLevel(block);
      const real_t dx = blockStorage.dx(level);
      const real_t dy = blockStorage.dy(level);
      const real_t dz = blockStorage.dz(level);
      Vector3<real_t> aabbExtensionByGhostLayers(real_c(numberOfGhostLayersToInclude) * dx,
                                                 real_c(numberOfGhostLayersToInclude) * dy,
                                                 real_c(numberOfGhostLayersToInclude) * dz);
      auto extendedBlockAABB = blockStorage.getAABB(block.getId()).getExtended( aabbExtensionByGhostLayers );

      // intersect the infinite (global) body with the block AABB, extended by its ghost layers
      // then determine the cell bounding box of the intersection
      blockStorage.getCellBBFromAABB( cellBB, body->getAABB().getIntersection( extendedBlockAABB ), level );

      // if infinite body does not intersect with the extended block AABB, return an empty interval
      if( cellBB.empty() ) return CellInterval();
   }

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

} // namespace pe_coupling
} // namespace walberla
