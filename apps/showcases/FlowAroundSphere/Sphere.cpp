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
//! \file Sphere.cpp
//! \author Markus Holzer <markus.holzer@fau.de>
//
//======================================================================================================================

#include "Sphere.h"

namespace walberla
{

bool Sphere::contains(const Vector3< real_t >& point) const
{
   return (point - center_).sqrLength() <= radius2_;
}

bool Sphere::contains(const AABB& aabb) const
{
   Vector3< real_t > p[8];
   p[0].set(aabb.xMin(), aabb.yMin(), aabb.zMin());
   p[1].set(aabb.xMax(), aabb.yMin(), aabb.zMin());
   p[2].set(aabb.xMin(), aabb.yMax(), aabb.zMin());
   p[3].set(aabb.xMax(), aabb.yMax(), aabb.zMin());
   p[4].set(aabb.xMin(), aabb.yMin(), aabb.zMax());
   p[5].set(aabb.xMax(), aabb.yMin(), aabb.zMax());
   p[6].set(aabb.xMin(), aabb.yMax(), aabb.zMax());
   p[7].set(aabb.xMax(), aabb.yMax(), aabb.zMax());
   return contains(p[0]) && contains(p[1]) && contains(p[2]) && contains(p[3]) && contains(p[4]) && contains(p[5]) &&
          contains(p[6]) && contains(p[7]);
}

real_t Sphere::delta(const Vector3< real_t >& fluid, const Vector3< real_t >& boundary) const
{
   WALBERLA_ASSERT(!contains(fluid))
   WALBERLA_ASSERT(contains(boundary))

   // http://devmag.org.za/2009/04/17/basic-collision-detection-in-2d-part-2/
   const Vector3< real_t > f = fluid - center_;
   const Vector3< real_t > d = (boundary - center_) - f;

   const real_t dDotd = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
   const real_t fDotf = f[0] * f[0] + f[1] * f[1] + f[2] * f[2];
   const real_t dDotf = d[0] * f[0] + d[1] * f[1] + d[2] * f[2];

   const real_t b = real_c(2.0) * dDotf;
   const real_t c = fDotf - radius2_;

   const real_t bb4ac = b * b - (real_c(4.0) * dDotd * c);
   WALBERLA_CHECK_GREATER_EQUAL(bb4ac, real_c(0.0))

   const real_t sqrtbb4ac = std::sqrt(bb4ac);
   const real_t alpha = std::min((-b + sqrtbb4ac) / (real_c(2.0) * dDotd), (-b - sqrtbb4ac) / (real_c(2.0) * dDotd));

   WALBERLA_CHECK_GREATER_EQUAL(alpha, real_c(0.0))
   WALBERLA_CHECK_LESS_EQUAL(alpha, real_c(1.0))

   return alpha;
}

void Sphere::setupBoundary(const std::shared_ptr< StructuredBlockForest >& sbfs, const BlockDataID flagFieldID)
{

   for (auto bIt = sbfs->begin(); bIt != sbfs->end(); ++bIt)
   {
      auto flagField       = bIt->getData< FlagField_T >(flagFieldID);
      const FlagField_T::flag_t inflowFlag  = flagField->registerFlag(setup_.inflowUID);
      const FlagField_T::flag_t outflowFlag = flagField->registerFlag(setup_.outflowUID);
      const FlagField_T::flag_t wallFlag    = flagField->registerFlag(setup_.wallUID);

      const cell_idx_t gls = cell_idx_c(flagField->nrOfGhostLayers()) - cell_idx_c(1);

      CellInterval blockBB(-1, -1, -1,
                           cell_idx_c(setup_.cellsPerBlock[0]), cell_idx_c(setup_.cellsPerBlock[1]), cell_idx_c(setup_.cellsPerBlock[2]));

      // inflow WEST
      if(sbfs->atDomainXMinBorder(*bIt)){
         CellInterval west(blockBB.xMin() - gls, blockBB.yMin() - gls, blockBB.zMin() - gls, blockBB.xMin(),
                           blockBB.yMax() + gls, blockBB.zMax() + gls);
         setBoundaryFromCellInterval(west, inflowFlag, flagField);
      }

      // outflow EAST
      if(sbfs->atDomainXMaxBorder(*bIt)){
         CellInterval east(blockBB.xMax(), blockBB.yMin() - gls, blockBB.zMin() - gls, blockBB.xMax() + gls,
                           blockBB.yMax() + gls, blockBB.zMax() + gls);
         setBoundaryFromCellInterval(east, outflowFlag, flagField);
      }

      // SOUTH
      if(sbfs->atDomainYMinBorder(*bIt))
      {
         CellInterval south(blockBB.xMin() - gls, blockBB.yMin() - gls, blockBB.zMin() - gls, blockBB.xMax() + gls,
                            blockBB.yMin(), blockBB.zMax() + gls);
         setBoundaryFromCellInterval(south, wallFlag, flagField);
      }

      // NORTH
      if(sbfs->atDomainYMaxBorder(*bIt)){
         CellInterval north( blockBB.xMin() - gls, blockBB.yMax(), blockBB.zMin() - gls, blockBB.xMax() + gls,
                             blockBB.yMax() + gls, blockBB.zMax() + gls );
         setBoundaryFromCellInterval(north, wallFlag, flagField);
      }

      // BOTTOM
      if(sbfs->atDomainZMinBorder(*bIt)){
         CellInterval bottom(blockBB.xMin() - gls, blockBB.yMin() - gls, blockBB.zMin() - gls, blockBB.xMax() + gls,
                             blockBB.yMax() + gls, blockBB.zMin());
         setBoundaryFromCellInterval(bottom, wallFlag, flagField);
      }

      // TOP
      if(sbfs->atDomainZMaxBorder(*bIt)){
         CellInterval top(blockBB.xMin() - gls, blockBB.yMin() - gls, blockBB.zMax(), blockBB.xMax() + gls,
                          blockBB.yMax() + gls, blockBB.zMax() + gls);
         setBoundaryFromCellInterval(top, wallFlag, flagField);
      }
   }

   checkConsistency(sbfs, flagFieldID);
   setupSphereBoundary(sbfs, flagFieldID);
}

void Sphere::setBoundaryFromCellInterval(CellInterval& cells, const FlagField_T::flag_t flag, FlagField_T* flagField)
{

   for (auto cell = cells.begin(); cell != cells.end(); ++cell){
      if(flagField->get(cell->x(), cell->y(), cell->z()) == FlagField_T::flag_t(0))
         flagField->addFlag(cell->x(), cell->y(), cell->z(), flag);
   }
}

void Sphere::checkConsistency(const std::shared_ptr< StructuredBlockForest >& sbfs, const BlockDataID flagFieldID){
   const uint_t depth = sbfs->getDepth();
   for (auto bIt = sbfs->begin(); bIt != sbfs->end(); ++bIt)
   {
      Block& b             = dynamic_cast< Block& >(*bIt);
      auto flagField       = b.getData< FlagField_T >(flagFieldID);
      if (sbfs->getLevel(b) < depth){
         for( auto it = flagField->beginWithGhostLayer(1); it != flagField->end(); ++it ){
            Vector3< real_t > cellCenter = sbfs->getBlockLocalCellCenter( b, it.cell() );
            sbfs->mapToPeriodicDomain(cellCenter);
            WALBERLA_CHECK(!contains(cellCenter), "The sphere must be completely located on the finest level")
         }
      }
   }
}

void Sphere::setupSphereBoundary(const std::shared_ptr< StructuredBlockForest >& sbfs, const BlockDataID flagFieldID){
   const uint_t depth = sbfs->getDepth();
   for (auto bIt = sbfs->begin(); bIt != sbfs->end(); ++bIt)
   {
      Block& b             = dynamic_cast< Block& >(*bIt);
      auto flagField       = b.getData< FlagField_T >(flagFieldID);
      uint8_t obstacleFlag = flagField->registerFlag(setup_.obstacleUID);

      if (sbfs->getLevel(b) == depth){
         for( auto it = flagField->beginWithGhostLayer(1); it != flagField->end(); ++it ){
            Vector3< real_t > cellCenter = sbfs->getBlockLocalCellCenter( b, it.cell() );
            sbfs->mapToPeriodicDomain(cellCenter);
            if (contains(cellCenter)) { flagField->addFlag(it.x(), it.y(), it.z(), obstacleFlag); }
         }
      }
   }
}

real_t wallDistance::operator()(const Cell& fluidCell, const Cell& boundaryCell,
                                const shared_ptr< StructuredBlockForest >& SbF, IBlock& block) const
{
   Vector3< real_t > boundary = SbF->getBlockLocalCellCenter( block, boundaryCell );
   Vector3< real_t > fluid = SbF->getBlockLocalCellCenter( block, fluidCell );
   SbF->mapToPeriodicDomain(boundary);
   SbF->mapToPeriodicDomain(fluid);

   WALBERLA_CHECK(!sphere_.contains(fluid), "fluid cell is in contained in sphere (" << fluid[0] << ", " << fluid[1] << ", " << fluid[2] << "). The block local cell is " << fluidCell)
   WALBERLA_CHECK(sphere_.contains(boundary), "boundary cell is not in contained in sphere (" << boundary[0] << ", " << boundary[1] << ", " << boundary[2] << "). The block local cell is " << boundaryCell)

   return sphere_.delta( fluid, boundary );
}
} // namespace walberla