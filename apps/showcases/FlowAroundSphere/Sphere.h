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
//! \file Sphere.h
//! \author Markus Holzer <markus.holzer@fau.de>
//
//======================================================================================================================
#pragma once

#include "blockforest/SetupBlock.h"
#include "blockforest/SetupBlockForest.h"
#include "blockforest/StructuredBlockForest.h"

#include "domain_decomposition/IBlock.h"

#include "field/FlagUID.h"

#include "core/DataTypes.h"
#include "core/math/AABB.h"
#include "core/math/Vector3.h"
#include "core/cell/Cell.h"

#include "stencil/D3Q7.h"
#include "stencil/D3Q27.h"

#include "Types.h"
#include "Setup.h"

namespace walberla
{

class Sphere
{
 public:
   Sphere(const Setup& setup) : setup_(setup)
   {
      const real_t px = setup_.sphereXPosition * setup_.dxC;
      const real_t py = setup_.sphereYPosition * setup_.dxC;
      const real_t pz = setup_.sphereZPosition * setup_.dxC;

      center_  = Vector3< real_t >(px, py, pz);
      radius_  = setup_.sphereRadius * setup_.dxC;
      radius2_ = radius_ * radius_;
   }

   bool operator()(const Vector3< real_t >& point) const { return contains(point); }

   bool contains(const Vector3< real_t >& point) const;
   bool contains(const AABB& aabb) const;

   real_t delta(const Vector3< real_t >& fluid, const Vector3< real_t >& boundary) const;
   Setup getSetup(){ return setup_; }
   void setupBoundary(const std::shared_ptr< StructuredBlockForest >& sbfs, const BlockDataID flagFieldID);
   void checkConsistency(const std::shared_ptr< StructuredBlockForest >& sbfs, const BlockDataID flagFieldID);
   void setupSphereBoundary(const std::shared_ptr< StructuredBlockForest >& sbfs, const BlockDataID flagFieldID);
   void setBoundaryFromCellInterval(CellInterval& cells, const FlagField_T::flag_t flag, FlagField_T* flagField);

 private:
   Setup setup_;
   Vector3< real_t > center_;
   real_t radius_;
   real_t radius2_;

}; // class Sphere

class SphereRefinementSelection
{
 public:
   SphereRefinementSelection(const Sphere& sphere, const uint_t level)
      : sphere_(sphere), level_(level)
   {
      auto setup = sphere_.getSetup();
      const real_t px = setup.sphereXPosition * setup.dxC;
      const real_t py = setup.sphereYPosition * setup.dxC;
      const real_t pz = setup.sphereZPosition * setup.dxC;

      center_  = Vector3< real_t >(px, py, pz);
      const real_t sphereRadius = setup.sphereRadius * setup.dxC;
      const real_t bufferDistance = setup.dxF;
      const real_t d = sphereRadius + bufferDistance;
      radius2_ = d * d;

      sphereBoundingBox1_ = AABB(center_[0], center_[1] - d, center_[2] - d,
                                center_[0] + (real_c(2.5) * sphereRadius), center_[1] + d, center_[2] + d);

      sphereBoundingBox2_ = AABB(center_[0], center_[1] - d, center_[2] - d,
                                 center_[0] + (real_c(5) * sphereRadius), center_[1] + d, center_[2] + d);
   }

   bool contains(const Vector3< real_t >& point) const
   {
      return (point - center_).sqrLength() <= radius2_;
   }

   bool intersects(const AABB& aabb) const
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
      return contains(p[0]) || contains(p[1]) || contains(p[2]) || contains(p[3]) || contains(p[4]) || contains(p[5]) ||
             contains(p[6]) || contains(p[7]);
   }

   void operator()(SetupBlockForest& forest)
   {
      if(level_ == 0)
         return;
      for (auto block = forest.begin(); block != forest.end(); ++block)
      {
         const AABB& aabb = block->getAABB();

         if (block->getLevel() < level_ && (intersects(aabb) || sphereBoundingBox1_.intersects(aabb)) )
            block->setMarker(true);

         if (block->getLevel() < (level_ - 1) && (intersects(aabb) || sphereBoundingBox2_.intersects(aabb)) )
            block->setMarker(true);
      }
   }

 private:
   Sphere sphere_;
   Vector3< real_t > center_;
   uint_t level_;
   AABB sphereBoundingBox1_;
   AABB sphereBoundingBox2_;
   real_t radius2_;

}; // class SphereRefinementSelection

class SphereBlockExclusion
{
 public:
   SphereBlockExclusion(const Sphere& sphere) : sphere_(sphere) {}

   bool operator()(const blockforest::SetupBlock& block)
   {
      const AABB aabb = block.getAABB();
      return static_cast< bool >(sphere_.contains(aabb));
   }

 private:
   Sphere sphere_;

}; // class SphereBlockExclusion


class wallDistance
{
 public:
   wallDistance(const Sphere& sphere) : sphere_(sphere) {}

   real_t operator()(const Cell& fluidCell, const Cell& boundaryCell, const shared_ptr< StructuredBlockForest >& SbF,
                     IBlock& block) const;

 private:
   Sphere sphere_;
}; // class wallDistance

} // namespace walberla