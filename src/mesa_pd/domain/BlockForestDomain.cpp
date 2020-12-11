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
//! \file BlockForestDomain.cpp
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include "BlockForestDomain.h"

#include <core/debug/CheckFunctions.h>
#include <core/math/AABB.h>
#include <core/mpi/MPIManager.h>

namespace walberla {
namespace mesa_pd {
namespace domain {


/// \post neighborSubdomains_ is sorted by rank
BlockForestDomain::BlockForestDomain(const std::shared_ptr<blockforest::BlockForest>& blockForest)
   : blockForest_(blockForest)
{
   refresh();
}

/// \post neighborSubdomains_ is sorted by rank
void BlockForestDomain::refresh()
{
   ownRank_ = mpi::MPIManager::instance()->rank();

   periodic_[0] = blockForest_->isPeriodic(0);
   periodic_[1] = blockForest_->isPeriodic(1);
   periodic_[2] = blockForest_->isPeriodic(2);

   localAABBs_.clear();
   neighborSubdomains_.clear();
   neighborProcesses_.clear();
   unionOfLocalAABBs_ = math::AABB(Vec3(real_t(0)), Vec3(real_t(0)));

   if (blockForest_->empty()) return;

   unionOfLocalAABBs_ = blockForest_->begin()->getAABB();
   for (auto& iBlk : *blockForest_)
   {
      const Block& blk = *static_cast<blockforest::Block*>(&iBlk);
      localAABBs_.push_back(blk.getAABB());
      unionOfLocalAABBs_.merge(blk.getAABB());
      for (uint_t nb = 0; nb < blk.getNeighborhoodSize(); ++nb)
      {
         if (int_c(blk.getNeighborProcess(nb)) == ownRank_) continue;

         //check if neighbor aabb is already present
         const BlockID& nbBlkId = blk.getNeighborId(nb);
         if (std::find_if(neighborSubdomains_.begin(),
                          neighborSubdomains_.end(),
                          [&nbBlkId](const auto& subdomain){return subdomain.blockID == nbBlkId;}) ==
             neighborSubdomains_.end())
         {
            neighborSubdomains_.emplace_back(int_c(blk.getNeighborProcess(nb)), nbBlkId, blk.getNeighborAABB(nb));
         }
      }
   }

   //sort by rank
   std::sort(neighborSubdomains_.begin(),
             neighborSubdomains_.end(),
             [](const auto& lhs, const auto& rhs){ return lhs.rank < rhs.rank; });

   //generate list of neighbor processes
   int prevRank = -1;
   for (auto& subdomain : neighborSubdomains_)
   {
      if ((prevRank != subdomain.rank) && (subdomain.rank != ownRank_))
      {
         neighborProcesses_.emplace_back(uint_c(subdomain.rank));
         prevRank = subdomain.rank;
      }
   }
}

bool BlockForestDomain::isContainedInProcessSubdomain(const uint_t rank, const Vec3& pt) const
{
   if (blockForest_->empty()) return false;

   if (uint_c(ownRank_) == rank)
   {
      // check if point is in local subdomain by checking all aabbs
      for (auto& aabb : localAABBs_)
      {
         if (aabb.contains(pt)) return true;
      }
   } else
   {
      WALBERLA_ASSERT(std::is_sorted(neighborSubdomains_.begin(),
                                     neighborSubdomains_.end(),
                                     [](const auto& lhs, const auto& rhs){ return lhs.rank < rhs.rank;}));

      auto begin = std::find_if(neighborSubdomains_.begin(),
                                neighborSubdomains_.end(),
                                [&rank](const auto& subdomain){ return subdomain.rank == int_c(rank); });
      if (begin == neighborSubdomains_.end()) return false; //no information for rank available
      auto end = std::find_if(begin,
                              neighborSubdomains_.end(),
                              [&rank](const auto& subdomain){ return subdomain.rank != int_c(rank); });

      for (auto it = begin; it != end; ++it)
      {
         if (it->aabb.contains(pt)) return true;
      }
   }
   return false;
}

bool   BlockForestDomain::isContainedInLocalSubdomain(const Vec3& pt,
                                                      const real_t& radius) const
{
   return std::any_of(localAABBs_.begin(),
                      localAABBs_.end(),
                      [&](auto& aabb)
                      {return isInsideAABB(pt, radius, aabb);});
}

bool BlockForestDomain::isContainedInProcessSubdomain(const Vec3& pt, const real_t& radius) const
{
   if (blockForest_->empty()) return false;

   //completely contained in local aabb?
   for (auto& aabb : localAABBs_)
   {
      if (aabb.contains(pt))
      {
         if (isInsideAABB(pt, radius, aabb))
         {
            return true;
         }

         break;
      }
   }

   //intersects one of the neighboring subdomains?
   return std::none_of(neighborSubdomains_.begin(),
                       neighborSubdomains_.end(),
                       [&](const auto &subdomain) {
                          return sqDistancePointToAABB(pt, subdomain.aabb) < radius * radius;
                       });
}

int BlockForestDomain::findContainingProcessRank(const Vec3& pt) const
{
   if (blockForest_->empty()) return -1;

   if (isContainedInProcessSubdomain(uint_c(ownRank_), pt)) return ownRank_;
   for( uint_t rank : getNeighborProcesses() )
   {
      if (isContainedInProcessSubdomain(rank, pt)) return int_c(rank);
   }
   return -1;
}
void BlockForestDomain::periodicallyMapToDomain(Vec3& pt) const
{
   blockForest_->mapToPeriodicDomain(pt);
}

std::vector<uint_t> BlockForestDomain::getNeighborProcesses() const
{
   return neighborProcesses_;
}

bool BlockForestDomain::intersectsWithProcessSubdomain(const uint_t rank, const Vec3& pt, const real_t& radius) const
{
   if (blockForest_->empty()) return false;

   if (uint_c(ownRank_) == rank)
   {
      //=====================
      // LOCAL DOMAIN
      if (isInsideGlobalDomain(pt, radius))
      {
         for (auto& aabb : localAABBs_)
         {
            if (sqDistancePointToAABB(pt, aabb) <= radius * radius) return true;
         }
      } else
      {
         for (auto& aabb : localAABBs_)
         {
            if (sqDistancePointToAABBPeriodic(pt, aabb, blockForest_->getDomain(), periodic_) <= radius * radius) return true;
         }
      }
   } else
   {
      //=====================
      // NEIGHBORING DOMAIN
      WALBERLA_ASSERT(std::is_sorted(neighborSubdomains_.begin(),
                                     neighborSubdomains_.end(),
                                     [](const auto& lhs, const auto& rhs){ return lhs.rank < rhs.rank;}));

      if (isInsideGlobalDomain(pt, radius))
      {
         size_t idx = 0;
         WALBERLA_ASSERT_LESS(idx, neighborSubdomains_.size());
         while (neighborSubdomains_[idx].rank != int_c(rank))
         {
            ++idx;
            WALBERLA_ASSERT_LESS(idx, neighborSubdomains_.size());
         }
         while (neighborSubdomains_[idx].rank == int_c(rank))
         {
            if (sqDistancePointToAABB(pt, neighborSubdomains_[idx].aabb) <= radius * radius) return true;
            ++idx;
            if (idx >= neighborSubdomains_.size()) break;
            WALBERLA_ASSERT_LESS(idx, neighborSubdomains_.size());
         }
      } else
      {
         size_t idx = 0;
         WALBERLA_ASSERT_LESS(idx, neighborSubdomains_.size());
         while (neighborSubdomains_[idx].rank != int_c(rank))
         {
            ++idx;
            WALBERLA_ASSERT_LESS(idx, neighborSubdomains_.size());
         }
         while (neighborSubdomains_[idx].rank == int_c(rank))
         {
            if (sqDistancePointToAABBPeriodic(pt, neighborSubdomains_[idx].aabb, blockForest_->getDomain(), periodic_) <= radius * radius) return true;
            ++idx;
            if (idx >= neighborSubdomains_.size()) break;
            WALBERLA_ASSERT_LESS(idx, neighborSubdomains_.size());
         }
      }
   }

   return false;
}

void BlockForestDomain::correctParticlePosition(Vec3& pt) const
{
   const Vec3 center = unionOfLocalAABBs_.center();
   const Vec3 dis = pt - center;

   const auto& domain = blockForest_->getDomain();

   if (periodic_[0] && (-domain.xSize() * 0.5 > dis[0])) pt[0] += domain.xSize();
   if (periodic_[0] && (+domain.xSize() * 0.5 < dis[0])) pt[0] -= domain.xSize();

   if (periodic_[1] && (-domain.ySize() * 0.5 > dis[1])) pt[1] += domain.ySize();
   if (periodic_[1] && (+domain.ySize() * 0.5 < dis[1])) pt[1] -= domain.ySize();

   if (periodic_[2] && (-domain.zSize() * 0.5 > dis[2])) pt[2] += domain.zSize();
   if (periodic_[2] && (+domain.zSize() * 0.5 < dis[2])) pt[2] -= domain.zSize();
}


bool BlockForestDomain::isInsideGlobalDomain(const Vec3& pt, const real_t& radius) const
{
   return isInsideAABB(pt, radius, blockForest_->getDomain());
}

} //namespace domain
} //namespace mesa_pd
} //namespace walberla
