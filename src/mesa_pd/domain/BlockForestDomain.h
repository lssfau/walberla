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
//! \file
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#pragma once

#include <mesa_pd/domain/IDomain.h>

#include <blockforest/BlockForest.h>

#include <memory>

namespace walberla {
namespace mesa_pd {
namespace domain {

class BlockForestDomain : public IDomain
{
public:
   BlockForestDomain(const std::shared_ptr<blockforest::BlockForest>& blockForest);

   /**
    * @brief If the BlockForest is changed this function has to be called in order to
    * update all interal caches!
    *
    * Updates the local caches for local and neighbor AABBs.
    */
   void refresh();

   bool   isContainedInProcessSubdomain(const uint_t rank, const Vec3& pt) const override;
   bool   isContainedInLocalSubdomain(const Vec3& pt, const real_t& radius) const override;
   /// Is the sphere defined by \p pt and \p radius completely inside the local subdomin?
   /// \attention Also take into account periodicity!
   /// \param pt center of the sphere
   /// \param radius radius of the sphere
   bool   isContainedInProcessSubdomain(const Vec3& pt, const real_t& radius) const;
   int    findContainingProcessRank(const Vec3& pt) const override;
   void   periodicallyMapToDomain(Vec3& pt) const override;
   std::vector<uint_t> getNeighborProcesses() const override;
   bool   intersectsWithProcessSubdomain(const uint_t rank, const Vec3& pt, const real_t& radius) const override;
   void   correctParticlePosition(Vec3& pt) const override;

   const math::AABB& getUnionOfLocalAABBs() const {return unionOfLocalAABBs_;}
   size_t getNumLocalAABBs() const {return localAABBs_.size();}
   size_t getNumNeighborSubdomains() const {return neighborSubdomains_.size();}
   size_t getNumNeighborProcesses() const {return neighborProcesses_.size();}
private:
   bool isInsideGlobalDomain(const Vec3& pt, const real_t& radius) const;

   std::shared_ptr<blockforest::BlockForest> blockForest_;

   struct Subdomain
   {
      Subdomain(const int r, const BlockID& id, const math::AABB& ab) : rank(r), blockID(id), aabb(ab) {}
      int rank;
      BlockID blockID;
      math::AABB aabb;
   };

   int ownRank_ = -1;
   std::array< bool, 3 > periodic_;

   std::vector<math::AABB> localAABBs_;
   math::AABB              unionOfLocalAABBs_;
   std::vector<Subdomain>  neighborSubdomains_;
   std::vector<uint_t>     neighborProcesses_;
};

} //namespace domain

inline real_t sqDistanceLineToPoint( const real_t& pt, const real_t& min, const real_t& max  )
{
   if (pt < min)
      return (min - pt) * (min - pt);
   if (pt > max)
      return (pt - max) * (pt - max);
   return real_t(0);
}

inline real_t sqDistancePointToAABB( const Vec3& pt, const math::AABB& aabb )
{
   real_t sq = 0.0;

   sq += sqDistanceLineToPoint( pt[0], aabb.xMin(), aabb.xMax() );
   sq += sqDistanceLineToPoint( pt[1], aabb.yMin(), aabb.yMax() );
   sq += sqDistanceLineToPoint( pt[2], aabb.zMin(), aabb.zMax() );

   return sq;
}

inline real_t sqDistancePointToAABBPeriodic( Vec3 pt,
                                             const math::AABB& aabb,
                                             const math::AABB& domain,
                                             const std::array< bool, 3 >& periodic )
{
   auto size = domain.sizes() * real_t(0.5);
   auto d = pt - aabb.center();

   if (periodic[0] && (d[0] < -size[0])) pt[0] += domain.sizes()[0];
   if (periodic[0] && (d[0] > +size[0])) pt[0] -= domain.sizes()[0];

   if (periodic[1] && (d[1] < -size[1])) pt[1] += domain.sizes()[1];
   if (periodic[1] && (d[1] > +size[1])) pt[1] -= domain.sizes()[1];

   if (periodic[2] && (d[2] < -size[2])) pt[2] += domain.sizes()[2];
   if (periodic[2] && (d[2] > +size[2])) pt[2] -= domain.sizes()[2];

   real_t sq = 0.0;

   sq += sqDistanceLineToPoint( pt[0], aabb.xMin(), aabb.xMax() );
   sq += sqDistanceLineToPoint( pt[1], aabb.yMin(), aabb.yMax() );
   sq += sqDistanceLineToPoint( pt[2], aabb.zMin(), aabb.zMax() );

   return sq;
}

inline bool isInsideAABB( const Vec3& pt,
                          const real_t radius,
                          const math::AABB& aabb)
{
   if (!aabb.contains(pt)) return false;
   if ((pt[0] - aabb.xMin()) < radius) return false;
   if ((aabb.xMax() - pt[0]) < radius) return false;
   if ((pt[1] - aabb.yMin()) < radius) return false;
   if ((aabb.yMax() - pt[1]) < radius) return false;
   if ((pt[2] - aabb.zMin()) < radius) return false;
   if ((aabb.zMax() - pt[2]) < radius) return false;
   return true;
}

} //namespace mesa_pd
} //namespace walberla
