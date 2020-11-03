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

#include <core/mpi/MPIManager.h>
#include <mesa_pd/domain/IDomain.h>

namespace walberla {
namespace mesa_pd {
namespace domain {

/**
 * Every process assumes the whole simulation space belongs to its subdomain.
 */
class InfiniteDomain : public IDomain
{
public:
   ///Everything belongs to the calling process.
   bool   isContainedInProcessSubdomain(const uint_t rank, const Vec3& /*pt*/) const override {return rank==rank_;}
   ///Everything belongs to the calling process.
   int    findContainingProcessRank(const Vec3& /*pt*/) const override {return static_cast<int>(rank_);}
   ///Nothing to do here since domain is infinite.
   void   periodicallyMapToDomain(Vec3& /*pt*/) const override {}
   ///If I own everything I do not have neighbors.
   std::vector<uint_t> getNeighborProcesses() const override {return {};}
   ///Everything belongs to my subdomain.
   bool   intersectsWithProcessSubdomain(const uint_t rank, const Vec3& /*pt*/, const real_t& /*radius*/) const override
   { return rank==rank_;}
   ///Nothing to do here.
   void   correctParticlePosition(Vec3& /*pt*/) const override {}
private:
   const uint_t rank_ = static_cast<uint_t>(MPIManager::instance()->rank());
};

} //namespace domain
} //namespace mesa_pd
} //namespace walberla
