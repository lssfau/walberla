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
//! \file InfiniteDomain.h
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#pragma once

#include <mesa_pd/domain/IDomain.h>

namespace walberla {
namespace mesa_pd {
namespace domain {

class InfiniteDomain : public IDomain
{
public:
   bool   isContainedInProcessSubdomain(const uint_t /*rank*/, const Vec3& /*pt*/) const override {return true;}
   int    findContainingProcessRank(const Vec3& /*pt*/) const override {return mpi::MPIManager::instance()->rank();}
   void   periodicallyMapToDomain(Vec3& /*pt*/) const override {}
   std::vector<uint_t> getNeighborProcesses() const override {return {};}
   bool   intersectsWithProcessSubdomain(const uint_t rank, const Vec3& /*pt*/, const real_t& /*radius*/) const override
   { return int_c(rank)==mpi::MPIManager::instance()->rank() ? true : false;}
   void   correctParticlePosition(Vec3& /*pt*/) const override {}
};

} //namespace domain
} //namespace mesa_pd
} //namespace walberla
