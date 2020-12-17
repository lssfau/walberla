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

#include <mesa_pd/data/DataTypes.h>

namespace walberla {
namespace mesa_pd {
namespace domain {

/**
 * Abstract base class for the local subdomain
 */
class IDomain
{
public:
   virtual ~IDomain() = default;

   /// Is the point \p pt located inside the subdomain of the process with rank \p rank?
   /// \attention This function is supposed to check the local subdomain and next neighbor subdomains (periodicity).
   /// No global check is required!
   /// \return If you have no information about the specified process return false.
   virtual bool isContainedInProcessSubdomain(const uint_t rank, const Vec3& pt) const = 0;

   /// Is the sphere located inside the local subdomain?
   /// \attention This function is used for an early out. Therefore returning true has to be correct, returning false
   /// can also be wrong.
   virtual bool isContainedInLocalSubdomain(const Vec3& /*pt*/, const real_t& /*radius*/) const { return false; }

   /// Find the process rank which is responsible for \p pt.
   /// \attention No global information is required, just check local process and adjacent processes (periodicity!).
   /// \pre \p pt will be always inside the global domain.
   /// \return Returns the process rank or -1 if cannot be found within next neighbors.
   virtual int findContainingProcessRank(const Vec3& pt) const = 0;

   /// Map the point \p pt periodically into the global domain.
   /// \post pt has to be located inside the global domain
   virtual void periodicallyMapToDomain(Vec3& pt) const = 0;

   /// Returns a vector of ranks from all neighboring processes.
   virtual std::vector<uint_t> getNeighborProcesses() const = 0;

   /// Does the sphere defined by \p pt and \p radius intersect with the subdomain
   /// of process \p rank?
   /// \param pt center of the sphere
   /// \param radius radius of the sphere
   /// \return If you have no information about the specified process return false.
   virtual bool intersectsWithProcessSubdomain(const uint_t rank, const Vec3& pt, const real_t& radius) const = 0;

   /// Correct the particle position in regard to the local subdomain.
   virtual void correctParticlePosition(Vec3& pt) const = 0;
};

} //namespace domain
} //namespace mesa_pd
} //namespace walberla
