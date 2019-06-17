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
//! \file HilbertCompareFunctor.h
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#pragma once

#include <core/math/AABB.h>
#include <core/math/Vector3.h>

#include <mesa_pd/data/DataTypes.h>
#include <mesa_pd/data/ParticleStorage.h>

#include <vector>

namespace walberla {
namespace mesa_pd {
namespace sorting {

/**
 * Compare functor which sorts particles according to the hilbert space filling curve.
 */
class HilbertCompareFunctor
{
public:
   /**
    * Subdivides a domain into cells and arranges the cells according to the sfc.
    * The position within the sfc gives the particles their weight.
    * \param domain domain which should be partitioned into cells
    * \param cells number of cells in every direction.
    * \attention cells has to be a power of 2!
    */
   HilbertCompareFunctor(const math::AABB& domain, const uint_t cells);
   bool operator()(const data::Particle bd1, const data::Particle bd2) const;
   double getWeight(const data::Particle p1) const;
private:
   uint_t discretize(const Vec3& pos) const;
   void initializeLookup();

   math::AABB  domain_;
   const uint_t cells_;
   Vec3 inverse_dx;
   std::vector<uint_t> hilbertLookup_;
};

inline
double HilbertCompareFunctor::getWeight(const data::Particle p1) const
{
   const auto hash1 = discretize(p1.getPosition() - domain_.minCorner());
   WALBERLA_ASSERT_LESS( hash1, cells_*cells_*cells_);

   return double_c(hilbertLookup_[hash1]);
}

} //namespace sorting
} //namespace mesa_pd
} //namespace walberla
