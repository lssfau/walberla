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
//! \file   sortParticleStorage.cpp
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include "sortParticleStorage.h"

namespace walberla {
namespace mesa_pd {

class RandomCompareFunctor
{
public:
   double getWeight(const data::Particle p1) const
   {
      math::seedRandomGenerator(static_cast<std::mt19937::result_type>(p1.getUid()));
      return math::realRandom<double>();
   }
};


void sortParticleStorage( data::ParticleStorage& ps,
                          const std::string& algorithm,
                          const math::AABB& domain,
                          const uint64_t cells )
{
   if( algorithm == "none" ) return;
   if( algorithm == "random" )
   {
      RandomCompareFunctor rand;
      ps.sort(rand);
      return;
   }
   if( algorithm == "linear" )
   {
      sorting::LinearizedCompareFunctor linear(domain, Vector3<uint_t>(cells));
      ps.sort(linear);
      return;
   }
   if( algorithm == "hilbert" )
   {
      auto tmp = uint_t(1) << math::uintMSBPosition(cells);
      if (tmp!=cells)
      {
         tmp <<= 1;
      }
      sorting::HilbertCompareFunctor hilbert(domain, tmp);
      ps.sort(hilbert);
      return;
   }
   WALBERLA_ABORT("unknown sorting algorithm: " << algorithm);
}

} // namespace mesa_pd
} // namespace walberla
