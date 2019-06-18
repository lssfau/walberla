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
//! \file   sortParticleStorage.h
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include <mesa_pd/data/ParticleStorage.h>
#include <mesa_pd/sorting/HilbertCompareFunctor.h>
#include <mesa_pd/sorting/LinearizedCompareFunctor.h>

#include <core/Abort.h>
#include <core/math/Random.h>
#include <core/logging/Logging.h>

#include <memory>
#include <string>

namespace walberla {
namespace mesa_pd {

void sortParticleStorage( data::ParticleStorage& ps,
                          const std::string& algorithm,
                          const math::AABB& domain,
                          const uint64_t cells );

} // namespace mesa_pd
} // namespace walberla
