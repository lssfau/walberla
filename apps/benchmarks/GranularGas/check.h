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
//! \file   check.h
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include <mesa_pd/data/ParticleStorage.h>

#include <blockforest/BlockForest.h>
#include <core/grid_generator/SCIterator.h>
#include <core/logging/Logging.h>

namespace walberla {
namespace mesa_pd {

void check( data::ParticleStorage& ps, blockforest::BlockForest& forest, real_t spacing, const Vec3& shift )
{
   WALBERLA_LOG_INFO_ON_ROOT("*** CHECKING RESULT - START ***");
   auto pIt = ps.begin();
   for (auto& iBlk : forest)
   {
      for (auto it = grid_generator::SCIterator(iBlk.getAABB(), Vector3<real_t>(spacing, spacing, spacing) * real_c(0.5) + shift, spacing);
           it != grid_generator::SCIterator();
           ++it, ++pIt)
      {
         WALBERLA_CHECK_UNEQUAL(pIt, ps.end());
         WALBERLA_CHECK_FLOAT_EQUAL(pIt->getPositionRef(), *it);
      }
   }
   WALBERLA_LOG_INFO_ON_ROOT("*** CHECKING RESULT - END ***");
}

} //namespace mesa_pd
} //namespace walberla
