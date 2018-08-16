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
//! \file GetInfo.cpp
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include "GetInfo.h"

#include "core/debug/Debug.h"
#include <core/mpi/Reduce.h>
#include <pe/rigidbody/BodyStorage.h>
#include <pe/Types.h>

namespace walberla {
namespace pe {

std::pair<int64_t, int64_t> getNumBodies( const domain_decomposition::BlockStorage& bs, const BlockDataID& storageID, const int rank )
{
   int64_t numParticles       = 0;
   int64_t numShadowParticles = 0;
   for (const auto& blk : bs)
   {
      Storage const * storage = blk.getData< Storage >( storageID );
      const BodyStorage& localStorage = (*storage)[0];
      const BodyStorage& shadowStorage = (*storage)[1];
      numParticles += localStorage.size();
      numShadowParticles += shadowStorage.size();
   }
   if (rank == -1)
   {
      mpi::allReduceInplace(numParticles, mpi::SUM);
      mpi::allReduceInplace(numShadowParticles, mpi::SUM);
   } else
   {
      mpi::reduceInplace(numParticles, mpi::SUM, rank);
      mpi::reduceInplace(numShadowParticles, mpi::SUM, rank);
   }
   return std::make_pair(numParticles, numShadowParticles);
}

}  // namespace pe
}  // namespace walberla
