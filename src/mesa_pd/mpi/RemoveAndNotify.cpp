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
//! \file RemoveAndNotify.cpp
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include "RemoveAndNotify.h"

#include <mesa_pd/data/DataTypes.h>
#include <mesa_pd/data/Flags.h>
#include <mesa_pd/mpi/notifications/PackNotification.h>
#include <mesa_pd/mpi/notifications/ParticleRemovalNotification.h>

namespace walberla {
namespace mesa_pd {
namespace mpi {

/**
 * Removes a particle from the local storage and informs ghost particle holders.
 *
 * This function removes the particle from the particle storage and generates deletion notifications.
 */
data::ParticleStorage::iterator removeAndNotify( walberla::mpi::BufferSystem& bs,
                                                 data::ParticleStorage& ps,
                                                 data::ParticleStorage::iterator& pIt )
{
   WALBERLA_ASSERT( !data::particle_flags::isSet( pIt->getFlags(), data::particle_flags::GHOST),
                    "Trying to remove ghost particle from the particle storage." );

   WALBERLA_ASSERT( !data::particle_flags::isSet( pIt->getFlags(), data::particle_flags::GLOBAL),
                    "Trying to remove a global particle from the particle storage." );

   if( !pIt->getGhostOwners().empty() )
   {
      // Notify registered processes (intersecting or interacting) of particle removal since they possess a shadow copy.
      for( auto ghostRank : pIt->getGhostOwnersRef() )
      {
         WALBERLA_LOG_DETAIL( "__Notify registered process " << ghostRank << " of deletion of particle " << pIt->getUid() );
         auto& sb = bs.sendBuffer(static_cast<walberla::mpi::MPIRank>(ghostRank));
         if (sb.isEmpty()) sb << walberla::uint8_c(0);
         packNotification(sb, ParticleRemovalNotification( *pIt ));
      }
   }

   pIt->getGhostOwnersRef().clear();
   return ps.erase( pIt );
}

}  // namespace mpi
}  // namespace mesa_pd
}  // namespace walberla
