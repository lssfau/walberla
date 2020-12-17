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

//======================================================================================================================
//
//  THIS FILE IS GENERATED - PLEASE CHANGE THE TEMPLATE !!!
//
//======================================================================================================================

#pragma once

#include <mesa_pd/data/DataTypes.h>
#include <mesa_pd/data/Flags.h>
#include <mesa_pd/data/ParticleStorage.h>
#include <mesa_pd/domain/IDomain.h>
#include <mesa_pd/mpi/notifications/PackNotification.h>
#include <mesa_pd/mpi/notifications/ParseMessage.h>
#include <mesa_pd/mpi/notifications/ParticleGhostCopyNotification.h>
#include <mesa_pd/mpi/notifications/ParticleMigrationNotification.h>
#include <mesa_pd/mpi/notifications/ParticleRemoteMigrationNotification.h>
#include <mesa_pd/mpi/notifications/ParticleRemovalNotification.h>
#include <mesa_pd/mpi/notifications/ParticleUpdateNotification.h>

#include <core/mpi/BufferSystem.h>
#include <core/logging/Logging.h>

namespace walberla {
namespace mesa_pd {
namespace mpi {

/**
 * Next neighbor synchronization kernel without ghost particles.
 *
 * This kernel conducts particle migrations if they move to another subdomain
 * but does not generate ghost particles on overlap. This can be useful for
 * particles without spatial extend like tracer particles.
 *
 * \ingroup mesa_pd_mpi
 */
class SyncNextNeighborsNoGhosts
{
public:
   void operator()(data::ParticleStorage& ps,
                   const domain::IDomain& domain) const;

   int64_t getBytesSent() const { return bytesSent_; }
   int64_t getBytesReceived() const { return bytesReceived_; }

   int64_t getNumberOfSends() const { return numberOfSends_; }
   int64_t getNumberOfReceives() const { return numberOfReceives_; }
private:
   void generateSynchronizationMessages(walberla::mpi::BufferSystem& bs,
                                        data::ParticleStorage& ps,
                                        const domain::IDomain& domain) const;
   mutable std::vector<uint_t> neighborRanks_; ///cache for neighbor ranks -> will be updated in operator()

   int numProcesses_ = walberla::mpi::MPIManager::instance()->numProcesses();
   int rank_         = walberla::mpi::MPIManager::instance()->rank();

   mutable int64_t bytesSent_ = 0;
   mutable int64_t bytesReceived_ = 0;
   mutable int64_t numberOfSends_ = 0;
   mutable int64_t numberOfReceives_ = 0;
};

}  // namespace mpi
}  // namespace mesa_pd
}  // namespace walberla