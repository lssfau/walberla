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

#include "SyncNextNeighborsBlockForest.h"
#include <mesa_pd/domain/BlockForestDomain.h>

#include <mesa_pd/mpi/RemoveAndNotify.h>

namespace walberla {
namespace mesa_pd {
namespace mpi {

void correctParticlePosition(Vec3& pt,
                             const Vec3& center,
                             const math::AABB& domain,
                             const std::array<bool, 3>& periodic)
{
   const Vec3 dis = pt - center;

   if (periodic[0] && (-domain.xSize() * 0.5 > dis[0])) pt[0] += domain.xSize();
   if (periodic[0] && (+domain.xSize() * 0.5 < dis[0])) pt[0] -= domain.xSize();

   if (periodic[1] && (-domain.ySize() * 0.5 > dis[1])) pt[1] += domain.ySize();
   if (periodic[1] && (+domain.ySize() * 0.5 < dis[1])) pt[1] -= domain.ySize();

   if (periodic[2] && (-domain.zSize() * 0.5 > dis[2])) pt[2] += domain.zSize();
   if (periodic[2] && (+domain.zSize() * 0.5 < dis[2])) pt[2] -= domain.zSize();
}

void SyncNextNeighborsBlockForest::operator()(data::ParticleStorage& ps,
                                              const std::shared_ptr<blockforest::BlockForest>& bf,
                                              const std::shared_ptr<domain::BlockForestDomain>& domain,
                                              const real_t dx) const
{
   if (numProcesses_ == 1) return;

   walberla::mpi::BufferSystem bs( walberla::mpi::MPIManager::instance()->comm() );

   WALBERLA_CHECK(!bs.isCommunicationRunning());
   WALBERLA_CHECK_EQUAL(bs.size(), 0);

   for (auto& blk : bf->getBlockMap())
   {
      for( uint_t i = uint_t(0); i != blk.second->getNeighborhoodSize(); ++i )
      {
         auto nbProcessRank = blk.second->getNeighborProcess(i);
         if (bs.sendBuffer(nbProcessRank).isEmpty())
         {
            // fill empty buffers with a dummy byte to force transmission
            bs.sendBuffer(nbProcessRank) << walberla::uint8_c(0);
         }
      }
   }

   generateSynchronizationMessages(bs, ps, bf, dx);

   // size of buffer is unknown and changes with each send
   bs.setReceiverInfoFromSendBufferState(false, true);
   bs.sendAll();

   // Receiving the updates for the remote rigid bodies from the connected processes
   WALBERLA_LOG_DETAIL( "Parsing of particle synchronization response starts..." );
   ParseMessage parseMessage;
   for( auto it = bs.begin(); it != bs.end(); ++it )
   {
      walberla::uint8_t tmp;
      it.buffer() >> tmp;
      while( !it.buffer().isEmpty() )
      {
         parseMessage(it.rank(), it.buffer(), ps, *domain);
      }
   }
   WALBERLA_LOG_DETAIL( "Parsing of particle synchronization response ended." );

   bytesSent_ = bs.getBytesSent();
   bytesReceived_ = bs.getBytesReceived();
   numberOfSends_ = bs.getNumberOfSends();
   numberOfReceives_ = bs.getNumberOfReceives();
}

void SyncNextNeighborsBlockForest::generateSynchronizationMessages(walberla::mpi::BufferSystem& bs,
                                                                   data::ParticleStorage& ps,
                                                                   const std::shared_ptr<blockforest::BlockForest>& bf,
                                                                   const real_t dx) const
{
   const uint_t ownRank = uint_c(rank_);
   std::array<bool, 3> periodic;
   periodic[0] = bf->isPeriodic(0);
   periodic[1] = bf->isPeriodic(1);
   periodic[2] = bf->isPeriodic(2);

   WALBERLA_LOG_DETAIL( "Assembling of particle synchronization message starts..." );

   // position update
   for( auto pIt = ps.begin(); pIt != ps.end(); )
   {
      WALBERLA_ASSERT_GREATER(pIt->getInteractionRadius(), 0_r, "Did you forget to set the interaction radius?");

      //skip all ghost particles
      if (data::particle_flags::isSet( pIt->getFlags(), data::particle_flags::GHOST))
      {
         ++pIt;
         continue;
      }

      //skip all particles that do not communicate (create ghost particles) on other processes
      if (data::particle_flags::isSet( pIt->getFlags(), data::particle_flags::NON_COMMUNICATING))
      {
         ++pIt;
         continue;
      }

      auto currentBlock = bf->getBlock(pIt->getCurrentBlock());
      WALBERLA_CHECK_NOT_NULLPTR(currentBlock, *pIt);
      if (isInsideAABB(pIt->getPosition(), pIt->getInteractionRadius() + dx, currentBlock->getAABB()))
      {
         //no sync needed
         //just delete ghost particles if there are any

         for (const auto& ghostOwner : pIt->getGhostOwners() )
         {
            auto& buffer( bs.sendBuffer(static_cast<walberla::mpi::MPIRank>(ghostOwner)) );

            WALBERLA_LOG_DETAIL( "Sending removal notification for particle " << pIt->getUid() << " to process " << ghostOwner );

            packNotification(buffer, ParticleRemovalNotification( *pIt ));
         }

         pIt->getGhostOwnersRef().clear();

         ++pIt;
         continue;
      }

      //correct position to make sure particle is always inside the domain!
      //everything is decided by the master particle therefore ghost particles are not touched
      if (!data::particle_flags::isSet( pIt->getFlags(), data::particle_flags::FIXED) &&
          !data::particle_flags::isSet( pIt->getFlags(), data::particle_flags::GHOST))
      {
         bf->mapToPeriodicDomain( pIt->getPositionRef() );
      }

      // Note: At this point we know that the particle was locally owned before the position update.
      WALBERLA_CHECK_EQUAL(pIt->getOwner(), ownRank);

      WALBERLA_LOG_DETAIL( "Processing local particle " << pIt->getUid() );

      // Update nearest neighbor processes.
      auto isInsideDomain = isInsideAABB(pIt->getPosition(), pIt->getInteractionRadius(), bf->getDomain());
      std::vector<int> ranksAlreadyTreated{int_c(ownRank)};
      for( uint_t nb = uint_t(0); nb < currentBlock->getNeighborhoodSize(); ++nb )
      {
         auto nbProcessRank = currentBlock->getNeighborProcess(nb);
         if (std::find(ranksAlreadyTreated.begin(), ranksAlreadyTreated.end(), int_c(nbProcessRank)) != ranksAlreadyTreated.end())
         {
            continue;
         }
         auto nbAABB = currentBlock->getNeighborAABB(nb);
         auto sqDistance = isInsideDomain
                           ? sqDistancePointToAABB(pIt->getPosition(), nbAABB)
                           : sqDistancePointToAABBPeriodic(pIt->getPosition(), nbAABB, bf->getDomain(), periodic);
         auto tmp = pIt->getInteractionRadius() + dx;
         if( sqDistance <  tmp*tmp)
         {
            ranksAlreadyTreated.emplace_back(int_c(nbProcessRank));
            auto ghostOwnerIt = std::find( pIt->getGhostOwners().begin(), pIt->getGhostOwners().end(), nbProcessRank );
            if( ghostOwnerIt != pIt->getGhostOwners().end() )
            {
               // already a ghost there -> update
               auto& buffer( bs.sendBuffer(nbProcessRank) );
               WALBERLA_LOG_DETAIL( "Sending update notification for particle " << pIt->getUid() << " to process " << (nbProcessRank) );
               packNotification(buffer, ParticleUpdateNotification( *pIt ));
            } else
            {
               // no ghost there -> create ghost
               auto& buffer( bs.sendBuffer(nbProcessRank) );
               WALBERLA_LOG_DETAIL( "Sending ghost copy notification for particle " << pIt->getUid() << " to process " << (nbProcessRank) );
               packNotification(buffer, ParticleGhostCopyNotification( *pIt ));
               pIt->getGhostOwnersRef().insert( int_c(nbProcessRank) );
            }
         }
      }
      for (auto ghostOwnerIt = pIt->getGhostOwnersRef().begin();
           ghostOwnerIt != pIt->getGhostOwnersRef().end();
           )
      {
         if (std::find(ranksAlreadyTreated.begin(),
                       ranksAlreadyTreated.end(),
                       int_c(*ghostOwnerIt)) == ranksAlreadyTreated.end())
         {
            // In case the rigid particle no longer intersects the remote process nor interacts with it but is registered,
            // send removal notification.
            auto& buffer( bs.sendBuffer(*ghostOwnerIt) );

            WALBERLA_LOG_DETAIL( "Sending removal notification for particle " << pIt->getUid() << " to process " << *ghostOwnerIt );

            packNotification(buffer, ParticleRemovalNotification( *pIt ));

            ghostOwnerIt = pIt->getGhostOwnersRef().erase(ghostOwnerIt);

            continue;
         }
         ++ghostOwnerIt;
      }

      //particle has left subdomain?
      if (currentBlock->getAABB().contains(pIt->getPosition()))
      {
         // particle still is locally owned after position update.
         WALBERLA_LOG_DETAIL( "Owner of particle " << pIt->getUid() << " is still process " << pIt->getOwner() );
      } else
      {
         //find new owner
         int ownerRank = -1;
         for( uint_t i = uint_t(0); i != currentBlock->getNeighborhoodSize(); ++i )
         {
            if (currentBlock->getNeighborAABB(i).contains(pIt->getPosition()))
            {
               ownerRank = int_c(currentBlock->getNeighborProcess(i));
            }
         }

         if( ownerRank != int_c(ownRank) )
         {
            WALBERLA_LOG_DETAIL( "Local particle " << pIt->getUid() << " is no longer on process " << ownRank << " but on process " << ownerRank );

            if( ownerRank < 0 )
            {
               // No owner found: Outflow condition.
               WALBERLA_LOG_DETAIL( "Sending deletion notifications for particle " << pIt->getUid() << " due to outflow." );

               // Registered processes receive removal notification in the remove() routine.
               pIt = removeAndNotify( bs, ps, pIt );

               continue;
            }

            WALBERLA_LOG_DETAIL( "Sending migration notification for particle " << pIt->getUid() << " to process " << ownerRank << "." );
            //WALBERLA_LOG_DETAIL( "Process registration list before migration: " << pIt->getGhostOwners() );

            // Set new owner and transform to ghost particle
            pIt->setOwner(ownerRank);
            data::particle_flags::set( pIt->getFlagsRef(), data::particle_flags::GHOST );

            // currently position is mapped to periodically to global domain,
            // this might not be the correct position for a ghost particle
            correctParticlePosition( pIt->getPositionRef(), currentBlock->getAABB().center(), bf->getDomain(), periodic );

            // Correct registration list (exclude new owner and us - the old owner) and
            // notify registered processes (except for new owner) of (remote) migration since they possess a ghost particle.
            auto ownerIt = std::find( pIt->getGhostOwners().begin(), pIt->getGhostOwners().end(), ownerRank );
            WALBERLA_CHECK_UNEQUAL(ownerIt, pIt->getGhostOwners().end(), "New owner has to be former ghost owner!" );

            pIt->getGhostOwnersRef().erase( ownerIt );

            for( auto ghostRank : pIt->getGhostOwners() )
            {
               auto& buffer( bs.sendBuffer(static_cast<walberla::mpi::MPIRank>(ghostRank)) );

               WALBERLA_LOG_DETAIL( "Sending remote migration notification for particle " << pIt->getUid() <<
                                    " to process " << ghostRank );

               packNotification(buffer, ParticleRemoteMigrationNotification( *pIt, ownerRank ));
            }

            pIt->getGhostOwnersRef().insert( int_c(ownRank) );

            // Send migration notification to new owner
            auto& buffer( bs.sendBuffer(ownerRank) );
            packNotification(buffer, ParticleMigrationNotification( *pIt ));

            pIt->getGhostOwnersRef().clear();

            continue;
         }
      }

      ++pIt;
   }

   WALBERLA_LOG_DETAIL( "Assembling of particle synchronization message ended." );
}

}  // namespace mpi
}  // namespace mesa_pd
}  // namespace walberla
