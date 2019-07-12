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
//! \file SyncGhostOwners.cpp
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

//======================================================================================================================
//
//  THIS FILE IS GENERATED - PLEASE CHANGE THE TEMPLATE !!!
//
//======================================================================================================================

#include "SyncGhostOwners.h"

#include <mesa_pd/mpi/RemoveAndNotify.h>

namespace walberla {
namespace mesa_pd {
namespace mpi {

void SyncGhostOwners::operator()( data::ParticleStorage& ps,
                                  const domain::IDomain& domain,
                                  const real_t dx,
                                  const bool syncNonCommunicatingBodies ) const
{
   if (numProcesses_ == 1) return;

   //==========================================================
   // STEP1: Update & Migrate
   //==========================================================
   updateAndMigrate( ps, domain, syncNonCommunicatingBodies );

   //==========================================================
   // STEP2: Check & Resolve
   //==========================================================
   checkAndResolveOverlap( ps, domain, dx, syncNonCommunicatingBodies );
}

void SyncGhostOwners::updateAndMigrate( data::ParticleStorage& ps,
                                        const domain::IDomain& domain,
                                        const bool syncNonCommunicatingBodies ) const
{
   using namespace walberla::mesa_pd::data::particle_flags;
   //==========================================================
   // STEP1: Update & Migrate
   //==========================================================

   WALBERLA_CHECK(!bs1.isCommunicationRunning());

   WALBERLA_LOG_DETAIL( "Assembling of Update&Migrate starts..." );
   std::set<walberla::mpi::MPIRank> recvRanks; // potential message senders
   for( auto pIt = ps.begin(); pIt != ps.end(); ++pIt)
   {
      if (isSet( pIt->getFlags(), GHOST))
      {
         if (!isSet( pIt->getFlags(), NON_COMMUNICATING) || syncNonCommunicatingBodies)
         {
            recvRanks.insert(pIt->getOwner());
         }
      }
   }

   for( auto pIt = ps.begin(); pIt != ps.end(); )
   {
      if (isSet( pIt->getFlags(), GHOST))
      {
         ++pIt;
         continue;
      }

      //==================
      // LOCAL

      //skip all particles that do not communicate (create ghost particles) on other processes
      if (isSet( pIt->getFlags(), NON_COMMUNICATING) && !syncNonCommunicatingBodies)
      {
         ++pIt;
         continue;
      }

      //correct position to make sure particle is always inside the domain!
      //everything is decided by the master particle therefore ghost particles are not touched
      if (!data::particle_flags::isSet( pIt->getFlags(), data::particle_flags::FIXED) &&
          !data::particle_flags::isSet( pIt->getFlags(), data::particle_flags::GHOST))
      {
         domain.periodicallyMapToDomain( pIt->getPositionRef() );
      }

      // Update
      for (auto ghostOwner : pIt->getGhostOwners())
      {
         WALBERLA_LOG_DETAIL( "Sending update notification for body " << pIt->getUid() << " to process " << ghostOwner );
         walberla::mpi::SendBuffer& sb = bs1.sendBuffer(static_cast<walberla::mpi::MPIRank>(ghostOwner));
         if (sb.isEmpty()) sb << walberla::uint8_c(0);
         packNotification(sb, ParticleUpdateNotification( *pIt ));
      }

      //particle has left subdomain?
      const auto newOwner = domain.findContainingProcessRank( pIt->getPosition() );
      if( newOwner != int_c(rank_) )
      {
         if ( newOwner < 0)
         {
            // No owner found: Outflow condition.
            WALBERLA_LOG_DETAIL( "Sending deletion notifications for body " << pIt->getUid() << " due to outflow." );

            //delete body
            pIt = removeAndNotify( bs1, ps, pIt );

            continue;
         }

         // Set new owner and transform to shadow copy
         pIt->setOwner( newOwner );
         set( pIt->getFlagsRef(), GHOST );

         // currently position is mapped to periodically to global domain,
         // this might not be the correct position for a ghost particle
         domain.correctParticlePosition( pIt->getPositionRef() );

         // Correct registration list (exclude new owner and us - the old owner) and
         // notify registered processes (except for new owner) of (remote) migration since they possess a ghost particle.
         auto ownerIt = std::find( pIt->getGhostOwners().begin(), pIt->getGhostOwners().end(), newOwner );
         WALBERLA_CHECK_UNEQUAL(ownerIt, pIt->getGhostOwners().end(), "New owner has to be former ghost owner!" );

         pIt->getGhostOwnersRef().erase( ownerIt );

         // Send remote migration notifications
         for( auto ghostRank : pIt->getGhostOwners() )
         {
            auto& buffer( bs1.sendBuffer(static_cast<walberla::mpi::MPIRank>(ghostRank)) );
            if (buffer.isEmpty()) buffer << walberla::uint8_c(0);

            WALBERLA_LOG_DETAIL( "Sending remote migration notification for particle " <<
                                 pIt->getUid() <<
                                 " to process " <<
                                 ghostRank );

            packNotification(buffer, ParticleRemoteMigrationNotification( *pIt, newOwner ));
         }

         pIt->getGhostOwnersRef().insert( int_c(rank_) );

         WALBERLA_LOG_DETAIL( "Sending migration notification for body " <<
                              pIt->getUid() <<
                              " to process " <<
                              (newOwner) );

         // Send migration notification to new owner
         auto& sb( bs1.sendBuffer(newOwner) );
         if (sb.isEmpty()) sb << walberla::uint8_c(0);
         packNotification(sb, ParticleMigrationNotification( *pIt ));

         pIt->getGhostOwnersRef().clear();

         continue;
      }
      ++pIt;
   }
   WALBERLA_LOG_DETAIL( "Assembling of Update&Migrate ended." );

   WALBERLA_LOG_DETAIL( "UM: number of recv " << recvRanks.size());
   bs1.setReceiverInfo(recvRanks, true);
   bs1.sendAll();
   WALBERLA_LOG_DETAIL( "UM: number of sends " << bs1.getNumberOfSends());

   // Receiving the updates for the remote rigid bodies from the connected processes
   WALBERLA_LOG_DETAIL( "Parsing of Update&Migrate starts..." );
   ParseMessage parseMessage;
   for( auto it = bs1.begin(); it != bs1.end(); ++it )
   {
      walberla::uint8_t tmp;
      it.buffer() >> tmp;
      while( !it.buffer().isEmpty() )
      {
         parseMessage(it.rank(), it.buffer(), ps, domain);
      }
   }
   WALBERLA_LOG_DETAIL( "Parsing of Update&Migrate ended." );
}

void SyncGhostOwners::checkAndResolveOverlap( data::ParticleStorage& ps,
                                              const domain::IDomain& domain,
                                              const real_t dx,
                                              const bool syncNonCommunicatingBodies ) const
{
   using namespace walberla::mesa_pd::data::particle_flags;
   //==========================================================
   // STEP2: Check&Resolve
   //==========================================================

   WALBERLA_CHECK(!bs2.isCommunicationRunning());

   //init buffers
   neighborRanks_ = domain.getNeighborProcesses();
   for( uint_t nbProcessRank : neighborRanks_ )
   {
      if (bs2.sendBuffer(nbProcessRank).isEmpty())
      {
         // fill empty buffers with a dummy byte to force transmission
         bs2.sendBuffer(nbProcessRank) << walberla::uint8_c(0);
      }
   }
   bs2.sendBuffer(int_c(rank_)) << walberla::uint8_c(0);

   WALBERLA_LOG_DETAIL( "Assembling of Check&Resolve starts..." );

   for( auto pIt = ps.begin(); pIt != ps.end(); )
   {
      //skip all particles that do not communicate (create ghost particles) on other processes
      if (isSet( pIt->getFlags(), NON_COMMUNICATING) && !syncNonCommunicatingBodies)
      {
          ++pIt;
          continue;
      }

      if (!isSet( pIt->getFlags(), GHOST))
      {
         //LOCAL

         walberla::mpi::SendBuffer& sbMaster = bs2.sendBuffer(pIt->getOwner());
         if (sbMaster.isEmpty()) sbMaster << walberla::uint8_c(0);

         // Update (nearest) neighbor processes.
         for( uint_t nbProcessRank : neighborRanks_ )
         {
            auto& sb = bs2.sendBuffer(nbProcessRank);
            if (sb.isEmpty()) sb << walberla::uint8_c(0);

            // dont send to owner!!
            if (pIt->getOwner() == int_c(nbProcessRank)) continue;
            // only send to neighbor which do not know this body
            if (pIt->getNeighborState().find( int_c(nbProcessRank) ) != pIt->getNeighborState().end()) continue;

            if( domain.intersectsWithProcessSubdomain( nbProcessRank, pIt->getPosition(), pIt->getInteractionRadius() + dx ) )
            {
               // no ghost there -> create ghost
               WALBERLA_LOG_DETAIL( "Sending copy notification for body " << pIt->getUid() << " to process " << (nbProcessRank) << "\n master: " << pIt->getOwner());
               packNotification(sb, ParticleCopyNotification( *pIt ));
               packNotification(sbMaster, NewGhostParticleNotification( *pIt, int_c(nbProcessRank) ));
               pIt->getNeighborStateRef().insert( int_c(nbProcessRank) );
            }
         }
      } else
      {
         //GHOST

         walberla::mpi::SendBuffer& sbMaster = bs2.sendBuffer(pIt->getOwner());
         if (sbMaster.isEmpty()) sbMaster << walberla::uint8_c(0);

         // Update (nearest) neighbor processes.
         for( uint_t nbProcessRank : neighborRanks_ )
         {
            auto& sb = bs2.sendBuffer(nbProcessRank);
            if (sb.isEmpty()) sb << walberla::uint8_c(0);

            if (pIt->getOwner() == int_c(nbProcessRank)) continue; // dont send to owner!!
            if (pIt->getNeighborState().find( int_c(nbProcessRank) ) != pIt->getNeighborState().end()) continue; // only send to neighbor which do not know this body

            if( domain.intersectsWithProcessSubdomain( nbProcessRank, pIt->getPosition(), pIt->getInteractionRadius() + dx ) )
            {
               // no ghost there -> create ghost
               WALBERLA_LOG_DETAIL( "Sending copy notification for body " << pIt->getUid() << " to process " << (nbProcessRank) << "\n master: " << pIt->getOwner());
               packNotification(sb, ParticleCopyNotification( *pIt ));
               packNotification(sbMaster, NewGhostParticleNotification( *pIt, int_c(nbProcessRank) ));
               pIt->getNeighborStateRef().insert( int_c(nbProcessRank) );
            }
         }

         if ( !domain.intersectsWithProcessSubdomain(uint_c(rank_), pIt->getPosition(), pIt->getInteractionRadius() + dx) )
         {
            // Delete
            // inform nearest neighbor processes.
            for( uint_t nbProcessRank : neighborRanks_ )
            {
               WALBERLA_LOG_DETAIL( "Sending removal information notification for body " << pIt->getUid() << " to process " << (nbProcessRank) );
               auto& sb = bs2.sendBuffer(nbProcessRank);
               if (sb.isEmpty()) sb << walberla::uint8_c(0);
               packNotification(sb, ParticleRemovalInformationNotification( *pIt ));
            }

            //notify owner
            WALBERLA_LOG_DETAIL( "Sending removal information notification for body " << pIt->getUid() << " to process " << (pIt->getOwner()) );
            auto& sb = bs2.sendBuffer(pIt->getOwner());
            if (sb.isEmpty()) sb << walberla::uint8_c(0);
            packNotification(sb, ParticleRemovalInformationNotification( *pIt ));

            pIt = ps.erase( pIt );
            continue;
         }
      }
      ++pIt;
   }

   std::set<walberla::mpi::MPIRank> recvRanks; // potential message senders
   // schedule receives
   for( auto pIt = ps.begin(); pIt != ps.end(); ++pIt)
   {
      if (isSet( pIt->getFlags(), GHOST)) continue;

      //skip all particles that do not communicate (create ghost particles) on other processes
      if (isSet( pIt->getFlags(), NON_COMMUNICATING) && !syncNonCommunicatingBodies) continue;

      for( auto ghostRank : pIt->getGhostOwners() )
      {
         recvRanks.insert(ghostRank);
      }
   }

   for( uint_t nbProcessRank : neighborRanks_ )
   {
      recvRanks.insert(int_c(nbProcessRank));
   }

   recvRanks.insert( int_c(rank_) );
   WALBERLA_LOG_DETAIL( "Assembling of Check&Resolve ended." );

   // size of buffer is unknown and changes with each send
   WALBERLA_LOG_DETAIL( "CR: number of recv " << recvRanks.size());
   bs2.setReceiverInfo(recvRanks, true);
   bs2.sendAll();
   WALBERLA_LOG_DETAIL( "CR: number of sends " << bs2.getNumberOfSends());

   // Receiving the updates for the remote rigid bodies from the connected processes
   WALBERLA_LOG_DETAIL( "Parsing of Check&Resolve starts..." );
   ParseMessage parseMessage;
   for( auto it = bs2.begin(); it != bs2.end(); ++it )
   {
      walberla::uint8_t tmp;
      it.buffer() >> tmp;
      while( !it.buffer().isEmpty() )
      {
         parseMessage(it.rank(), it.buffer(), ps, domain);
      }
   }
   WALBERLA_LOG_DETAIL( "Parsing of Check&Resolve ended." );
}

}  // namespace mpi
}  // namespace mesa_pd
}  // namespace walberla
