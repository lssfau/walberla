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
//! \file BlockForestDataHandling.cpp
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include <mesa_pd/domain/BlockForestDataHandling.h>

#include <mesa_pd/data/DataTypes.h>
#include <mesa_pd/data/ParticleStorage.h>

#include <mesa_pd/mpi/notifications/ParticleCopyNotification.h>

#include "blockforest/BlockDataHandling.h"
#include "blockforest/BlockForest.h"
#include "core/Abort.h"

namespace walberla{
namespace mesa_pd{
namespace domain {

namespace internal {
ParticleDeleter::ParticleDeleter(const std::shared_ptr<data::ParticleStorage>& ps,
                                 const math::AABB& aabb) :
   ps_(ps),
   aabb_(aabb)
{}

ParticleDeleter::~ParticleDeleter()
{
   for (auto pIt = ps_->begin(); pIt != ps_->end(); )
   {
      if (aabb_.contains(pIt->getPosition()))
      {
         pIt = ps_->erase(pIt);
         continue;
      }
      ++pIt;
   }
}

bool operator==(const ParticleDeleter& lhs, const ParticleDeleter& rhs)
{
   return ((lhs.ps_ == rhs.ps_) && (lhs.aabb_ == rhs.aabb_));
}
} //namespace internal

internal::ParticleDeleter* BlockForestDataHandling::initialize( IBlock * const block )
{
   return new internal::ParticleDeleter(ps_, block->getAABB());
}

void BlockForestDataHandling::serialize( IBlock * const block, const BlockDataID & /*id*/, mpi::SendBuffer & buffer )
{
   decltype(ps_->size()) numOfParticles = 0;

   //allocate bytes to store the number of particles which are sent
   //this number is only known at the end -> this value is written at the end
   auto ptr = buffer.allocate<decltype(ps_->size())>();

   for( auto pIt = ps_->begin(); pIt != ps_->end(); ++pIt)
   {
      //since ps_ is process local we have to skip particles which are not located on this block
      if ( !block->getAABB().contains( pIt->getPosition()) ) continue;
      //skip ghosts
      if (data::particle_flags::isSet( pIt->getFlags(), data::particle_flags::GHOST)) continue;
      //skip globals
      if (data::particle_flags::isSet( pIt->getFlags(), data::particle_flags::GLOBAL)) continue;

      buffer << ParticleCopyNotification( *pIt );
      ++numOfParticles;
   }

   *ptr = numOfParticles;
}

internal::ParticleDeleter* BlockForestDataHandling::deserialize( IBlock * const block )
{
   return initialize(block);
}

void BlockForestDataHandling::deserialize( IBlock * const block, const BlockDataID & id, mpi::RecvBuffer & buffer )
{
   deserializeImpl( block, id, buffer);
}

void BlockForestDataHandling::serializeCoarseToFine( Block * const block, const BlockDataID & /*id*/, mpi::SendBuffer & buffer, const uint_t child )
{
   // get child aabb
   const auto childID   = BlockID(block->getId(), child);
   const auto childAABB = block->getForest().getAABBFromBlockId(childID);
   //WALBERLA_LOG_DEVEL( (child & uint_t(1)) << (child & uint_t(2)) << (child & uint_t(4)) << "\naabb: " << aabb << "\nchild: " << childAABB );

   decltype(ps_->size()) numOfParticles = 0;

   //allocate bytes to store the number of particles which are sent
   //this number is only known at the end -> this value is written at the end
   auto ptr = buffer.allocate<decltype(ps_->size())>();

   for (auto pIt = ps_->begin(); pIt != ps_->end(); ++pIt)
   {
      //since ps_ is process local we have to skip particles which are not located on this block
      if ( !block->getAABB().contains( pIt->getPosition()) ) continue;
      //skip ghosts
      if (data::particle_flags::isSet( pIt->getFlags(), data::particle_flags::GHOST)) continue;
      //skip globals
      if (data::particle_flags::isSet( pIt->getFlags(), data::particle_flags::GLOBAL)) continue;

      if( childAABB.contains( pIt->getPosition()) )
      {
         buffer << ParticleCopyNotification( *pIt );
         ++numOfParticles;
      }
   }
   *ptr = numOfParticles;
}

void BlockForestDataHandling::serializeFineToCoarse( Block * const block, const BlockDataID & /*id*/, mpi::SendBuffer & buffer )
{
   decltype(ps_->size()) numOfParticles = 0;

   //allocate bytes to store the number of particles which are sent
   //this number is only known at the end -> this value is written at the end
   auto ptr = buffer.allocate<decltype(ps_->size())>();

   for (auto pIt = ps_->begin(); pIt != ps_->end(); ++pIt)
   {
      //since ps_ is process local we have to skip particles which are not located on this block
      if ( !block->getAABB().contains( pIt->getPosition()) ) continue;
      //skip ghosts
      if (data::particle_flags::isSet( pIt->getFlags(), data::particle_flags::GHOST)) continue;
      //skip globals
      if (data::particle_flags::isSet( pIt->getFlags(), data::particle_flags::GLOBAL)) continue;

      buffer << ParticleCopyNotification( *pIt );
      ++numOfParticles;
   }

   *ptr = numOfParticles;
}

internal::ParticleDeleter* BlockForestDataHandling::deserializeCoarseToFine( Block * const block )
{
   return initialize(block);
}

internal::ParticleDeleter* BlockForestDataHandling::deserializeFineToCoarse( Block * const block )
{
   return initialize(block);
}

void BlockForestDataHandling::deserializeCoarseToFine( Block * const block, const BlockDataID & id, mpi::RecvBuffer & buffer )
{
   deserializeImpl( block, id, buffer);
}


void BlockForestDataHandling::deserializeFineToCoarse( Block * const block, const BlockDataID & id, mpi::RecvBuffer & buffer, const uint_t /*child*/ )
{
   deserializeImpl( block, id, buffer);
}


void BlockForestDataHandling::deserializeImpl( IBlock * const block, const BlockDataID & /*id*/, mpi::RecvBuffer & buffer )
{
   decltype(ps_->size()) numBodies = 0;
   buffer >> numBodies;

   while( numBodies > 0 )
   {
      typename ParticleCopyNotification::Parameters objparam;
      buffer >> objparam;

      auto pIt = createNewParticle(*ps_, objparam);
      WALBERLA_CHECK(!data::particle_flags::isSet(pIt->getFlags(), data::particle_flags::GHOST));
      pIt->setOwner( MPIManager::instance()->rank() );

      if ( !block->getAABB().contains( pIt->getPosition()) )
      {
         WALBERLA_ABORT("Loaded particle not contained within block!\n" << "aabb: " << block->getAABB() << "\nparticle:" << *pIt );
      }

      --numBodies;
   }
}

std::shared_ptr<BlockForestDataHandling> createBlockForestDataHandling(const std::shared_ptr<data::ParticleStorage>& ps)
{
   return make_shared<BlockForestDataHandling >( ps );
}

} //namespace domain
} //namespace mesa_pd
} //namespace walberla
