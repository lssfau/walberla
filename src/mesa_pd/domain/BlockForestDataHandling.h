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
//! \file BlockForestDataHandling.h
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#pragma once

#include "blockforest/BlockDataHandling.h"
#include "blockforest/BlockForest.h"
#include "core/Abort.h"

namespace walberla{
namespace mesa_pd{

namespace data { class ParticleStorage; }

namespace domain {

namespace internal {
/**
 * Class which gets initialized by BlockForestDataHandling
 *
 * This class is needed to delete particles if blocks get destroyed.
 * Since the data of MESA_PD is isolated from the blocks one has to
 * take care of correct removal of particles.
 *
 * If a block gets moved by the load balancing the remaining block is
 * deleted. Therefore all particles have to be destroyed.
 * However, if the DataHandling is used for checkpointing the particles
 * should not be deleted after serialization. Coupling the deletion
 * to the lifetime of the block solves both problems.
 *
 * \attention Has to be defined in the header file since data management
 * of blocks needs complete types.
 */
class ParticleDeleter
{
   friend bool operator==(const ParticleDeleter& lhs, const ParticleDeleter& rhs);
public:
   ParticleDeleter(const std::shared_ptr<data::ParticleStorage>& ps,
                   const math::AABB& aabb);

   ~ParticleDeleter();
private:
   std::shared_ptr<data::ParticleStorage> ps_;
   math::AABB aabb_; ///< AABB of the associated block.
};

bool operator==(const ParticleDeleter& lhs, const ParticleDeleter& rhs);
}

class BlockForestDataHandling: public blockforest::BlockDataHandling<internal::ParticleDeleter>
{
public:
   BlockForestDataHandling(const std::shared_ptr<data::ParticleStorage>& ps) : ps_(ps) {}
   virtual ~BlockForestDataHandling() {}

   virtual internal::ParticleDeleter* initialize( IBlock * const block ) override;

   virtual void serialize( IBlock * const block, const BlockDataID & id, mpi::SendBuffer & buffer ) override;
   virtual internal::ParticleDeleter* deserialize( IBlock * const block ) override;
   virtual void deserialize( IBlock * const block, const BlockDataID & id, mpi::RecvBuffer & buffer ) override;

   virtual void serializeCoarseToFine( Block * const block, const BlockDataID & id, mpi::SendBuffer & buffer, const uint_t child ) override;
   virtual void serializeFineToCoarse( Block * const block, const BlockDataID & id, mpi::SendBuffer & buffer ) override;

   virtual internal::ParticleDeleter* deserializeCoarseToFine( Block * const block ) override;
   virtual internal::ParticleDeleter* deserializeFineToCoarse( Block * const block ) override;

   virtual void deserializeCoarseToFine( Block * const block, const BlockDataID & id, mpi::RecvBuffer & buffer ) override;
   virtual void deserializeFineToCoarse( Block * const block, const BlockDataID & id, mpi::RecvBuffer & buffer, const uint_t child ) override;

private:
   void deserializeImpl( IBlock * const block, const BlockDataID & id, mpi::RecvBuffer & buffer );

   std::shared_ptr<data::ParticleStorage> ps_;
};

std::shared_ptr<BlockForestDataHandling> createBlockForestDataHandling(const std::shared_ptr<data::ParticleStorage>& ps);

} //namespace domain
} //namespace mesa_pd
} //namespace walberla
