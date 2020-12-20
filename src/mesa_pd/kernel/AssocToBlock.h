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

#pragma once

#include <mesa_pd/data/DataTypes.h>
#include <mesa_pd/data/IAccessor.h>

#include <blockforest/BlockForest.h>

namespace walberla {
namespace mesa_pd {
namespace kernel {

/**
 * Kernel which updates the currentBlock property of all local properties.
 * All particles are checked against the blocks in the BlockForest and the property
 * is set accordingly.
 *
 * \attention This kernel must only be run on local particles. Ghost particles do not have
 * a corresponding block!
 * \post currentBlock property of all local particles is up-to-date.
 * \ingroup mesa_pd_kernel
 */
class AssocToBlock
{
public:
   explicit AssocToBlock(const std::shared_ptr<BlockForest>& bf) : bf_(bf) {}

   template <typename Accessor>
   void operator()(const size_t i, Accessor& ac) const;
private:
   std::shared_ptr<BlockForest> bf_ = nullptr;
};

template <typename Accessor>
inline void AssocToBlock::operator()(const size_t idx,
                                     Accessor& ac) const
{
   blockforest::Block* currentBlock = bf_->getBlock(ac.getCurrentBlock(idx));

   if (currentBlock != nullptr)
   {
      if (currentBlock->getAABB().contains(ac.getPosition(idx)))
      {
         return;
      } else
      {
         currentBlock = nullptr;
      }
   }

   for (auto& blk : bf_->getBlockMap())
   {
      if (blk.second->getAABB().contains(ac.getPosition(idx)))
      {
         ac.setCurrentBlock(idx, blk.second->getId());
         currentBlock = blk.second.get();
         return;
      }
   }

   //cannot happen if called only for local particles!
   //no "owning" block was found within the BlockForest...
   if  (currentBlock == nullptr)
   {
      WALBERLA_LOG_DEVEL( ac.getPosition(idx) );
      for (auto& blk : bf_->getBlockMap())
      {
         WALBERLA_LOG_DEVEL(blk.second->getAABB());
      }
      using namespace walberla::mesa_pd::data::particle_flags;
      WALBERLA_LOG_DEVEL_VAR(isSet(ac.getFlags(idx), GHOST));
   }
   WALBERLA_CHECK_NOT_NULLPTR(currentBlock, ac.getPosition(idx));
}

} //namespace kernel
} //namespace mesa_pd
} //namespace walberla
