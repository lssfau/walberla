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
//! \file AssocToBlock.h
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
#include <mesa_pd/data/IAccessor.h>

#include <blockforest/BlockForest.h>

namespace walberla {
namespace mesa_pd {
namespace kernel {

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
   blockforest::Block*& currentBlock = ac.getCurrentBlockRef(idx);

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
         currentBlock = blk.second.get();
         return;
      }
   }

   WALBERLA_CHECK_NOT_NULLPTR(currentBlock, ac.getPosition(idx) << "\n" << bf_->begin()->getAABB());
}

} //namespace kernel
} //namespace mesa_pd
} //namespace walberla
