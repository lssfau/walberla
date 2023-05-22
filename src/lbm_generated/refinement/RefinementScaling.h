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
//! \file RefinementScaling.h
//! \author Markus Holzer <markus.holzer@fau.de>
//
//======================================================================================================================

#pragma once

#include "blockforest/BlockDataHandling.h"

#include "domain_decomposition/IBlock.h"
#include "domain_decomposition/StructuredBlockStorage.h"

namespace walberla
{
namespace lbm_generated
{

class DefaultRefinementScaling : public blockforest::AlwaysInitializeBlockDataHandling< real_t >
{
 public:
   DefaultRefinementScaling(const weak_ptr< StructuredBlockStorage >& blocks, const real_t parameter)
      : blocks_(blocks), parameter_(parameter){};

   real_t* initialize(IBlock* const block) override
   {
      WALBERLA_ASSERT_NOT_NULLPTR(block)
      auto blocks = blocks_.lock();
      WALBERLA_CHECK_NOT_NULLPTR(blocks)

      level_ = block->getBlockStorage().getLevel(*block);

      const real_t level_scale_factor = real_c(uint_t(1) << level_);
      const real_t one                = real_c(1.0);
      const real_t half               = real_c(0.5);

      return new real_t(parameter_ / (level_scale_factor * (-parameter_ * half + one) + parameter_ * half));
   }
   bool operator==(const DefaultRefinementScaling& other) const { return level_ == other.level_; }

 private:
   const weak_ptr< StructuredBlockStorage > blocks_;
   const real_t parameter_;

   uint_t level_;
};

} // namespace lbm_generated
} // namespace walberla