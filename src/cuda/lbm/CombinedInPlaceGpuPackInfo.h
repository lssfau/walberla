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
//! \file CombinedInPlaceGpuPackInfo.h
//! \author Frederik Hennig <frederik.hennig@fau.de>
//
//======================================================================================================================

#pragma once

#define IS_EVEN(x) ((x & 1) ^ 1)

#include "cuda/communication/GeneratedGPUPackInfo.h"

#include "lbm/inplace_streaming/TimestepTracker.h"

namespace walberla {
namespace lbm {

template< typename EvenPackInfo, typename OddPackInfo >
class CombinedInPlaceGpuPackInfo : public cuda::GeneratedGPUPackInfo
{
 public:
   template< typename... Args >
   CombinedInPlaceGpuPackInfo(std::shared_ptr< lbm::TimestepTracker >& tracker, Args&&... args)
      : tracker_(tracker), evenPackInfo_(std::forward< Args >(args)...), oddPackInfo_(std::forward< Args >(args)...)
   {}

   virtual ~CombinedInPlaceGpuPackInfo() = default;

   void pack(stencil::Direction dir, unsigned char* buffer, IBlock* block, cudaStream_t stream) override
   {
      if (IS_EVEN(tracker_->getCounter()))
      {
         evenPackInfo_.pack(dir, buffer, block, stream);
      }
      else
      {
         oddPackInfo_.pack(dir, buffer, block, stream);
      }
   }

   void unpack(stencil::Direction dir, unsigned char* buffer, IBlock* block, cudaStream_t stream) override {
      if (IS_EVEN(tracker_->getCounter()))
      {
         evenPackInfo_.unpack(dir, buffer, block, stream);
      }
      else
      {
         oddPackInfo_.unpack(dir, buffer, block, stream);
      }
   }

   uint_t size(stencil::Direction dir, IBlock* block) override {
      if (IS_EVEN(tracker_->getCounter()))
      {
         return evenPackInfo_.size(dir, block);
      }
      else
      {
         return oddPackInfo_.size(dir, block);
      }
   }

 private:
   const std::shared_ptr< lbm::TimestepTracker >& tracker_;
   EvenPackInfo evenPackInfo_;
   OddPackInfo oddPackInfo_;
};

} // namespace lbm
} // namespace walberla


