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
//! \file CudaRAII.h
//! \ingroup gpu
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================
#pragma once

#include "ErrorChecking.h"

namespace walberla
{
namespace gpu
{

class StreamRAII
{
 public:
   ~StreamRAII()
   {
      if (stream_ != nullptr) { WALBERLA_GPU_CHECK(gpuStreamDestroy(stream_)) }
   }

   StreamRAII(StreamRAII&& other) noexcept
   {
      stream_       = other.stream_;
      other.stream_ = nullptr;
   }

   StreamRAII(const StreamRAII&) = delete;

   void operator=(const StreamRAII&) = delete;

   operator gpuStream_t() const { return stream_; }

   static StreamRAII defaultStream()
   {
      StreamRAII result;
      result.stream_ = nullptr;
      return result;
   }

   static StreamRAII newPriorityStream(int priority)
   {
      StreamRAII result;
      WALBERLA_GPU_CHECK(gpuStreamCreateWithPriority(&result.stream_, gpuStreamDefault, priority))
      return result;
   }

   static StreamRAII newStream()
   {
      StreamRAII result;
      WALBERLA_GPU_CHECK(gpuStreamCreate(&result.stream_))
      return result;
   }

 private:
   StreamRAII() = default;

   gpuStream_t stream_;
};

class EventRAII
{
 public:
   explicit EventRAII()
   {
      event = gpuEvent_t();
      WALBERLA_GPU_CHECK(gpuEventCreate(&event))
   }

   ~EventRAII()
   {
      if (event != gpuEvent_t()) { WALBERLA_GPU_CHECK(gpuEventDestroy(event)) }
   }

   EventRAII(const EventRAII&) = delete;

   void operator=(const EventRAII&) = delete;

   EventRAII(EventRAII&& other) noexcept
   {
      event       = other.event;
      other.event = gpuEvent_t();
   }

   operator gpuEvent_t() const { return event; }

 private:
   gpuEvent_t event;
};

} // namespace gpu
} // namespace walberla