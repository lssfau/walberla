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
//! \ingroup cuda
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================
#pragma once

#include "ErrorChecking.h"

namespace walberla {
namespace cuda {


   class StreamRAII
   {
   public:
      ~StreamRAII()
      {
         if( stream_ != 0 ) {
            WALBERLA_CUDA_CHECK( cudaStreamDestroy( stream_ ));
         }
      }

      StreamRAII( StreamRAII &&other )
      {
         stream_ = other.stream_;
         other.stream_ = 0;
      }

      StreamRAII( const StreamRAII & ) = delete;

      void operator=( const StreamRAII & ) = delete;

      operator cudaStream_t() const { return stream_; }


      static StreamRAII defaultStream()
      {
         StreamRAII result;
         result.stream_ = 0;
         return result;
      }

      static StreamRAII newPriorityStream( int priority )
      {
         StreamRAII result;
         WALBERLA_CUDA_CHECK( cudaStreamCreateWithPriority( &result.stream_, cudaStreamDefault, priority ));
         return result;
      }

      static StreamRAII newStream()
      {
         StreamRAII result;
         WALBERLA_CUDA_CHECK( cudaStreamCreate( &result.stream_ ));
         return result;
      }

   private:
      StreamRAII() {}

      cudaStream_t stream_;
   };


   class EventRAII
   {
   public:
      explicit EventRAII()
      {
         event = cudaEvent_t();
         WALBERLA_CUDA_CHECK( cudaEventCreate( &event ));
      }

      ~EventRAII()
      {
         if( event != cudaEvent_t() )
         {
            WALBERLA_CUDA_CHECK( cudaEventDestroy( event ));
         }
      }

      EventRAII( const EventRAII & ) = delete;

      void operator=( const EventRAII & ) = delete;

      EventRAII( EventRAII &&other )
      {
         event = other.event;
         other.event = cudaEvent_t();
      }

      operator cudaEvent_t() const { return event; }

   private:
      cudaEvent_t event;
   };


} // namespace cuda
} // namespace walberla