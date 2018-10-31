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
//! \file ParallelStreams.cpp
//! \ingroup cuda
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================


#include "cuda/ParallelStreams.h"

namespace walberla {
namespace cuda {


   ParallelSection::ParallelSection(ParallelStreams * parent, cudaStream_t mainStream)
     : parent_( parent ), mainStream_( mainStream ), counter_( 0 )
   {
      WALBERLA_CUDA_CHECK( cudaEventCreate(&startEvent_) );
      WALBERLA_CUDA_CHECK( cudaEventRecord( startEvent_, mainStream_ ) );
   }

   ParallelSection::~ParallelSection()
   {
      synchronize();
      WALBERLA_CUDA_CHECK( cudaEventDestroy(startEvent_) );
   }

   void ParallelSection::next()
   {
      if( counter_ > 0 ) {
         WALBERLA_CUDA_CHECK( cudaEventRecord( parent_->events_[counter_ - 1], parent_->sideStreams_[counter_ - 1] ) );
      }
      else {
         WALBERLA_CUDA_CHECK( cudaEventRecord( parent_->mainEvent_, mainStream_ ) );
      }
      ++counter_;

      parent_->ensureSize( counter_ );

      WALBERLA_CUDA_CHECK( cudaStreamWaitEvent( stream(), startEvent_, 0 ));
   }

   void ParallelSection::run(const std::function<void( cudaStream_t)> & f)
   {
      f( stream() );
      next();
   }

   void ParallelSection::synchronize()
   {
      for( uint_t i=0; i < counter_; ++i )
         for( uint_t j=0; j < counter_; ++j )
         {
            if( i == j )
               continue;

            auto & event  = i == 0 ? parent_->mainEvent_ : parent_->events_[i - 1];
            cudaStream_t stream = j == 0 ? mainStream_ : parent_->sideStreams_[j - 1];
            WALBERLA_CUDA_CHECK( cudaStreamWaitEvent( stream, event, 0 ));
         }

      WALBERLA_CUDA_CHECK( cudaEventRecord( startEvent_, mainStream_ ) );
   }

   cudaStream_t ParallelSection::stream()
   {
      return counter_ == 0 ? mainStream_ : parent_->sideStreams_[counter_ - 1];
   }



   ParallelStreams::ParallelStreams( int priority )
           : streamPriority_( priority )
   {
   }

   ParallelSection ParallelStreams::parallelSection( cudaStream_t stream ) {
      return ParallelSection(this, stream);
   }

   void ParallelStreams::ensureSize( uint_t size ) {
      for( uint_t i = sideStreams_.size(); i < size; ++i )
      {
         sideStreams_.emplace_back( StreamRAII::newPriorityStream(streamPriority_));
         events_.emplace_back( EventRAII() );
      }
   }

   void ParallelStreams::setStreamPriority( int priority )
   {
      streamPriority_ = priority;
      sideStreams_.clear();
      events_.clear();
   }



} // namespace cuda
} // namespace walberla