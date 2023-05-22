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
//! \file ParallelStreams.h
//! \ingroup gpu
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================
#pragma once
#include <vector>

#include "gpu/ErrorChecking.h"
#include "gpu/GPURAII.h"

namespace walberla {
namespace gpu
{

   class ParallelStreams;

   class ParallelSection
   {
   public:
      ~ParallelSection();
      void run( const std::function<void( gpuStream_t )> &f );

      gpuStream_t stream();
      void next();

   private:
      friend class ParallelStreams;

      ParallelSection( ParallelStreams *parent, gpuStream_t mainStream );
      void synchronize();

      ParallelStreams * parent_;
      gpuStream_t mainStream_;
      gpuEvent_t startEvent_;
      uint_t counter_;
   };


   //*******************************************************************************************************************
   /*!
    * Helper class to run CUDA/HIP operations on parallel streams
    *
    * This class introduces "side streams" that overlap with one "main stream". In a parallel section, multiple
    * kernels (or other CUDA/HIP operations) are scheduled to the streams. The first "run" is scheduled on the main stream
    * all subsequent operations on the side streams. The passed priority affects only the side streams. When
    * the parallel section goes out of scope the side streams are synchronized to the main stream via CUDA/HIP events.
    *
    * Example:
    *
    * \code
    * ParallelStreams streams;
    * {
    *   // new scope for the parallel section
    *   ParallelSection sec = streams.parallelSection( mainGPUStream );
    *   sec.run([&] ( gpuStream_t sideStream ) {
    *       // run something on the side stream
    *   });
    *   // after the parallel section goes out of scope the side streams are synchronized to the main stream
    * }
    *
    * \endcode
    *
    */
   //*******************************************************************************************************************
   class ParallelStreams
   {
   public:
      ParallelStreams( int priority = 0 );
      ParallelSection parallelSection( gpuStream_t stream );
      void setStreamPriority( int priority );

   private:
      friend class ParallelSection;

      void ensureSize( uint_t size );

      std::vector<StreamRAII> sideStreams_;
      std::vector<EventRAII> events_;
      EventRAII mainEvent_;
      int streamPriority_;
   };


} // namespace gpu
} // namespace walberla