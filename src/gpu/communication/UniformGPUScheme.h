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
//! \file UniformGPUScheme.h
//! \ingroup gpu
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#pragma once

#include "blockforest/StructuredBlockForest.h"

#include "core/mpi/BufferSystem.h"
#include "core/mpi/MPIWrapper.h"

#include "domain_decomposition/IBlock.h"

#include "stencil/Directions.h"

#include <thread>

#include "gpu/GPURAII.h"
#include "gpu/GPUWrapper.h"
#include "gpu/ParallelStreams.h"
#include "gpu/communication/CustomMemoryBuffer.h"
#include "gpu/communication/GeneratedGPUPackInfo.h"

namespace walberla {
namespace gpu
{
namespace communication {



   template<typename Stencil>
   class UniformGPUScheme
   {
   public:
       explicit UniformGPUScheme( weak_ptr<StructuredBlockForest> bf,
                                  bool sendDirectlyFromGPU = false,
                                  bool useLocalCommunication = true,
                                  const int tag = 5432 );

       explicit UniformGPUScheme( weak_ptr<StructuredBlockForest> bf,
                                 const Set<SUID> & requiredBlockSelectors,
                                 const Set<SUID> & incompatibleBlockSelectors,
                                 bool sendDirectlyFromGPU = false,
                                 bool useLocalCommunication = true,
                                 const int tag = 5432 );

       void addPackInfo( const shared_ptr<GeneratedGPUPackInfo> &pi );

       void startCommunication( gpuStream_t stream = nullptr);
       void wait( gpuStream_t stream = nullptr);

       void operator()( gpuStream_t stream = nullptr )         { communicate( stream ); }
       inline void communicate( gpuStream_t stream = nullptr ) { startCommunication(stream); wait(stream); }

       std::function<void()> getCommunicateFunctor( gpuStream_t stream = nullptr );
       std::function<void()> getStartCommunicateFunctor( gpuStream_t stream = nullptr );
       std::function<void()> getWaitFunctor( gpuStream_t stream = nullptr );

   private:
       void setupCommunication();

       weak_ptr<StructuredBlockForest> blockForest_;
       uint_t forestModificationStamp_;

       bool setupBeforeNextCommunication_;
       bool communicationInProgress_;
       bool sendFromGPU_;
       bool useLocalCommunication_;

       using CpuBuffer_T = gpu::communication::PinnedMemoryBuffer;
       using GpuBuffer_T = gpu::communication::GPUMemoryBuffer;

       mpi::GenericBufferSystem<CpuBuffer_T, CpuBuffer_T> bufferSystemCPU_;
       mpi::GenericBufferSystem<GpuBuffer_T, GpuBuffer_T> bufferSystemGPU_;

       std::vector<shared_ptr<GeneratedGPUPackInfo> > packInfos_;

       ParallelStreams parallelSectionManager_;

       struct Header
       {
           BlockID blockId;
           stencil::Direction dir;
       };
       std::map<mpi::MPIRank, std::vector<Header> > headers_;

       Set<SUID> requiredBlockSelectors_;
       Set<SUID> incompatibleBlockSelectors_;
   };


} // namespace communication
} // namespace gpu
} // namespace walberla

#include "UniformGPUScheme.impl.h"
