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

#include "gpu/GPUWrapper.h"
#include "gpu/communication/CustomMemoryBuffer.h"
#include "gpu/communication/GeneratedGPUPackInfo.h"

namespace walberla {
namespace gpu
{
namespace communication {



/**
 * \brief Communication scheme for buffered communication in uniform block grids.
 *
 * Synchronize a set of \ref GPUField between GPU devices.
 * Communication between fields on the same process: use direct copy
 * via \ref GeneratedGPUPackInfo::communicateLocal.
 * Communication between different processes: use a buffered communication scheme;
 * when multiple fields have been changed they can be synchronized at once,
 * using one MPI message per communication partner.
 *
 *   \code
 *      UniformGPUScheme<stencil::D3Q19> scheme;  // the stencil defines the communication neighbors
 *      scheme.addPackInfo( make_shared<gpu::communication::MemcpyPackInfo<FieldType> >( idOfFirstField ) );
 *      scheme.addPackInfo( make_shared<gpu::communication::MemcpyPackInfo<FieldType> >( idOfSecondField ) );
 *
 *      // either synchronous communication...
 *      scheme();
 *
 *      // .. or asynchronous:
 *      scheme.startCommunication();
 *      functionWhichDoesNotNeedCommunicatedValues();
 *      scheme.wait();
 *   \endcode
 *
 * This scheme sends one message per communication step and neighbor device.
 * Therefore all contents that have to be sent are packed into a single buffer.
 * Multiple \ref GeneratedGPUPackInfo can be registered to send their contents in a single step.
 *
 * When running multiple \ref UniformGPUScheme concurrently, different MPI tags
 * have to be used for the schemes: the tag can be passed in the constructor.
 */
   template<typename Stencil>
   class UniformGPUScheme
   {
   public:
       explicit UniformGPUScheme( const weak_ptr< StructuredBlockForest >& bf,
                                  const bool sendDirectlyFromGPU = false,
                                  const bool useLocalCommunication = true,
                                  const int tag = 5432 );

       explicit UniformGPUScheme( const weak_ptr< StructuredBlockForest >& bf,
                                  const Set<SUID> & requiredBlockSelectors,
                                  const Set<SUID> & incompatibleBlockSelectors,
                                  const bool sendDirectlyFromGPU = false,
                                  const bool useLocalCommunication = true,
                                  const int tag = 5432 );
       ~UniformGPUScheme()
       {
          for (uint_t i = 0; i < Stencil::Q; ++i)
             WALBERLA_GPU_CHECK(gpuStreamDestroy(streams_[i]))
       }

       void addPackInfo( const shared_ptr<GeneratedGPUPackInfo> &pi );

       void startCommunication();
       void wait();

       void operator()()         { communicate( ); }
       inline void communicate() { startCommunication(); wait(); }

       std::function<void()> getCommunicateFunctor();
       std::function<void()> getStartCommunicateFunctor();
       std::function<void()> getWaitFunctor();

   private:
       void setupCommunication();

       weak_ptr<StructuredBlockForest> blockForest_;
       uint_t forestModificationStamp_;

       bool setupBeforeNextCommunication_;
       bool communicationInProgress_;
       const bool sendFromGPU_;
       const bool useLocalCommunication_;

       using CpuBuffer_T = gpu::communication::PinnedMemoryBuffer;
       using GpuBuffer_T = gpu::communication::GPUMemoryBuffer;

       mpi::GenericBufferSystem<CpuBuffer_T, CpuBuffer_T> bufferSystemCPU_;
       mpi::GenericBufferSystem<GpuBuffer_T, GpuBuffer_T> bufferSystemGPU_;

       std::vector<shared_ptr<GeneratedGPUPackInfo> > packInfos_;

       struct Header
       {
           BlockID blockId;
           stencil::Direction dir;
       };
       std::map<mpi::MPIRank, std::vector<Header> > headers_;

       Set<SUID> requiredBlockSelectors_;
       Set<SUID> incompatibleBlockSelectors_;

       gpuStream_t streams_[Stencil::Q];
   };


} // namespace communication
} // namespace gpu
} // namespace walberla

#include "UniformGPUScheme.impl.h"
