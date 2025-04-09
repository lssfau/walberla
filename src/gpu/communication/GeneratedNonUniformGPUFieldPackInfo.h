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
//! \file GeneratedNonUniformGPUFieldPackInfo.h
//! \ingroup gpu
//! \author Philipp Suffa <philipp.suffa@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/DataTypes.h"
#include "core/cell/CellInterval.h"
#include "gpu/communication/GeneratedNonUniformGPUPackInfo.h"
#include "field/refinement/PackInfoHelper.h"

namespace walberla {

/*******************************************************************************************************************
 * \brief  Class for packing non-uniform fields, which are not PDF fields, such as velocity or pressure fields, on GPU.
 *
 * It overrides all needed function for non-uniform packing from the walberla::gpu::GeneratedNonUniformGPUPackInfo.
 *
 * The functions are gathering the relevant fields and defining the pack and unpack intervals for the packing kernels.
 * The actual packing kernels are generated and imported as template parameter PackingKernels_T.
 */
//*******************************************************************************************************************

template< typename Field_T, typename PackingKernels_T >
class GeneratedNonUniformGPUFieldPackInfo : public walberla::gpu::GeneratedNonUniformGPUPackInfo
{
 public:
   using value_type = typename Field_T::value_type;

   GeneratedNonUniformGPUFieldPackInfo( const BlockDataID fieldID )
      : fieldID_(fieldID) {};

   ~GeneratedNonUniformGPUFieldPackInfo() = default;

   bool constantDataExchange() const override { return true; };
   bool threadsafeReceiving() const override { return false; };

   void packDataEqualLevelImpl(const Block* sender, stencil::Direction dir, GpuBuffer_T & byte_buffer, gpuStream_t stream) const override;
   void unpackDataEqualLevel( Block * receiver, stencil::Direction dir, GpuBuffer_T & byte_buffer, gpuStream_t stream) override;
   void communicateLocalEqualLevel( const Block * sender, Block * receiver, stencil::Direction dir, gpuStream_t stream) override;

   void addForLocalEqualLevelComm(const Block* /*sender*/, Block* /*receiver*/, stencil::Direction /*dir*/) override {};
   void communicateLocalEqualLevel(uint64_t /*level*/, uint8_t /*timestep*/, gpuStream_t /*stream*/) override {};

   void packDataCoarseToFineImpl(const Block* coarseSender, const BlockID& fineReceiver, stencil::Direction dir, GpuBuffer_T & byte_buffer, gpuStream_t stream) const override;
   void unpackDataCoarseToFine      (       Block * fineReceiver, const BlockID & coarseSender, stencil::Direction dir, GpuBuffer_T & byte_buffer, gpuStream_t stream) override;
   void communicateLocalCoarseToFine( const Block * /*coarseSender*/, Block * /*fineReceiver*/, stencil::Direction /*dir*/, gpuStream_t /*stream*/) override {};
   void communicateLocalCoarseToFine( const Block * coarseSender, Block * fineReceiver, stencil::Direction dir, GpuBuffer_T & byte_buffer, gpuStream_t stream) override;

   void packDataFineToCoarseImpl(const Block* fineSender, const BlockID& coarseReceiver, stencil::Direction dir, GpuBuffer_T & byte_buffer, gpuStream_t stream) const override;
   void unpackDataFineToCoarse      (       Block * coarseReceiver, const BlockID & fineSender,     stencil::Direction dir, GpuBuffer_T & byte_buffer, gpuStream_t stream) override;
   void communicateLocalFineToCoarse( const Block * /*fineSender*/, Block * /*coarseReceiver*/, stencil::Direction /*dir*/, gpuStream_t /*stream*/) override {};
   void communicateLocalFineToCoarse( const Block * fineSender, Block * coarseReceiver, stencil::Direction dir, GpuBuffer_T & byte_buffer, gpuStream_t stream) override;

   uint_t sizeEqualLevelSend( const Block * sender, stencil::Direction dir) const override;
   uint_t sizeCoarseToFineSend ( const Block * coarseSender, const BlockID & fineReceiver, stencil::Direction dir) const override;
   uint_t sizeCoarseToFineReceive ( Block* fineReceiver, stencil::Direction dir) const override;
   uint_t sizeFineToCoarseSend ( const Block * fineSender, stencil::Direction dir) const override;

   void sync() override {};

 private:

   const BlockDataID fieldID_;
   PackingKernels_T kernels_;

};
} //namespace walberla

#include "gpu/communication/GeneratedNonUniformGPUFieldPackInfo.impl.h"
