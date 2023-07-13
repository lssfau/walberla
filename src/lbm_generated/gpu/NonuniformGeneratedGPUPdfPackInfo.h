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
//! \file NonuniformGeneratedGPUPdfPackInfo.h
//! \author Markus Holzer <markus.holzer@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/DataTypes.h"

#include "gpu/GPUWrapper.h"
#include "gpu/communication/GeneratedNonUniformGPUPackInfo.h"

#include "lbm_generated/gpu/NonuniformGPUCommData.h"
#include "lbm_generated/field/PdfField.h"

namespace walberla::lbm_generated
{
using stencil::Direction;

namespace internal
{
/*
 * Base Template for Packing Kernels Wrapper. This wrapper is required for passing the time step to
 * kernels generated for in-place streaming patterns. The generated code should not be templated.
 */
template< typename PdfField_T, bool inplace >
class NonuniformGPUPackingKernelsWrapper
{
 public:
   void packAll(PdfField_T* srcField, CellInterval ci, unsigned char* outBuffer, gpuStream_t stream ) const  = 0;
   void unpackAll(PdfField_T* dstField, CellInterval ci, unsigned char* inBuffer, gpuStream_t stream ) const = 0;
   void localCopyAll(PdfField_T* srcField, CellInterval srcInterval, PdfField_T* dstField,
                     CellInterval dstInterval, gpuStream_t stream ) const                                    = 0;

   void packDirection(PdfField_T* srcField, CellInterval ci, unsigned char* outBuffer, Direction dir, gpuStream_t stream ) const  = 0;
   void unpackDirection(PdfField_T* dstField, CellInterval ci, unsigned char* inBuffer, Direction dir, gpuStream_t stream ) const = 0;
   void localCopyDirection(PdfField_T* srcField, CellInterval srcInterval, PdfField_T* dstField,
                           CellInterval dstInterval, Direction dir, gpuStream_t stream) const               = 0;

   void unpackRedistribute(PdfField_T* dstField, CellInterval ci, unsigned char* inBuffer,
                           stencil::Direction dir, gpuStream_t stream ) const = 0;

   void packPartialCoalescence(PdfField_T* srcField, PartialCoalescenceMaskFieldGPU* maskField, CellInterval ci,
                               unsigned char* outBuffer, Direction dir, gpuStream_t stream ) const                                   = 0;
   void zeroCoalescenceRegion(PdfField_T* dstField, CellInterval ci, Direction dir) const                      = 0;
   void unpackCoalescence(PdfField_T* dstField, CellInterval ci, unsigned char* inBuffer, Direction dir, gpuStream_t stream ) const = 0;

   uint_t size(CellInterval ci, Direction dir) const                   = 0;
   uint_t size(CellInterval ci) const                                  = 0;
   uint_t redistributeSize(CellInterval ci) const                      = 0;
   uint_t partialCoalescenceSize(CellInterval ci, Direction dir) const = 0;
};

/*
 * Template Specialization for two-fields patterns, with trivial method wrappers.
 */
template< typename PdfField_T >
class NonuniformGPUPackingKernelsWrapper< PdfField_T, false >
{
 public:
   using LatticeStorageSpecification_T = typename PdfField_T::LatticeStorageSpecification;
   using PackingKernels_T              = typename LatticeStorageSpecification_T::PackKernels;

   void packAll(PdfField_T* srcField, CellInterval ci, unsigned char* outBuffer, gpuStream_t stream = nullptr) const
   {
      kernels_.packAll(srcField, ci, outBuffer, stream);
   }

   void unpackAll(PdfField_T* dstField, CellInterval ci, unsigned char* inBuffer, gpuStream_t stream = nullptr) const
   {
      kernels_.unpackAll(dstField, ci, inBuffer, stream);
   }

   void localCopyAll(PdfField_T* srcField, CellInterval srcInterval, PdfField_T* dstField,
                     CellInterval dstInterval, gpuStream_t stream = nullptr) const
   {
      kernels_.localCopyAll(srcField, srcInterval, dstField, dstInterval, stream);
   }

   void packDirection(PdfField_T* srcField, CellInterval ci, unsigned char* outBuffer, Direction dir, gpuStream_t stream = nullptr) const
   {
      kernels_.packDirection(srcField, ci, outBuffer, dir, stream);
   }

   void unpackDirection(PdfField_T* dstField, CellInterval ci, unsigned char* inBuffer, Direction dir, gpuStream_t stream = nullptr) const
   {
      kernels_.unpackDirection(dstField, ci, inBuffer, dir, stream);
   }

   void localCopyDirection(PdfField_T* srcField, CellInterval srcInterval, PdfField_T* dstField,
                           CellInterval dstInterval, Direction dir, gpuStream_t stream) const
   {
      kernels_.localCopyDirection(srcField, srcInterval, dstField, dstInterval, dir, stream);
   }

   void unpackRedistribute(PdfField_T* dstField, CellInterval ci, unsigned char* inBuffer,
                           stencil::Direction dir, gpuStream_t stream = nullptr) const
   {
      kernels_.unpackRedistribute(dstField, ci, inBuffer, dir, stream);
   }

   void packPartialCoalescence(PdfField_T* srcField, PartialCoalescenceMaskFieldGPU* maskField, CellInterval ci,
                               unsigned char* outBuffer, Direction dir, gpuStream_t stream = nullptr) const
   {
      kernels_.packPartialCoalescence(srcField, maskField, ci, outBuffer, dir, stream);
   }

   void unpackCoalescence(PdfField_T* dstField, CellInterval ci, unsigned char* inBuffer, Direction dir, gpuStream_t stream = nullptr) const
   {
      kernels_.unpackCoalescence(dstField, ci, inBuffer, dir, stream);
   }

   void zeroCoalescenceRegion(PdfField_T* dstField, CellInterval ci, Direction dir, gpuStream_t stream = nullptr) const
   {
      kernels_.zeroCoalescenceRegion(dstField, ci, dir, stream);
   }

   uint_t size(CellInterval ci, Direction dir) const { return kernels_.size(ci, dir); }
   uint_t size(CellInterval ci) const { return kernels_.size(ci); }
   uint_t redistributeSize(CellInterval ci) const { return kernels_.redistributeSize(ci); }
   uint_t partialCoalescenceSize(CellInterval ci, Direction dir) const
   {
      return kernels_.partialCoalescenceSize(ci, dir);
   }

 private:
   PackingKernels_T kernels_;
};

/*
 * Template Specialization for in-place patterns, extracting the timestep from the lattice model.
 */
template< typename PdfField_T >
class NonuniformGPUPackingKernelsWrapper< PdfField_T, true >
{
 public:
   using LatticeStorageSpecification_T = typename PdfField_T::LatticeStorageSpecification;
   using PackingKernels_T              = typename LatticeStorageSpecification_T::PackKernels;

   void packAll(PdfField_T* srcField, CellInterval ci, unsigned char* outBuffer, gpuStream_t stream = nullptr) const
   {
      uint8_t timestep = srcField->getTimestep();
      kernels_.packAll(srcField, ci, outBuffer, timestep, stream);
   }

   void unpackAll(PdfField_T* dstField, CellInterval ci, unsigned char* inBuffer, gpuStream_t stream = nullptr) const
   {
      uint8_t timestep = dstField->getTimestep();
      kernels_.unpackAll(dstField, ci, inBuffer, timestep, stream);
   }

   void localCopyAll(PdfField_T* srcField, CellInterval srcInterval, PdfField_T* dstField,
                     CellInterval dstInterval, gpuStream_t stream = nullptr) const
   {
      uint8_t timestep = srcField->getTimestep();
      WALBERLA_ASSERT_EQUAL(timestep, dstField->getTimestep())
      kernels_.localCopyAll(srcField, srcInterval, dstField, dstInterval, timestep, stream);
   }

   void packDirection(PdfField_T* srcField, CellInterval ci, unsigned char* outBuffer, Direction dir, gpuStream_t stream = nullptr) const
   {
      uint8_t timestep = srcField->getTimestep();
      kernels_.packDirection(srcField, ci, outBuffer, dir, timestep, stream);
   }

   void unpackDirection(PdfField_T* dstField, CellInterval ci, unsigned char* inBuffer, Direction dir, gpuStream_t stream = nullptr) const
   {
      uint8_t timestep = dstField->getTimestep();
      kernels_.unpackDirection(dstField, ci, inBuffer, dir, timestep, stream);
   }

   void localCopyDirection(PdfField_T* srcField, CellInterval srcInterval, PdfField_T* dstField,
                           CellInterval dstInterval, Direction dir, gpuStream_t stream) const
   {
      uint8_t timestep = srcField->getTimestep();
      WALBERLA_ASSERT_EQUAL(timestep, dstField->getTimestep())
      kernels_.localCopyDirection(srcField, srcInterval, dstField, dstInterval, dir, timestep, stream);
   }

   void unpackRedistribute(PdfField_T* dstField, CellInterval ci, unsigned char* inBuffer,
                           stencil::Direction dir, gpuStream_t stream = nullptr) const
   {
      uint8_t timestep = dstField->getTimestep();
      kernels_.unpackRedistribute(dstField, ci, inBuffer, dir, timestep, stream);
   }

   void packPartialCoalescence(PdfField_T* srcField, PartialCoalescenceMaskFieldGPU* maskField, CellInterval ci,
                               unsigned char* outBuffer, Direction dir, gpuStream_t stream = nullptr) const
   {
      uint8_t timestep = srcField->getTimestep();
      kernels_.packPartialCoalescence(srcField, maskField, ci, outBuffer, dir, timestep, stream);
   }

   void zeroCoalescenceRegion(PdfField_T* dstField, CellInterval ci, Direction dir, gpuStream_t stream = nullptr) const
   {
      uint8_t timestep = dstField->getTimestep();
      kernels_.zeroCoalescenceRegion(dstField, ci, dir, timestep, stream);
   }

   void unpackCoalescence(PdfField_T* dstField, CellInterval ci, unsigned char* inBuffer, Direction dir, gpuStream_t stream = nullptr) const
   {
      uint8_t timestep = dstField->getTimestep();
      kernels_.unpackCoalescence(dstField, ci, inBuffer, dir, timestep, stream);
   }

   uint_t size(CellInterval ci, Direction dir) const { return kernels_.size(ci, dir); }
   uint_t size(CellInterval ci) const { return kernels_.size(ci); }
   uint_t redistributeSize(CellInterval ci) const { return kernels_.redistributeSize(ci); }
   uint_t partialCoalescenceSize(CellInterval ci, Direction dir) const
   {
      return kernels_.partialCoalescenceSize(ci, dir);
   }

 private:
   PackingKernels_T kernels_;
};
} // namespace internal

/***********************************************************************************************************************
 *                                                  Class Declaration                                                  *
 **********************************************************************************************************************/

template< typename PdfField_T >
class NonuniformGeneratedGPUPdfPackInfo : public walberla::gpu::GeneratedNonUniformGPUPackInfo
{
 public:
   using LatticeStorageSpecification_T = typename PdfField_T::LatticeStorageSpecification;
   using Stencil                       = typename LatticeStorageSpecification_T::Stencil;
   using CommunicationStencil          = typename LatticeStorageSpecification_T::CommunicationStencil;
   using CommData_T                    = NonuniformGPUCommData< LatticeStorageSpecification_T >;

   NonuniformGeneratedGPUPdfPackInfo(const BlockDataID pdfFieldID, const BlockDataID commDataID)
      : pdfFieldID_(pdfFieldID), commDataID_(commDataID){};

   bool constantDataExchange() const override { return true; };
   bool threadsafeReceiving() const override { return false; };

   /// Equal Level
   void unpackDataEqualLevel(Block* receiver, Direction dir, GpuBuffer_T& buffer, gpuStream_t stream) override;
   void communicateLocalEqualLevel(const Block* sender, Block* receiver, stencil::Direction dir,
                                   gpuStream_t stream) override;

   /// Coarse to Fine
   void unpackDataCoarseToFine(Block* fineReceiver, const BlockID& coarseSender, stencil::Direction dir,
                               GpuBuffer_T& buffer, gpuStream_t stream) override;
   void communicateLocalCoarseToFine(const Block* coarseSender, Block* fineReceiver, stencil::Direction dir, gpuStream_t stream) override;
   void communicateLocalCoarseToFine(const Block* coarseSender, Block* fineReceiver, stencil::Direction dir,
                                     GpuBuffer_T& buffer, gpuStream_t stream) override;

   /// Fine to Coarse
   void prepareCoalescence(Block* coarseReceiver, gpuStream_t gpuStream = nullptr);
   void unpackDataFineToCoarse(Block* coarseReceiver, const BlockID& fineSender, stencil::Direction dir,
                               GpuBuffer_T& buffer, gpuStream_t stream) override;

   void communicateLocalFineToCoarse(const Block* fineSender, Block* coarseReceiver, stencil::Direction dir, gpuStream_t stream) override;
   void communicateLocalFineToCoarse(const Block* fineSender, Block* coarseReceiver, stencil::Direction dir,
                                     GpuBuffer_T& buffer, gpuStream_t stream) override;

   uint_t sizeEqualLevelSend(const Block* sender, stencil::Direction dir) override;
   uint_t sizeCoarseToFineSend(const Block* coarseSender, const BlockID& fineReceiver, stencil::Direction dir) override;
   uint_t sizeCoarseToFineReceive ( Block* fineReceiver, stencil::Direction dir) override;
   uint_t sizeFineToCoarseSend(const Block* fineSender, stencil::Direction dir) override;

 protected:
   void packDataEqualLevelImpl(const Block* sender, stencil::Direction dir, GpuBuffer_T& buffer, gpuStream_t stream) const override;

   void packDataCoarseToFineImpl(const Block* coarseSender, const BlockID& fineReceiver, stencil::Direction dir,
                                 GpuBuffer_T& buffer, gpuStream_t stream) const override;
   void packDataFineToCoarseImpl(const Block* fineSender, const BlockID& coarseReceiver, stencil::Direction dir,
                                 GpuBuffer_T& buffer, gpuStream_t stream) const override;

 private:
   /// Helper Functions
   /// As in PdfFieldPackInfo.h
   Vector3< cell_idx_t > getNeighborShift(const BlockID& fineBlock, stencil::Direction dir) const;
   bool areNeighborsInDirection(const Block* block, const BlockID& neighborID,
                                Vector3< cell_idx_t > dirVec) const;

   CellInterval intervalHullInDirection(const CellInterval& ci, Vector3< cell_idx_t > tangentialDir,
                                        cell_idx_t width) const;
   bool skipsThroughCoarseBlock(const Block* block, Direction dir) const;

   void getCoarseBlockCommIntervals(const BlockID& fineBlockID, Direction dir, const PdfField_T* field,
                                    std::vector< std::pair< Direction, CellInterval > >& intervals) const;
   void getFineBlockCommIntervals(const BlockID& fineBlockID, Direction dir, const PdfField_T* field,
                                  std::vector< std::pair< Direction, CellInterval > >& intervals) const;

   CellInterval getCoarseBlockCoalescenceInterval(const Block* coarseBlock, const BlockID& fineBlockID, Direction dir,
                                                  const PdfField_T* field) const;

   const BlockDataID pdfFieldID_;
   internal::NonuniformGPUPackingKernelsWrapper< PdfField_T, LatticeStorageSpecification_T::inplace > kernels_;

 public:
   const BlockDataID commDataID_;
};

/***********************************************************************************************************************
 *                                                  Factory Functions                                                  *
 **********************************************************************************************************************/

template< typename PdfField_T >
std::shared_ptr< NonuniformGeneratedGPUPdfPackInfo< PdfField_T > >
   setupNonuniformGPUPdfCommunication(const std::weak_ptr< StructuredBlockForest >& blocks,
                                      BlockDataID pdfFieldID,
                                      const std::string& dataIdentifier = "NonuniformGPUCommData");

} // namespace walberla::lbm_generated

#include "lbm_generated/gpu/NonuniformGeneratedGPUPdfPackInfo.impl.h"
