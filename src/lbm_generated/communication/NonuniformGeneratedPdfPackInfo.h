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
//! \file NonuniformGeneratedPdfPackInfo.h
//! \author Frederik Hennig <frederik.hennig@fau.de>
//! \author Markus Holzer <markus.holzer@fau.de>
//
//======================================================================================================================

#pragma once

#include "blockforest/communication/NonUniformPackInfo.h"

#include "core/DataTypes.h"
#include "core/mpi/RecvBuffer.h"
#include "core/mpi/SendBuffer.h"

#include "lbm_generated/communication/NonuniformCommData.h"
#include "lbm_generated/field/PdfField.h"

namespace walberla::lbm_generated {
using stencil::Direction;

namespace internal
{
/*
 * Base Template for Packing Kernels Wrapper. This wrapper is required for passing the time step to
 * kernels generated for in-place streaming patterns. The generated code should not be templated.
 */
template< typename PdfField_T, bool inplace >
class NonuniformPackingKernelsWrapper
{
 public:

   void packAll(PdfField_T* srcField, CellInterval& ci, unsigned char* outBuffer) const  = 0;
   void unpackAll(PdfField_T* dstField, CellInterval& ci, unsigned char* inBuffer) const = 0;
   void localCopyAll(PdfField_T* srcField, CellInterval& srcInterval, PdfField_T* dstField,
                     CellInterval& dstInterval) const                                    = 0;

   void packDirection(PdfField_T* srcField, CellInterval& ci, unsigned char* outBuffer, Direction dir) const  = 0;
   void unpackDirection(PdfField_T* dstField, CellInterval& ci, unsigned char* inBuffer, Direction dir) const = 0;
   void localCopyDirection(PdfField_T* srcField, CellInterval& srcInterval, PdfField_T* dstField,
                           CellInterval& dstInterval, Direction dir) const                                    = 0;

   void unpackRedistribute(PdfField_T* dstField, CellInterval& ci, unsigned char* inBuffer,
                           stencil::Direction dir) const = 0;

   void packPartialCoalescence(PdfField_T* srcField, PartialCoalescenceMaskField* maskField, CellInterval& ci,
                               unsigned char* outBuffer, Direction dir) const = 0;
   void zeroCoalescenceRegion(PdfField_T* dstField, CellInterval& ci, Direction dir) const = 0;
   void unpackCoalescence(PdfField_T* dstField, CellInterval& ci, unsigned char* inBuffer, Direction dir) const = 0;

   uint_t size(CellInterval& ci, Direction dir) const                   = 0;
   uint_t size(CellInterval& ci) const                                  = 0;
   uint_t redistributeSize(CellInterval& ci) const                      = 0;
   uint_t partialCoalescenceSize(CellInterval& ci, Direction dir) const = 0;
};

/*
 * Template Specialization for two-fields patterns, with trivial method wrappers.
 */
template< typename PdfField_T >
class NonuniformPackingKernelsWrapper< PdfField_T, false >
{
 public:
   using LatticeStorageSpecification_T = typename PdfField_T::LatticeStorageSpecification;
   using PackingKernels_T = typename LatticeStorageSpecification_T::PackKernels;

   void packAll(PdfField_T* srcField, CellInterval& ci, unsigned char* outBuffer) const
   {
      kernels_.packAll(srcField, ci, outBuffer);
   }

   void unpackAll(PdfField_T* dstField, CellInterval& ci, unsigned char* inBuffer) const
   {
      kernels_.unpackAll(dstField, ci, inBuffer);
   }

   void localCopyAll(PdfField_T* srcField, CellInterval& srcInterval, PdfField_T* dstField,
                     CellInterval& dstInterval) const
   {
      kernels_.localCopyAll(srcField, srcInterval, dstField, dstInterval);
   }

   void packDirection(PdfField_T* srcField, CellInterval& ci, unsigned char* outBuffer, Direction dir) const
   {
      kernels_.packDirection(srcField, ci, outBuffer, dir);
   }

   void unpackDirection(PdfField_T* dstField, CellInterval& ci, unsigned char* inBuffer, Direction dir) const
   {
      kernels_.unpackDirection(dstField, ci, inBuffer, dir);
   }

   void localCopyDirection(PdfField_T* srcField, CellInterval& srcInterval, PdfField_T* dstField,
                           CellInterval& dstInterval, Direction dir) const
   {
      kernels_.localCopyDirection(srcField, srcInterval, dstField, dstInterval, dir);
   }

   void unpackRedistribute(PdfField_T* dstField, CellInterval& ci, unsigned char* inBuffer,
                           stencil::Direction dir) const
   {
      kernels_.unpackRedistribute(dstField, ci, inBuffer, dir);
   }

   void packPartialCoalescence(PdfField_T* srcField, PartialCoalescenceMaskField* maskField, CellInterval& ci,
                               unsigned char* outBuffer, Direction dir) const
   {
      kernels_.packPartialCoalescence(srcField, maskField, ci, outBuffer, dir);
   }

   void unpackCoalescence(PdfField_T* dstField, CellInterval& ci, unsigned char* inBuffer, Direction dir) const
   {
      kernels_.unpackCoalescence(dstField, ci, inBuffer, dir);
   }

   void zeroCoalescenceRegion(PdfField_T* dstField, CellInterval& ci, Direction dir) const
   {
      kernels_.zeroCoalescenceRegion(dstField, ci, dir);
   }

   uint_t size(CellInterval& ci, Direction dir) const { return kernels_.size(ci, dir); }
   uint_t size(CellInterval& ci) const { return kernels_.size(ci); }
   uint_t redistributeSize(CellInterval& ci) const { return kernels_.redistributeSize(ci); }
   uint_t partialCoalescenceSize(CellInterval& ci, Direction dir) const
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
class NonuniformPackingKernelsWrapper< PdfField_T, true >
{
 public:
   using LatticeStorageSpecification_T = typename PdfField_T::LatticeStorageSpecification;
   using PackingKernels_T = typename LatticeStorageSpecification_T::PackKernels;

   void packAll(PdfField_T* srcField, CellInterval& ci, unsigned char* outBuffer) const
   {
      uint8_t timestep = srcField->getTimestep();
      kernels_.packAll(srcField, ci, outBuffer, timestep);
   }

   void unpackAll(PdfField_T* dstField, CellInterval& ci, unsigned char* inBuffer) const
   {
      uint8_t timestep = dstField->getTimestep();
      kernels_.unpackAll(dstField, ci, inBuffer, timestep);
   }

   void localCopyAll(PdfField_T* srcField, CellInterval& srcInterval, PdfField_T* dstField,
                     CellInterval& dstInterval) const
   {
      uint8_t timestep = srcField->getTimestep();
      WALBERLA_ASSERT_EQUAL(timestep, dstField->getTimestep())
      kernels_.localCopyAll(srcField, srcInterval, dstField, dstInterval, timestep);
   }

   void packDirection(PdfField_T* srcField, CellInterval& ci, unsigned char* outBuffer, Direction dir) const
   {
      uint8_t timestep = srcField->getTimestep();
      kernels_.packDirection(srcField, ci, outBuffer, dir, timestep);
   }

   void unpackDirection(PdfField_T* dstField, CellInterval& ci, unsigned char* inBuffer, Direction dir) const
   {
      uint8_t timestep = dstField->getTimestep();
      kernels_.unpackDirection(dstField, ci, inBuffer, dir, timestep);
   }

   void localCopyDirection(PdfField_T* srcField, CellInterval& srcInterval, PdfField_T* dstField,
                           CellInterval& dstInterval, Direction dir) const
   {
      uint8_t timestep = srcField->getTimestep();
      WALBERLA_ASSERT_EQUAL(timestep, dstField->getTimestep())
      kernels_.localCopyDirection(srcField, srcInterval, dstField, dstInterval, dir, timestep);
   }

   void unpackRedistribute(PdfField_T* dstField, CellInterval& ci, unsigned char* inBuffer,
                           stencil::Direction dir) const
   {
      uint8_t timestep = dstField->getTimestep();
      kernels_.unpackRedistribute(dstField, ci, inBuffer, dir, timestep);
   }

   void packPartialCoalescence(PdfField_T* srcField, PartialCoalescenceMaskField* maskField, CellInterval& ci,
                               unsigned char* outBuffer, Direction dir) const
   {
      uint8_t timestep = srcField->getTimestep();
      kernels_.packPartialCoalescence(srcField, maskField, ci, outBuffer, dir, timestep);
   }

   void zeroCoalescenceRegion(PdfField_T* dstField, CellInterval& ci, Direction dir) const
   {
      uint8_t timestep = dstField->getTimestep();
      kernels_.zeroCoalescenceRegion(dstField, ci, dir, timestep);
   }

   void unpackCoalescence(PdfField_T* dstField, CellInterval& ci, unsigned char* inBuffer, Direction dir) const
   {
      uint8_t timestep = dstField->getTimestep();
      kernels_.unpackCoalescence(dstField, ci, inBuffer, dir, timestep);
   }

   uint_t size(CellInterval& ci, Direction dir) const { return kernels_.size(ci, dir); }
   uint_t size(CellInterval& ci) const { return kernels_.size(ci); }
   uint_t redistributeSize(CellInterval& ci) const { return kernels_.redistributeSize(ci); }
   uint_t partialCoalescenceSize(CellInterval& ci, Direction dir) const
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
class NonuniformGeneratedPdfPackInfo : public blockforest::communication::NonUniformPackInfo
{
 public:
   using LatticeStorageSpecification_T = typename PdfField_T::LatticeStorageSpecification;
   using PackingKernels_T = typename LatticeStorageSpecification_T::PackKernels;
   using Stencil      = typename LatticeStorageSpecification_T::Stencil;
   using CommunicationStencil = typename LatticeStorageSpecification_T::CommunicationStencil;
   using CommData_T           = NonuniformCommData< LatticeStorageSpecification_T >;


   NonuniformGeneratedPdfPackInfo(const BlockDataID pdfFieldID, const BlockDataID commDataID)
      : pdfFieldID_(pdfFieldID), commDataID_(commDataID){};

   bool constantDataExchange() const override { return true; };
   bool threadsafeReceiving() const override { return false; };

   /// Equal Level
   void unpackDataEqualLevel(Block* receiver, Direction dir, mpi::RecvBuffer& buffer) override;
   void communicateLocalEqualLevel(const Block* sender, Block* receiver, stencil::Direction dir) override;

   /// Coarse to Fine
   void unpackDataCoarseToFine(Block* fineReceiver, const BlockID& coarseSender, stencil::Direction dir,
                               mpi::RecvBuffer& buffer) override;
   void communicateLocalCoarseToFine(const Block* coarseSender, Block* fineReceiver, stencil::Direction dir) override;

   /// Fine to Coarse
   void prepareCoalescence(Block* coarseReceiver);
   void unpackDataFineToCoarse(Block* coarseReceiver, const BlockID& fineSender, stencil::Direction dir,
                               mpi::RecvBuffer& buffer) override;

   void communicateLocalFineToCoarse(const Block* fineSender, Block* coarseReceiver, stencil::Direction dir) override;

 protected:
   void packDataEqualLevelImpl(const Block* sender, stencil::Direction dir, mpi::SendBuffer& buffer) const override;

   void packDataCoarseToFineImpl(const Block* coarseSender, const BlockID& fineReceiver, stencil::Direction dir,
                                 mpi::SendBuffer& buffer) const override;
   void packDataFineToCoarseImpl(const Block* fineSender, const BlockID& coarseReceiver, stencil::Direction dir,
                                 mpi::SendBuffer& buffer) const override;

 private:
   /// Helper Functions
   /// As in PdfFieldPackInfo.h
   Vector3< cell_idx_t > getNeighborShift(const BlockID& fineBlock, stencil::Direction dir) const;
   bool areNeighborsInDirection(const Block * block, const BlockID & neighborID, const Vector3< cell_idx_t> dirVec) const;

   CellInterval intervalHullInDirection(const CellInterval& ci, const Vector3< cell_idx_t > tangentialDir,
                                        cell_idx_t width) const;
   bool skipsThroughCoarseBlock(const Block* block, const Direction dir) const;

   void getCoarseBlockCommIntervals(const BlockID& fineBlockID, const Direction dir, const PdfField_T* field,
                                    std::vector< std::pair< Direction, CellInterval > >& intervals) const;
   void getFineBlockCommIntervals(const BlockID& fineBlockID, const Direction dir, const PdfField_T* field,
                                  std::vector< std::pair< Direction, CellInterval > >& intervals) const;

   CellInterval getCoarseBlockCoalescenceInterval(const Block * coarseBlock, const BlockID & fineBlockID,
                                                  Direction dir, const PdfField_T * field) const;

   const BlockDataID pdfFieldID_;
   internal::NonuniformPackingKernelsWrapper< PdfField_T, LatticeStorageSpecification_T::inplace > kernels_;

 public:
   const BlockDataID commDataID_;
};

/***********************************************************************************************************************
 *                                                  Factory Functions                                                  *
 **********************************************************************************************************************/

template< typename PdfField_T>
std::shared_ptr< NonuniformGeneratedPdfPackInfo< PdfField_T > >
   setupNonuniformPdfCommunication(const std::weak_ptr< StructuredBlockForest >& blocks, const BlockDataID pdfFieldID,
                                   const std::string& dataIdentifier = "NonuniformCommData");

} // walberla::lbm_generated

#include "lbm_generated/communication/NonuniformGeneratedPdfPackInfo.impl.h"
