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
//! \file UniformGeneratedGPUPdfPackInfo.h
//! \ingroup lbm
//! \author Markus Holzer <markus.holzer@fau.de>
//! \brief Class Template for Lattice Boltzmann PDF Pack Infos using code-generated kernels
//
//======================================================================================================================

#pragma once

#include "core/DataTypes.h"
#include "core/cell/CellInterval.h"

#include "gpu/GPUWrapper.h"
#include "gpu/communication/GeneratedGPUPackInfo.h"

#include "lbm/field/PdfField.h"

#include "stencil/Directions.h"

namespace walberla
{
using gpu::GeneratedGPUPackInfo;

namespace lbm_generated
{
using stencil::Direction;

namespace internal
{
/*
 * Base Template for Packing Kernels Wrapper. This wrapper is required for passing the time step to
 * kernels generated for in-place streaming patterns. The generated code should not be templated.
 */
template< typename PdfField_T, bool inplace >
class UniformPackingGPUKernelsWrapper
{
 public:
   void packAll(PdfField_T* srcField, CellInterval& ci, unsigned char* outBuffer, gpuStream_t stream) const  = 0;
   void unpackAll(PdfField_T* dstField, CellInterval& ci, unsigned char* inBuffer, gpuStream_t stream) const = 0;
   void localCopyAll(PdfField_T* srcField, CellInterval& srcInterval, PdfField_T* dstField, CellInterval& dstInterval,
                     gpuStream_t stream) const                                                               = 0;

   void packDirection(PdfField_T* srcField, CellInterval& ci, unsigned char* outBuffer, Direction dir,
                      gpuStream_t stream) const                                                = 0;
   void unpackDirection(PdfField_T* dstField, CellInterval& ci, unsigned char* inBuffer, Direction dir,
                        gpuStream_t stream) const                                              = 0;
   void localCopyDirection(PdfField_T* srcField, CellInterval& srcInterval, PdfField_T* dstField,
                           CellInterval& dstInterval, Direction dir, gpuStream_t stream) const = 0;

   uint_t size(CellInterval& ci, Direction dir) const = 0;
   uint_t size(CellInterval& ci) const                = 0;
};

/*
 * Template Specialization for two-fields patterns, with trivial method wrappers.
 */
template< typename PdfField_T >
class UniformPackingGPUKernelsWrapper< PdfField_T, false >
{
 public:
   using LatticeStorageSpecification_T = typename PdfField_T::LatticeStorageSpecification;
   using PackingKernels_T              = typename LatticeStorageSpecification_T::PackKernels;

   void packAll(PdfField_T* srcField, CellInterval& ci, unsigned char* outBuffer, gpuStream_t stream) const
   {
      kernels_.packAll(srcField, ci, outBuffer, stream);
   }

   void unpackAll(PdfField_T* dstField, CellInterval& ci, unsigned char* inBuffer, gpuStream_t stream) const
   {
      kernels_.unpackAll(dstField, ci, inBuffer, stream);
   }

   void localCopyAll(PdfField_T* srcField, CellInterval& srcInterval, PdfField_T* dstField, CellInterval& dstInterval,
                     gpuStream_t stream) const
   {
      kernels_.localCopyAll(srcField, srcInterval, dstField, dstInterval, stream);
   }

   void packDirection(PdfField_T* srcField, CellInterval& ci, unsigned char* outBuffer, Direction dir,
                      gpuStream_t stream) const
   {
      kernels_.packDirection(srcField, ci, outBuffer, dir, stream);
   }

   void unpackDirection(PdfField_T* dstField, CellInterval& ci, unsigned char* inBuffer, Direction dir,
                        gpuStream_t stream) const
   {
      kernels_.unpackDirection(dstField, ci, inBuffer, dir, stream);
   }

   void localCopyDirection(PdfField_T* srcField, CellInterval& srcInterval, PdfField_T* dstField,
                           CellInterval& dstInterval, Direction dir, gpuStream_t stream) const
   {
      kernels_.localCopyDirection(srcField, srcInterval, dstField, dstInterval, dir, stream);
   }

   uint_t size(CellInterval& ci, Direction dir) const { return kernels_.size(ci, dir); }
   uint_t size(CellInterval& ci) const { return kernels_.size(ci); }

 private:
   PackingKernels_T kernels_;
};

/*
 * Template Specialization for in-place patterns, extracting the timestep from the lattice model.
 */
template< typename PdfField_T >
class UniformPackingGPUKernelsWrapper< PdfField_T, true >
{
 public:
   using LatticeStorageSpecification_T = typename PdfField_T::LatticeStorageSpecification;
   using PackingKernels_T              = typename LatticeStorageSpecification_T::PackKernels;

   void packAll(PdfField_T* srcField, CellInterval& ci, unsigned char* outBuffer, gpuStream_t stream) const
   {
      uint8_t timestep = srcField->getTimestep();
      kernels_.packAll(srcField, ci, outBuffer, timestep, stream);
   }

   void unpackAll(PdfField_T* dstField, CellInterval& ci, unsigned char* inBuffer, gpuStream_t stream) const
   {
      uint8_t timestep = dstField->getTimestep();
      kernels_.unpackAll(dstField, ci, inBuffer, timestep, stream);
   }

   void localCopyAll(PdfField_T* srcField, CellInterval& srcInterval, PdfField_T* dstField, CellInterval& dstInterval,
                     gpuStream_t stream) const
   {
      uint8_t timestep = srcField->getTimestep();
      WALBERLA_ASSERT_EQUAL(timestep, dstField->getTimestep())
      kernels_.localCopyAll(srcField, srcInterval, dstField, dstInterval, timestep, stream);
   }

   void packDirection(PdfField_T* srcField, CellInterval& ci, unsigned char* outBuffer, Direction dir,
                      gpuStream_t stream) const
   {
      uint8_t timestep = srcField->getTimestep();
      kernels_.packDirection(srcField, ci, outBuffer, dir, timestep, stream);
   }

   void unpackDirection(PdfField_T* dstField, CellInterval& ci, unsigned char* inBuffer, Direction dir,
                        gpuStream_t stream) const
   {
      uint8_t timestep = dstField->getTimestep();
      kernels_.unpackDirection(dstField, ci, inBuffer, dir, timestep, stream);
   }

   void localCopyDirection(PdfField_T* srcField, CellInterval& srcInterval, PdfField_T* dstField,
                           CellInterval& dstInterval, Direction dir, gpuStream_t stream) const
   {
      uint8_t timestep = srcField->getTimestep();
      WALBERLA_ASSERT_EQUAL(timestep, dstField->getTimestep())
      kernels_.localCopyDirection(srcField, srcInterval, dstField, dstInterval, dir, timestep, stream);
   }

   uint_t size(CellInterval& ci, Direction dir) const { return kernels_.size(ci, dir); }
   uint_t size(CellInterval& ci) const { return kernels_.size(ci); }

 private:
   PackingKernels_T kernels_;
};
} // namespace internal

/**
 * Pack Info class template for lattice Boltzmann PDF fields. Relies on a code-generated
 * class providing kernel implementations for packing, unpacking and local copying of data.
 *
 * This template relies on a PackingKernels implementation generated by lbmpy_walberla.packing_kernels.
 * The code generated part provides the kernels for transferring data between communication buffers
 * and fields. The iteration slices are constructed by this class.
 *
 * The code-generated substructure enables the usage of arbitrary, in particular in-place streaming
 * patterns.
 *
 * @tparam  PackingKernels_T Type of a PackingKernels implementation generated using
 *          `lbmpy_walberla.generate_packing_kernels`.
 *
 * \ingroup lbm
 */
template< typename PdfField_T >
class UniformGeneratedGPUPdfPackInfo : public GeneratedGPUPackInfo
{
 public:
   using LatticeStorageSpecification_T = typename PdfField_T::LatticeStorageSpecification;
   using PackingKernels_T              = typename LatticeStorageSpecification_T::PackKernels;
   using Stencil                       = typename LatticeStorageSpecification_T::Stencil;

   UniformGeneratedGPUPdfPackInfo(const BlockDataID pdfFieldID, cell_idx_t cellLayersToSend = 1, bool sendAll = false)
      : pdfFieldID_(pdfFieldID), ghostLayersToSend_(cellLayersToSend), sendAll_(sendAll)
   {}

   void pack(stencil::Direction dir, unsigned char* buffer, IBlock* block, gpuStream_t stream) override;
   void communicateLocal(stencil::Direction dir, const IBlock* sender, IBlock* receiver, gpuStream_t stream) override;
   void unpack(stencil::Direction dir, unsigned char* buffer, IBlock* block, gpuStream_t stream) override;
   uint_t size(stencil::Direction dir, IBlock* block) override;

 private:
   const BlockDataID pdfFieldID_;
   internal::UniformPackingGPUKernelsWrapper< PdfField_T, LatticeStorageSpecification_T::inplace > kernels_;
   cell_idx_t ghostLayersToSend_;
   bool sendAll_;
};

template< typename PdfField_T >
void UniformGeneratedGPUPdfPackInfo< PdfField_T >::unpack(stencil::Direction dir, unsigned char* buffer, IBlock* block,
                                                          gpuStream_t stream)
{
   auto field = block->getData< PdfField_T >(pdfFieldID_);
   CellInterval ci;
   field->getGhostRegion(dir, ci, ghostLayersToSend_, false);

   if (sendAll_) { kernels_.unpackAll(field, ci, buffer, stream); }
   else { kernels_.unpackDirection(field, ci, buffer, dir, stream); }
}

template< typename PdfField_T >
void UniformGeneratedGPUPdfPackInfo< PdfField_T >::pack(stencil::Direction dir, unsigned char* buffer, IBlock* block,
                                                        gpuStream_t stream)
{
   auto field = const_cast< IBlock* >(block)->getData< PdfField_T >(pdfFieldID_);
   CellInterval ci;
   field->getSliceBeforeGhostLayer(dir, ci, ghostLayersToSend_, false);

   if (sendAll_) { kernels_.packAll(field, ci, buffer, stream); }
   else { kernels_.packDirection(field, ci, buffer, dir, stream); }
}

template< typename PdfField_T >
void UniformGeneratedGPUPdfPackInfo< PdfField_T >::communicateLocal(stencil::Direction dir, const IBlock* sender,
                                                                    IBlock* receiver, gpuStream_t stream)
{
   auto srcField = const_cast< IBlock* >(sender)->getData< PdfField_T >(pdfFieldID_);
   auto dstField = receiver->getData< PdfField_T >(pdfFieldID_);

   CellInterval srcRegion;
   CellInterval dstRegion;
   srcField->getSliceBeforeGhostLayer(dir, srcRegion, ghostLayersToSend_, false);
   dstField->getGhostRegion(stencil::inverseDir[dir], dstRegion, ghostLayersToSend_, false);

   if (sendAll_) { kernels_.localCopyAll(srcField, srcRegion, dstField, dstRegion, stream); }
   else { kernels_.localCopyDirection(srcField, srcRegion, dstField, dstRegion, dir, stream); }
}

template< typename PdfField_T >
uint_t UniformGeneratedGPUPdfPackInfo< PdfField_T >::size(stencil::Direction dir, IBlock* block)
{
   auto field = block->getData< PdfField_T >(pdfFieldID_);
   CellInterval ci;
   field->getGhostRegion(dir, ci, 1, false);

   uint_t elementsPerCell = kernels_.size(ci, dir);
   return elementsPerCell;
}

} // namespace lbm_generated
} // namespace walberla