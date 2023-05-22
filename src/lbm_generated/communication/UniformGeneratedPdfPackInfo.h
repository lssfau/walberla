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
//! \file UniformGeneratedPdfPackInfo.h
//! \ingroup lbm
//! \author Frederik Hennig <frederik.hennig@fau.de>
//! \brief Class Template for Lattice Boltzmann PDF Pack Infos using code-generated kernels
//
//======================================================================================================================

#pragma once

#include "communication/UniformPackInfo.h"

#include "core/DataTypes.h"
#include "core/cell/CellInterval.h"

#include "lbm/field/PdfField.h"

#include "stencil/Directions.h"

namespace walberla
{
using communication::UniformPackInfo;

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
class UniformPackingKernelsWrapper
{
 public:

   void packAll(PdfField_T* srcField, CellInterval& ci, unsigned char* outBuffer) const  = 0;
   void unpackAll(PdfField_T* dstField, CellInterval& ci, unsigned char* inBuffer) const = 0;
   void localCopyAll(PdfField_T* srcField, CellInterval& srcInterval, PdfField_T* dstField,
                     CellInterval& dstInterval) const                                    = 0;

   void packDirection(PdfField_T* srcField, CellInterval& ci, unsigned char* outBuffer, const Direction dir) const  = 0;
   void unpackDirection(PdfField_T* dstField, CellInterval& ci, unsigned char* inBuffer, const Direction dir) const = 0;
   void localCopyDirection(PdfField_T* srcField, CellInterval& srcInterval, PdfField_T* dstField,
                           CellInterval& dstInterval, const Direction dir) const                                    = 0;

   uint_t size(CellInterval& ci, const Direction dir) const = 0;
   uint_t size(CellInterval& ci) const                = 0;
};

/*
 * Template Specialization for two-fields patterns, with trivial method wrappers.
 */
template< typename PdfField_T >
class UniformPackingKernelsWrapper< PdfField_T, false >
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

   void packDirection(PdfField_T* srcField, CellInterval& ci, unsigned char* outBuffer, const Direction dir) const
   {
      kernels_.packDirection(srcField, ci, outBuffer, dir);
   }

   void unpackDirection(PdfField_T* dstField, CellInterval& ci, unsigned char* inBuffer, const Direction dir) const
   {
      kernels_.unpackDirection(dstField, ci, inBuffer, dir);
   }

   void localCopyDirection(PdfField_T* srcField, CellInterval& srcInterval, PdfField_T* dstField,
                           CellInterval& dstInterval, const Direction dir) const
   {
      kernels_.localCopyDirection(srcField, srcInterval, dstField, dstInterval, dir);
   }

   uint_t size(CellInterval& ci, const Direction dir) const { return kernels_.size(ci, dir); }
   uint_t size(CellInterval& ci) const { return kernels_.size(ci); }

 private:
   PackingKernels_T kernels_;
};

/*
 * Template Specialization for in-place patterns, extracting the timestep from the lattice model.
 */
template< typename PdfField_T >
class UniformPackingKernelsWrapper< PdfField_T, true >
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

   void packDirection(PdfField_T* srcField, CellInterval& ci, unsigned char* outBuffer, const Direction dir) const
   {
      uint8_t timestep = srcField->getTimestep();
      kernels_.packDirection(srcField, ci, outBuffer, dir, timestep);
   }

   void unpackDirection(PdfField_T* dstField, CellInterval& ci, unsigned char* inBuffer, const Direction dir) const
   {
      uint8_t timestep = dstField->getTimestep();
      kernels_.unpackDirection(dstField, ci, inBuffer, dir, timestep);
   }

   void localCopyDirection(PdfField_T* srcField, CellInterval& srcInterval, PdfField_T* dstField,
                           CellInterval& dstInterval, const Direction dir) const
   {
      uint8_t timestep = srcField->getTimestep();
      WALBERLA_ASSERT_EQUAL(timestep, dstField->getTimestep())
      kernels_.localCopyDirection(srcField, srcInterval, dstField, dstInterval, dir, timestep);
   }

   uint_t size(CellInterval& ci, const Direction dir) const { return kernels_.size(ci, dir); }
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
class UniformGeneratedPdfPackInfo : public UniformPackInfo
{
 public:
   using LatticeStorageSpecification_T = typename PdfField_T::LatticeStorageSpecification;
   using PackingKernels_T = typename LatticeStorageSpecification_T::PackKernels;
   using Stencil      = typename LatticeStorageSpecification_T::Stencil;

   /**
    * Constructor.
    *
    * @param pdfFieldID ID of the associated walberla::lbm::PdfField
    * @param cellLayersToSend The amount of cell layers that should be communicated
    * @param sendAll If true, instead of only those populations streaming in subdirections of the communication
    *                direction, all populations will always be communicated.
    *                \warning Be careful when using this option with any streaming pattern other than
    *                the pull pattern. Other patterns store at least some of their post-collision
    *                populations in neighbouring cells. This might lead to out-of-bounds errors when
    *                copying to the outermost ghost layer! Solve this by adding an additional ghost layer
    *                as a safety margin.
    */
   UniformGeneratedPdfPackInfo(const BlockDataID pdfFieldID, cell_idx_t cellLayersToSend = 1, bool sendAll = false)
      : pdfFieldID_(pdfFieldID), ghostLayersToSend_(cellLayersToSend), sendAll_(sendAll)
   {}

   bool constantDataExchange() const override { return true; }
   bool threadsafeReceiving() const override { return true; }

   void unpackData(IBlock * receiver, Direction dir, mpi::RecvBuffer & buffer) override;
   void communicateLocal(const IBlock * sender, IBlock * receiver, Direction dir) override;

 protected:
   void packDataImpl(const IBlock * sender, Direction dir, mpi::SendBuffer & buffer) const override;

 private:
   const BlockDataID pdfFieldID_;
   internal::UniformPackingKernelsWrapper< PdfField_T, LatticeStorageSpecification_T::inplace > kernels_;
   cell_idx_t ghostLayersToSend_;
   bool sendAll_;
};

template< typename PdfField_T >
void UniformGeneratedPdfPackInfo< PdfField_T >::unpackData( IBlock * receiver, Direction dir, mpi::RecvBuffer& buffer)
{
   auto field = receiver->getData< PdfField_T >(pdfFieldID_);
   CellInterval ci;
   field->getGhostRegion(dir, ci, ghostLayersToSend_, false);

   if (sendAll_)
   {
      unsigned char* bufferPtr = buffer.skip(kernels_.size(ci));
      kernels_.unpackAll(field, ci, bufferPtr);
   }
   else
   {
      uint_t size              = kernels_.size(ci, dir);
      unsigned char* bufferPtr = buffer.skip(size);
      kernels_.unpackDirection(field, ci, bufferPtr, dir);
   }
}

template< typename PdfField_T >
void UniformGeneratedPdfPackInfo< PdfField_T >::communicateLocal(const IBlock* sender, IBlock* receiver, Direction dir)
{
   auto srcField = const_cast< IBlock* >(sender)->getData< PdfField_T >(pdfFieldID_);
   auto dstField = receiver->getData< PdfField_T >(pdfFieldID_);

   CellInterval srcRegion;
   CellInterval dstRegion;
   srcField->getSliceBeforeGhostLayer(dir, srcRegion, ghostLayersToSend_, false);
   dstField->getGhostRegion(stencil::inverseDir[dir], dstRegion, ghostLayersToSend_, false);

   if (sendAll_) {
      kernels_.localCopyAll(srcField, srcRegion, dstField, dstRegion);
   }
   else
   {
      kernels_.localCopyDirection(srcField, srcRegion, dstField, dstRegion, dir);
   }
}

template< typename PdfField_T>
void UniformGeneratedPdfPackInfo< PdfField_T >:: packDataImpl(const IBlock* sender, Direction dir, mpi::SendBuffer& buffer) const
{
   auto field = const_cast< IBlock* >(sender)->getData< PdfField_T >(pdfFieldID_);
   CellInterval ci;
   field->getSliceBeforeGhostLayer(dir, ci, ghostLayersToSend_, false);

   if (sendAll_)
   {
      unsigned char* bufferPtr = buffer.forward(kernels_.size(ci));
      kernels_.packAll(field, ci, bufferPtr);
   }
   else
   {
      unsigned char* bufferPtr = buffer.forward(kernels_.size(ci, dir));
      kernels_.packDirection(field, ci, bufferPtr, dir);
   }
}

} // namespace lbm
} // namespace walberla