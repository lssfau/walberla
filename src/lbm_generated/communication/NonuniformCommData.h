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
//! \file NonuniformCommData.h
//! \author Frederik Hennig <frederik.hennig@fau.de>
//! \author Markus Holzer <markus.holzer@fau.de>
//
//======================================================================================================================

#pragma once

#include "blockforest/StructuredBlockForest.h"
#include "blockforest/BlockDataHandling.h"

#include "domain_decomposition/IBlock.h"

#include "field/FlagField.h"

#include "lbm_generated/field/PdfField.h"

#include "stencil/Directions.h"

#define USE_CELL_INTERVALS

namespace walberla::lbm_generated {

using PartialCoalescenceMaskField = FlagField< uint32_t >;

namespace util {
   void forEachSubdirection(const Vector3< cell_idx_t > mainDirection, const std::function< void(Vector3< cell_idx_t >) >& func);
   bool forEachSubdirectionCancel(const Vector3< cell_idx_t > mainDirection,
                               const std::function< bool(Vector3< cell_idx_t >) >& func);
   void getSubdirections(const Vector3< cell_idx_t > mainDirection, std::vector< Vector3< cell_idx_t > > subdirs);

   template< typename Stencil_T >
   void forEachOrthogonalDirection(Vector3<cell_idx_t> d, std::function< void(Vector3< cell_idx_t >) > func);
} // namespace util

template< typename LatticeStorageSpecification_T >
class NonuniformCommData
{
 private:
   void registerFlags();
   void computeBitMask();

 public:
   using Stencil              = typename LatticeStorageSpecification_T::Stencil;
   using CommunicationStencil = typename LatticeStorageSpecification_T::CommunicationStencil;

#if defined(USE_CELL_INTERVALS)
   NonuniformCommData(IBlock* const block, uint_t xSize, uint_t ySize, uint_t zSize)
      : block_(block), maskField_(xSize, ySize, zSize, 2),
        interiorInterval(0, 0, 0, cell_idx_c(xSize) - 1, cell_idx_c(ySize) - 1, cell_idx_c(zSize) - 1)
   {
      registerFlags();
      computeBitMask();
   };
#else
   NonuniformCommData(IBlock* const block, const BlockDataID pdfFieldID, uint_t xSize, uint_t ySize, uint_t zSize)
      : block_(block), pdfFieldID_(pdfFieldID), maskField_(xSize, ySize, zSize, 2)
   {
      registerFlags();
      computeBitMask();
   };
#endif

   bool operator==(const NonuniformCommData& other) { return this == &other; }
   bool operator!=(const NonuniformCommData& other) { return this != &other; }

   PartialCoalescenceMaskField& getMaskField() { return maskField_; }
   const PartialCoalescenceMaskField& getMaskField() const { return maskField_; }

 private:
#if defined(USE_CELL_INTERVALS)
   void prepareIntervals();
   void setFlagOnInterval(const CellInterval & ci, const uint_t fIdx);
#else
   void prepareFlags();
   void resetCornerSkippingOriginFlags();
#endif

   void setupCornerSkippingOrigins(stencil::Direction commDir);
   void setupBitMaskSlice(stencil::Direction commDir, stencil::Direction streamDir);

   bool haveSmallestIdInIntersection(Vector3<cell_idx_t> cornerDir);

   const IBlock* const block_;
   PartialCoalescenceMaskField maskField_;

#if defined(USE_CELL_INTERVALS)
   const CellInterval interiorInterval;
   std::vector< CellInterval > passThroughIntervals_;
   std::vector< CellInterval > cornerSkippingOriginIntervals_;
#endif
};


template< typename LatticeStorageSpecification_T >
class NonuniformCommDataHandling
   : public blockforest::AlwaysInitializeBlockDataHandling< NonuniformCommData< LatticeStorageSpecification_T > >
{
 public:
   using CommmData_T = NonuniformCommData< LatticeStorageSpecification_T >;

   NonuniformCommDataHandling(const weak_ptr< StructuredBlockForest >& blocks)
      : blocks_(blocks){};

   CommmData_T* initialize(IBlock* const block) override
   {
      WALBERLA_ASSERT_NOT_NULLPTR(block)
      auto blocks = blocks_.lock();
      WALBERLA_CHECK_NOT_NULLPTR(blocks)

      return new CommmData_T(block, blocks->getNumberOfXCells(*block), blocks->getNumberOfYCells(*block),
                             blocks->getNumberOfZCells(*block));
   }

 private:
   const weak_ptr< StructuredBlockStorage > blocks_;
};

} // walberla::lbm_generated

#include "lbm_generated/communication/NonuniformCommData.impl.h"
