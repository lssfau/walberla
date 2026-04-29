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
//! \file Field.hpp
//! \author Frederik Hennig <frederik.hennig@fau.de>
//
//======================================================================================================================

#pragma once

#include "blockforest/StructuredBlockForest.h"

#include "core/Stdlib.hpp"

#include "field/GhostRegions.h"
#include "field/Layout.h"
#include WALBERLA_STDLIB(array)
#include WALBERLA_STDLIB(span)

#include <concepts>
#include <cstdint>

#include "./BufferSystem.hpp"
#include "./Initializers.hpp"
#include "./MemoryTags.hpp"

namespace walberla::v8
{
namespace memory
{

/**
 * @defgroup v8core-fields Fields
 * @ingroup v8core-memory
 * @brief Data structures for distributed numerical fields
 */

/**
 * @brief Iteration slices and cell intervals over a field
 * @ingroup v8core-fields
 */
class FieldSlices
{
 public:
   FieldSlices(std::array< size_t, 3 > fieldShape, size_t numGhostLayers)
      : shape_{ fieldShape }, numGhostLayers_{ numGhostLayers }
   {}

   /**
    * @brief Get the interior cell interval for each of a fields' per-block views
    */
   CellInterval interior() const
   {
      return { { 0, 0, 0 },
               {
                  cell_idx_c(shape_[0] - 1),
                  cell_idx_c(shape_[1] - 1),
                  cell_idx_c(shape_[2] - 1),
               } };
   }

   /**
    * @brief Get the full cell interval for each of a fields' per-block views, including the interior and all ghost
    * layers
    */
   CellInterval withGhostLayers() const { return withGhostLayers(numGhostLayers_); }

   /**
    * @brief Get an expanded cell interval for each of a fields' per-block views, including the interior and the given
    * number of ghost layers
    */
   CellInterval withGhostLayers(size_t gls) const
   {
      WALBERLA_DEBUG_SECTION()
      {
         if (gls > numGhostLayers_)
         {
            throw std::range_error{ "Given number of ghost layers exceeds ghost layers of field." };
         }
      }

      return { { -cell_idx_c(gls), -cell_idx_c(gls), -cell_idx_c(gls) },
               {
                  cell_idx_c(shape_[0] + gls - 1),
                  cell_idx_c(shape_[1] + gls - 1),
                  cell_idx_c(shape_[2] + gls - 1),
               } };
   }

   /**
    * @brief Get a cell interval for a ghost layer slice in a certain direction
    * 
    * @param dir The direction in which the ghost-layer slice should be taken
    * @param thickness The number of ghost layers to include
    * @param fullSlice If `true`, fully expand the slice in normal directions to include all ghost layers
    */
   CellInterval ghostSlice(stencil::Direction dir, size_t thickness = 1, bool fullSlice = false) const
   {
      std::array< cell_idx_t, 3 > fieldShape{ cell_idx_c(shape_[0]), cell_idx_c(shape_[1]), cell_idx_c(shape_[2]) };
      return field::detail::getGhostRegion(fieldShape, numGhostLayers_, dir, cell_idx_c(thickness), fullSlice);
   }

   /**
    * @brief Get a cell interval for a border slice in a certain direction
    * 
    * @param dir The direction in which the border slice should be taken
    * @param thickness The number of cell layers to include
    * @param fullSlice If `true`, fully expand the slice in normal directions to include all ghost layers
    */
   CellInterval borderSlice(stencil::Direction dir, size_t thickness = 1, bool fullSlice = false) const
   {
      std::array< cell_idx_t, 3 > fieldShape{ cell_idx_c(shape_[0]), cell_idx_c(shape_[1]), cell_idx_c(shape_[2]) };
      return field::detail::getSliceBeforeGhostLayer(fieldShape, numGhostLayers_, dir, cell_idx_c(thickness),
                                                     fullSlice);
   }

 private:
   std::array< size_t, 3 > shape_;
   size_t numGhostLayers_;
};

/**
 * @brief Primary distributed 3D numerical field data structure.
 * @ingroup v8core-fields
 */
template< typename TElement, size_t fSize, MemTag TMemTag >
   requires(std::semiregular< std::remove_cv_t< TElement > >)
class Field
{
 public:
   using Self = Field< TElement, fSize, TMemTag >;

   using element_type = TElement;
   using value_type   = std::remove_cv_t< element_type >;

   using memory_tag = TMemTag;

   using BufferSystemType = BufferSystem< element_type, 4, memory_tag >;

   static constexpr size_t F_SIZE = fSize;
   static constexpr bool IS_CONST = std::is_const_v< element_type >;

   // Constructors
   template< typename = void >
      requires(!IS_CONST)
   Field(StructuredBlockForest& blocks, uint_t numGhostLayers = 1, value_type initValue = value_type{},
         field::Layout layout = field::Layout::fzyx)
      : numGhostLayers_{ numGhostLayers }, layout_{ layout }, bufferSystem_{ createBufferSystem(blocks) }
   {
      auto init = initializers::BufferInitializer_T< memory_tag, value_type >{ initValue };

      for (auto& b : blocks)
      {
         auto bView = bufferSystem_.view(b);
         init(bView);
      }
   }

   template< typename = void >
      requires(IS_CONST)
   Field(const Field< value_type, F_SIZE, memory_tag >& other)
      : numGhostLayers_{ other.numGhostLayers() }, layout_{ other.layout() }, bufferSystem_{ other.bufferSystem() }
   {}

   size_t numGhostLayers() const { return numGhostLayers_; }

   field::Layout layout() const { return layout_; }

   const BufferSystemType& bufferSystem() const { return bufferSystem_; }

   FieldSlices slices() const
   {
      const auto& indexing = bufferSystem_.indexing();
      std::array< size_t, 3 > fieldShape{ indexing.shape()[0], indexing.shape()[1], indexing.shape()[2] };
      return { fieldShape, numGhostLayers_ };
   }

 private:
   size_t numGhostLayers_;
   field::Layout layout_;
   BufferSystemType bufferSystem_;

   BufferSystemType createBufferSystem(StructuredBlockForest& blocks);
};

template< typename TField >
concept IField =
   std::same_as< TField, Field< typename TField::element_type, TField::F_SIZE, typename TField::memory_tag > >;

template< typename TElement, size_t fSize, MemTag TMemTag >
   requires(std::semiregular< std::remove_cv_t< TElement > >)
inline auto Field< TElement, fSize, TMemTag >::createBufferSystem(StructuredBlockForest& blocks) -> BufferSystemType
{
   constexpr size_t RANK{ 4 };

   const stdlib::array< size_t, RANK > shape{
      blocks.getNumberOfXCellsPerBlock(),
      blocks.getNumberOfYCellsPerBlock(),
      blocks.getNumberOfZCellsPerBlock(),
      F_SIZE,
   };

   const stdlib::array< size_t, RANK > allocShape{
      shape[0] + 2 * numGhostLayers_,
      shape[1] + 2 * numGhostLayers_,
      shape[2] + 2 * numGhostLayers_,
      F_SIZE,
   };

   const stdlib::array< size_t, RANK > offsets{
      numGhostLayers_,
      numGhostLayers_,
      numGhostLayers_,
      0,
   };

   stdlib::array< size_t, RANK > linearizationOrder{};
   bool addPaddingForAlignment{ false };

   switch (layout_)
   {
   case field::fzyx: {
      linearizationOrder     = { 0, 1, 2, 3 };
      addPaddingForAlignment = true;
      break;
   }

   case field::zyxf: {
      linearizationOrder     = { 3, 0, 1, 2 };
      addPaddingForAlignment = false;
      break;
   }
   }

   return BufferSystemType::create(blocks, shape)
      .allocShape(allocShape)
      .allocOffsets(offsets)
      .linearizationOrder(linearizationOrder)
      .addLinePaddingForAlignment(addPaddingForAlignment)
      .build();
}

/**
 * @brief View into a field's block-local portion
 * @ingroup v8core-fields
 */
template< typename TElement, size_t fSize, MemTag TMemTag >
   requires(std::semiregular< std::remove_cv_t< TElement > >)
class FieldView
{
 public:
   using element_type   = TElement;
   using value_type     = std::remove_cv_t< element_type >;
   using pointer_type   = element_type*;
   using reference_type = element_type&;
   using memory_tag     = TMemTag;

   static constexpr bool IS_CONST = std::is_const_v< element_type >;
   static constexpr size_t F_SIZE = fSize;

   using BufferViewType = BufferView< element_type, 4, memory_tag >;

   using NonConstFieldType = Field< value_type, F_SIZE, memory_tag >;
   using ConstFieldType    = Field< const value_type, F_SIZE, memory_tag >;

   WALBERLA_HOST FieldView(const NonConstFieldType& field, IBlock& block)
      : FieldView(field.bufferSystem().view(block), field.numGhostLayers(), field.layout())
   {}

   template< typename = void >
      requires(IS_CONST)
   WALBERLA_HOST FieldView(const NonConstFieldType& field, const IBlock& block)
      : FieldView(field.bufferSystem().view(block), field.numGhostLayers(), field.layout())
   {}

   template< typename = void >
      requires(IS_CONST)
   WALBERLA_HOST FieldView(const ConstFieldType& field, const IBlock& block)
      : FieldView(field.bufferSystem().view(block), field.numGhostLayers(), field.layout())
   {}

   /**
    * Convert non-const to const field view
    */
   template< typename = void >
      requires(IS_CONST)
   WALBERLA_HOST_DEVICE FieldView(const FieldView< value_type, F_SIZE, memory_tag >& other)
      : bufferView_{ other.bufferView() }, numGhostLayers_{ other.numGhostLayers() }, layout_{ other.layout() }
   {}

   WALBERLA_HOST_DEVICE size_t numGhostLayers() const { return numGhostLayers_; }

   WALBERLA_HOST_DEVICE field::Layout layout() const { return layout_; }

   WALBERLA_HOST_DEVICE const BufferViewType& bufferView() const { return bufferView_; }

   WALBERLA_HOST_DEVICE stdlib::span< const size_t, BufferViewType::RANK > shape() const { return bufferView_.shape(); }

   WALBERLA_HOST_DEVICE stdlib::span< const size_t, BufferViewType::RANK > allocShape() const
   {
      return bufferView_.allocShape();
   }

   WALBERLA_HOST_DEVICE stdlib::span< const size_t, BufferViewType::RANK > allocOffsets() const
   {
      return bufferView_.allocOffsets();
   }

   WALBERLA_HOST_DEVICE stdlib::span< const size_t, BufferViewType::RANK > strides() const
   {
      return bufferView_.strides();
   }

   WALBERLA_HOST FieldSlices slices() const
   {
      const auto& indexing = bufferView_.indexing();
      std::array< size_t, 3 > fieldShape{ indexing.shape()[0], indexing.shape()[1], indexing.shape()[2] };
      return { fieldShape, numGhostLayers_ };
   }

   WALBERLA_HOST_DEVICE reference_type operator()(Cell c, cell_idx_t f) const
   {
      return bufferView_(c.x(), c.y(), c.z(), f);
   }

   WALBERLA_HOST_DEVICE reference_type operator()(cell_idx_t x, cell_idx_t y, cell_idx_t z, cell_idx_t f) const
   {
      return bufferView_(x, y, z, f);
   }

   template< typename = void >
      requires(F_SIZE == 1)
   WALBERLA_HOST_DEVICE reference_type operator()(Cell c) const
   {
      return bufferView_(c.x(), c.y(), c.z(), 0);
   }

   template< typename = void >
      requires(F_SIZE == 1)
   WALBERLA_HOST_DEVICE reference_type operator()(cell_idx_t x, cell_idx_t y, cell_idx_t z) const
   {
      return bufferView_(x, y, z, 0);
   }

 private:
   FieldView(BufferViewType bufferView, size_t numGhostLayers, field::Layout layout)
      : bufferView_{ std::move(bufferView) }, numGhostLayers_{ numGhostLayers }, layout_{ layout }
   {}

   BufferViewType bufferView_;
   size_t numGhostLayers_;
   field::Layout layout_;
};

template< typename TView >
concept IFieldView =
   std::same_as< TView, FieldView< typename TView::element_type, TView::F_SIZE, typename TView::memory_tag > >;

template< std::semiregular TValue, size_t fSize, MemTag TMemTag >
using ConstFieldView = FieldView< const TValue, fSize, TMemTag >;

// View type aliases

template< IField TField >
using FieldViewType = FieldView< typename TField::element_type, TField::F_SIZE, typename TField::memory_tag >;

template< IField TField >
using ConstFieldViewType = FieldView< const typename TField::value_type, TField::F_SIZE, typename TField::memory_tag >;

// Deduction Guides for FieldView

template< typename FieldType >
FieldView(const FieldType&, IBlock&)
   -> FieldView< typename FieldType::element_type, FieldType::F_SIZE, typename FieldType::memory_tag >;

template< typename FieldType >
FieldView(const FieldType&, const IBlock&)
   -> FieldView< const typename FieldType::value_type, FieldType::F_SIZE, typename FieldType::memory_tag >;

} // namespace memory

using memory::ConstFieldView;
using memory::Field;
using memory::FieldView;
using memory::IField;
using memory::IFieldView;

} // namespace walberla::v8