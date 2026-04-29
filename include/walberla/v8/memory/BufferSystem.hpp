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
//! \file BufferSystem.hpp
//! \author Frederik Hennig <frederik.hennig@fau.de>
//
//======================================================================================================================

#pragma once

#include "blockforest/Block.h"
#include "blockforest/BlockDataHandling.h"
#include "blockforest/BlockID.h"
#include "blockforest/StructuredBlockForest.h"

#include "core/Stdlib.hpp"
#include WALBERLA_STDLIB(array)
#include WALBERLA_STDLIB(span)

#include "domain_decomposition/BlockDataID.h"

#include <algorithm>
#include <concepts>
#include <cstdint>
#include <ranges>
#include <span>
#include <utility>

#include "./MemoryTags.hpp"
#include "./MemoryTraits.hpp"

namespace walberla::v8::memory
{

/**
 * @defgroup v8core-buffersystem Buffer Systems
 * @ingroup v8core-memory
 * @brief Low-level memory management; buffer systems, allocators and initializers.
 */

namespace buffer_system
{

template< std::semiregular TValue, typename TAllocator >
class LinearBuffer
{
 public:
   using Self = LinearBuffer< TValue, TAllocator >;

   using value_type     = TValue;
   using allocator_type = TAllocator;

   LinearBuffer(size_t linearAllocSize, size_t linearAllocOffset, const allocator_type& allocator)
      : alloc_{ allocator }, linearAllocSize_{ linearAllocSize }, allocData_{ alloc_.allocate(linearAllocSize) },
        data_{ allocData_ + linearAllocOffset }
   {}

   ~LinearBuffer() { alloc_.deallocate(allocData_, linearAllocSize_); }

   LinearBuffer(const LinearBuffer&)            = delete;
   LinearBuffer& operator=(const LinearBuffer&) = delete;

   Self cloneUninitialized() const
   {
      const size_t allocOffset{ data_ - allocData_ };
      return { linearAllocSize_, allocOffset, alloc_ };
   }

   [[nodiscard]] size_t linearAllocSize() const { return linearAllocSize_; }

   [[nodiscard]] value_type* data() { return data_; }
   [[nodiscard]] const value_type* data() const { return data_; }

   [[nodiscard]] value_type* allocData() { return allocData_; }
   [[nodiscard]] const value_type* allocData() const { return allocData_; }

   [[nodiscard]] stdlib::span< value_type > linearView() { return { data_, allocData_ + linearAllocSize_ }; }
   [[nodiscard]] stdlib::span< const value_type > linearView() const { return { data_, allocData_ + linearAllocSize_ }; }

   [[nodiscard]] stdlib::span< value_type > linearAllocView() { return { allocData_, linearAllocSize_ }; }
   [[nodiscard]] stdlib::span< const value_type > linearAllocView() const { return { allocData_, linearAllocSize_ }; }

   // To make waLBerla data manager happy
   bool operator==(const Self& other) const { return this == &other; }

   void swap(Self& other)
   {
      std::swap(alloc_, other.alloc_);
      std::swap(linearAllocSize_, other.linearAllocSize_);
      std::swap(allocData_, other.allocData_);
      std::swap(data_, other.data_);
   }

 private:
   allocator_type alloc_;
   size_t linearAllocSize_;

   value_type* allocData_;
   value_type* data_;
};

template< std::semiregular TValue, typename TAllocator >
class LinearBufferDataHandling
   : public blockforest::AlwaysInitializeBlockDataHandling< LinearBuffer< TValue, TAllocator > >
{
 public:
   using buffer_value_type = TValue;
   using allocator_type    = TAllocator;
   using buffer_type       = LinearBuffer< buffer_value_type, allocator_type >;

   LinearBufferDataHandling(size_t linearAllocSize, size_t linearAllocOffset, const allocator_type& allocator)
      : linearAllocSize_{ linearAllocSize }, linearAllocOffset_{ linearAllocOffset }, allocator_{ allocator }
   {}

   buffer_type* initialize(IBlock* /* block */) override
   {
      return new buffer_type(linearAllocSize_, linearAllocOffset_, allocator_);
   }

 private:
   size_t linearAllocSize_;
   size_t linearAllocOffset_;

   allocator_type allocator_;
};

template< size_t TRank >
class BufferIndexing
{
 public:
   using index_type             = std::ptrdiff_t;
   static constexpr size_t RANK = TRank;

   WALBERLA_HOST_DEVICE BufferIndexing(
      stdlib::array< size_t, RANK > shape, 
      stdlib::array< size_t, RANK > allocShape,
      stdlib::array< size_t, RANK > allocOffsets, 
      stdlib::array< size_t, RANK > strides)
      : shape_{ shape }, allocShape_{ allocShape }, allocOffsets_{ allocOffsets }, strides_{ strides }
   {}

   WALBERLA_HOST_DEVICE stdlib::span< const size_t, RANK > shape() const { return shape_; }

   WALBERLA_HOST_DEVICE stdlib::span< const size_t, RANK > allocShape() const { return allocShape_; }

   WALBERLA_HOST_DEVICE stdlib::span< const size_t, RANK > allocOffsets() const { return allocOffsets_; }

   WALBERLA_HOST_DEVICE stdlib::span< const size_t, RANK > strides() const { return strides_; }

   WALBERLA_HOST_DEVICE bool operator==(const BufferIndexing& other) const
   {
      return shape_        == other.shape_ &&
             allocShape_   == other.allocShape_ &&
             allocOffsets_ == other.allocOffsets_ &&
             strides_      == other.strides_;
   }
   
   WALBERLA_HOST_DEVICE bool operator!=(const BufferIndexing& other) const { return !(*this == other); }

   WALBERLA_HOST_DEVICE size_t linearBufferSize() const
   {
      size_t sz{ 1 };
      for (size_t i = 0; i < RANK; ++i)
      {
         sz *= shape_[i];
      }
      return sz;
   }

   WALBERLA_HOST_DEVICE size_t linearBufferAllocSize() const
   {
      size_t sz{ 1 };
      for (size_t i = 0; i < RANK; ++i)
      {
         sz *= allocShape_[i];
      }
      return sz;
   }

   WALBERLA_HOST_DEVICE size_t linearBufferOffset() const
   {
      size_t offset{ 0 };
      for (size_t i = 0; i < RANK; ++i)
      {
         offset += allocOffsets_[i] * strides_[i];
      }
      return offset;
   }

   WALBERLA_HOST_DEVICE index_type linearIndex(std::convertible_to< index_type > auto... indices) const
   {
      return linearIndexImpl(std::make_index_sequence< RANK >{}, indices...);
   }

 private:
   template< size_t... Coords >
   WALBERLA_HOST_DEVICE index_type linearIndexImpl(std::index_sequence< Coords... >, std::convertible_to< index_type > auto... indices) const
   {
      return (index_type(0) + ... + (static_cast< index_type >(indices) * static_cast< index_type >(strides_[Coords])));
   }

   stdlib::array< size_t, RANK > shape_;
   stdlib::array< size_t, RANK > allocShape_;
   stdlib::array< size_t, RANK > allocOffsets_;
   stdlib::array< size_t, RANK > strides_;
};

template< typename TElement, size_t TRank, MemTag TMemTag >
   requires(std::semiregular< std::remove_cv_t< TElement > >)
class BufferView
{
 public:
   using element_type   = TElement;
   using value_type     = std::remove_cv_t< element_type >;
   using pointer_type   = element_type*;
   using reference_type = element_type&;

   static constexpr bool IS_CONST = std::is_const_v< element_type >;

   static constexpr size_t RANK = TRank;

   using memory_tag = TMemTag;

   using IndexingType   = BufferIndexing< TRank >;
   using index_type     = IndexingType::index_type;

   WALBERLA_HOST_DEVICE BufferView(element_type* data, element_type* allocData, const IndexingType& indexing)
      : data_{ data }, allocData_{ allocData }, indexing_{ indexing }
   {}

   /**
    * Convert non-const to const buffer view
    */
   template< typename = void >
      requires(IS_CONST)
   WALBERLA_HOST_DEVICE BufferView(const BufferView< value_type, RANK, memory_tag >& other)
      : data_{ other.data() }, allocData_{ other.allocData() }, indexing_{ other.indexing() }
   {}

   [[nodiscard]] WALBERLA_HOST_DEVICE reference_type operator()(std::convertible_to< index_type > auto... indices) const
   {
      return data_[indexing_.linearIndex(indices...)];
   }

   WALBERLA_HOST_DEVICE pointer_type data() const { return data_; }

   WALBERLA_HOST_DEVICE pointer_type allocData() const { return allocData_; }

   WALBERLA_HOST_DEVICE pointer_type dataAt(std::convertible_to< index_type > auto... indices) const
   {
      return &((*this)(indices...));
   }

   WALBERLA_HOST_DEVICE const IndexingType& indexing() const { return indexing_; }

   WALBERLA_HOST_DEVICE stdlib::span< const size_t, RANK > shape() const { return indexing_.shape(); }

   WALBERLA_HOST_DEVICE stdlib::span< const size_t, RANK > allocShape() const { return indexing_.allocShape(); }

   WALBERLA_HOST_DEVICE stdlib::span< const size_t, RANK > allocOffsets() const { return indexing_.allocOffsets(); }

   WALBERLA_HOST_DEVICE stdlib::span< const size_t, RANK > strides() const { return indexing_.strides(); }

   [[nodiscard]] WALBERLA_HOST_DEVICE stdlib::span< element_type > linearView() const
   {
      return { data_, allocData_ + indexing_.linearBufferAllocSize() };
   }

   [[nodiscard]] WALBERLA_HOST_DEVICE stdlib::span< element_type > linearAllocView() const
   {
      return { allocData_, indexing_.linearBufferAllocSize() };
   }

 private:
   pointer_type data_;
   pointer_type allocData_;
   IndexingType indexing_;
};

template< typename TView >
concept IBufferView =
   std::same_as< TView, BufferView< typename TView::element_type, TView::RANK, typename TView::memory_tag > >;

template< std::semiregular TValue, size_t TRank, MemTag TMemTag >
class BufferSystemBuilder;

template< typename TElement, size_t TRank, MemTag TMemTag >
   requires(std::semiregular< std::remove_cv_t< TElement > >)
class BufferSystem
{
 public:
   using Self = BufferSystem< TElement, TRank, TMemTag >;

   using element_type           = TElement;
   using value_type             = std::remove_cv_t< element_type >;
   static constexpr size_t RANK = TRank;

   using MT            = MemoryTraits< TMemTag, value_type >;
   using AllocatorType = MT::AllocatorType;

   static constexpr bool IS_CONST = std::is_const_v< element_type >;
   using BufferType               = std::conditional_t< IS_CONST, const LinearBuffer< value_type, AllocatorType >,
                                                        LinearBuffer< value_type, AllocatorType > >;
   using ConstBufferType          = const LinearBuffer< value_type, AllocatorType >;
   using IndexingType             = BufferIndexing< RANK >;

   using BlockDataIDType = std::conditional_t< IS_CONST, ConstBlockDataID, BlockDataID >;

   using memory_tag = TMemTag;

   using ViewType      = BufferView< element_type, RANK, memory_tag >;
   using ConstViewType = BufferView< const element_type, RANK, memory_tag >;

   static BufferSystemBuilder< value_type, RANK, memory_tag > create(StructuredBlockForest& blocks,
                                                                     const stdlib::array< size_t, RANK >& shape);

   template< typename = void >
      requires(!IS_CONST)
   BufferSystem(StructuredBlockForest& blocks, IndexingType indexing, AllocatorType allocator) : indexing_{ indexing }
   {
      using DataHandling_T = LinearBufferDataHandling< value_type, AllocatorType >;
      auto dh = std::make_shared< DataHandling_T >(indexing_.linearBufferAllocSize(), indexing_.linearBufferOffset(),
                                                   allocator);

      buffersId_ = blocks.addBlockData(dh);
   }

   // Convert non-const to const buffer system
   template< typename = void >
      requires(IS_CONST)
   BufferSystem(const BufferSystem< value_type, RANK, memory_tag >& other)
      : indexing_{ other.indexing() }, buffersId_{ other.buffersId() }
   {}

   const IndexingType& indexing() const { return indexing_; }

   BlockDataIDType buffersId() const { return buffersId_; }

   ViewType view(IBlock& block) const
   {
      auto& buf = buffer(block);
      return { buf.data(), buf.allocData(), indexing_ };
   }

   ConstViewType view(const IBlock& block) const
   {
      auto& buf = buffer(block);
      return { buf.data(), buf.allocData(), indexing_ };
   }

   BufferType& buffer(IBlock& block) const { return *block.uncheckedFastGetData< BufferType >(buffersId_); }

   ConstBufferType& buffer(const IBlock& block) const { return *block.uncheckedFastGetData< BufferType >(buffersId_); }

   /**
    * @brief Create a new empty buffer associated with this buffer system.
    *
    */
   std::unique_ptr< BufferType > createEmptyBuffer() const
   {
      AllocatorType allocator{ MT::getAllocator(sizeof(value_type) * indexing_.linearBufferOffset()) };
      return std::make_unique< BufferType >(indexing_.linearBufferAllocSize(), indexing_.linearBufferOffset(),
                                            allocator);
   }

 private:
   IndexingType indexing_;
   BlockDataIDType buffersId_;
};

template< std::semiregular TValue, size_t TRank, MemTag TMemTag >
class BufferSystemBuilder
{
 public:
   using value_type             = TValue;
   static constexpr size_t RANK = TRank;
   using memory_tag             = TMemTag;

   using Self             = BufferSystemBuilder< value_type, RANK, memory_tag >;
   using MT               = MemoryTraits< memory_tag, value_type >;
   using IndexingType     = BufferIndexing< RANK >;
   using AllocatorType    = MT::AllocatorType;
   using BufferSystemType = BufferSystem< value_type, RANK, memory_tag >;

   BufferSystemBuilder(StructuredBlockForest& blocks, const stdlib::array< size_t, RANK >& shape)
      : blocks_{ &blocks }, shape_{ shape }, allocShape_{ shape }
   {}

   Self& allocShape(const stdlib::array< size_t, RANK >& allocShape)
   {
      allocShape_ = allocShape;
      return *this;
   }

   Self& allocOffsets(const stdlib::array< size_t, RANK >& allocOffsets)
   {
      allocOffsets_ = allocOffsets;
      return *this;
   }

   Self& linearizationOrder(const stdlib::array< size_t, RANK >& linearizationOrder)
   {
      linearizationOrder_ = linearizationOrder;
      return *this;
   }

   Self& addLinePaddingForAlignment(bool addPadding = true)
   {
      withPaddingForAlignment_ = addPadding;
      return *this;
   }

   Self& alignment(size_t alignmentBytes)
   {
      alignmentBytes_ = alignmentBytes;
      return *this;
   }

   BufferSystemType build()
   {
      if (withPaddingForAlignment_) { addPaddingToAllocShape(); }

      auto indexing = createIndexing();
      AllocatorType allocator{ MT::getAllocator(sizeof(value_type) * indexing.linearBufferOffset()) };

      return { *blocks_, indexing, allocator };
   }

 private:
   void addPaddingToAllocShape()
   {
      const size_t alignmentElements{ alignmentBytes_ / sizeof(value_type) };
      const stdlib::array< size_t, RANK >& minimumAllocShape{ allocShape_ };
      const size_t fastestCoord{ linearizationOrder_[0] };
      const size_t fastestDimensionMinimum{ minimumAllocShape[fastestCoord] };
      const size_t requiredPadding{ (alignmentElements - (fastestDimensionMinimum % alignmentElements)) %
                                    alignmentElements };

      allocShape_[fastestCoord] += requiredPadding;
   }

   IndexingType createIndexing()
   {
      stdlib::array< size_t, RANK > strides;
      size_t nextStride{ 1 };
      for (auto coord : linearizationOrder_)
      {
         strides[coord] = nextStride;
         nextStride *= allocShape_[coord];
      }

      return { shape_, allocShape_, allocOffsets_, strides };
   }

   static constexpr stdlib::array< size_t, RANK > zeroOffsets()
   {
      stdlib::array< size_t, RANK > zeros;
      for (size_t i = 0; i < RANK; ++i)
      {
         zeros[i] = 0;
      }
      return zeros;
   }

   static constexpr stdlib::array< size_t, RANK > rightmostFirst()
   {
      stdlib::array< size_t, RANK > linearizationOrder;
      for (size_t i = 0; i < RANK; ++i)
      {
         linearizationOrder[i] = RANK - i - 1;
      }
      return linearizationOrder;
   }

   StructuredBlockForest* blocks_;

   stdlib::array< size_t, RANK > shape_;
   stdlib::array< size_t, RANK > allocShape_;
   stdlib::array< size_t, RANK > allocOffsets_{ zeroOffsets() };
   stdlib::array< size_t, RANK > linearizationOrder_{ rightmostFirst() };

   bool withPaddingForAlignment_{ false };
   size_t alignmentBytes_{ MT::alignment() };
};

template< typename TElement, size_t TRank, MemTag TMemTag >
   requires(std::semiregular< std::remove_cv_t< TElement > >)
inline auto BufferSystem< TElement, TRank, TMemTag >::create(StructuredBlockForest& blocks,
                                                             const stdlib::array< size_t, RANK >& shape)
   -> BufferSystemBuilder< value_type, RANK, memory_tag >
{
   return { blocks, shape };
}

} // namespace buffer_system

using buffer_system::BufferSystem;
using buffer_system::BufferView;
using buffer_system::IBufferView;

} // namespace walberla::v8::memory