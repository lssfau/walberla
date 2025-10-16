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
//! \file SparseIndexList.hpp
//! \author Frederik Hennig <frederik.hennig@fau.de>
//
//======================================================================================================================

#pragma once

#include "blockforest/StructuredBlockForest.h"

#include "core/DataTypes.h"
#include "core/cell/Cell.h"

#include "domain_decomposition/BlockDataID.h"

#include <type_traits>
#include <vector>
#include <span>

#include "../memory/MemoryTags.hpp"

namespace walberla::experimental::sweep
{
namespace internal
{
/**
 * @brief Index vectors for sparse sweeps
 */
template< typename IndexStruct, typename Allocator_T >
class IndexListBuffer
{
 public:
   using Vector_T = std::vector< IndexStruct, Allocator_T >;

   IndexListBuffer() = default;
   bool operator==(const IndexListBuffer& other ) const { return &other.vec_ == &vec_; }

   Vector_T& vector() { return vec_; }

   IndexStruct* data() { return vec_.data(); }
   const IndexStruct* data() const { return vec_.data(); }
   size_t size() const { return vec_.size(); }

 private:
   Vector_T vec_;
};
} // namespace internal

struct CellIdx
{
   int64_t x;
   int64_t y;
   int64_t z;

   CellIdx(const Cell& c) : x{ int64_c(c.x()) }, y{ int64_c(c.y()) }, z{ int64_c(c.z()) } {}

   bool operator==(const CellIdx& other) const { return std::tie(x, y, z) == std::tie(other.x, other.y, other.z); }
};

template< typename IndexStruct = CellIdx, memory::MemTag memtag_t = memory::memtag::host >
class SparseIndexList
{
 public:
   using Buffer_T = internal::IndexListBuffer< IndexStruct, memory::StandardAllocator_T< memtag_t, IndexStruct > >;

   explicit SparseIndexList(StructuredBlockForest& sbfs) : bufferId_{ createBuffers(sbfs) } {}

   Buffer_T::Vector_T& getVector(IBlock& block)
   {
      Buffer_T* buf = block.template getData< Buffer_T >(bufferId_);
      return buf->vector();
   }

   std::span< IndexStruct > view(IBlock& block){
    Buffer_T& buf = *block.template getData< Buffer_T >(bufferId_);
    return { buf.data(), buf.size() };
   }

   std::span< const IndexStruct > view(IBlock& block) const {
    const Buffer_T& buf = *block.template getData< Buffer_T >(bufferId_);
    return { buf.data(), buf.size() };
   }

   const BlockDataID& bufferId() const { return bufferId_; }

 private:
   static BlockDataID createBuffers(StructuredBlockForest& sbfs)
   {
      auto createIdxVector = [](IBlock* const, StructuredBlockStorage* const) { return new Buffer_T(); };
      return sbfs.addStructuredBlockData< Buffer_T >(createIdxVector);
   }

   BlockDataID bufferId_;
};
} // namespace walberla::experimental::sweep
