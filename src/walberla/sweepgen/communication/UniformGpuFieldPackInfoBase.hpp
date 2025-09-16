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
//! \file UniformGpuFieldPackInfoBase.hpp
//! \author Frederik Hennig <frederik.hennig@fau.de>
//
//======================================================================================================================

#pragma once

#if defined(WALBERLA_BUILD_WITH_CUDA) || defined(WALBERLA_BUILD_WITH_HIP)

#include "core/all.h"

#include "gpu/communication/GeneratedGPUPackInfo.h"

#include <concepts>
#include <span>

namespace walberla::sweepgen
{

namespace detail
{
template< typename T >
concept UniformGpuFieldPackInfoImpl = requires(T impl,                                              //
                                               typename T::Field_T& field,                          //
                                               std::span< typename T::Field_T::value_type > buffer, //
                                               stencil::Direction dir,                              //
                                               CellInterval& ci,                                    //
                                               gpuStream_t stream                                   //
) {
   typename T::Field_T;

   { impl.doPack(field, buffer, dir, ci, stream) } -> std::same_as< void >;

   { impl.doUnpack(field, buffer, dir, ci, stream) } -> std::same_as< void >;

   { impl.doLocalCopy(field, ci, field, ci, dir, stream) } -> std::same_as< void >;

   { impl.elementsPerCell(dir) } -> std::same_as< uint_t >;
};
} // namespace detail

template< typename Impl >
class UniformGpuFieldPackInfoBase : public gpu::GeneratedGPUPackInfo
{
 public:
   // static_assert( detail::UniformGpuFieldPackInfoImpl< Impl >, "Impl does not satisfy contraints.");

   UniformGpuFieldPackInfoBase(BlockDataID fieldId, uint_t sliceWidth = 1)
      : fieldId_{ fieldId }, sliceWidth_{ cell_idx_c(sliceWidth) }
   {}

   void pack(stencil::Direction dir, unsigned char* rawBuffer, IBlock* block, gpuStream_t stream) override
   {
      using Field_T    = typename Impl::Field_T;
      using value_type = typename Field_T::value_type;
      Field_T * field   = block->getData< Field_T >(fieldId_);
      CellInterval ci;
      field->getSliceBeforeGhostLayer(dir, ci, sliceWidth_, false);
      std::span< value_type > buffer{ ( value_type* ) rawBuffer, this->size(dir, block) };
      impl().doPack(field, buffer, dir, ci, stream);
   }

   void unpack(stencil::Direction dir, unsigned char* rawBuffer, IBlock* block, gpuStream_t stream) override
   {
      using Field_T    = typename Impl::Field_T;
      using value_type = typename Field_T::value_type;
      Field_T * field   = block->getData< Field_T >(fieldId_);
      CellInterval ci;
      field->getGhostRegion(dir, ci, sliceWidth_, false);
      std::span< value_type > buffer{ (value_type*) rawBuffer, this->size(dir, block) };
      stencil::Direction commDir{ stencil::inverseDir[ dir ] };
      impl().doUnpack(field, buffer, commDir, ci, stream);
   }

   void communicateLocal  ( stencil::Direction dir, const IBlock *sender, IBlock *receiver, gpuStream_t stream ) override {
      using Field_T    = typename Impl::Field_T;

      Field_T * srcField = const_cast< IBlock * >(sender)->getData< Field_T >(fieldId_);
      Field_T * dstField = receiver->getData< Field_T >(fieldId_);

      CellInterval srcRegion;
      CellInterval dstRegion;
      srcField->getSliceBeforeGhostLayer(dir, srcRegion, sliceWidth_, false);
      dstField->getGhostRegion(stencil::inverseDir[dir], dstRegion, sliceWidth_, false);

      impl().doLocalCopy(srcField, srcRegion, dstField, dstRegion, dir, stream);
   }

   uint_t size(stencil::Direction dir, IBlock* block) override
   {
      using Field_T = typename Impl::Field_T;
      using value_type = typename Field_T::value_type;

      const Field_T * field = block->getData< Field_T >(fieldId_);
      CellInterval ci;
      field->getGhostRegion(dir, ci, sliceWidth_, false);

      uint_t elementsPerCell{ impl().elementsPerCell(dir) };
      return elementsPerCell * ci.numCells() * sizeof( value_type );
   }

 protected:
   BlockDataID fieldId_;
   cell_idx_t sliceWidth_;

 private:
   Impl& impl() { return static_cast< Impl& >(*this); }
};

} // namespace walberla::experimental::communication

#endif
