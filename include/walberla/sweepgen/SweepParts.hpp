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
//! \file SweepParts.hpp
//! \author Frederik Hennig <frederik.hennig@fau.de>
//
//======================================================================================================================

#pragma once

#include "domain_decomposition/IBlock.h"

#include <memory>
#include <stdexcept>

#include "walberla/v8/Memory.hpp"

namespace walberla::sweepgen::sweep_parts
{

/**
 * @brief Manages a single shadow buffer for a field to facilitate data swapping
 *
 * Internal component of SweepGen that realizes field-swapping logic
 * by managing a shadow buffer of a field.
 */
template< v8::memory::IField TField >
class ShadowBufferCache
{
 public:
   using FieldType      = TField;
   using BufferType     = FieldType::BufferSystemType::BufferType;
   using BufferViewType = FieldType::BufferSystemType::ViewType;

   ShadowBufferCache() = default;

   BufferViewType view(const FieldType& field)
   {
      if (shadowBuffer_ == nullptr)
      {
         shadowBuffer_ = field.bufferSystem().createEmptyBuffer();
      }

      return { shadowBuffer_->data(), shadowBuffer_->allocData(), field.bufferSystem().indexing() };
   }

   void swapBuffers(const FieldType& field, IBlock& block)
   {
      if (shadowBuffer_ == nullptr)
      {
         throw std::logic_error{ "`swapBuffers` called before shadow buffer was initialized" };
      }

      auto& blockBuffer = field.bufferSystem().buffer(block);
      shadowBuffer_->swap(blockBuffer);
   }

 private:
   std::unique_ptr< BufferType > shadowBuffer_;
};

} // namespace walberla::sweepgen::sweep_parts