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
//! \file FieldVtkExport.hpp
//! \author Frederik Hennig <frederik.hennig@fau.de>
//
//======================================================================================================================

#pragma once

#include "vtk/BlockCellDataWriter.h"

#include <optional>

#include "./Field.hpp"

namespace walberla::v8::memory
{

template< IField TField, typename TOutput >
   requires(TField::memory_tag::isHostAccessible)
class FieldVtkWriter : public vtk::BlockCellDataWriter< TOutput, TField::F_SIZE >
{
 public:
   using Base       = vtk::BlockCellDataWriter< TOutput, TField::F_SIZE >;
   using FieldType  = TField;
   using ViewType   = memory::ConstFieldViewType< FieldType >;
   using OutputType = TOutput;

   FieldVtkWriter(const FieldType& field, const std::string& outputId) : Base(outputId), field_{ field } {}

 protected:
   void configure() override { view_ = ViewType(field_, *this->block_); }

   using Base::evaluate;

   OutputType evaluate(cell_idx_t x, cell_idx_t y, cell_idx_t z, cell_idx_t f) override
   {
      WALBERLA_ASSERT(view_.has_value());
      return numeric_cast< OutputType >((*view_)(x, y, z, f)); // NOLINT(bugprone-unchecked-optional-access)
   }

 private:
   FieldType field_;
   std::optional< ViewType > view_;
};

} // namespace walberla::v8::memory