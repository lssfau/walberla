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
//! \file StencilRanges.hpp
//! \author Frederik Hennig <frederik.hennig@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/Stdlib.hpp"

#include "core/Macros.h"

#include "core/math/Vector3.h"

#include "stencil/Iterator.h"
#include "stencil/Directions.h"

#include WALBERLA_STDLIB(array)
#include WALBERLA_STDLIB(span)

#include <concepts>
#include <span>

namespace walberla::v8::stencil_ranges
{

using stencil::Direction;

class Velocities {
public:
   constexpr Velocities() : cs_( stencil::c ) {}

   constexpr Vector3< cell_idx_t > operator[] (Direction dir) const {
      return Vector3< cell_idx_t >{
         cs_[0][dir],
         cs_[1][dir],
         cs_[2][dir],
      };
   }

   constexpr cell_idx_t x(Direction dir) const { return cs_[0][dir]; }
   constexpr cell_idx_t y(Direction dir) const { return cs_[1][dir]; }
   constexpr cell_idx_t z(Direction dir) const { return cs_[2][dir]; }

private:
   std::array< std::array< int, stencil::NR_OF_DIRECTIONS >, 3 > cs_;
};

template< typename TStencil >
class StencilRange
{
 public:
   using Stencil = TStencil;

   constexpr StencilRange() : dirs_(TStencil::dir) {}
   constexpr StencilRange(size_t start, size_t end) : dirs_(TStencil::dir), start_{ start }, end_{ end } {}

   constexpr const Direction* begin() const { return &dirs_[start_]; }
   constexpr const Direction* end() const { return dirs_.data() + end_; }

 private:
   std::array< Direction, Stencil::Q > dirs_;
   size_t start_{ 0 };
   size_t end_{ Stencil::Q };
};

template< typename TStencil >
constexpr auto all()
{
   return StencilRange< TStencil >();
}

template< typename TStencil >
constexpr auto noCenter()
{
   return StencilRange< TStencil >(TStencil::noCenterFirstIdx, TStencil::Q);
}

template< typename TStencil >
class Subdirections
{
 public:
   static constexpr size_t Q = TStencil::Q;
   static constexpr size_t N = Q / 2;

   constexpr WALBERLA_HOST_DEVICE Subdirections()
      : grid_(TStencil::d_per_d), sizes_(TStencil::d_per_d_length) {}

   constexpr WALBERLA_HOST_DEVICE
   stdlib::span< const Direction > of (Direction primary) const {
      return { grid_[primary].data(), sizes_[primary] };
   }

   constexpr WALBERLA_HOST_DEVICE
   size_t size(Direction dir) const {
      return sizes_[dir];
   }

 private:
   std::array< std::array< Direction, N >, stencil::NR_OF_DIRECTIONS > grid_;
   std::array< size_t, Q > sizes_;
};

template< typename TStencil >
class SubdirectionsMatrix {
public:
   constexpr SubdirectionsMatrix()
      : subdirs_( TStencil::subdirs ), subdirsPtr_( TStencil::subdirsPtr ) {}

   WALBERLA_HOST_DEVICE constexpr
   std::span< const Direction > operator() (Direction dir) const {
      const size_t dirIdx(dir);
      return {
         subdirs_.begin() + subdirsPtr_[dirIdx],
         subdirs_.begin() + subdirsPtr_[dirIdx + 1]
      };
   }
private:
   std::array< Direction, TStencil::subdirsLength > subdirs_;
   std::array< size_t, stencil::NR_OF_DIRECTIONS + 1 > subdirsPtr_;
};

template< typename TStencil >
class SubdirectionsRange {
public:
   constexpr SubdirectionsRange(Direction primary) : matrix_(), primary_{primary} {}

   constexpr auto begin() const {
      return matrix_(primary_).begin();
   }

   constexpr auto end() const {
      return matrix_(primary_).end();
   }

   constexpr size_t size() const {
      return matrix_(primary_).size();
   }

   constexpr Direction operator[] (size_t idx) const {
      return matrix_(primary_)[idx];
   }

private:
   SubdirectionsMatrix< TStencil > matrix_;
   Direction primary_;
};


template< typename TStencil >
WALBERLA_HOST_DEVICE
constexpr auto subdirections(Direction dir) {
   return SubdirectionsRange< TStencil >(dir);
}


} // namespace walberla::v8::stencil_ranges