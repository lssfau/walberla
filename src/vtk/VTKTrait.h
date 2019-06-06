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
//! \file VTKTrait.h
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#pragma once

#include <core/DataTypes.h>
#include <core/math/Vector3.h>

namespace walberla {
namespace vtk {

template <typename T>
struct VTKTrait {};

template <>
struct VTKTrait<int8_t>
{
   using type = int8_t;
   constexpr static char const * const type_string = "Int8";
   constexpr static const uint_t components = 1;
};

template <>
struct VTKTrait<int16_t>
{
   using type = int16_t;
   constexpr static char const * const type_string = "Int16";
   constexpr static const uint_t components = 1;
};

template <>
struct VTKTrait<int32_t>
{
   using type = int32_t;
   constexpr static char const * const type_string = "Int32";
   constexpr static const uint_t components = 1;
};

template <>
struct VTKTrait<int64_t>
{
   using type = int64_t;
   constexpr static char const * const type_string = "Int64";
   constexpr static const uint_t components = 1;
};

template <>
struct VTKTrait<uint8_t>
{
   using type = uint8_t;
   constexpr static char const * const type_string = "UInt8";
   constexpr static const uint_t components = 1;
};

template <>
struct VTKTrait<uint16_t>
{
   using type = uint16_t;
   constexpr static char const * const type_string = "UInt16";
   constexpr static const uint_t components = 1;
};

template <>
struct VTKTrait<uint32_t>
{
   using type = uint32_t;
   constexpr static char const * const type_string = "UInt32";
   constexpr static const uint_t components = 1;
};

template <>
struct VTKTrait<uint64_t>
{
   using type = uint64_t;
   constexpr static char const * const type_string = "UInt64";
   constexpr static const uint_t components = 1;
};

template <>
struct VTKTrait<float>
{
   using type = float;
   constexpr static char const * const type_string = "Float32";
   constexpr static const uint_t components = 1;
};

template <>
struct VTKTrait<double>
{
   using type = double;
   constexpr static char const * const type_string = "Float64";
   constexpr static const uint_t components = 1;
};

template <typename T>
struct VTKTrait<math::Vector3<T>>
{
   using type = typename VTKTrait<T>::type;
   constexpr static char const * const type_string = VTKTrait<T>::type_string;
   constexpr static const uint_t components = 3;
};

} // namespace vtk
} // namespace walberla
