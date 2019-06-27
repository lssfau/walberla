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
//! \file WriteOutput.h
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#pragma once

#include <mesa_pd/data/DataTypes.h>
#include <mesa_pd/vtk/TensorGlyph.h>

#include <core/debug/CheckFunctions.h>
#include <vtk/Base64Writer.h>
#include <vtk/UtilityFunctions.h>

#include <ostream>

namespace walberla {
namespace mesa_pd {
namespace vtk {

template <typename T>
inline
void writeOutput(std::ostream& os, const T& data, const uint_t component)
{
   static_assert (std::is_arithmetic<T>::value, "this function only supports arithmetic data types" );
   WALBERLA_ASSERT_EQUAL(component, 0);
   WALBERLA_UNUSED(component);
   walberla::vtk::toStream(os, data);
}

template <>
inline
void writeOutput(std::ostream& os, const TensorGlyph& data, const uint_t component)
{
   WALBERLA_ASSERT_LESS(component, 6);
   walberla::vtk::toStream(os, data[component]);
}

template <>
inline
void writeOutput(std::ostream& os, const Vec3& data, const uint_t component)
{
   WALBERLA_ASSERT_LESS(component, 3);
   walberla::vtk::toStream(os, data[component]);
}

template <typename T>
inline
void writeOutput(walberla::vtk::Base64Writer& b64, const T& data, const uint_t component)
{
   static_assert (std::is_arithmetic<T>::value, "this function only supports arithmetic data types" );
   WALBERLA_ASSERT_EQUAL(component, 0);
   WALBERLA_UNUSED(component);
   b64 << data;
}

template <>
inline
void writeOutput(walberla::vtk::Base64Writer& b64, const TensorGlyph& data, const uint_t component)
{
   WALBERLA_ASSERT_LESS(component, 6);
   b64 << data[component];
}

template <>
inline
void writeOutput(walberla::vtk::Base64Writer& b64, const Vec3& data, const uint_t component)
{
   WALBERLA_ASSERT_LESS(component, 3);
   b64 << data[component];
}

} // namespace vtk
} // namespace pe
} // namespace walberla
