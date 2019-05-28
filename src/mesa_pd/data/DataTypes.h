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
//! \file DataTypes.h
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#pragma once

#include <core/DataTypes.h>
#include <core/math/Matrix3.h>
#include <core/math/Quaternion.h>
#include <core/math/Rot3.h>
#include <core/math/Vector3.h>

#include <mesa_pd/data/Flags.h>

namespace walberla {
namespace mesa_pd {

using Mat3 = math::Matrix3<real_t>;
using Rot3 = math::Rot3<real_t>;
using Quat = math::Quaternion<real_t>;
using Vec3 = math::Vector3<real_t>;

}
}
