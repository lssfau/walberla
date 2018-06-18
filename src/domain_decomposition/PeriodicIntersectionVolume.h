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
//! \file PeriodicIntersectionVolume.h
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/DataTypes.h"
#include "core/math/AABB.h"
#include "core/math/Vector3.h"

#include <array>

namespace walberla {
namespace domain_decomposition {

real_t periodicIntersectionVolume( const std::array< bool, 3 > & periodic, const math::AABB & domain, const math::AABB & box1, const math::AABB & box2 );
real_t periodicIntersectionVolume( const std::array< bool, 3 > & periodic, const math::AABB & domain, const math::AABB& box1, const math::AABB & box2, const real_t dx );

} // namespace domain_decomposition
} // namespace walberla
