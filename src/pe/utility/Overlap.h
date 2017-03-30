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
//! \file Overlap.h
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#pragma once

#include <core/DataTypes.h>

namespace walberla {
namespace pe {

/*!\brief Calculates overlap between two spheres.
 * \param d distance between two spheres
 * \param r1 radius of sphere 1
 * \param r2 radius of sphere 2
 * \return overlap volume
 */
real_t getSphereSphereOverlap(const real_t d, const real_t r1, const real_t r2);

}  // namespace pe
}  // namespace walberla
