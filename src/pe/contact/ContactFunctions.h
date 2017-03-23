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
//! \file ContactFunctions.h
//! \author Klaus Iglberger
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#pragma once

#include "pe/contact/Contact.h"
#include "pe/Materials.h"

#include "core/math/AABB.h"

namespace walberla {
namespace pe {

inline bool shouldContactBeTreated( ContactID c, const math::AABB& blkAABB );
inline real_t         getRestitution( ConstContactID c );
inline real_t         getStiffness(ConstContactID c);
inline real_t         getDampingN(ConstContactID c);
inline real_t         getDampingT(ConstContactID c);
//   inline real_t         getCorParam(ConstContactID c);
inline real_t         getFriction(ConstContactID c);

} // namespace pe
} // namespace walberla

#include "ContactFunctions.impl.h"
