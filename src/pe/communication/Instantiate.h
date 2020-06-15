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
//! \file Instantiate.h
//! \author Tobias Preclik
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//! \brief Marshalling of objects for data transmission or storage.
//
//======================================================================================================================

#pragma once


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include "core/Abort.h"
#include "core/debug/demangle.h"
#include "core/math/AABB.h"
#include "core/math/Vector3.h"

#include "pe/Types.h"

namespace walberla {
namespace pe {
namespace communication {

inline
void correctBodyPosition(const math::AABB& domain, const Vec3& center, Vec3& pos)
{
   Vec3 dis = pos - center;

   if (-domain.xSize() * 0.5 > dis[0]) pos[0] += domain.xSize();
   if (+domain.xSize() * 0.5 < dis[0]) pos[0] -= domain.xSize();

   if (-domain.ySize() * 0.5 > dis[1]) pos[1] += domain.ySize();
   if (+domain.ySize() * 0.5 < dis[1]) pos[1] -= domain.ySize();

   if (-domain.zSize() * 0.5 > dis[2]) pos[2] += domain.zSize();
   if (+domain.zSize() * 0.5 < dis[2]) pos[2] -= domain.zSize();

   WALBERLA_ASSERT(dis.sqrLength() >= (pos - center).sqrLength());
}

template < class BodyT >
std::unique_ptr<BodyT> instantiate( mpi::RecvBuffer& /*buffer*/, const math::AABB& /*domain*/, const math::AABB& /*block*/, BodyT*& /*newBody*/ )
{
   WALBERLA_ABORT( "Body instantiation not implemented! (" << debug::demangle(typeid(BodyT).name()) << ")" );
}

}  // namespace communication
}  // namespace pe
}  // namespace walberla
