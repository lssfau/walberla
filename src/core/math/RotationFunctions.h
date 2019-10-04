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
//! \file RotationFunctions.h
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#pragma once

#include <core/math/Rot3.h>
#include <core/math/Vector3.h>

namespace walberla {
namespace math {

template< typename Type >
inline void rotateAroundPoint( Vector3<Type>& pos,
                               Rot3<Type>& rot,
                               const Vector3<Type>& rotateAround,
                               const Rot3<Type>& dRot )
{
   const auto dp = pos - rotateAround;
   pos = rotateAround + dRot.getMatrix() * dp;
   rot.rotate(dRot);
}

} // mpi
} // walberla
