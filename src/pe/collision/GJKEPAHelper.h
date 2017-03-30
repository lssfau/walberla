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
//! \file GJKHelper.h
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#pragma once

//*************************************************************************************************
// Includes
//*************************************************************************************************

#include "pe/contact/Contact.h"
#include "pe/Types.h"

#include "core/math/Vector3.h"

namespace walberla {
namespace pe {

static unsigned long maxGJKIterations = 100ul;
static real_t epaTolerance = real_t(0.000001);

inline
void setMaxGJKIterations(const unsigned long& iter)
{
   maxGJKIterations = iter;
}
inline
unsigned long getMaxGJKIterations()
{
   return maxGJKIterations;
}

inline
void setEPATolerance(const real_t& tol)
{
   epaTolerance = tol;
}
inline
real_t getEPATolerance()
{
   return epaTolerance;
}

bool collideGJK( ConstGeomID bd1,
                 ConstGeomID bd2,
                 Vec3& contactPoint,
                 Vec3& contactNormal,
                 real_t& penetrationDepth,
                 const unsigned long numIterations = maxGJKIterations,
                 const real_t epaTol = epaTolerance );

template <typename Container>
bool collideGJK( GeomID bd1,
                 GeomID bd2,
                 Container& container,
                 const unsigned long numIterations = maxGJKIterations,
                 const real_t epaTol = epaTolerance )
{
   real_t penetrationDepth;
   Vec3   contactPoint;
   Vec3   contactNormal;
   bool retVal = collideGJK(bd1, bd2, contactPoint, contactNormal, penetrationDepth, numIterations, epaTol);
   if (retVal)
   {
      container.push_back( Contact(bd1, bd2, contactPoint, contactNormal, penetrationDepth) );
   }
   return retVal;
}

} // namespace pe
} // namespace walberla
