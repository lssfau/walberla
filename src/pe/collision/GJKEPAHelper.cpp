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
//! \file GJKHelper.cpp
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

//*************************************************************************************************
// Includes
//*************************************************************************************************

#include "GJKEPAHelper.h"

#include "pe/rigidbody/GeomPrimitive.h"

extern "C" {
   #include "pe/extern/libccd/ccd/ccd.h"
   #include "pe/extern/libccd/ccd/quat.h"
}

#include "core/logging/Logging.h"

namespace walberla {
namespace pe {

Vec3       convertVec3(const ccd_vec3_t& vec) { return Vec3(real_c(vec.v[0]), real_c(vec.v[1]), real_c(vec.v[2])); }
ccd_vec3_t convertVec3(const Vec3& vec)       { ccd_vec3_t ret; ret.v[0] = vec[0]; ret.v[1] = vec[1]; ret.v[2] = vec[2]; return ret; }

void support(const void *obj, const ccd_vec3_t *dir, ccd_vec3_t *vec)
{
    ConstGeomID bd = reinterpret_cast<ConstGeomID> (obj);
    Vec3 d = convertVec3(*dir);
    Vec3 sup = bd->support( d );
    *vec = convertVec3(sup);
}

bool collideGJK( ConstGeomID bd1,
                 ConstGeomID bd2,
                 Vec3& contactPoint,
                 Vec3& contactNormal,
                 real_t& penetrationDepth,
                 const unsigned long numIterations,
                 const real_t epaTol )
{
    ccd_t ccd;
    CCD_INIT(&ccd); // initialize ccd_t struct

    // set up ccd_t struct
    ccd.support1       = support;       // support function for first object
    ccd.support2       = support;       // support function for second object
    ccd.max_iterations = numIterations; // maximal number of iterations
    ccd.epa_tolerance  = epaTol;

    ccd_vec3_t dir, pos;
    ccd_real_t penetrationDepthCCD;
    int intersect = ccdGJKPenetration(reinterpret_cast<const void*> (bd1), reinterpret_cast<const void*> (bd2), &ccd, &penetrationDepthCCD, &dir, &pos);
    penetrationDepth = real_c(penetrationDepthCCD);
    contactPoint  = convertVec3(pos);
    contactNormal = -convertVec3(dir);
    penetrationDepth *= -1;

    return (intersect == 0);
}

} // namespace pe
} // namespace walberla
