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
//! \file PeIntersectionRatio.h
//! \ingroup pe_coupling
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//! \brief Specialization of the intersectionRatio function for certain pe bodies
//
//======================================================================================================================

#include "pe_coupling/geometry/PeIntersectionRatio.h"

namespace walberla {
namespace lbm {

real_t intersectionRatioSpherePe( const pe::Sphere & sphere,
                                  const Vector3<real_t> & fluidPoint,
                                  const Vector3<real_t> & direction )
{
   WALBERLA_ASSERT( !walberla::geometry::contains( sphere, fluidPoint ), "fluidPoint: " << fluidPoint );
   WALBERLA_ASSERT(  walberla::geometry::contains( sphere, fluidPoint + direction ), "fluidPoint + direction: " << fluidPoint + direction );

   // get the physical sphere center
   Vector3< real_t > sphereCenter( sphere.getPosition() );

   real_t dirLength = direction.length();
   Vector3< real_t > l = sphereCenter - fluidPoint;
   real_t s = l * direction / dirLength;
   real_t lsqr = l.sqrLength();
   real_t msqr = lsqr - s * s;
   real_t rsqr = real_c( sphere.getRadius() * sphere.getRadius() );
   real_t q = std::sqrt( rsqr - msqr );
   real_t delta = ( lsqr > rsqr ) ? s - q : s + q;
   delta /= dirLength;

   WALBERLA_ASSERT_GREATER_EQUAL( delta, real_t( 0 ) );
   WALBERLA_ASSERT_LESS_EQUAL( delta, real_t( 1 ) );

   return delta;

}


} // namespace lbm
} // namespace walberla
