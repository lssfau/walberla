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

#pragma once

#include "lbm/geometry/IntersectionRatio.h"
#include "pe_coupling/geometry/PeBodyOverlapFunctions.h"

#include "pe/rigidbody/RigidBody.h"
#include "pe/rigidbody/Ellipsoid.h"
#include "pe/rigidbody/Plane.h"
#include "pe/rigidbody/Sphere.h"
#include "pe/rigidbody/Squirmer.h"
#include "pe/Types.h"


namespace walberla {
namespace lbm {

real_t intersectionRatioSpherePe( const pe::Sphere & sphere,
                                  const Vector3<real_t> & fluidPoint,
                                  const Vector3<real_t> & direction );

real_t intersectionRatioPlanePe( const pe::Plane & plane,
                                 const Vector3<real_t> & fluidPoint,
                                 const Vector3<real_t> & direction );

real_t intersectionRatioEllipsoidPe( const pe::Ellipsoid & ellipsoid,
                                     const Vector3<real_t> & fluidPoint,
                                     const Vector3<real_t> & direction );

inline real_t intersectionRatio( const pe::RigidBody & peRigidBody,
                                 const Vector3<real_t> & fluidPoint,
                                 const Vector3<real_t> & direction,
                                 const real_t epsilon )
{
   if( peRigidBody.getTypeID() == pe::Sphere::getStaticTypeID() || peRigidBody.getTypeID() == pe::Squirmer::getStaticTypeID() )
   {
      const pe::Sphere & sphere = static_cast< const pe::Sphere & >( peRigidBody );
      const real_t ratio = intersectionRatioSpherePe( sphere, fluidPoint, direction );
      WALBERLA_ASSERT_LESS_EQUAL( std::fabs( ( fluidPoint + ratio * direction - sphere.getPosition() ).length() - sphere.getRadius() ), epsilon );
      return ratio;
   }
   else if ( peRigidBody.getTypeID() == pe::Plane::getStaticTypeID() )
   {
      const pe::Plane & plane = static_cast< const pe::Plane & >( peRigidBody );
      const real_t ratio = intersectionRatioPlanePe( plane, fluidPoint, direction );
      WALBERLA_ASSERT_FLOAT_EQUAL( ( fluidPoint + ratio * direction - plane.getPosition() ) * plane.getNormal(), real_t(0) );
      return ratio;
   }
   else if ( peRigidBody.getTypeID() == pe::Ellipsoid::getStaticTypeID() )
   {
      const pe::Ellipsoid & ellipsoid = static_cast< const pe::Ellipsoid & >( peRigidBody );
      const real_t ratio = intersectionRatioEllipsoidPe( ellipsoid, fluidPoint, direction );
      return ratio;
   }
   // Add more pe bodies here if specific intersectionRatio(...) function is available
   else
   {
      return lbm::intersectionRatio< pe::RigidBody >( peRigidBody, fluidPoint, direction, epsilon );
   }
}


} // namespace lbm
} // namespace walberla
