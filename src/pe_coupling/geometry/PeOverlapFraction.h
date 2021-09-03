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
//! \file PeOverlapFraction.h
//! \ingroup pe_coupling
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//
//======================================================================================================================

#pragma once

#include "geometry/bodies/BodyOverlapFunctions.h"

#include "pe/rigidbody/RigidBody.h"
#include "pe/rigidbody/Sphere.h"
#include "pe/rigidbody/Squirmer.h"
#include "pe/rigidbody/Plane.h"

#include "PeBodyOverlapFunctions.h"

namespace walberla {
namespace pe_coupling{


inline real_t overlapFractionPe( const pe::RigidBody & peRigidBody, const Vector3<real_t> & cellMidpoint,
                                 const Vector3<real_t> & dx, uint_t maxDepth=4 )
{
   if( peRigidBody.getTypeID() == pe::Sphere::getStaticTypeID() || peRigidBody.getTypeID() == pe::Squirmer::getStaticTypeID() )
   {
      const pe::Sphere & sphere = static_cast< const pe::Sphere & >( peRigidBody );
      return geometry::overlapFraction( sphere, cellMidpoint, dx, maxDepth );
   }
   if( peRigidBody.getTypeID() == pe::Plane::getStaticTypeID() )
   {
      const pe::Plane & plane = static_cast< const pe::Plane & >( peRigidBody );
      return geometry::overlapFraction( plane, cellMidpoint, dx, maxDepth );
   }

   // Add more pe bodies here if specific fastOverlapCheck(...) and contains(...) function is available
   // else: fallback to default implementation
   return geometry::overlapFraction( peRigidBody, cellMidpoint, dx, maxDepth );

}


} // namespace pe_coupling
} // namespace walberla
