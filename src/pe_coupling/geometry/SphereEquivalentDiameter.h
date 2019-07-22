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
//! \file SphereEquivalentDiameter.h
//! \ingroup pe_coupling
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//
//======================================================================================================================

#pragma once


#include "pe/rigidbody/RigidBody.h"
#include "pe/rigidbody/Sphere.h"
#include "pe/rigidbody/Squirmer.h"

namespace walberla {
namespace pe_coupling {


// calculates sphere-equivalent diameter (diameter of a sphere with same volume as given body)
real_t getSphereEquivalentDiameter( pe::RigidBody & body )
{
   if( body.getTypeID() == pe::Sphere::getStaticTypeID() || body.getTypeID() == pe::Squirmer::getStaticTypeID() )
   {
      pe::Sphere & sphere = static_cast<pe::Sphere &>( body );
      real_t radius = sphere.getRadius();
      return real_t(2) * radius;
   } else {
      const real_t preFac = real_t(6) / math::pi;
      return std::cbrt( body.getVolume() * preFac );
   }
}

} // namespace pe_coupling
} // namespace walberla
