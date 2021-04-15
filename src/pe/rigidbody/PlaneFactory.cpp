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
//! \file PlaneFactory.cpp
//! \author Klaus Iglberger
//! \author Tobias Preclik
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include "PlaneFactory.h"

#include "pe/Materials.h"
#include "pe/rigidbody/BodyStorage.h"
#include "pe/rigidbody/Plane.h"

#include <core/logging/Logging.h>
#include <core/UniqueID.h>

namespace walberla {
namespace pe {

PlaneID createPlane( BodyStorage& globalStorage, id_t uid, Vec3 normal, const Vec3& gpos, MaterialID material)
{
   WALBERLA_ASSERT_UNEQUAL( Plane::getStaticTypeID(), std::numeric_limits<id_t>::max(), "Plane TypeID not initialized!");

   // Checking the normal of the plane
   if( floatIsEqual(normal.sqrLength(), real_c(0) ) )
      throw std::invalid_argument( "Invalid plane normal!" );

   // Normalizing the plane normal
   normal = normal.getNormalized();

   const id_t sid( UniqueID<RigidBody>::createGlobal() );

   PlanePtr pl = std::make_unique<Plane>( sid, uid, gpos, normal, normal*gpos, material );
   PlaneID plane = static_cast<PlaneID>(&globalStorage.add(std::move(pl)));

   // Logging the successful creation of the plane
   WALBERLA_LOG_DETAIL( "Created plane " << sid << "\n"
          << "   User-ID      = " << uid << "\n"
          << "   Anchor point = " << gpos << "\n"
          << "   Normal       = " << normal << "\n"
          << "   Displacement = " << plane->getDisplacement() << "\n"
          << "   Material     = " << Material::getName( material ));

   return plane;
}

}  // namespace pe
}  // namespace walberla
