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
//! \file CylindricalBoundaryFactory.cpp
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include "CylindricalBoundaryFactory.h"

#include "pe/Materials.h"
#include "pe/rigidbody/BodyStorage.h"
#include "pe/rigidbody/CylindricalBoundary.h"

#include <core/logging/Logging.h>
#include <core/UniqueID.h>

namespace walberla {
namespace pe {

CylindricalBoundaryID createCylindricalBoundary( BodyStorage& globalStorage,
                                                 id_t uid, const Vec3& gpos, const real_t radius,
                                                 MaterialID material)
{
   WALBERLA_ASSERT_UNEQUAL( CylindricalBoundary::getStaticTypeID(), std::numeric_limits<id_t>::max(), "CylindricalBoundary TypeID not initialized!");

   const id_t sid( UniqueID<RigidBody>::createGlobal() );

   CylindricalBoundaryPtr cbPtr = std::make_unique<CylindricalBoundary>( sid, uid, gpos, radius, material );

   CylindricalBoundaryID cb = static_cast<CylindricalBoundaryID>(&globalStorage.add(std::move(cbPtr)));

   // Logging the successful creation of the plane
   WALBERLA_LOG_DETAIL( "Created cylindrical boundary " << sid << "\n" << cb);

   return cb;
}

}  // namespace pe
}  // namespace walberla
