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
//! \file CapsuleFactory.cpp
//! \author Klaus Iglberger
//! \author Tobias Preclik
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include "BoxFactory.h"

#include "pe/Materials.h"
#include "pe/rigidbody/BodyStorage.h"
#include "pe/rigidbody/Capsule.h"

#include <core/logging/Logging.h>
#include <core/UniqueID.h>

namespace walberla {
namespace pe {

CapsuleID createCapsule(   BodyStorage& globalStorage, BlockStorage& blocks, BlockDataID storageID,
                           id_t uid, const Vec3& gpos, const real_t radius, const real_t length,
                           MaterialID material,
                           bool global, bool communicating, bool infiniteMass )
{
   WALBERLA_ASSERT_UNEQUAL( Capsule::getStaticTypeID(), std::numeric_limits<id_t>::max(), "Capsule TypeID not initialized!");

   // Checking the radius and the length
   WALBERLA_ASSERT_GREATER( radius, real_t(0), "Invalid capsule radius" );
   WALBERLA_ASSERT_GREATER( length, real_t(0), "Invalid capsule length" );

   CapsuleID capsule = nullptr;

   if (global)
   {
      const id_t sid = UniqueID<RigidBody>::createGlobal();
      WALBERLA_ASSERT_EQUAL(communicating, false);
      WALBERLA_ASSERT_EQUAL(infiniteMass, true);
      CapsulePtr cp = std::make_unique<Capsule>(sid, uid, gpos, Quat(), radius, length, material, global, false, true);
      capsule = static_cast<CapsuleID>(&globalStorage.add(std::move(cp)));
   } else
   {
      for (auto& block : blocks){
         if (block.getAABB().contains(gpos))
         {
            const id_t sid( UniqueID<RigidBody>::create() );

            BodyStorage& bs = (*block.getData<Storage>(storageID))[0];
            CapsulePtr cp = std::make_unique<Capsule>(sid, uid, gpos, Quat(), radius, length, material, global, communicating, infiniteMass);
            cp->MPITrait.setOwner(Owner(MPIManager::instance()->rank(), block.getId().getID()));
            capsule = static_cast<CapsuleID>(&bs.add( std::move(cp) ));
         }
      }
   }

   if (capsule != nullptr)
   {
      // Logging the successful creation of the box
      WALBERLA_LOG_DETAIL(
               "Created capsule " << capsule->getSystemID() << "\n"
               << "   User-ID         = " << uid << "\n"
               << "   Global position = " << gpos << "\n"
               << "   Radius          = " << radius << "\n"
               << "   Length          = " << length << "\n"
               << "   Material        = " << Material::getName( material )
               );
   }

   return capsule;
}

}  // namespace pe
}  // namespace walberla
