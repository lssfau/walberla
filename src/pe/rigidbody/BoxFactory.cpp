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
//! \file BoxFactory.cpp
//! \author Klaus Iglberger
//! \author Tobias Preclik
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include "BoxFactory.h"

#include "pe/Materials.h"
#include "pe/rigidbody/BodyStorage.h"
#include "pe/rigidbody/Box.h"

#include <core/logging/Logging.h>
#include <core/UniqueID.h>

namespace walberla {
namespace pe {

BoxID createBox(       BodyStorage& globalStorage, BlockStorage& blocks, BlockDataID storageID,
                       id_t uid, const Vec3& gpos, const Vec3& lengths,
                       MaterialID material,
                       bool global, bool communicating, bool infiniteMass )
{
   WALBERLA_ASSERT_UNEQUAL( Box::getStaticTypeID(), std::numeric_limits<id_t>::max(), "Box TypeID not initalized!");

   // Checking the side lengths
   if( lengths[0] <= real_t(0) || lengths[1] <= real_t(0) || lengths[2] <= real_t(0) )
      throw std::invalid_argument( "Invalid side length" );

   BoxID box = nullptr;

   if (global)
   {
      const id_t sid = UniqueID<RigidBody>::createGlobal();
      WALBERLA_ASSERT_EQUAL(communicating, false);
      WALBERLA_ASSERT_EQUAL(infiniteMass, true);
      BoxPtr bx = std::make_unique<Box>(sid, uid, gpos, Quat(), lengths, material, global, false, true);
      box = static_cast<BoxID>(&globalStorage.add(std::move(bx)));
   } else
   {
      for (auto& block : blocks){
         if (block.getAABB().contains(gpos))
         {
            const id_t sid( UniqueID<RigidBody>::create() );

            BodyStorage& bs = (*block.getData<Storage>(storageID))[0];
            BoxPtr bx = std::make_unique<Box>(sid, uid, gpos, Quat(), lengths, material, global, communicating, infiniteMass);
            bx->MPITrait.setOwner(Owner(MPIManager::instance()->rank(), block.getId().getID()));
            box = static_cast<BoxID>(&bs.add(std::move(bx)));
         }
      }
   }

   if (box != nullptr)
   {
      // Logging the successful creation of the box
      WALBERLA_LOG_DETAIL(
                "Created box " << box->getSystemID() << "\n"
             << "   User-ID         = " << uid << "\n"
             << "   Global position = " << gpos << "\n"
             << "   side length     = " << lengths << "\n"
             << "   LinVel          = " << box->getLinearVel() << "\n"
             << "   Material        = " << Material::getName( material )
               );
   }

   return box;
}

}  // namespace pe
}  // namespace walberla
