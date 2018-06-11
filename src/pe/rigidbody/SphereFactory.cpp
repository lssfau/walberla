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
//! \file SphereFactory.cpp
//! \author Klaus Iglberger
//! \author Tobias Preclik
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include "SphereFactory.h"

#include "pe/Materials.h"
#include "pe/rigidbody/BodyStorage.h"
#include "pe/rigidbody/Sphere.h"

#include <core/logging/Logging.h>
#include <core/UniqueID.h>

namespace walberla {
namespace pe {

SphereID createSphere( BodyStorage& globalStorage, BlockStorage& blocks, BlockDataID storageID,
                       id_t uid, const Vec3& gpos, real_t radius,
                       MaterialID material,
                       bool global, bool communicating, bool infiniteMass )
{
   WALBERLA_ASSERT_UNEQUAL( Sphere::getStaticTypeID(), std::numeric_limits<id_t>::max(), "Sphere TypeID not initalized!");
   // Checking the radius
   if( radius <= real_c(0) )
      throw std::invalid_argument( "Invalid sphere radius" );

   SphereID sphere = NULL;

   if (global)
   {
      const id_t sid = UniqueID<RigidBody>::createGlobal();
      WALBERLA_ASSERT_EQUAL(communicating, false);
      WALBERLA_ASSERT_EQUAL(infiniteMass, true);
      sphere = new Sphere(sid, uid, gpos, Vec3(0,0,0), Quat(), radius, material, global, false, true);
      globalStorage.add(sphere);
   } else
   {
      for (auto it = blocks.begin(); it != blocks.end(); ++it){
         IBlock* block = (&(*it));
         if (block->getAABB().contains(gpos))
         {
            const id_t sid( UniqueID<RigidBody>::create() );

            Storage* bs = block->getData<Storage>(storageID);
            sphere = new Sphere(sid, uid, gpos, Vec3(0,0,0), Quat(), radius, material, global, communicating, infiniteMass);
            sphere->MPITrait.setOwner(Owner(MPIManager::instance()->rank(), block->getId().getID()));
            (*bs)[0].add(sphere);
         }
      }
   }

   if (sphere != NULL)
   {
      // Logging the successful creation of the sphere
      WALBERLA_LOG_DETAIL(
                "Created sphere " << sphere->getSystemID() << "\n"
             << "   User-ID         = " << uid << "\n"
             << "   Global position = " << gpos << "\n"
             << "   Radius          = " << radius << "\n"
             << "   LinVel          = " << sphere->getLinearVel() << "\n"
             << "   Material        = " << Material::getName( material )
               );
   }

   return sphere;
}

}  // namespace pe
}  // namespace walberla
