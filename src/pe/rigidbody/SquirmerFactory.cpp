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
//! \file SquirmerFactory.cpp
//! \author Michael Kuron <mkuron@icp.uni-stuttgart.de>
//
//======================================================================================================================

#include "SquirmerFactory.h"

#include "pe/Materials.h"
#include "pe/rigidbody/BodyStorage.h"
#include "pe/rigidbody/Squirmer.h"

#include <core/logging/Logging.h>
#include <core/UniqueID.h>

namespace walberla {
namespace pe {

SquirmerID createSquirmer( BodyStorage& globalStorage, BlockStorage& blocks, BlockDataID storageID,
                        id_t uid, const Vec3& gpos, real_t radius,
                        real_t squirmerVelocity, real_t squirmerBeta,
                        MaterialID material,
                        bool global, bool communicating, bool infiniteMass )
{
   WALBERLA_ASSERT_UNEQUAL( Squirmer::getStaticTypeID(), std::numeric_limits<id_t>::max(), "Squirmer TypeID not initalized!");
   // Checking the radius
   if( radius <= real_c(0) )
      throw std::invalid_argument( "Invalid squirmer radius" );

   SquirmerID squirmer = nullptr;

   if (global)
   {
      const id_t sid = UniqueID<RigidBody>::createGlobal();
      WALBERLA_ASSERT_EQUAL(communicating, false);
      WALBERLA_ASSERT_EQUAL(infiniteMass, true);
      SquirmerPtr sq = std::make_unique<Squirmer>(sid, uid, gpos, Quat(), radius, squirmerVelocity, squirmerBeta, material, global, false, true);
      squirmer = static_cast<SquirmerID>(&globalStorage.add(std::move(sq)));
   } else
   {
      for (auto& block : blocks){
         if (block.getAABB().contains(gpos))
         {
            const id_t sid( UniqueID<RigidBody>::create() );

            BodyStorage& bs = (*block.getData<Storage>(storageID))[0];
            SquirmerPtr sq = std::make_unique<Squirmer>(sid, uid, gpos, Quat(), radius, squirmerVelocity, squirmerBeta, material, global, communicating, infiniteMass);
            sq->MPITrait.setOwner(Owner(MPIManager::instance()->rank(), block.getId().getID()));
            squirmer = static_cast<SquirmerID>(&bs.add( std::move(sq) ));
         }
      }
   }

   if (squirmer != nullptr)
   {
      // Logging the successful creation of the squirmer
      WALBERLA_LOG_DETAIL(
                "Created squirmer " << squirmer->getSystemID() << "\n"
             << "   User-ID         = " << uid << "\n"
             << "   Global position = " << gpos << "\n"
             << "   Radius          = " << radius << "\n"
             << "   LinVel          = " << squirmer->getLinearVel() << "\n"
             << "   Material        = " << Material::getName( material )
               );
   }

   return squirmer;
}

}  // namespace pe
}  // namespace walberla

