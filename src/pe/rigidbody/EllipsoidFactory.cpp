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
//! \file EllipsoidFactory.cpp
//! \author Tobias Leemann <tobias.leemann@fau.de>
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include "EllipsoidFactory.h"

#include "pe/Materials.h"
#include "pe/rigidbody/BodyStorage.h"
#include "pe/rigidbody/Ellipsoid.h"

#include <core/logging/Logging.h>
#include <core/UniqueID.h>

namespace walberla {
namespace pe {

EllipsoidID createEllipsoid( BodyStorage& globalStorage, BlockStorage& blocks, BlockDataID storageID,
                             id_t uid, const Vec3& gpos, const Vec3& semiAxes,
                             MaterialID material,
                             bool global, bool communicating, bool infiniteMass )
{
   WALBERLA_ASSERT_UNEQUAL( Ellipsoid::getStaticTypeID(), std::numeric_limits<id_t>::max(), "Ellipsoid TypeID not initalized!");
   // Checking the semiAxes
   if( semiAxes[0] <= real_c(0) || semiAxes[1] <= real_c(0) || semiAxes[2] <= real_c(0) )
      throw std::invalid_argument( "Invalid Ellipsoid semi-axes" );

   EllipsoidID ellipsoid = nullptr;

   if (global)
   {
      const id_t sid = UniqueID<RigidBody>::createGlobal();
      WALBERLA_ASSERT_EQUAL(communicating, false);
      WALBERLA_ASSERT_EQUAL(infiniteMass, true);
      EllipsoidPtr el = std::make_unique<Ellipsoid>(sid, uid, gpos, Quat(), semiAxes, material, global, false, true);
      ellipsoid = static_cast<EllipsoidID>(&globalStorage.add(std::move(el)));
   } else
   {
      for (auto& block : blocks){
         if (block.getAABB().contains(gpos))
         {
            const id_t sid( UniqueID<RigidBody>::create() );

            BodyStorage& bs = (*block.getData<Storage>(storageID))[0];
            EllipsoidPtr el = std::make_unique<Ellipsoid>(sid, uid, gpos, Quat(), semiAxes, material, global, communicating, infiniteMass);
            el->MPITrait.setOwner(Owner(MPIManager::instance()->rank(), block.getId().getID()));
            ellipsoid = static_cast<EllipsoidID>(&bs.add(std::move(el)));
         }
      }
   }

   if (ellipsoid != nullptr)
   {
      // Logging the successful creation of the Ellipsoid
      WALBERLA_LOG_DETAIL(
               "Created Ellipsoid " << ellipsoid->getSystemID() << "\n"
               << "   User-ID         = " << uid << "\n"
               << "   Global position = " << gpos << "\n"
               << "   Semi-axes       = " << semiAxes << "\n"
               << "   LinVel          = " << ellipsoid->getLinearVel() << "\n"
               << "   Material        = " << Material::getName( material )
               );
   }

   return ellipsoid;
}

}  // namespace pe
}  // namespace walberla
