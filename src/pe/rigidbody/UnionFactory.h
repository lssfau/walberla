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
//! \file UnionFactory.h
//! \author Klaus Iglberger
//! \author Tobias Preclik
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#pragma once

#include "pe/Materials.h"
#include "pe/rigidbody/BodyStorage.h"
#include "pe/rigidbody/Box.h"
#include "pe/rigidbody/Capsule.h"
#include "pe/rigidbody/Plane.h"
#include "pe/rigidbody/Sphere.h"
#include "pe/rigidbody/Union.h"
#include "pe/Types.h"

#include "blockforest/BlockForest.h"
#include "core/debug/Debug.h"

namespace walberla {
namespace pe {

//=================================================================================================
//
//  BOX SETUP FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Setup of a new Union.
 * \ingroup rigidbody
 *
 * \param globalStorage process local global storage
 * \param blocks storage of all the blocks on this process
 * \param storageID BlockDataID of the BlockStorage block datum
 * \param uid The user-specific ID of the union.
 * \param gpos The global position of the center of the union.
 * \param global specifies if the union should be created in the global storage
 * \param communicating specifies if the union should take part in synchronization (syncNextNeighbour, syncShadowOwner)
 * \param infiniteMass specifies if the union has infinite mass and will be treated as an obstacle
 * \return Handle for the new union.
 * \exception std::invalid_argument Invalid box radius.
 * \exception std::invalid_argument Invalid global box position.
 */
template <typename BodyTypeTuple>
Union<BodyTypeTuple>* createUnion(   BodyStorage& globalStorage, BlockStorage& blocks, BlockDataID storageID,
                                     id_t uid, const Vec3& gpos,
                                     bool global = false, bool communicating = true, bool infiniteMass = false )
{
   WALBERLA_ASSERT_UNEQUAL( Union<BodyTypeTuple>::getStaticTypeID(), std::numeric_limits<id_t>::max(), "Union TypeID not initalized!");

   Union<BodyTypeTuple>* bd = NULL;

   if (global)
   {
      const id_t sid = UniqueID<RigidBody>::createGlobal();
      WALBERLA_ASSERT_EQUAL(communicating, false);
      WALBERLA_ASSERT_EQUAL(infiniteMass, true);
      bd = new Union<BodyTypeTuple>(sid, uid, gpos, Vec3(0,0,0), Quat(), global, false, true);
      globalStorage.add(bd);
   } else
   {
      for (auto it = blocks.begin(); it != blocks.end(); ++it){
         IBlock* block = (&(*it));
         if (block->getAABB().contains(gpos))
         {
            const id_t sid( UniqueID<RigidBody>::create() );

            Storage* bs = block->getData<Storage>(storageID);
            bd = new Union<BodyTypeTuple>(sid, uid, gpos, Vec3(0,0,0), Quat(), global, communicating, infiniteMass);
            bd->MPITrait.setOwner(Owner(MPIManager::instance()->rank(), block->getId().getID()));
            (*bs)[0].add(bd);
         }
      }
   }

   if (bd != NULL)
   {
      // Logging the successful creation of the box
      WALBERLA_LOG_DETAIL(
               "Created union " << bd->getSystemID() << "\n"
            << "   User-ID  = " << bd->getID()
               );
   }

   return bd;
}

template <typename BodyTypeTuple>
BoxID createBox( Union<BodyTypeTuple>* un,
                 id_t uid, const Vec3& gpos, const Vec3& lengths,
                 MaterialID material,
                 bool global, bool communicating, bool infiniteMass )
{
   if (Box::getStaticTypeID() == std::numeric_limits<id_t>::max())
      throw std::runtime_error("Box TypeID not initalized!");

   // union not on this process/block -> terminate creation
   if (un == NULL)
      throw std::invalid_argument( "createSphere: Union argument is NULL" );

   // main union not on this process/block -> terminate creation
   if ( un->isRemote() )
      throw std::logic_error( "createSphere: Union is remote" );

   // Checking the side lengths
   if( lengths[0] <= real_t(0) || lengths[1] <= real_t(0) || lengths[2] <= real_t(0) )
      throw std::invalid_argument( "Invalid side length" );

   BoxID box = NULL;
   id_t  sid = 0;

   if (global)
   {
      sid = UniqueID<RigidBody>::createGlobal();
      WALBERLA_ASSERT_EQUAL(communicating, false);
      WALBERLA_ASSERT_EQUAL(infiniteMass, true);
   } else
   {
      sid = UniqueID<RigidBody>::create();
   }

   box = new Box(sid, uid, gpos, Vec3(0,0,0), Quat(), lengths, material, global, communicating, infiniteMass);
   box->MPITrait.setOwner( un->MPITrait.getOwner() );
   un->add(box);

   if (box != NULL)
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

template <typename BodyTypeTuple>
CapsuleID createCapsule( Union<BodyTypeTuple>* un,
                         id_t uid, const Vec3& gpos, const real_t radius, const real_t length,
                         MaterialID material,
                         bool global, bool communicating, bool infiniteMass )
{
   if (Capsule::getStaticTypeID() == std::numeric_limits<id_t>::max())
      throw std::runtime_error("Capsule TypeID not initalized!");

   // union not on this process/block -> terminate creation
   if (un == NULL)
      throw std::invalid_argument( "createSphere: Union argument is NULL" );

   // main union not on this process/block -> terminate creation
   if ( un->isRemote() )
      throw std::logic_error( "createSphere: Union is remote" );

   // Checking the radius
   if( radius <= real_c(0) )
      throw std::invalid_argument( "Invalid capsule radius" );

   // Checking the length
   if( length <= real_c(0) )
      throw std::invalid_argument( "Invalid capsule length" );

   CapsuleID capsule = NULL;
   id_t      sid     = 0;

   if (global)
   {
      sid = UniqueID<RigidBody>::createGlobal();
      WALBERLA_ASSERT_EQUAL(communicating, false);
      WALBERLA_ASSERT_EQUAL(infiniteMass, true);
   } else
   {
      sid = UniqueID<RigidBody>::create();
   }

   capsule = new Capsule(sid, uid, gpos, Vec3(0,0,0), Quat(), radius, length, material, global, communicating, infiniteMass);
   capsule->MPITrait.setOwner( un->MPITrait.getOwner() );
   un->add(capsule);

   if (capsule != NULL)
   {
      WALBERLA_LOG_DETAIL("Created capsule " << capsule->getSystemID() << "\n" << capsule);
   }

   return capsule;
}

template <typename BodyTypeTuple>
SphereID createSphere( Union<BodyTypeTuple>* un,
                       id_t uid, const Vec3& gpos, real_t radius,
                       MaterialID material = Material::find("iron"),
                       bool global = false, bool communicating = true, bool infiniteMass = false )
{
   if (Sphere::getStaticTypeID() == std::numeric_limits<id_t>::max())
      throw std::runtime_error("Sphere TypeID not initalized!");

   // union not on this process/block -> terminate creation
   if (un == NULL)
      throw std::invalid_argument( "createSphere: Union argument is NULL" );

   // main union not on this process/block -> terminate creation
   if ( un->isRemote() )
      throw std::logic_error( "createSphere: Union is remote" );

   // Checking the radius
   if( radius <= real_c(0) )
      throw std::invalid_argument( "Invalid sphere radius" );

   id_t sid(0);
   SphereID sphere = NULL;

   if (global)
   {
      sid = UniqueID<RigidBody>::createGlobal();
      WALBERLA_ASSERT_EQUAL(communicating, false);
      WALBERLA_ASSERT_EQUAL(infiniteMass, true);

   } else
   {
      sid = UniqueID<RigidBody>::create();
   }

   sphere = new Sphere(sid, uid, gpos, Vec3(0,0,0), Quat(), radius, material, global, communicating, infiniteMass);
   sphere->MPITrait.setOwner( un->MPITrait.getOwner() );
   un->add( sphere );

   if (sphere != NULL)
   {
      // Logging the successful creation of the sphere
      WALBERLA_LOG_DETAIL(
                "Created sphere " << sphere->getSystemID() << " as part of union " << un->getSystemID() << "\n"
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
