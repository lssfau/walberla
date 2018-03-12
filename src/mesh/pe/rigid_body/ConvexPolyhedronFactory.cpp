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
//! \file ConvexPolyhedronFactory.cpp
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//
//======================================================================================================================

#include "ConvexPolyhedronFactory.h"

//*************************************************************************************************
// Includes
//*************************************************************************************************

#include "core/debug/Debug.h"

#include "mesh/pe/rigid_body/ConvexPolyhedron.h"
#include "mesh/MeshOperations.h"
#include "mesh/QHull.h"
#include "mesh/TriangleMeshes.h"


#include <vector>

namespace walberla {
namespace mesh {
namespace pe {

ConvexPolyhedronID createConvexPolyhedron( BodyStorage& globalStorage, BlockStorage& blocks, BlockDataID storageID,
                                           id_t uid, Vec3 gpos, const std::vector< Vec3 > & pointCloud,
                                           MaterialID material,
                                           bool global, bool communicating, bool infiniteMass )
{
   WALBERLA_ASSERT_UNEQUAL( ConvexPolyhedron::getStaticTypeID(), std::numeric_limits<id_t>::max(), "ConvexPolyhedron TypeID not initalized!");

   // Checking the side lengths
   if( pointCloud.size() < size_t(4) )
      WALBERLA_ABORT( "Polyhedron needs at leat 4 points!" );
   
   shared_ptr< TriangleMesh > mesh = make_shared<TriangleMesh>();
   mesh::QHull<TriangleMesh> qhull( pointCloud, mesh );
   qhull.run();

   return createConvexPolyhedron( globalStorage, blocks, storageID, uid, gpos, *mesh, material, global, communicating, infiniteMass );
}



ConvexPolyhedronID createConvexPolyhedron( BodyStorage& globalStorage, BlockStorage& blocks, BlockDataID storageID,
                                           id_t uid, Vec3 gpos, TriangleMesh mesh,
                                           MaterialID material,
                                           bool global, bool communicating, bool infiniteMass )
{
   WALBERLA_ASSERT_UNEQUAL( ConvexPolyhedron::getStaticTypeID(), std::numeric_limits<id_t>::max(), "ConvexPolyhedron TypeID not initalized!");

   ConvexPolyhedronID poly = NULL;

   Vec3 centroid = toWalberla( computeCentroid( mesh ) );
   translate( mesh, -centroid );

   gpos += centroid;

   if (global)
   {
      const id_t sid = UniqueID<RigidBody>::createGlobal();
      WALBERLA_CHECK_EQUAL(communicating, false, "Global bodies can not be communicating!" );
      WALBERLA_CHECK_EQUAL(infiniteMass, true, "Global bodies must have infinite mass!" );

      poly = new ConvexPolyhedron(sid, uid, gpos, Vec3(0,0,0), Quat(), mesh, material, global, false, true);
      globalStorage.add(poly);
   } else
   {
      for (auto it = blocks.begin(); it != blocks.end(); ++it){
         IBlock* block = (&(*it));
         if (block->getAABB().contains(gpos))
         {
            const id_t sid( UniqueID<RigidBody>::create() );

            Storage* bs = block->getData<Storage>(storageID);
            poly = new ConvexPolyhedron(sid, uid, gpos, Vec3(0,0,0), Quat(), mesh, material, global, communicating, infiniteMass);
            poly->MPITrait.setOwner(Owner(MPIManager::instance()->rank(), block->getId().getID()));
            (*bs)[0].add(poly);
         }
      }
   }

   if (poly != NULL)
   {
      // Logging the successful creation of the box
      WALBERLA_LOG_DETAIL(
                "Created ConvexPolyhedron " << poly->getSystemID() << "\n"
             << "   User-ID         = " << uid << "\n"
             << "   Global position = " << gpos << "\n"
             << "   Material        = " << Material::getName( material )
               );
   }

   return poly;
}

} // namespace pe
} // namespace mesh
} // namespace walberla
