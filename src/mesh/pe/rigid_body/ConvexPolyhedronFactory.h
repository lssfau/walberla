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
//! \file ConvexPolyhedronFactory.h
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//
//======================================================================================================================

#pragma once

//*************************************************************************************************
// Includes
//*************************************************************************************************

#include "ConvexPolyhedron.h"

#include "domain_decomposition/BlockStorage.h"

#include "pe/rigidbody/BodyStorage.h"
#include "pe/Materials.h"

namespace walberla {
namespace mesh {
namespace pe {

//=================================================================================================
//
//  CONVEXPOLYHEDRON SETUP FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/**
 * \ingroup pe
 * \brief Setup of a new ConvexPolyhedron.
 *
 * \param globalStorage process local global storage
 * \param blocks storage of all the blocks on this process
 * \param storageID BlockDataID of the BlockStorage block datum
 * \param uid The user-specific ID of the box.
 * \param gpos The global position of the center of the box.
 * \param pointCloud A point cloud which convex hull defines the polyhedron
 * \param material The material of the box.
 * \param global specifies if the box should be created in the global storage
 * \param communicating specifies if the box should take part in synchronization (syncNextNeighbour, syncShadowOwner)
 * \param infiniteMass specifies if the box has infinite mass and will be treated as an obstacle
 * \return Handle for the new box.
 * \exception std::invalid_argument Invalid box radius.
 * \exception std::invalid_argument Invalid global box position.
 *
 * This function creates a box primitive in the \b pe simulation system. The box with
 * user-specific ID \a uid is placed at the global position \a gpos, has the side lengths \a lengths,
 * and consists of the material \a material.
 *
 * The following code example illustrates the setup of a box:
 * \snippet PeDocumentationSnippets.cpp Create a Box
 *
 */
ConvexPolyhedronID createConvexPolyhedron( BodyStorage& globalStorage, BlockStorage& blocks, BlockDataID storageID,
                                           id_t uid, const Vec3& gpos, const std::vector< Vec3 > & pointCloud,
                                           MaterialID material = Material::find("iron"),
                                           bool global = false, bool communicating = true, bool infiniteMass = false );
//*************************************************************************************************


//*************************************************************************************************
/**
 * \ingroup pe
 * \brief Setup of a new ConvexPolyhedron.
 *
 * \param globalStorage process local global storage
 * \param blocks storage of all the blocks on this process
 * \param storageID BlockDataID of the BlockStorage block datum
 * \param uid The user-specific ID of the box.
 * \param gpos The global position of the center of the box.
 * \param mesh Surface mesh of convex polyhedron
 * \param material The material of the box.
 * \param global specifies if the box should be created in the global storage
 * \param communicating specifies if the box should take part in synchronization (syncNextNeighbour, syncShadowOwner)
 * \param infiniteMass specifies if the box has infinite mass and will be treated as an obstacle
 * \return Handle for the new box.
 * \exception std::invalid_argument Invalid box radius.
 * \exception std::invalid_argument Invalid global box position.
 *
 * This function creates a box primitive in the \b pe simulation system. The box with
 * user-specific ID \a uid is placed at the global position \a gpos, has the side lengths \a lengths,
 * and consists of the material \a material.
 *
 * The following code example illustrates the setup of a box:
 * \snippet PeDocumentationSnippets.cpp Create a Box
 *
 */
ConvexPolyhedronID createConvexPolyhedron( BodyStorage& globalStorage, BlockStorage& blocks, BlockDataID storageID,
                                           id_t uid, Vec3 gpos, TriangleMesh mesh,
                                           MaterialID material = Material::find("iron"),
                                           bool global = false, bool communicating = true, bool infiniteMass = false );
//*************************************************************************************************


} // namespace pe
} // namespace mesh
} // namespace walberla
