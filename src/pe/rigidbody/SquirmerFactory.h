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
//! \file SquirmerFactory.h
//! \author Michael Kuron <mkuron@icp.uni-stuttgart.de>
//
//======================================================================================================================

#pragma once

#include "pe/Materials.h"
#include "pe/Types.h"

#include "blockforest/BlockForest.h"
#include "core/debug/Debug.h"

namespace walberla {
namespace pe {

//=================================================================================================
//
//  SQUIRMER SETUP FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/**
 * \ingroup pe
 * \brief Setup of a new Squirmer.
 *
 * \param globalStorage process local global storage
 * \param blocks storage of all the blocks on this process
 * \param storageID BlockDataID of the BlockStorage block datum
 * \param uid The user-specific ID of the sphere.
 * \param gpos The global position of the center of the sphere.
 * \param radius The radius of the sphere \f$ (0..\infty) \f$.
 * \param squirmerVelocity The velocity of the squirmer.
 * \param squirmerBeta The dipolar characteristic of the squirmer.
 * \param material The material of the sphere.
 * \param global specifies if the sphere should be created in the global storage
 * \param communicating specifies if the sphere should take part in synchronization (syncNextNeighbour, syncShadowOwner)
 * \param infiniteMass specifies if the sphere has infinite mass and will be treated as an obstacle
 * \return Handle for the new sphere.
 * \exception std::invalid_argument Invalid sphere radius.
 * \exception std::invalid_argument Invalid global sphere position.
 *
 * This function creates a squirmer primitive in the \b pe simulation system. The squirmer with
 * user-specific ID \a uid is placed at the global position \a gpos, has the radius \a radius,
 * and consists of the material \a material. Its self-propulsion velocity is \a squirmerVelocity
 * and its dipolar characteristic is \a squirmerBeta.
 *
 */
SquirmerID createSquirmer( BodyStorage& globalStorage, BlockStorage& blocks, BlockDataID storageID,
                           id_t uid, const Vec3& gpos, real_t radius,
                           real_t squirmerVelocity, real_t squirmerBeta,
                           MaterialID material = Material::find("iron"),
                           bool global = false, bool communicating = true, bool infiniteMass = false );
//*************************************************************************************************

}  // namespace pe
}  // namespace walberla

