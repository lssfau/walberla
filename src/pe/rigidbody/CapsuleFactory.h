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
//! \file CapsuleFactory.h
//! \author Klaus Iglberger
//! \author Tobias Preclik
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
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
//  CAPSULE SETUP FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/**
 * \ingroup pe
 * \brief Setup of a new Capsule.
 *
 * \param globalStorage process local global storage
 * \param blocks storage of all the blocks on this process
 * \param storageID BlockDataID of the BlockStorage block datum
 * \param uid The user-specific ID of the capsule.
 * \param gpos The global position of the center of the capsule.
 * \param radius The radius of the cylinder part and the end caps \f$ (0..\infty) \f$.
 * \param length The length of the cylinder part of the capsule \f$ (0..\infty) \f$.
 * \param material The material of the capsule.
 * \param global specifies if the capsule should be created in the global storage
 * \param communicating specifies if the capsule should take part in synchronization (syncNextNeighbour, syncShadowOwner)
 * \param infiniteMass specifies if the capsule has infinite mass and will be treated as an obstacle
 * \return Handle for the new capsule.
 *
 * The following code example illustrates the setup of a capsule:

   \snippet PeDocumentationSnippets.cpp Create a Capsule
 */
CapsuleID createCapsule(   BodyStorage& globalStorage, BlockStorage& blocks, BlockDataID storageID,
                           id_t uid, const Vec3& gpos, const real_t radius, const real_t length,
                           MaterialID material = Material::find("iron"),
                           bool global = false, bool communicating = true, bool infiniteMass = false );
//*************************************************************************************************

}  // namespace pe
}  // namespace walberla
