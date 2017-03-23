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
//! \file BoxFactory.h
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
//  BOX SETUP FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/**
 * \ingroup pe
 * \brief Setup of a new Box.
 *
 * \param globalStorage process local global storage
 * \param blocks storage of all the blocks on this process
 * \param storageID BlockDataID of the BlockStorage block datum
 * \param uid The user-specific ID of the box.
 * \param gpos The global position of the center of the box.
 * \param lengths The side length of the box \f$ (0..\infty) \f$.
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
BoxID createBox(       BodyStorage& globalStorage, BlockStorage& blocks, BlockDataID storageID,
                       id_t uid, const Vec3& gpos, const Vec3& lengths,
                       MaterialID material = Material::find("iron"),
                       bool global = false, bool communicating = true, bool infiniteMass = false );
//*************************************************************************************************

}  // namespace pe
}  // namespace walberla
