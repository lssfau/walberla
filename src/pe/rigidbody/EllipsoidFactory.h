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
//! \file EllipsoidFactory.h
//! \author Tobias Leemann
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
//  Ellipsoid SETUP FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/**
 * \ingroup pe
 * \brief Setup of a new Ellipsoid.
 *
 * \param globalStorage process local global storage
 * \param blocks storage of all the blocks on this process
 * \param storageID BlockDataID of the BlockStorage block datum
 * \param uid The user-specific ID of the Ellipsoid.
 * \param gpos The global position of the center of the Ellipsoid.
 * \param semiAxes The semiAxes of the Ellipsoid \f$ (0..\infty) \f$.
 * \param material The material of the Ellipsoid.
 * \param global specifies if the Ellipsoid should be created in the global storage
 * \param communicating specifies if the Ellipsoid should take part in synchronization (syncNextNeighbour, syncShadowOwner)
 * \param infiniteMass specifies if the Ellipsoid has infinite mass and will be treated as an obstacle
 * \return Handle for the new Ellipsoid.
 * \exception std::invalid_argument Invalid Ellipsoid semi-axes.
 * \exception std::invalid_argument Invalid global Ellipsoid position.
 *
 * This function creates a Ellipsoid primitive in the \b pe simulation system. The Ellipsoid with
 * user-specific ID \a uid is placed at the global position \a gpos, has the semi-axes \a semiAxes,
 * and consists of the material \a material.
 */
EllipsoidID createEllipsoid( BodyStorage& globalStorage, BlockStorage& blocks, BlockDataID storageID,
                       id_t uid, const Vec3& gpos, const Vec3& semiAxes,
                       MaterialID material = Material::find("iron"),
                       bool global = false, bool communicating = true, bool infiniteMass = false );
//*************************************************************************************************

}  // namespace pe
}  // namespace walberla
