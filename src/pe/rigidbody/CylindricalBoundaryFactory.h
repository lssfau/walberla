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
//! \file CylindricalBoundaryFactory.h
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#pragma once

#include "pe/rigidbody/BodyStorage.h"
#include "pe/Materials.h"
#include "pe/Types.h"

#include "core/debug/Debug.h"

namespace walberla {
namespace pe {

//=================================================================================================
//
//  CYLINDRICAL BOUNDARY SETUP FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/**
 * \ingroup pe
 * \brief Setup of a new Cylindrical Boundary.
 *
 * \param globalStorage process local global storage
 * \param uid The user-specific ID.
 * \param gpos One point located on the central axis.
 * \param radius radius of the cylinder
 * \param material The material of the boundary.
 * \return Handle for the new boundary.
 */
CylindricalBoundaryID createCylindricalBoundary( BodyStorage& globalStorage,
                     id_t uid, const Vec3& gpos, const real_t radius,
                     MaterialID material = Material::find("iron"));
//*************************************************************************************************

}  // namespace pe
}  // namespace walberla
