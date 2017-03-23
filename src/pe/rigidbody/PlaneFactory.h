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
//! \file PlaneFactory.h
//! \author Klaus Iglberger
//! \author Tobias Preclik
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
//  PLANE SETUP FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/**
 * \ingroup pe
 * \brief Setup of a new Plane.
 *
 * \param globalStorage process local global storage
 * \param uid The user-specific ID of the plane.
 * \param normal Normal of the plane.
 * \param gpos One point located on the plane.
 * \param material The material of the plane.
 * \return Handle for the new plane.
 *
 * This function creates a plane primitive in the simulation system. The plane with
 * user-specific ID \a uid is placed at the global position \a gpos, oriented with \a normal
 * and consists of the material \a material
 *
 * The following code example illustrates the setup of a plane:
 *
 * \snippet PeDocumentationSnippets.cpp Create a Plane
 *
 */
PlaneID createPlane( BodyStorage& globalStorage,
                     id_t uid, Vec3 normal, const Vec3& gpos,
                     MaterialID material = Material::find("iron"));
//*************************************************************************************************

}  // namespace pe
}  // namespace walberla
