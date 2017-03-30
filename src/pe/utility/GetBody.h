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
//! \file GetBody.h
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#pragma once

#include "pe/BlockFunctions.h"
#include "pe/rigidbody/BodyStorage.h"
#include "pe/Types.h"

#include "domain_decomposition/BlockStorage.h"

namespace walberla {
namespace pe {

//*************************************************************************************************
/*!\brief Tries to locate a rigid body and to return a pointer.
 *
 * Tries to locate a rigidy body by its system ID (sid). If the body is found its BodyID is returned,
 * otherwise NULL. The search space can be defined by storageType.
 *
 * \return a pointer to the rigid body if it was found, NULL otherwise
 */
BodyID getBody(BodyStorage& globalStorage, BlockStorage& blocks, BlockDataID storageID, walberla::id_t sid, int storageSelect = (StorageSelect::LOCAL | StorageSelect::SHADOW | StorageSelect::GLOBAL) );

}
}
