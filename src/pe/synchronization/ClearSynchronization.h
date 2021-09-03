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
//! \file ClearSynchronization.h
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#pragma once

//*************************************************************************************************
// Includes
//*************************************************************************************************

#include "pe/rigidbody/BodyStorage.h"

#include "domain_decomposition/BlockStorage.h"

namespace walberla {
namespace pe {

//*************************************************************************************************
/*!\brief Removes all synchronization information.
 *
 * \param blocks BlockStorage of the simulation.
 * \param storageID BlockDataID of the pe storage.
 */
inline
void clearSynchronization( BlockStorage& blocks, BlockDataID storageID )
{
   for (auto it = blocks.begin(); it != blocks.end(); ++it)
   {
      IBlock & currentBlock      = *it;
      Storage * storage          = currentBlock.getData< Storage >( storageID );
      BodyStorage& localStorage  = (*storage)[0];
      BodyStorage& shadowStorage = (*storage)[1];

      for( auto bodyIt = localStorage.begin(); bodyIt != localStorage.end(); ++bodyIt )
      {
         bodyIt->MPITrait.clearShadowOwners();
         bodyIt->MPITrait.clearBlockStates();
      }

      shadowStorage.clear();
   }
}
//*************************************************************************************************

}  // namespace pe
}  // namespace walberla
