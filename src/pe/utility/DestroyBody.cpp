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
//! \file DestroyBody.cpp
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include "DestroyBody.h"

#include "pe/rigidbody/RigidBody.h"

namespace walberla {
namespace pe {

void destroyBodyBySID(BodyStorage& globalStorage, BlockStorage& blocks, BlockDataID storageID, const walberla::id_t sid)
{
   {
      auto bodyIt = globalStorage.find(sid);
      if (bodyIt != globalStorage.end())
      {
         globalStorage.remove( bodyIt );
      }
   }

   for (auto it = blocks.begin(); it != blocks.end(); ++it)
   {
      IBlock& block = *it;
      Storage& storage  = *(block.getData<Storage>(storageID));

      {
         BodyStorage& localStorage = storage[0];
         auto bodyIt = localStorage.find(sid);
         if (bodyIt != localStorage.end())
         {
            localStorage.remove(bodyIt);
         }
      }

      {
         BodyStorage& shadowStorage = storage[1];
         auto bodyIt = shadowStorage.find(sid);
         if (bodyIt != shadowStorage.end())
         {
            shadowStorage.remove(bodyIt);
         }
      }
   }
}

void destroyBodyByUID(BodyStorage& globalStorage, BlockStorage& blocks, BlockDataID storageID, const walberla::id_t uid)
{
   destroyBody(globalStorage, blocks, storageID, internal::CheckUID(uid) );
}

}  // namespace pe
}  // namespace walberla
