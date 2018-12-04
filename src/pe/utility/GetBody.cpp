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
//! \file GetBody.cpp
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include "GetBody.h"

namespace walberla {
namespace pe {

BodyID getBody(BodyStorage& globalStorage, BlockStorage& blocks, BlockDataID storageID, walberla::id_t sid, int storageSelect)
{
   if (storageSelect & StorageSelect::GLOBAL)
   {
      auto bodyIt = globalStorage.find(sid);
      if (bodyIt != globalStorage.end())
      {
         return bodyIt.getBodyID();
      }
   }

   for (auto it = blocks.begin(); it != blocks.end(); ++it)
   {
      IBlock& block = *it;
      Storage& storage  = *(block.getData<Storage>(storageID));
      if (storageSelect & StorageSelect::LOCAL)
      {
         BodyStorage& localStorage = storage[0];
         auto bodyIt = localStorage.find(sid);
         if (bodyIt != localStorage.end())
         {
            return bodyIt.getBodyID();
         }
      }
      if (storageSelect & StorageSelect::SHADOW)
      {
         BodyStorage& shadowStorage = storage[1];
         auto bodyIt = shadowStorage.find(sid);
         if (bodyIt != shadowStorage.end())
         {
            return bodyIt.getBodyID();
         }
      }
   }

   return nullptr;
}

}  // namespace pe
}  // namespace walberla
