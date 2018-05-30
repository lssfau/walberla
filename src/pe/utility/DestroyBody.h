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
//! \file DestroyBody.h
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

namespace internal {
/// UnaryPredicate to check of rigid body has given uid.
class CheckUID
{
public:
   inline explicit CheckUID(const walberla::id_t uid) : uid_(uid) {}
   inline bool operator()( ConstBodyID bd) const { return (bd->getID() == uid_); }
private:
   const walberla::id_t uid_;
};
}

///Removes all rigid bodies which match the unary predicate p.
/// \attention Has to be called the same way on all processes to work correctly!
template < class UnaryPredicate >
void destroyBody(BodyStorage& globalStorage, BlockStorage& blocks, BlockDataID storageID, const UnaryPredicate& p)
{
   {
      for (auto bodyIt = globalStorage.begin(); bodyIt != globalStorage.end(); )
      {
         if ( p(bodyIt.getBodyID()) )
         {
            bodyIt = globalStorage.remove( bodyIt );
         } else
         {
            ++bodyIt;
         }
      }
   }

   for (auto it = blocks.begin(); it != blocks.end(); ++it)
   {
      IBlock& block = *it;
      Storage& storage  = *(block.getData<Storage>(storageID));

      {
         BodyStorage& localStorage = storage[StorageType::LOCAL];
         for (auto bodyIt = localStorage.begin(); bodyIt != localStorage.end(); )
         {
            if ( p(bodyIt.getBodyID()) )
            {
               bodyIt = localStorage.remove( bodyIt );
            } else
            {
               ++bodyIt;
            }
         }
      }

      {
         BodyStorage& shadowStorage = storage[StorageType::SHADOW];
         for (auto bodyIt = shadowStorage.begin(); bodyIt != shadowStorage.end(); )
         {
            if ( p(bodyIt.getBodyID()) )
            {
               bodyIt = shadowStorage.remove( bodyIt );
            } else
            {
               ++bodyIt;
            }
         }
      }
   }
}

/// removes the rigid body with the given sid
/// \attention Has to be called the same way on all processes to work correctly!
void destroyBodyBySID(BodyStorage& globalStorage, BlockStorage& blocks, BlockDataID storageID, const walberla::id_t sid );
/// removes all rigid bodies matching the given uid
/// \attention Has to be called the same way on all processes to work correctly!
void destroyBodyByUID(BodyStorage& globalStorage, BlockStorage& blocks, BlockDataID storageID, const walberla::id_t uid );

}
}
