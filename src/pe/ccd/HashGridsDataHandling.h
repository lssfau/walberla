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
//! \file HashGridsDataHandling.h
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#pragma once

#include "HashGrids.h"

#include "pe/rigidbody/BodyStorage.h"

#include "blockforest/BlockDataHandling.h"

namespace walberla{
namespace pe{
namespace ccd {

class HashGridsDataHandling : public blockforest::AlwaysInitializeBlockDataHandling<HashGrids>{
public:
   HashGridsDataHandling(const shared_ptr<BodyStorage>& globalStorage, const BlockDataID& storageID) : globalStorage_(globalStorage), storageID_(storageID) {}
   HashGrids * initialize( IBlock * const block ) override
   {
      Storage* storage = block->getData< Storage >( storageID_ );
      return new HashGrids(*globalStorage_, (*storage)[0], (*storage)[1]);
   }
private:
   shared_ptr<BodyStorage>  globalStorage_;
   BlockDataID              storageID_;
};

inline
shared_ptr<HashGridsDataHandling> createHashGridsDataHandling(const shared_ptr<BodyStorage>& globalStorage,const BlockDataID& storageID)
{
   return make_shared<HashGridsDataHandling>( globalStorage, storageID );
}

}
}
}
