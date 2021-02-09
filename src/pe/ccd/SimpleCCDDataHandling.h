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
//! \file SimpleCCDDataHandling.h
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#pragma once

#include "SimpleCCD.h"

#include "pe/rigidbody/BodyStorage.h"

#include "blockforest/BlockDataHandling.h"

namespace walberla{
namespace pe{
namespace ccd {

class SimpleCCDDataHandling : public blockforest::AlwaysInitializeBlockDataHandling<SimpleCCD>{
public:
   SimpleCCDDataHandling(const shared_ptr<BodyStorage>& globalStorage, const BlockDataID& storageID) : globalStorage_(globalStorage), storageID_(storageID) {}
   SimpleCCD * initialize( IBlock * const block ) override
   {
      Storage* storage = block->getData< Storage >( storageID_ );
      return new SimpleCCD(*globalStorage_, *storage);
   }
private:
   shared_ptr<BodyStorage> globalStorage_;
   BlockDataID storageID_;
};

inline
shared_ptr<SimpleCCDDataHandling> createSimpleCCDDataHandling(const shared_ptr<BodyStorage>& globalStorage,const BlockDataID& storageID)
{
   return make_shared<SimpleCCDDataHandling>( globalStorage, storageID );
}

}
}
}
