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
//! \file BodyStorageDataHandling.h
//! \ingroup pe
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#pragma once

#include "BodyIterators.h"
#include "BodyStorage.h"
#include "pe/Types.h"

#include "pe/communication/PackNotification.h"
#include "pe/communication/RigidBodyCopyNotification.h"
#include "pe/communication/DynamicMarshalling.h"

#include "blockforest/BlockDataHandling.h"

namespace walberla{
namespace pe{


template<typename BodyTuple>
class StorageDataHandling : public domain_decomposition::BlockDataHandling<Storage>{
public:
    Storage * initialize( IBlock * const /*block*/ ) {return new Storage();}

    /// must be thread-safe !
    virtual inline void serialize( IBlock * const block, const BlockDataID & id, mpi::SendBuffer & buffer );
    /// must be thread-safe !
    virtual inline  Storage * deserialize( IBlock * const block );
    /// must be thread-safe !
    virtual inline void deserialize( IBlock * const block, const BlockDataID & id, mpi::RecvBuffer & buffer );
};

template<typename BodyTuple>
inline void StorageDataHandling<BodyTuple>::serialize( IBlock * const block, const BlockDataID & id, mpi::SendBuffer & buffer )
{
   using namespace walberla::pe::communication;
   BodyStorage& localBodyStorage = (*(block->getData< Storage >( id )))[0];
   buffer << localBodyStorage.size();
   for (auto bodyIt = localBodyStorage.begin(); bodyIt != localBodyStorage.end(); ++bodyIt)
   {
      marshal( buffer, RigidBodyCopyNotification( **bodyIt ) );
      MarshalDynamically<BodyTuple>::execute( buffer, **bodyIt );
   }
}

template<typename BodyTuple>
inline  Storage * StorageDataHandling<BodyTuple>::deserialize( IBlock * const block )
{
   return initialize(block);
}

template<typename BodyTuple>
inline void StorageDataHandling<BodyTuple>::deserialize( IBlock * const block, const BlockDataID & id, mpi::RecvBuffer & buffer )
{
   using namespace walberla::pe::communication;

   BodyStorage& localBodyStorage = (*(block->getData< Storage >( id )))[0];
   auto numBodies = localBodyStorage.size(); // just to get type right!!! can be replaced by decltype (C++11)
   buffer >> numBodies;

   while( numBodies > 0 )
   {
      typename RigidBodyCopyNotification::Parameters objparam;
      unmarshal( buffer, objparam );

      const auto inf = math::Limits<real_t>::inf();
      BodyID bd = UnmarshalDynamically<BodyTuple>::execute(buffer, objparam.geomType_, math::AABB(-inf, -inf, -inf, inf, inf, inf), block->getAABB());
      bd->setRemote( false );

      WALBERLA_ASSERT_EQUAL(localBodyStorage.find( bd->getSystemID() ), localBodyStorage.end());
      localBodyStorage.add(bd);

      --numBodies;
   }
}

template <typename BodyTuple>
shared_ptr<StorageDataHandling<BodyTuple> > createStorageDataHandling()
{
   return make_shared<StorageDataHandling<BodyTuple> >( );
}

}
}
