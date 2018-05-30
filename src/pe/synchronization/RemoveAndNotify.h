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
//! \file RemoveAndNotify.h
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#pragma once


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include "pe/communication/PackNotification.h"
#include "pe/communication/RigidBodyDeletionNotification.h"
#include "pe/rigidbody/BodyStorage.h"

#include "core/mpi/BufferSystem.h"


namespace walberla {
namespace pe {

//*************************************************************************************************
/*!\brief Removes a rigid body from the local body storage and informs shadow copy holders.
 *
 * \param body An iterator pointing to the rigid body to be removed from the collision system.
 * \return An iterator pointing to the rigid body after the erase body.
 *
 * This function removes the rigid body from the body storage and generates deletion notifications.
 */
inline
BodyStorage::iterator removeAndNotify( Owner me, mpi::BufferSystem& bs, BodyStorage& localStorage, BodyStorage::iterator body )
{
   using namespace walberla::pe::communication;

   WALBERLA_ASSERT( !body->isRemote(), "Trying to remove remote body from central body storage." );

   if( !body->isGlobal() && !body->isRemote() && body->MPITrait.sizeShadowOwners() != 0 )
   {
      // Notify registered processes (intersecting or interacting) of body removal since they possess a shadow copy.
      for( auto it = body->MPITrait.beginShadowOwners(); it != body->MPITrait.endShadowOwners(); ++it ) {
         WALBERLA_LOG_DETAIL( "__Notify registered process " << (*it) << " of deletion of body " << body->getSystemID() );
         mpi::SendBuffer& sb = bs.sendBuffer(it->rank_);
         if (sb.isEmpty()) sb << walberla::uint8_c(0);
         packNotification(me, *it, sb, RigidBodyDeletionNotification( *body ));
      }
   }

   body->MPITrait.clearShadowOwners();
   return localStorage.remove( body );
}
//*************************************************************************************************

}  // namespace pe
}  // namespace walberla
