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
//! \file RigidBodyRemoteMigrationNotification.h
//! \author Tobias Preclik
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//! \brief Header file for the RigidBodyRemoteMigrationNotification class
//
//======================================================================================================================

#pragma once

//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <pe/rigidbody/RigidBody.h>
#include "NotificationType.h"


namespace walberla {
namespace pe {
namespace communication {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Wrapper class for rigid body remote migration notifications.
 *
 * The RigidBodyRemoteMigrationNotification class is a wrapper class for marshaling and unmarshaling
 * rigid body remote migration notifications. When receiving a remote migration notification this
 * indicates that a body of which the local process has a shadow copy migrates from the sending
 * process to another one. Remote migration notices may only be sent by the previous owner of the
 * body.
 */
class RigidBodyRemoteMigrationNotification {
public:
   struct Parameters {
      id_t sid_;
      Owner to_;
   };

   inline RigidBodyRemoteMigrationNotification( const RigidBody& b, const Owner& to ) : body_(b), to_(to) {}
   const RigidBody& body_;
   const Owner to_;
};
//*************************************************************************************************



//*************************************************************************************************
/*!\brief Marshaling rigid body remote migration notifications.
 *
 * \param buffer The buffer to be filled.
 * \param obj The object to be marshaled.
 * \return void
 *
 * The remote migration notifications just consists of the system ID of the body migrating and the
 * rank of the process it is migrating to.
 */
template< typename Buffer >
inline void marshal( Buffer& buffer, const RigidBodyRemoteMigrationNotification& obj ) {
   // TODO fixed sizes

   buffer << obj.body_.getSystemID();
   buffer << obj.to_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Unmarshaling rigid body remote migration notifications.
 *
 * \param buffer The buffer from where to read.
 * \param objparam The object to be reconstructed.
 * \return void
 *
 * The remote migration notifications just consists of the system ID of the body migrating and the
 * rank of the process it is migrating to.
 */
template< typename Buffer >
inline void unmarshal( Buffer& buffer, RigidBodyRemoteMigrationNotification::Parameters& objparam ) {
   // TODO fixed sizes

   buffer >> objparam.sid_;
   buffer >> objparam.to_;
}
//*************************************************************************************************

//*************************************************************************************************
/*!\brief Returns the notification type of a rigid body migration.
 * \return The notification type of a rigid body migration.
 */
template<>
inline NotificationType notificationType<RigidBodyRemoteMigrationNotification>() {
   return rigidBodyRemoteMigrationNotification;
}
//*************************************************************************************************


}  // namespace communication
}  // namespace pe
}  // namespace walberla
