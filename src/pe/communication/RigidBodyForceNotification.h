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
//! \file RigidBodyForceNotification.h
//! \author Tobias Preclik
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//! \brief Header file for the RigidBodyForceNotification class
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
/*!\brief Wrapper class for rigid body force and torque contribution notifications.
 *
 * The RigidBodyForceNotification class is a wrapper class for marshaling and unmarshaling rigid body
 * force and torque contribution notifications. They may only be sent by processes registered
 * to have a shadow copy. They may only be sent to the owner of the rigid body.
 */
class RigidBodyForceNotification {
public:
   struct Parameters {
      id_t sid_;
      Vec3 f_, tau_;
   };

   inline explicit RigidBodyForceNotification( const RigidBody& b ) : body_(b) {}
   const RigidBody& body_;
};
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Marshaling rigid body force notifications.
 *
 * \param buffer The buffer to be filled.
 * \param obj The object to be marshaled.
 * \return void
 *
 * The force notification consists of the force and torque contribution to a body.
 */
template< typename Buffer >
inline void marshal( Buffer& buffer, const RigidBodyForceNotification& obj ) {
   buffer << obj.body_.getSystemID();
   buffer << obj.body_.getForce();
   buffer << obj.body_.getTorque();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Unmarshaling rigid body updates.
 *
 * \param buffer The buffer from where to read.
 * \param objparam The object to be reconstructed.
 * \return void
 *
 * The update consists of the position, orientation and linear and angular velocities.
 */
template< typename Buffer >
inline void unmarshal( Buffer& buffer, RigidBodyForceNotification::Parameters& objparam ) {
   buffer >> objparam.sid_;
   buffer >> objparam.f_;
   buffer >> objparam.tau_;
}
//*************************************************************************************************

//*************************************************************************************************
/*!\brief Returns the notification type of a rigid body force notification.
 * \return The notification type of a rigid body notification.
 */
template<>
inline NotificationType notificationType<RigidBodyForceNotification>() {
   return rigidBodyForceNotification;
}
//*************************************************************************************************

}  // namespace communication
}  // namespace pe
}  // namespace walberla
