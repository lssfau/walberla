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
//! \file RigidBodyNewShadowCopyNotification.h
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//! \brief Header file for the RigidBodyNewShadowCopyNotification class
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
/*!\brief Wrapper class for rigid body new shadow copy notifications.
 *
 * The RigidBodyNewShadowCopyNotification class is a wrapper class for marshaling and unmarshaling
 * of rigid body new shadow copy notifications. This notification is send to the owner of a rigid
 * body to signal that a new shadow copy exists and the shadow copy holder list should be updated.
 */
class RigidBodyNewShadowCopyNotification {
public:
   struct Parameters {
      id_t sid_;
      Owner newOwner_;
   };

   inline RigidBodyNewShadowCopyNotification( const RigidBody& b, const Owner& newOwner ) : body_(b), newOwner_(newOwner) {}
   const RigidBody& body_;
   const Owner& newOwner_;
};
//*************************************************************************************************



//*************************************************************************************************
/*!\brief Marshaling rigid body new shadow copy notifications.
 *
 * \param buffer The buffer to be filled.
 * \param obj The object to be marshaled.
 * \return void
 *
 * The new shadow copy notifications just consists of the system ID of the body to remove.
 */
template< typename Buffer >
inline void marshal( Buffer& buffer, const RigidBodyNewShadowCopyNotification& obj ) {
   buffer << obj.body_.getSystemID();
   buffer << obj.newOwner_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Unmarshaling rigid body new shadow copy notifications.
 *
 * \param buffer The buffer from where to read.
 * \param objparam The object to be reconstructed.
 * \return void
 *
 * The new shadow copy notifications just consists of the system ID of the body to remove.
 */
template< typename Buffer >
inline void unmarshal( Buffer& buffer, RigidBodyNewShadowCopyNotification::Parameters& objparam ) {
   buffer >> objparam.sid_;
   buffer >> objparam.newOwner_;
}
//*************************************************************************************************

//*************************************************************************************************
/*!\brief Returns the notification type of a rigid body new shadow copy notification.
 * \return The notification type of a rigid body new shadow copy notification.
 */
template<>
inline NotificationType notificationType<RigidBodyNewShadowCopyNotification>() {
   return rigidBodyNewShadowCopyNotification;
}
//*************************************************************************************************


}  // namespace communication
}  // namespace pe
}  // namespace walberla
