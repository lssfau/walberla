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
//! \file RigidBodyRemovalInformationNotification.h
//! \author Tobias Preclik
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//! \brief Header file for the RigidBodyRemovalNotification class
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
/*!\brief Wrapper class for rigid body removal notifications.
 *
 * The RigidBodyRemovalInformationNotification class is used to signal other processes that a
 * shadow copy was destroyed.
 */
class RigidBodyRemovalInformationNotification {
public:
   struct Parameters {
      id_t sid_;
      Owner owner_;
   };

   inline explicit RigidBodyRemovalInformationNotification( const RigidBody& b ) : body_(b) {}
   const RigidBody& body_;
};
//*************************************************************************************************



//*************************************************************************************************
/*!\brief Marshaling rigid body removal information notifications.
 *
 * \param buffer The buffer to be filled.
 * \param obj The object to be marshaled.
 * \return void
 *
 * The removal information notifications just consists of the system ID of the body to remove.
 */
template< typename Buffer >
inline void marshal( Buffer& buffer, const RigidBodyRemovalInformationNotification& obj ) {
   buffer << obj.body_.getSystemID();
   buffer << obj.body_.MPITrait.getOwner();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Unmarshaling rigid body removal information notifications.
 *
 * \param buffer The buffer from where to read.
 * \param objparam The object to be reconstructed.
 * \return void
 *
 * The removal information notifications just consists of the system ID of the body to remove.
 */
template< typename Buffer >
inline void unmarshal( Buffer& buffer, RigidBodyRemovalInformationNotification::Parameters& objparam ) {
   buffer >> objparam.sid_;
   buffer >> objparam.owner_;
}
//*************************************************************************************************

//*************************************************************************************************
/*!\brief Returns the notification type of a rigid body removal information.
 * \return The notification type of a rigid body removal information.
 */
template<>
inline NotificationType notificationType<RigidBodyRemovalInformationNotification>() {
   return rigidBodyRemovalInformationNotification;
}
//*************************************************************************************************


}  // namespace communication
}  // namespace pe
}  // namespace walberla
