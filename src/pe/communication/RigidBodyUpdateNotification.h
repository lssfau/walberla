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
//! \file RigidBodyUpdateNotification.h
//! \author Tobias Preclik
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//! \brief Header file for the RigidBodyUpdateNotification class
//
//======================================================================================================================

#pragma once

//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <pe/rigidbody/RigidBody.h>
#include "NotificationType.h"
#include "Marshalling.h"


namespace walberla {
namespace pe {
namespace communication {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Wrapper class for rigid body updates.
 *
 * The RigidBodyUpdateNotification class is a wrapper class for marshalling and unmarshalling rigid
 * body updates. The class template should be the
 * collision system configuration and only serves the purpose to prevent compiler errors for
 * specific collision system choices. The marshalling of the list of shadow copy holders can only
 * be performed if the MPIRigidBodyTrait class is part of the rigid body class hierarchy, which is
 * not the case for all collision systems. However, if it is not part of the class hierarchy then
 * the function is also never called and the compiler need not compile the code and the template
 * parameter prevents it.
 */
class RigidBodyUpdateNotification {
public:
   struct Parameters {
      id_t sid_;
      Vec3 gpos_, v_, w_;
      Quat q_;
   };

   inline explicit RigidBodyUpdateNotification( const RigidBody& b ) : body_(b) {}
   const RigidBody& body_;
};
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Marshalling rigid body updates.
 *
 * \param buffer The buffer to be filled.
 * \param obj The object to be marshalled.
 * \return void
 *
 * The update consists of the position, orientation and linear and angular velocities and the
 * list of shadow copy holders.
 */
template< typename Buffer >
inline void marshal( Buffer& buffer, const RigidBodyUpdateNotification& obj ) {

   buffer << obj.body_.getSystemID();
   buffer << obj.body_.getPosition();
   buffer << obj.body_.getQuaternion();
   buffer << obj.body_.getLinearVel();
   buffer << obj.body_.getAngularVel();

}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Unmarshalling rigid body updates.
 *
 * \param buffer The buffer from where to read.
 * \param objparam The object to be reconstructed.
 * \return void
 *
 * The update consists of the position, orientation and linear and angular velocities and the
 * list of shadow copy holders.
 */
template< typename Buffer >
inline void unmarshal( Buffer& buffer, typename RigidBodyUpdateNotification::Parameters& objparam ) {
   buffer >> objparam.sid_;
   buffer >> objparam.gpos_;
   buffer >> objparam.q_;
   buffer >> objparam.v_;
   buffer >> objparam.w_;

}
//*************************************************************************************************

//*************************************************************************************************
/*!\brief Returns the notification type of a rigid body update.
 * \return The notification type of a rigid body update.
 */
template<>
inline NotificationType notificationType< RigidBodyUpdateNotification >() {
   return rigidBodyUpdateNotification;
}
//*************************************************************************************************

}  // namespace communication
}  // namespace pe
}  // namespace walberla
