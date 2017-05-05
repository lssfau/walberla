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
//! \file RigidBodyVelocityCorrectionNotification.h
//! \author Tobias Preclik
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//! \brief Header file for the RigidBodyVelocityCorrectionNotification class
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
/*!\brief Wrapper class for rigid body velocity corrections.
 *
 * The RigidBodyVelocityCorrectionNotification class is a wrapper class for marshalling and unmarshalling rigid body
 * velocity correction. It includes the system ID of the body the correction applies to and the linear and angular
 * velocity corrections.
 */
class RigidBodyVelocityCorrectionNotification {
public:
   struct Parameters {
      id_t sid_;
      Vec3 dv_, dw_;
   };

   inline RigidBodyVelocityCorrectionNotification( const RigidBody& b, const Vec3& dv, const Vec3& dw ) : body_(b), dv_(dv), dw_(dw) {}
   const RigidBody& body_;
   const Vec3& dv_, dw_;
};
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Marshalling rigid body updates.
 *
 * \param buffer The buffer to be filled.
 * \param obj The object to be marshalled.
 * \return void
 *
 * The update consists of the linear and angular velocity corrections.
 */
template< typename Buffer >
inline void marshal( Buffer& buffer, const RigidBodyVelocityCorrectionNotification& obj ) {
   buffer << obj.body_.getSystemID();
   buffer << obj.dv_;
   buffer << obj.dw_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Unmarshalling rigid body updates.
 *
 * \param buffer The buffer from where to read.
 * \param objparam The object to be reconstructed.
 * \return void
 *
 * The update consists of the linear and angular velocity corrections.
 */
template< typename Buffer >
inline void unmarshal( Buffer& buffer, RigidBodyVelocityCorrectionNotification::Parameters& objparam ) {
   buffer >> objparam.sid_;
   buffer >> objparam.dv_;
   buffer >> objparam.dw_;
}
//*************************************************************************************************

//*************************************************************************************************
/*!\brief Returns the notification type of a rigid body migration.
 * \return The notification type of a rigid body velocity correction.
 */
template<>
inline NotificationType notificationType<RigidBodyVelocityCorrectionNotification>() {
   return rigidBodyVelocityCorrectionNotification;
}
//*************************************************************************************************

}  // namespace communication
}  // namespace pe
}  // namespace walberla
