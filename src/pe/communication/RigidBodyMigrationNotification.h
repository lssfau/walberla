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
//! \file RigidBodyMigrationNotification.h
//! \author Tobias Preclik
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//! \brief Header file for the RigidBodyMigrationNotification class
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
/*!\brief Wrapper class for rigid body migration notifications.
 *
 * The RigidBodyMigrationNotification class is a wrapper class for marshaling and unmarshaling
 * rigid body migration notifications. When receiving a migration notification this indicates that
 * a body migrates from the sending neighbor to the local process and the local process takes over
 * ownership of the body. Migration notices may only be sent if the new owner already obtained a
 * shadow copy previously. They may also only be sent by a neighboring process.
 */
class RigidBodyMigrationNotification {
public:
   struct Parameters {
      id_t sid_;
      std::vector<Owner> reglist_;
   };

   inline explicit RigidBodyMigrationNotification( const RigidBody& b ) : body_(b) {}
   const RigidBody& body_;
};
//*************************************************************************************************



//*************************************************************************************************
/*!\brief Marshaling rigid body migration notifications.
 *
 * \param buffer The buffer to be filled.
 * \param obj The object to be marshaled.
 * \return void
 *
 * The migration notifications just consists of the system ID of the body migrating.
 */
template< typename Buffer >
inline void marshal( Buffer& buffer, const RigidBodyMigrationNotification& obj ) {
   buffer << obj.body_.getSystemID();

   buffer << obj.body_.MPITrait.sizeShadowOwners();
   for( auto it = obj.body_.MPITrait.beginShadowOwners(); it != obj.body_.MPITrait.endShadowOwners(); ++it )
      buffer << (*it);
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Unmarshaling rigid body migration notifications.
 *
 * \param buffer The buffer from where to read.
 * \param objparam The object to be reconstructed.
 * \return void
 *
 * The migration notifications just consists of the system ID of the body migrating.
 */
template< typename Buffer >
inline void unmarshal( Buffer& buffer, RigidBodyMigrationNotification::Parameters& objparam ) {
   buffer >> objparam.sid_;

   size_t n;
   buffer >> n;
   objparam.reglist_.reserve( n );
   for( size_t i = 0; i < n; ++i ) {
      Owner p;
      buffer >> p;
      objparam.reglist_.push_back( p );
   }
}
//*************************************************************************************************

//*************************************************************************************************
/*!\brief Returns the notification type of a rigid body migration.
 * \return The notification type of a rigid body migration.
 */
template<>
inline NotificationType notificationType<RigidBodyMigrationNotification>() {
   return rigidBodyMigrationNotification;
}
//*************************************************************************************************


}  // namespace communication
}  // namespace pe
}  // namespace walberla
