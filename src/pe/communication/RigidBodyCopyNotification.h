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
//! \file RigidBodyCopyNotification.h
//! \author Tobias Preclik
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//! \brief Header file for the RigidBodyCopyNotification class
//
//======================================================================================================================

#pragma once


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include "NotificationType.h"
#include "DynamicMarshalling.h"
#include <pe/rigidbody/RigidBody.h>
#include "core/DataTypes.h"

#include <vector>


namespace walberla {
namespace pe {
namespace communication {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Wrapper class for rigid body copies.
 *
 * The RigidBodyCopyNotification class is a wrapper class for marshaling and unmarshaling rigid body
 * copies. When receiving a copy notification a new shadow copy of the encapsulated body should be
 * created. Copy notification must originate from the owner. The class template should be the
 * collision system configuration and only serves the purpose to prevent compiler errors for
 * specific collision system choices. The marshalling of the list of shadow copy holders can only
 * be performed if the MPIRigidBodyTrait class is part of the rigid body class hierarchy, which is
 * not the case for all collision systems. However, if it is not part of the class hierarchy then
 * the function is also never called and the compiler need not compile the code and the template
 * parameter prevents it.
 */
class RigidBodyCopyNotification {
public:
   struct Parameters {
      id_t geomType_;
   };

   inline explicit RigidBodyCopyNotification( const RigidBody& b ) : body_(b) {}
   const RigidBody& body_;
};
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Marshaling rigid body copies.
 *
 * \param buffer The buffer to be filled.
 * \param obj The object to be marshaled.
 * \return void
 *
 * TODO
 */
template< typename Buffer >
inline void marshal( Buffer& buffer, const RigidBodyCopyNotification& obj ) {
   buffer << obj.body_.getTypeID();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Unmarshaling rigid body copies.
 *
 * \param buffer The buffer from where to read.
 * \param objparam The object to be reconstructed.
 * \return void
 *
 * TODO
 */
template< typename Buffer >
inline void unmarshal( Buffer& buffer, RigidBodyCopyNotification::Parameters& objparam ) {
   buffer >> objparam.geomType_;
}
//*************************************************************************************************

//*************************************************************************************************
/*!\brief Returns the notification type of a rigid body copy.
 * \return The notification type of a rigid body copy.
 */
template<>
inline NotificationType notificationType< RigidBodyCopyNotification >() {
   return rigidBodyCopyNotification;
}
//*************************************************************************************************

}  // namespace communication
}  // namespace pe
}  // namespace walberla
