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
//! \file NotificationType.h
//! \author Tobias Preclik
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//! \brief Header file for the notification types
//
//======================================================================================================================

#pragma once


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include "core/DataTypes.h"

namespace walberla {
namespace pe {
namespace communication {

//=================================================================================================
//
//  NOTIFICATION TYPES
//
//=================================================================================================

//*************************************************************************************************
//! Associate a unique number to notifications for identifying/tagging them.
enum NotificationType {
   rigidBodyDeletionNotification = 1,
   rigidBodyRemovalNotification,
   rigidBodyCopyNotification,
   rigidBodyForceNotification,
   rigidBodyUpdateNotification,
   rigidBodyMigrationNotification,
   rigidBodyRemoteMigrationNotification,
   rigidBodyVelocityUpdateNotification,
   rigidBodyVelocityCorrectionNotification,
   rigidBodyNewShadowCopyNotification,
   rigidBodyRemovalInformationNotification,
};
//*************************************************************************************************




//=================================================================================================
//
//  NOTIFICATION UTILITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\name Notification utility functions */
//@{

/*!\brief Returns the notification type of the template parameter.
 * \return The notification type of the template parameter.
 */
template<typename T>
NotificationType notificationType();
//@}
//*************************************************************************************************







//=================================================================================================
//
//  NOTIFICATION TYPE MARSHALLING FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\name Notification type marshalling functions */
//@{
template< typename Buffer > inline void marshal( Buffer& buffer, const NotificationType& type );
template< typename Buffer > inline void unmarshal( Buffer& buffer, NotificationType& type );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief TODO
 * \return void
 */
template< typename Buffer >
inline void marshal( Buffer& buffer, const NotificationType& type ) {
   buffer << static_cast<uint8_t>( type );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief TODO
 * \return void
 */
template< typename Buffer >
inline void unmarshal( Buffer& buffer, NotificationType& type ) {
   uint8_t tmp;
   buffer >> tmp;
   type = static_cast<NotificationType>( tmp );
}
//*************************************************************************************************


}  // namespace communication
}  // namespace pe
}  // namespace walberla
