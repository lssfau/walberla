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
//! \file
//! \author Tobias Preclik
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//! \brief Header file for the notification types
//
//======================================================================================================================

#pragma once

#include <core/DataTypes.h>

namespace walberla {
namespace mesa_pd {

//=================================================================================================
//
//  NOTIFICATION TYPES
//
//=================================================================================================

//*************************************************************************************************
//! Associate a unique number to notifications for identifying/tagging them.
enum NotificationType : uint8_t
{
   PARTICLE_DELETION_NOTIFICATION = 1,
   PARTICLE_REMOVAL_NOTIFICATION,
   PARTICLE_COPY_NOTIFICATION,
   PARTICLE_GHOST_COPY_NOTIFICATION,
   PARTICLE_FORCE_NOTIFICATION,
   PARTICLE_UPDATE_NOTIFICATION,
   PARTICLE_MIGRATION_NOTIFICATION,
   PARTICLE_REMOTE_MIGRATION_NOTIFICATION,
   PARTICLE_VELOCITY_UPDATE_NOTIFICATION,
   PARTICLE_VELOCITY_CORRECTION_NOTIFICATION,
   NEW_GHOST_PARTICLE_NOTIFICATION,
   PARTICLE_REMOVAL_INFORMATION_NOTIFICATION
};
//*************************************************************************************************




//=================================================================================================
//
//  NOTIFICATION UTILITY FUNCTIONS
//
//=================================================================================================

/**
 * Returns the notification type of the template parameter.
 *
 * To be specialized by the notifications.
 * \return The notification type of the template parameter.
 */
template< typename V >
struct NotificationTrait
{};

}  // namespace mesa_pd
}  // namespace walberla
