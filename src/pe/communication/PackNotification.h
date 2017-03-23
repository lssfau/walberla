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
//! \file PackNotification.h
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//! \brief Header file for packNotification
//
//======================================================================================================================

#pragma once


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include "NotificationType.h"

#include "pe/rigidbody/Owner.h"

#include "core/mpi/SendBuffer.h"


namespace walberla {
namespace pe {
namespace communication {

template <class T>
inline void packNotification(const Owner& sender, const Owner& receiver, mpi::SendBuffer& sb, const T& notification)
{
   sb << sender.blockID_;
   sb << receiver.blockID_;
   marshal( sb, notificationType< T >() );
   marshal( sb, notification );
}

template <class T>
inline void packNotificationWithoutSender(const Owner& receiver, mpi::SendBuffer& sb, const T& notification)
{
   sb << receiver.blockID_;
   marshal( sb, notificationType< T >() );
   marshal( sb, notification );
}

}  // namespace communication
}  // namespace pe
}  // namespace walberla
