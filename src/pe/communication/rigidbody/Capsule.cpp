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
//! \file Capsule.cpp
//! \author Tobias Preclik
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//! \brief Marshalling of objects for data transmission or storage.
//
//======================================================================================================================

#include "Capsule.h"

namespace walberla {
namespace pe {
namespace communication {

//*************************************************************************************************
void marshal( mpi::SendBuffer& buffer, const Capsule& obj )
{
   marshal( buffer, static_cast<const GeomPrimitive&>( obj ) );
   buffer << obj.getRadius();
   buffer << obj.getLength();
}
//*************************************************************************************************


//*************************************************************************************************
void unmarshal( mpi::RecvBuffer& buffer, CapsuleParameters& objparam )
{
   unmarshal( buffer, static_cast<GeomPrimitiveParameters&>( objparam ) );
   buffer >> objparam.radius_;
   buffer >> objparam.length_;
}
//*************************************************************************************************

}  // namespace communication
}  // namespace pe
}  // namespace walberla
