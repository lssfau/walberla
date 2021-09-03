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
//! \file GeomPrimitive.cpp
//! \author Klaus Iglberger
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//! \brief Source file for the GeomPrimitive class
//
//======================================================================================================================


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <pe/rigidbody/GeomPrimitive.h>


namespace walberla {
namespace pe {

//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Constructor for the GeomPrimitive class.
 *
 * \param typeID The geometry type of the geometric primitive.
 * \param sid The unique system-specific ID of the geometric primitive.
 * \param uid The user-specific ID of the geometric primitive.
 * \param material The material of the geometric primitive.
 */
GeomPrimitive::GeomPrimitive( id_t const typeID, id_t sid, id_t uid, MaterialID material )
   : RigidBody(typeID, sid, uid)  // Initializing the base object
   , material_(material)                     // The material of the geometric primitive
{}
//*************************************************************************************************




//=================================================================================================
//
//  DESTRUCTOR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Destructor for the Primitive class.
 */
GeomPrimitive::~GeomPrimitive()
= default;
//*************************************************************************************************

}  // namespace pe
}  // namespace walberla
