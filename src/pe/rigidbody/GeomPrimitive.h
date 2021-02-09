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
//! \file GeomPrimitive.h
//! \author Klaus Iglberger
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//! \brief Header file for the GeomPrimitive class
//
//======================================================================================================================

#pragma once


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <pe/rigidbody/RigidBody.h>
#include <pe/Types.h>
#include <core/DataTypes.h>

namespace walberla {
namespace pe {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Base class for all primitive geometries.
 *
 * The GeomPrimitive class is the abstract basis for all geometric primitives of the rigid body
 * physics module. The class provides the common data members and the common functionality for
 * all geometry classes.
 */
class GeomPrimitive : public RigidBody
{
public:
   //**Get functions*******************************************************************************
   /*!\name Get functions */
   //@{
   inline MaterialID getMaterial() const;
   //@}
   //**********************************************************************************************

protected:
   //**Constructor*********************************************************************************
   /*!\name Constructor */
   //@{
   explicit GeomPrimitive( id_t const typeID, id_t sid, id_t uid, MaterialID material );
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   /*!\name Destructor */
   //@{
   ~GeomPrimitive() override = 0;
   //@}
   //**********************************************************************************************

   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   MaterialID material_;  //!< The material of the geometric primitive.
   //@}
   //**********************************************************************************************
};
//*************************************************************************************************




//=================================================================================================
//
//  GET FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Returns the material of the geometric primitive.
 *
 * \return The coefficient of restitution of the rigid body.
 */
inline MaterialID GeomPrimitive::getMaterial() const
{
   return material_;
}
//*************************************************************************************************

} // namespace pe
}  // namespace walberla
