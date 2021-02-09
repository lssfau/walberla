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
//! \file CylindricalBoundary.h
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#pragma once


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <pe/rigidbody/GeomPrimitive.h>
#include <pe/Types.h>
#include <core/math/Constants.h>
#include <core/math/Matrix3.h>
#include <core/math/Vector3.h>
#include <core/DataTypes.h>
#include <core/logging/Logging.h>
#include <pe/Config.h>


namespace walberla {
namespace pe {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/**
 * \ingroup pe
 */
class CylindricalBoundary : public GeomPrimitive
{
public:
   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit CylindricalBoundary( id_t sid, id_t uid, const Vec3& gpos, const real_t radius,
                                 MaterialID material );
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   /*!\name Destructor */
   //@{
   ~CylindricalBoundary() override;
   //@}
   //**********************************************************************************************

public:
   //**Get functions*******************************************************************************
   /*!\name Get functions */
   //@{
   inline real_t  getRadius() const;
   inline Vec3    getAxis() const;
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   static inline id_t getStaticTypeID();
   //@}
   //**********************************************************************************************

   //**Output functions****************************************************************************
   /*!\name Output functions */
   //@{
   void print( std::ostream& os, const char* tab ) const override;
   //@}
   //**********************************************************************************************

protected:
   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   bool containsRelPointImpl ( real_t px, real_t py, real_t pz ) const override;
   bool isSurfaceRelPointImpl( real_t px, real_t py, real_t pz ) const override;

   void calcBoundingBox() override;  // Calculation of the axis-aligned bounding box
   //@}
   //**********************************************************************************************

   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   real_t  radius_;  //!< The radius of the cylinder.
   //@}
   //**********************************************************************************************

private:
   static id_t staticTypeID_;  //< type id of sphere, will be set by SetBodyTypeIDs
   static void setStaticTypeID(id_t typeID) {staticTypeID_ = typeID;}

   //** friend declaration
   /// needed to be able to set static type ids with setStaticTypeID
   template <class T, int N>
   friend struct SetBodyTypeIDs;
};
//*************************************************************************************************



//=================================================================================================
//
//  GET FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Returns the radius of the cylinder.
 *
 * \return The radius of the cylinder.
 */
inline real_t  CylindricalBoundary::getRadius() const
{
   return radius_;
}
//*************************************************************************************************

//*************************************************************************************************
/*!\brief Returns the main axis of the cylinder.
 *
 * \return The main axis of the cylinder.
 */
inline Vec3 CylindricalBoundary::getAxis() const
{
   return getRotation() * Vec3(1,0,0);
}
//*************************************************************************************************


//=================================================================================================
//
//  UTILITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Returns unique type id of this type
 *
 * \return geometry specific type id
 */
inline id_t CylindricalBoundary::getStaticTypeID()
{
   return staticTypeID_;
}




//=================================================================================================
//
//  GLOBAL OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\name Cylinder operators */
//@{
std::ostream& operator<<( std::ostream& os, const CylindricalBoundary& cb );
std::ostream& operator<<( std::ostream& os, CylindricalBoundaryID cb );
//@}
//*************************************************************************************************

} // namespace pe
} // namespace walberla
