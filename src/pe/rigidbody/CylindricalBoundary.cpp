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
//! \file CylindricalBoundary.cpp
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include "CylindricalBoundary.h"

//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <iomanip>
#include <iostream>
#include <stdexcept>
#include "core/math/Shims.h"
#include <pe/Materials.h>
#include <core/math/Matrix3.h>
#include <core/math/Limits.h>
#include <core/debug/Debug.h>

namespace walberla {
namespace pe {

//=================================================================================================
//
//  CONSTRUCTOR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Constructor for the CylindricalBoundary class.
 *
 * \param sid Unique system-specific ID for the cylindrical boundary.
 * \param uid User-specific ID for the cylindrical boundary.
 * \param gpos Global geometric center of the cylindrical boundary.
 * \param radius The radius of the cylinder.
 * \param material The material of the cylindrical boundary.
 *
 * The cylindrical boundary is created lying along the x-axis.
 */
CylindricalBoundary::CylindricalBoundary( id_t sid, id_t uid, const Vec3& gpos, const real_t radius,
                                          MaterialID material )
   : GeomPrimitive( getStaticTypeID(), sid, uid, material )           // Initializing the base object
   , radius_(radius)                                                  // Radius of the cylinder
{
   //boundaries are always considered locally and have infinite mass
   setGlobal( true );
   setCommunicating( false );
   setMassAndInertiaToInfinity();
   setFinite( false );

   // Checking the radius
   // Since the constructor is never directly called but only used in a small number
   // of functions that already check the cylinder arguments, only asserts are used here to
   // double check the arguments.
   WALBERLA_ASSERT_GREATER( radius, real_t(0), "Invalid cylinder radius"  );

   // Initializing the instantiated cylinder
   setPosition(gpos);
   setOrientation(Quat());

   // Setting the axis-aligned bounding box
   CylindricalBoundary::calcBoundingBox();
}
//*************************************************************************************************




//*************************************************************************************************
/*!\brief Destructor for the cylindrical boundary class.
 */
CylindricalBoundary::~CylindricalBoundary()
{
   WALBERLA_LOG_DETAIL( "Destroyed CylindricalBoundary " << sid_ );
}
//*************************************************************************************************




//=================================================================================================
//
//  UTILITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Checks, whether a point in body relative coordinates lies inside the cylindrical boundary.
 *
 * \param px The x-component of the relative coordinate.
 * \param py The y-component of the relative coordinate.
 * \param pz The z-component of the relative coordinate.
 * \return \a true if the point lies inside the cylindrical boundary, \a false if not.
 */
bool CylindricalBoundary::containsRelPointImpl( real_t /*px*/, real_t py, real_t pz ) const
{
   return ( math::sq(py) + math::sq(pz) ) >= ( radius_ * radius_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checks, whether a point in body relative coordinates lies on the surface of the cylindrical boundary.
 *
 * \param px The x-component of the relative coordinate.
 * \param py The y-component of the relative coordinate.
 * \param pz The z-component of the relative coordinate.
 * \return \a true if the point lies on the surface of the cylindrical boundary, \a false if not.
 *
 * The tolerance level of the check is pe::surfaceThreshold.
 */
bool CylindricalBoundary::isSurfaceRelPointImpl( real_t /*px*/, real_t py, real_t pz ) const
{
   return ( std::fabs( ( math::sq(py) + math::sq(pz) ) - ( radius_ * radius_ ) ) <= surfaceThreshold*surfaceThreshold );
}
//*************************************************************************************************

//*************************************************************************************************
/*!\brief Calculation of the bounding box of the cylindrical boundary.
 *
 * \return (-inf, -inf, -inf, +inf, +inf, +inf )
 */
void CylindricalBoundary::calcBoundingBox()
{
   aabb_ = math::AABB(
         -math::Limits<real_t>::inf(),
         -math::Limits<real_t>::inf(),
         -math::Limits<real_t>::inf(),
         +math::Limits<real_t>::inf(),
         +math::Limits<real_t>::inf(),
         +math::Limits<real_t>::inf());

   WALBERLA_ASSERT( aabb_.contains( getPosition() ), "Invalid bounding box detected" );
}
//*************************************************************************************************


//=================================================================================================
//
//  OUTPUT FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Output of the current state of a cylindrical boundary.
 *
 * \param os Reference to the output stream.
 * \param tab Indentation in front of every line of the cylindrical boundary output.
 * \return void
 */
void CylindricalBoundary::print( std::ostream& os, const char* tab ) const
{
   using std::setw;

   os << tab << " Cylindrical Boundary " << getID() << " with radius " << getRadius() << "\n";

   os << tab << "   Fixed: " << isFixed() << " , sleeping: " << !isAwake() << "\n";

   os << tab << "   System ID         = " << getSystemID() << "\n"
      << tab << "   Total mass        = " << getMass() << "\n"
      << tab << "   Material          = " << Material::getName( material_ ) << "\n"
      << tab << "   Global position   = " << getPosition() << "\n"
      << tab << "   Linear velocity   = " << getLinearVel() << "\n"
      << tab << "   Angular velocity  = " << getAngularVel() << "\n";
   Mat3 R = getRotation();
   os << tab << "   Bounding box      = " << getAABB() << "\n"
      << tab << "   Quaternion        = " << getQuaternion() << "\n"
      << tab << "   Rotation matrix   = ( " << setw(9) << R[0] << " , " << setw(9) << R[1] << " , " << setw(9) << R[2] << " )\n"
      << tab << "                       ( " << setw(9) << R[3] << " , " << setw(9) << R[4] << " , " << setw(9) << R[5] << " )\n"
      << tab << "                       ( " << setw(9) << R[6] << " , " << setw(9) << R[7] << " , " << setw(9) << R[8] << " )\n";
}
//*************************************************************************************************

//=================================================================================================
//
//  GLOBAL OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Global output operator for cylindrical boundaryies.
 *
 * \param os Reference to the output stream.
 * \param c Reference to a cylindrical boundary object.
 * \return Reference to the output stream.
 */
std::ostream& operator<<( std::ostream& os, const CylindricalBoundary& c )
{
   os << "--" << "CYLINDRICAL BOUNDARY PARAMETERS"
      << "------------------------------------------------------------\n";
   c.print( os, "" );
   os << "--------------------------------------------------------------------------------\n"
      << std::endl;
   return os;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Global output operator for cylindrical boundaryies handles.
 *
 * \param os Reference to the output stream.
 * \param c Constant cylindrical boundary handle.
 * \return Reference to the output stream.
 */
std::ostream& operator<<( std::ostream& os, CylindricalBoundaryID c )
{
   os << "--" << "CYLINDRICAL BOUNDARY PARAMETERS"
      << "------------------------------------------------------------\n";
   c->print( os, "" );
   os << "--------------------------------------------------------------------------------\n"
      << std::endl;
   return os;
}
//*************************************************************************************************

id_t CylindricalBoundary::staticTypeID_ = std::numeric_limits<id_t>::max();

} // namespace pe
} // namespace walberla


