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
//! \file Ellipsoid.cpp
//! \author Klaus Iglberger
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include "Ellipsoid.h"

//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <cmath>
#include <pe/Materials.h>
#include <core/math/Matrix3.h>
#include <core/debug/Debug.h>

namespace walberla {
namespace pe {

//=================================================================================================
//
//  CONSTRUCTOR
//
//=================================================================================================


//*************************************************************************************************
//*************************************************************************************************
/*!\brief Constructor for the Ellipsoid class.
 *
 * \param sid Unique system-specific ID for the Ellipsoid.
 * \param uid User-specific ID for the Ellipsoid.
 * \param gpos Global geometric center of the Ellipsoid.
 * \param rpos The relative position within the body frame of a superordinate body.
 * \param q The orientation of the Ellipsoid's body frame in the global world frame.
 * \param radius The radius of the Ellipsoid \f$ (0..\infty) \f$.
 * \param material The material of the Ellipsoid.
 * \param global specifies if the Ellipsoid should be created in the global storage
 * \param communicating specifies if the Ellipsoid should take part in synchronization (syncNextNeighbour, syncShadowOwner)
 * \param infiniteMass specifies if the Ellipsoid has infinite mass and will be treated as an obstacle
 */
Ellipsoid::Ellipsoid( id_t sid, id_t uid, const Vec3& gpos, const Quat& q,
                const Vec3& semiAxes, MaterialID material,
                const bool global, const bool communicating, const bool infiniteMass )
   : Ellipsoid::Ellipsoid( getStaticTypeID(), sid, uid, gpos, q, semiAxes, material, global, communicating, infiniteMass )
{}
Ellipsoid::Ellipsoid( id_t const typeId, id_t sid, id_t uid, const Vec3& gpos, const Quat& q,
                const Vec3& semiAxes, MaterialID material,
                const bool global, const bool communicating, const bool infiniteMass )
   : GeomPrimitive( typeId, sid, uid, material )  // Initialization of the parent class
   , semiAxes_(semiAxes)
{
   // Checking the radius
   // Since the Ellipsoid constructor is never directly called but only used in a small number
   // of functions that already check the Ellipsoid arguments, only asserts are used here to
   // double check the arguments.
   WALBERLA_ASSERT( semiAxes_[0] > real_c(0), "Invalid Ellipsoid radius" );
   WALBERLA_ASSERT( semiAxes_[1] > real_c(0), "Invalid Ellipsoid radius" );
   WALBERLA_ASSERT( semiAxes_[2] > real_c(0), "Invalid Ellipsoid radius" );

   // Setting the center of the Ellipsoid
   setPosition(gpos);
   setOrientation(q);

   setGlobal( global );
   if (infiniteMass)
   {
      setMassAndInertiaToInfinity();
   } else
   {
      auto mass = calcMass( semiAxes, Material::getDensity( material ) );
      setMassAndInertia( mass, calcInertia( mass, semiAxes ) );
   }
   setCommunicating( communicating );
   setFinite( true );

   // Setting the axis-aligned bounding box
   Ellipsoid::calcBoundingBox();
}
//*************************************************************************************************




//=================================================================================================
//
//  DESTRUCTOR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Destructor for the Ellipsoid class.
 */
Ellipsoid::~Ellipsoid()
{
   // Logging the destruction of the Ellipsoid
   WALBERLA_LOG_DETAIL( "Destroyed Ellipsoid " << sid_ );
}
//*************************************************************************************************




//=================================================================================================
//
//  UTILITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Checks, whether a point in body relative coordinates lies inside the Ellipsoid.
 *
 * \param px The x-component of the relative coordinate.
 * \param py The y-component of the relative coordinate.
 * \param pz The z-component of the relative coordinate.
 * \return \a true if the point lies inside the Ellipsoid, \a false if not.
 */
bool Ellipsoid::containsRelPointImpl( real_t px, real_t py, real_t pz ) const
{
return ( (px * px)/(semiAxes_[0] * semiAxes_[0]) + (py * py)/(semiAxes_[1] * semiAxes_[1]) 
		+ (pz * pz)/(semiAxes_[2] * semiAxes_[2]) <= real_t(1.0) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checks, whether a point in body relative coordinates lies on the surface of the Ellipsoid.
 *
 * \param px The x-component of the relative coordinate.
 * \param py The y-component of the relative coordinate.
 * \param pz The z-component of the relative coordinate.
 * \return \a true if the point lies on the surface of the Ellipsoid, \a false if not.
 *
 * The (relative) tolerance level of the check is pe::surfaceThreshold.
 */
bool Ellipsoid::isSurfaceRelPointImpl( real_t px, real_t py, real_t pz ) const
{
   return floatIsEqual( (px * px)/(semiAxes_[0] * semiAxes_[0]) + (py * py)/(semiAxes_[1] * semiAxes_[1]) 
		+ (pz * pz)/(semiAxes_[2] * semiAxes_[2]), real_t(1.0), pe::surfaceThreshold);
}
//*************************************************************************************************




//=================================================================================================
//
//  OUTPUT FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Output of the current state of a Ellipsoid.
 *
 * \param os Reference to the output stream.
 * \param tab Indentation in front of every line of the Ellipsoid output.
 * \return void
 */
void Ellipsoid::print( std::ostream& os, const char* tab ) const
{
   using std::setw;

   Mat3 R = getRotation();

   os << tab << " Ellipsoid " << uid_ << " with semi-axis " << semiAxes_ << "\n";

   os << tab << "   Fixed: " << isFixed() << " , sleeping: " << !isAwake() << "\n";

   os << tab << "   System ID         = " << getSystemID() << "\n"
      << tab << "   Total mass        = " << getMass() << "\n"
      << tab << "   Material          = " << Material::getName( material_ ) << "\n"
      << tab << "   Owner             = " << MPITrait.getOwner() << "\n"
      << tab << "   Global position   = " << getPosition() << "\n"
      << tab << "   Relative position = " << getRelPosition() << "\n"
      << tab << "   Linear velocity   = " << getLinearVel() << "\n"
      << tab << "   Angular velocity  = " << getAngularVel() << "\n";

//   if( verboseMode )
//   {
      os << tab << "   Bounding box      = " << getAABB() << "\n"
         << tab << "   Quaternion        = " << getQuaternion() << "\n"
         << tab << "   Rotation matrix   = ( " << setw(9) << R[1] << " , " << setw(9) << R[1] << " , " << setw(9) << R[2] << " )\n"
         << tab << "                       ( " << setw(9) << R[3] << " , " << setw(9) << R[4] << " , " << setw(9) << R[5] << " )\n"
         << tab << "                       ( " << setw(9) << R[6] << " , " << setw(9) << R[7] << " , " << setw(9) << R[8] << " )\n";

      os << std::setiosflags(std::ios::right)
         << tab << "   Moment of inertia = ( " << setw(9) << I_[0] << " , " << setw(9) << I_[1] << " , " << setw(9) << I_[2] << " )\n"
         << tab << "                       ( " << setw(9) << I_[3] << " , " << setw(9) << I_[4] << " , " << setw(9) << I_[5] << " )\n"
         << tab << "                       ( " << setw(9) << I_[6] << " , " << setw(9) << I_[7] << " , " << setw(9) << I_[8] << " )\n"
         << std::resetiosflags(std::ios::right);
//   }
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Global output operator for Ellipsoids.
 *
 * \param os Reference to the output stream.
 * \param s Reference to a constant Ellipsoid object.
 * \return Reference to the output stream.
 */
std::ostream& operator<<( std::ostream& os, const Ellipsoid& s )
{
   os << "--" << "Ellipsoid PARAMETERS"
      << "-------------------------------------------------------------\n";
   s.print( os, "" );
   os << "--------------------------------------------------------------------------------\n"
      << std::endl;
   return os;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Global output operator for Ellipsoid handles.
 *
 * \param os Reference to the output stream.
 * \param s Constant Ellipsoid handle.
 * \return Reference to the output stream.
 */
std::ostream& operator<<( std::ostream& os, ConstEllipsoidID s )
{
   os << "--" << "Ellipsoid PARAMETERS"
      << "-------------------------------------------------------------\n";
   s->print( os, "" );
   os << "--------------------------------------------------------------------------------\n"
      << std::endl;
   return os;
}
//*************************************************************************************************

id_t Ellipsoid::staticTypeID_ = std::numeric_limits<id_t>::max();

} // namespace pe
}  // namespace walberla
