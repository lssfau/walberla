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
//! \file Sphere.cpp
//! \author Klaus Iglberger
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include "Sphere.h"

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
/*!\brief Constructor for the Sphere class.
 *
 * \param sid Unique system-specific ID for the sphere.
 * \param uid User-specific ID for the sphere.
 * \param gpos Global geometric center of the sphere.
 * \param q The orientation of the sphere's body frame in the global world frame.
 * \param radius The radius of the sphere \f$ (0..\infty) \f$.
 * \param material The material of the sphere.
 * \param global specifies if the sphere should be created in the global storage
 * \param communicating specifies if the sphere should take part in synchronization (syncNextNeighbour, syncShadowOwner)
 * \param infiniteMass specifies if the sphere has infinite mass and will be treated as an obstacle
 */
Sphere::Sphere( id_t sid, id_t uid, const Vec3& gpos, const Quat& q,
                real_t radius, MaterialID material,
                const bool global, const bool communicating, const bool infiniteMass )
   : Sphere::Sphere( getStaticTypeID(), sid, uid, gpos, q, radius, material, global, communicating, infiniteMass )
{}

Sphere::Sphere( id_t const typeId, id_t sid, id_t uid, const Vec3& gpos, const Quat& q,
                real_t radius, MaterialID material,
                const bool global, const bool communicating, const bool infiniteMass )
   : GeomPrimitive( typeId, sid, uid, material )  // Initialization of the parent class
   , radius_(radius)
{
   // Checking the radius
   // Since the sphere constructor is never directly called but only used in a small number
   // of functions that already check the sphere arguments, only asserts are used here to
   // double check the arguments.
   WALBERLA_ASSERT( radius > real_c(0), "Invalid sphere radius" );

   // Setting the center of the sphere
   setPosition(gpos);
   setOrientation(q);

   setGlobal( global );
   if (infiniteMass)
   {
      setMassAndInertia( std::numeric_limits<real_t>::infinity(), Mat3(real_t(0)) );
   } else
   {
      auto mass = calcMass( radius, Material::getDensity( material ) );
      setMassAndInertia( mass, calcInertia( mass, radius ) );
   }
   setCommunicating( communicating );
   setFinite( true );

   // Setting the axis-aligned bounding box
   Sphere::calcBoundingBox();
}
//*************************************************************************************************




//=================================================================================================
//
//  DESTRUCTOR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Destructor for the Sphere class.
 */
Sphere::~Sphere()
{
   // Logging the destruction of the sphere
   WALBERLA_LOG_DETAIL( "Destroyed sphere " << sid_ );
}
//*************************************************************************************************




//=================================================================================================
//
//  UTILITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Checks, whether a point in body relative coordinates lies inside the sphere.
 *
 * \param px The x-component of the relative coordinate.
 * \param py The y-component of the relative coordinate.
 * \param pz The z-component of the relative coordinate.
 * \return \a true if the point lies inside the sphere, \a false if not.
 */
bool Sphere::containsRelPointImpl( real_t px, real_t py, real_t pz ) const
{
   const Vec3 rpos( px, py, pz );
   return ( rpos.sqrLength() <= ( radius_*radius_ ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checks, whether a point in body relative coordinates lies on the surface of the sphere.
 *
 * \param px The x-component of the relative coordinate.
 * \param py The y-component of the relative coordinate.
 * \param pz The z-component of the relative coordinate.
 * \return \a true if the point lies on the surface of the sphere, \a false if not.
 *
 * The tolerance level of the check is pe::surfaceThreshold.
 */
bool Sphere::isSurfaceRelPointImpl( real_t px, real_t py, real_t pz ) const
{
   const Vec3 rpos( px, py, pz );
   return floatIsEqual( rpos.sqrLength(), radius_*radius_, pe::surfaceThreshold );
}
//*************************************************************************************************




//=================================================================================================
//
//  OUTPUT FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Output of the current state of a sphere.
 *
 * \param os Reference to the output stream.
 * \param tab Indentation in front of every line of the sphere output.
 * \return void
 */
void Sphere::print( std::ostream& os, const char* tab ) const
{
   using std::setw;

   os << tab << " Sphere " << uid_ << " with radius " << radius_ << "\n";

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
   Mat3 R = getRotation();
   os << tab << "   Bounding box      = " << getAABB() << "\n"
      << tab << "   Quaternion        = " << getQuaternion() << "\n"
      << tab << "   Rotation matrix   = ( " << setw(9) << R[0] << " , " << setw(9) << R[1] << " , " << setw(9) << R[2] << " )\n"
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
/*!\brief Global output operator for spheres.
 *
 * \param os Reference to the output stream.
 * \param s Reference to a constant sphere object.
 * \return Reference to the output stream.
 */
std::ostream& operator<<( std::ostream& os, const Sphere& s )
{
   os << "--" << "SPHERE PARAMETERS"
      << "-------------------------------------------------------------\n";
   s.print( os, "" );
   os << "--------------------------------------------------------------------------------\n"
      << std::endl;
   return os;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Global output operator for sphere handles.
 *
 * \param os Reference to the output stream.
 * \param s Constant sphere handle.
 * \return Reference to the output stream.
 */
std::ostream& operator<<( std::ostream& os, ConstSphereID s )
{
   os << "--" << "SPHERE PARAMETERS"
      << "-------------------------------------------------------------\n";
   s->print( os, "" );
   os << "--------------------------------------------------------------------------------\n"
      << std::endl;
   return os;
}
//*************************************************************************************************

id_t Sphere::staticTypeID_ = std::numeric_limits<id_t>::max();

} // namespace pe
}  // namespace walberla
