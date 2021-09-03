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
//! \author Klaus Iglberger
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include "Capsule.h"

//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <iomanip>
#include <iostream>
#include <stdexcept>
#include "core/math/Shims.h"
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
/*!\brief Constructor for the CapsuleBase class.
 *
 * \param sid Unique system-specific ID for the capsule.
 * \param uid User-specific ID for the capsule.
 * \param gpos Global geometric center of the capsule.
 * \param q The orientation of the capsule's body frame in the global world frame.
 * \param radius The radius of the cylinder part and the end caps \f$ (0..\infty) \f$.
 * \param length The length of the cylinder part \f$ (0..\infty) \f$.
 * \param material The material of the capsule.
 * \param global specifies if the capsule should be created in the global storage
 * \param communicating specifies if the capsule should take part in synchronization (syncNextNeighbour, syncShadowOwner)
 * \param infiniteMass specifies if the capsule has infinite mass and will be treated as an obstacle
 *
 * The capsule is created lying along the x-axis.
 */
Capsule::Capsule( id_t sid, id_t uid, const Vec3& gpos, const Quat& q,
                  real_t  radius, real_t  length, MaterialID material,
                  const bool global, const bool communicating, const bool infiniteMass )
   : GeomPrimitive( getStaticTypeID(), sid, uid, material )           // Initializing the base object
   , radius_(radius)                                                  // Radius of the capsule
   , length_(length)                                                  // Length of the capsule
{
   // Checking the radius and the length
   // Since the capsule constructor is never directly called but only used in a small number
   // of functions that already check the capsule arguments, only asserts are used here to
   // double check the arguments.
   WALBERLA_ASSERT_GREATER( radius, real_t(0), "Invalid capsule radius"  );
   WALBERLA_ASSERT_GREATER( length, real_t(0), "Invalid capsule length"  );

   // Initializing the instantiated capsule
   setPosition(gpos);
   setOrientation(q);

   setGlobal( global );
   if (infiniteMass)
   {
      setMassAndInertiaToInfinity();
   } else
   {
      auto mass = calcMass( radius, length, Material::getDensity( material ) );
      setMassAndInertia( mass, calcInertia( radius, length, Material::getDensity( material )  ) );
   }
   setCommunicating( communicating );
   setFinite( true );

   // Setting the axis-aligned bounding box
   Capsule::calcBoundingBox();
}
//*************************************************************************************************




//*************************************************************************************************
/*!\brief Destructor for the CapsuleBase class.
 */
Capsule::~Capsule()
{
   WALBERLA_LOG_DETAIL( "Destroyed capsule " << sid_ );
}
//*************************************************************************************************




//=================================================================================================
//
//  UTILITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Checks, whether a point in body relative coordinates lies inside the capsule.
 *
 * \param px The x-component of the relative coordinate.
 * \param py The y-component of the relative coordinate.
 * \param pz The z-component of the relative coordinate.
 * \return \a true if the point lies inside the capsule, \a false if not.
 */
bool Capsule::containsRelPointImpl( real_t px, real_t py, real_t pz ) const
{
   const real_t xabs( std::fabs( px ) );         // Absolute x-distance
   const real_t hlength( real_t(0.5) * length_ );  // Capsule half length

   if( xabs > hlength ) {
      return ( ( math::sq(xabs-hlength) + math::sq(py) + math::sq(pz) ) <= ( radius_ * radius_ ) );
   }
   else {
      return ( math::sq(py) + math::sq(pz) ) <= ( radius_ * radius_ );
   }
}
//*************************************************************************************************

//*************************************************************************************************
/*!\brief Checks, whether a point in body relative coordinates lies on the surface of the capsule.
 *
 * \param px The x-component of the relative coordinate.
 * \param py The y-component of the relative coordinate.
 * \param pz The z-component of the relative coordinate.
 * \return \a true if the point lies on the surface of the capsule, \a false if not.
 *
 * The tolerance level of the check is pe::surfaceThreshold.
 */
bool Capsule::isSurfaceRelPointImpl( real_t px, real_t py, real_t  pz ) const
{
   const real_t  xabs( std::fabs( px ) );         // Absolute x-distance
   const real_t  hlength( real_t(0.5) * length_ );  // Capsule half length

   if( xabs > hlength ) {
      return ( std::fabs( math::sq(xabs-hlength) + math::sq(py) + math::sq(pz) - radius_*radius_ ) <= surfaceThreshold*surfaceThreshold );
   }
   else {
      return ( std::fabs( ( math::sq(py) + math::sq(pz) ) - ( radius_ * radius_ ) ) <= surfaceThreshold*surfaceThreshold );
   }
}
//*************************************************************************************************

//*************************************************************************************************
/*!\brief Calculation of the bounding box of the capsule.
 *
 * \return void
 *
 * This function updates the axis-aligned bounding box of the capsule primitive according to the
 * current position and orientation of the capsule. Note that the bounding box is increased in
 * all dimensions by pe::contactThreshold to guarantee that rigid bodies in close proximity of
 * the capsule are also considered during the collision detection process.
 */
void Capsule::calcBoundingBox()
{
   Mat3 R = getRotation();
   Vec3 gpos = getPosition();
   const real_t  xlength( std::fabs( R[0]*length_ )*real_t (0.5) + radius_ + contactThreshold );
   const real_t  ylength( std::fabs( R[3]*length_ )*real_t (0.5) + radius_ + contactThreshold );
   const real_t  zlength( std::fabs( R[6]*length_ )*real_t (0.5) + radius_ + contactThreshold );
   aabb_ = math::AABB(
           gpos[0] - xlength,
         gpos[1] - ylength,
         gpos[2] - zlength,
         gpos[0] + xlength,
         gpos[1] + ylength,
         gpos[2] + zlength);

   WALBERLA_ASSERT( aabb_.contains( getPosition() ), "Invalid bounding box detected" );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Calculation of the moment of inertia in reference to the body frame of the capsule.
 *
 * \return void
 */
Mat3 Capsule::calcInertia( const real_t radius, const real_t length, const real_t density)
{
   const real_t  sphereMass( real_t (4)/real_t (3) * math::pi * radius*radius*radius * density );
   const real_t  cylinderMass( math::pi * radius*radius * length * density );

   // 'Ia' represent the moment of inertia along the x-axis. 'Ia' contains the following two parts:
   //  - cylinder :  I = (1/2)*mass*radius^2
   //  - sphere   :  I = (2/5)*mass*radius^2
   const real_t  Ia( radius*radius * ( real_t (0.5)*cylinderMass + real_t (0.4)*sphereMass ) );

   // 'Ib' represent the moment of inertia along the y- and z-axis. 'Ib' contains the following two parts,
   // where full_length is the full length of the cylinder part and half_length is (1/2)*full_length:
   //  - cylinder :  I = mass*( (1/4)*radius^2 + (1/12)*full_length^2 )
   //  - sphere   :  I = mass*( (2/5)*radius^2 + half_length^2 + (3/4)*half_length*radius )
   const real_t  Ib( cylinderMass*( real_t (0.25)*radius*radius + real_t (1)/real_t (12)*length*length ) +
                     sphereMass*( real_t (0.4)*radius*radius + real_t (0.375)*radius*length + real_t (0.25)*length*length ) );

   // Setting the moment of inertia (capsule is aligned along the x-axis)
   return Mat3::makeDiagonalMatrix(Ia, Ib, Ib);
}
//*************************************************************************************************

//=================================================================================================
//
//  OUTPUT FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Output of the current state of a capsule.
 *
 * \param os Reference to the output stream.
 * \param tab Indentation in front of every line of the capsule output.
 * \return void
 */
void Capsule::print( std::ostream& os, const char* tab ) const
{
   using std::setw;

   os << tab << " Capsule " << getID() << " with radius " << getRadius() << " and length " << getLength() << "\n";

   os << tab << "   Fixed: " << isFixed() << " , sleeping: " << !isAwake() << "\n";

   os << tab << "   System ID         = " << getSystemID() << "\n"
      << tab << "   Total mass        = " << getMass() << "\n"
      << tab << "   Material          = " << Material::getName( material_ ) << "\n"
      << tab << "   Global position   = " << getPosition() << "\n"
      << tab << "   Relative position = " << getRelPosition() << "\n"
      << tab << "   Linear velocity   = " << getLinearVel() << "\n"
      << tab << "   Angular velocity  = " << getAngularVel() << "\n";

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
}
//*************************************************************************************************

//=================================================================================================
//
//  GLOBAL OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Global output operator for capsules.
 *
 * \param os Reference to the output stream.
 * \param c Reference to a constant capsule object.
 * \return Reference to the output stream.
 */
std::ostream& operator<<( std::ostream& os, const Capsule& c )
{
   os << "--" << "CAPSULE PARAMETERS"
      << "------------------------------------------------------------\n";
   c.print( os, "" );
   os << "--------------------------------------------------------------------------------\n"
      << std::endl;
   return os;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Global output operator for capsule handles.
 *
 * \param os Reference to the output stream.
 * \param c Constant capsule handle.
 * \return Reference to the output stream.
 */
std::ostream& operator<<( std::ostream& os, ConstCapsuleID c )
{
   os << "--" << "CAPSULE PARAMETERS"
      << "------------------------------------------------------------\n";
   c->print( os, "" );
   os << "--------------------------------------------------------------------------------\n"
      << std::endl;
   return os;
}
//*************************************************************************************************

id_t Capsule::staticTypeID_ = std::numeric_limits<id_t>::max();

} // namespace pe
} // namespace walberla


