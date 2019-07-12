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
//! \file Box.cpp
//! \author Klaus Iglberger
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include "Box.h"

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
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Instantiation constructor for the Box class.
 *
 * \param sid Unique system-specific ID for the box.
 * \param uid User-specific ID for the box.
 * \param gpos Global geometric center of the box.
 * \param rpos The relative position within the body frame of a superordinate body.
 * \param q The orientation of the box's body frame in the global world frame.
 * \param lengths Side lengths of the box \f$ (0..\infty) \f$.
 * \param material The material of the box.
 * \param visible Specifies if the box is visible in a visualization.
 * \param fixed \a true to fix the box, \a false to unfix it.
 */
Box::Box( id_t sid, id_t uid, const Vec3& gpos, const Quat& q,
          const Vec3& lengths, MaterialID material,
          const bool global, const bool communicating, const bool infiniteMass )
   : GeomPrimitive( getStaticTypeID(), sid, uid, material )  // Initialization of the parent class
   , lengths_(lengths)
{
   // Checking the side lengths
   // Since the box constructor is never directly called but only used in a small number
   // of functions that already check the box arguments, only asserts are used here to
   // double check the arguments.
   WALBERLA_ASSERT_GREATER( lengths[0], real_t(0), "Invalid side length in x-dimension" );
   WALBERLA_ASSERT_GREATER( lengths[1], real_t(0), "Invalid side length in y-dimension" );
   WALBERLA_ASSERT_GREATER( lengths[2], real_t(0), "Invalid side length in z-dimension" );

   // Initializing the instantiated box
   setPosition(gpos);
   setOrientation(q);                      // Setting the orientation

   setGlobal( global );
   if (infiniteMass)
   {
      setMassAndInertiaToInfinity();
   } else
   {
      auto mass = calcMass( lengths, Material::getDensity( material ) );
      setMassAndInertia( mass, calcInertia( lengths, mass ) );
   }
   setCommunicating( communicating );
   setFinite( true );

   // Setting the axis-aligned bounding box
   Box::calcBoundingBox();
}
//*************************************************************************************************




//=================================================================================================
//
//  DESTRUCTOR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Destructor for the Box class.
 */
Box::~Box()
{
   // Logging the destruction of the box
   WALBERLA_LOG_DETAIL( "Destroyed box " << sid_ );
}
//*************************************************************************************************



//=================================================================================================
//
//  UTILITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Checks, whether a point in body relative coordinates lies inside the box.
 *
 * \param px The x-component of the relative coordinate.
 * \param py The y-component of the relative coordinate.
 * \param pz The z-component of the relative coordinate.
 * \return \a true if the point lies inside the box, \a false if not.
 */
bool Box::containsRelPointImpl( real_t px, real_t py, real_t pz ) const
{
   return std::fabs(px) <= real_t(0.5)*lengths_[0] &&
         std::fabs(py) <= real_t(0.5)*lengths_[1] &&
         std::fabs(pz) <= real_t(0.5)*lengths_[2];
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checks, whether a point in body relative coordinates lies on the surface of the box.
 *
 * \param px The x-component of the relative coordinate.
 * \param py The y-component of the relative coordinate.
 * \param pz The z-component of the relative coordinate.
 * \return \a true if the point lies on the surface of the box, \a false if not.
 *
 * The tolerance level of the check is pe::surfaceThreshold.
 */
bool Box::isSurfaceRelPointImpl( real_t px, real_t py, real_t pz ) const
{
   // Checking if the body relative point lies on one of the x-faces
   if( std::fabs( real_t(0.5)*lengths_[0] - std::fabs(px) ) <= surfaceThreshold &&
       std::fabs(py) < real_t(0.5)*lengths_[1] + surfaceThreshold &&
       std::fabs(pz) < real_t(0.5)*lengths_[2] + surfaceThreshold ) {
      return true;
   }
   // Checking if the body relative point lies on one of the y-faces
   else if( std::fabs( real_t(0.5)*lengths_[1] - std::fabs(py) ) <= surfaceThreshold &&
            std::fabs(pz) < real_t(0.5)*lengths_[2] + surfaceThreshold &&
            std::fabs(px) < real_t(0.5)*lengths_[0] + surfaceThreshold ) {
      return true;
   }
   // Checking if the body relative point lies on one of the z-faces
   else if( std::fabs( real_t(0.5)*lengths_[2] - std::fabs(pz) ) <= surfaceThreshold &&
            std::fabs(px) < real_t(0.5)*lengths_[0] + surfaceThreshold &&
            std::fabs(py) < real_t(0.5)*lengths_[1] + surfaceThreshold ) {
      return true;
   }
   else return false;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Calculates the depth of a point relative to the box's frame of reference.
 *
 * \param px The x-component of the relative coordinate.
 * \param py The y-component of the relative coordinate.
 * \param pz The z-component of the relative coordinate.
 * \return Depth of the relative point.
 *
 * Returns a positive value, if the point lies inside the box and a negative value, if the point
 * lies outside the box. The returned depth is calculated relative to the closest side of the box.
 */
real_t Box::getRelDepth( real_t px, real_t py, real_t pz ) const
{
   const real_t xdepth( std::fabs(px)-real_t(0.5)*lengths_[0] );
   const real_t ydepth( std::fabs(py)-real_t(0.5)*lengths_[1] );
   const real_t zdepth( std::fabs(pz)-real_t(0.5)*lengths_[2] );

   // Calculating the depth for relative points outside the box
   if( xdepth >= real_t(0) ) {
      if( ydepth >= real_t(0) ) {
         if( zdepth >= real_t(0) ) {
            return -std::sqrt( math::sq(xdepth) + math::sq(ydepth) + math::sq(zdepth) );
         }
         else return -std::sqrt( math::sq(xdepth) + math::sq(ydepth) );
      }
      else if( zdepth >= real_t(0) ) {
         return -std::sqrt( math::sq(xdepth) + math::sq(zdepth) );
      }
      else return -xdepth;
   }
   else if( ydepth >= real_t(0) ) {
      if( zdepth >= real_t(0) ) {
         return -std::sqrt( math::sq(ydepth) + math::sq(zdepth) );
      }
      else return -ydepth;
   }
   else if( zdepth >= real_t(0) ) {
      return -zdepth;
   }

   // Relative point is inside the box => calculating the
   // depth depending on the distance to the nearest face
   else return -std::max( xdepth, std::max(ydepth, zdepth) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Calculates the depth of a point relative to the box's frame of reference.
 *
 * \param rpos The relative coordinate.
 * \return Depth of the relative point.
 *
 * Returns a positive value, if the point lies inside the box and a negative value, if the point
 * lies outside the box. The returned depth is calculated relative to the closest side of the box.
 */
real_t Box::getRelDepth( const Vec3& rpos ) const
{
   return getRelDepth(rpos[0], rpos[1], rpos[2]);
}
//*************************************************************************************************





//=================================================================================================
//
//  OUTPUT FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Output of the current state of a box.
 *
 * \param os Reference to the output stream.
 * \param tab Indentation in front of every line of the box output.
 * \return void
 */
void Box::print( std::ostream& os, const char* tab ) const
{
   using std::setw;

   os << tab << " Box " << getID() << " with side lengths " << getLengths() << "\n";

   os << tab << "   Fixed: " << isFixed() << " , sleeping: " << !isAwake() << "\n";

   os << tab << "   System ID         = " << getSystemID() << "\n"
      << tab << "   Total mass        = " << getMass() << "\n"
      << tab << "   Global position   = " << getPosition() << "\n"
      << tab << "   Relative position = " << getRelPosition() << "\n"
      << tab << "   Linear velocity   = " << getLinearVel() << "\n"
      << tab << "   Angular velocity  = " << getAngularVel() << "\n";

   //if( verboseMode )
   {
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
}
//*************************************************************************************************


//=================================================================================================
//
//  UTILITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Calculation of the bounding box of the box.
 *
 * \return void
 *
 * This function updates the axis-aligned bounding box of the box primitive according to the
 * current position and orientation of the box. Note that the bounding box is increased in
 * all dimensions by pe::contactThreshold to guarantee that rigid bodies in close proximity
 * of the box are also considered during the collision detection process.
 */
void Box::calcBoundingBox()
{
   using std::fabs;
   Mat3 R = getRotation();
   const real_t xlength( real_t(0.5) * ( fabs(R[0]*lengths_[0]) + fabs(R[1]*lengths_[1]) + fabs(R[2]*lengths_[2]) ) + contactThreshold );
   const real_t ylength( real_t(0.5) * ( fabs(R[3]*lengths_[0]) + fabs(R[4]*lengths_[1]) + fabs(R[5]*lengths_[2]) ) + contactThreshold );
   const real_t zlength( real_t(0.5) * ( fabs(R[6]*lengths_[0]) + fabs(R[7]*lengths_[1]) + fabs(R[8]*lengths_[2]) ) + contactThreshold );
   aabb_ = math::AABB(
            getPosition()[0] - xlength,
         getPosition()[1] - ylength,
         getPosition()[2] - zlength,
         getPosition()[0] + xlength,
         getPosition()[1] + ylength,
         getPosition()[2] + zlength
         );

   //WALBERLA_ASSERT( aabb_.isValid()        , "Invalid bounding box detected" );
   WALBERLA_ASSERT( aabb_.contains( getPosition() ), "Invalid bounding box detected" );
}

//*************************************************************************************************


//*************************************************************************************************
/*!\brief Calculation of the moment of inertia in reference to the body frame of the box.
 *
 * \return void
 */
Mat3 Box::calcInertia(const Vec3& length, const real_t mass)
{
   return Mat3::makeDiagonalMatrix(
            mass/static_cast<real_t>( 12 ) * ( length[1]*length[1] + length[2]*length[2] ),
         mass/static_cast<real_t>( 12 ) * ( length[0]*length[0] + length[2]*length[2] ),
         mass/static_cast<real_t>( 12 ) * ( length[0]*length[0] + length[1]*length[1] ));
}
//*************************************************************************************************


//=================================================================================================
//
//  GLOBAL OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Global output operator for boxes.
 *
 * \param os Reference to the output stream.
 * \param b Reference to a constant box object.
 * \return Reference to the output stream.
 */
std::ostream& operator<<( std::ostream& os, const Box& b )
{
   os << "--" << "BOX PARAMETERS"
      << "----------------------------------------------------------------\n";
   b.print( os, "" );
   os << "--------------------------------------------------------------------------------\n"
      << std::endl;
   return os;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Global output operator for box handles.
 *
 * \param os Reference to the output stream.
 * \param b Constant box handle.
 * \return Reference to the output stream.
 */
std::ostream& operator<<( std::ostream& os, ConstBoxID b )
{
   os << "--" << "BOX PARAMETERS"
      << "----------------------------------------------------------------\n";
   b->print( os, "" );
   os << "--------------------------------------------------------------------------------\n"
      << std::endl;
   return os;
}
//*************************************************************************************************

id_t Box::staticTypeID_ = std::numeric_limits<id_t>::max();

} // namespace pe
} // namespace walberla


