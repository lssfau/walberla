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
//! \file Plane.cpp
//! \author Klaus Iglberger
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include "Plane.h"

#include <cmath>

#include <pe/Thresholds.h>
#include "pe/Materials.h"

#include <core/math/Matrix3.h>
#include <core/math/Quaternion.h>
#include "core/math/Limits.h"
#include "core/debug/Debug.h"
#include "core/logging/Logging.h"


namespace walberla {
namespace pe {

//=================================================================================================
//
//  CONSTRUCTOR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Constructor for the Plane class.
 *
 * \param sid Unique system-specific ID for the plane.
 * \param uid User-specific ID for the plane.
 * \param gpos The global position (anchor point) of the plane.
 * \param normal The plane's normal in reference to the global world frame, \f$ |n| = 1 \f$.
 * \param d The displacement of the plane.
 * \param material The material of the plane.
 *
 * The plane equation is: \f$ ax + by + cz = d \f$.\n
 * \a a, \a b and \a c are the x, y and z coordinate of the normal vector and \a d is the distance
 * from the origin to the plane.
 */
Plane::Plane( id_t sid, id_t uid,
              const Vec3& gpos, const Vec3& normal, real_t d,
              MaterialID material )
   : GeomPrimitive( getStaticTypeID(), sid, uid, material )  // Initialization of the parent class
   , normal_( normal )                                                       // Normal of the plane
   , d_( d )
{
   //planes are always considered locally and have infinite mass
   setGlobal( true );
   setCommunicating( false );
   setMassAndInertiaToInfinity();
   setFinite( false );

   // Checking the mass properties of the plane
   WALBERLA_ASSERT            ( std::isinf(mass_), "Mass of plane is not infinity" );
   WALBERLA_ASSERT_FLOAT_EQUAL( invMass_ , real_c(0), "Inverse mass of plane is not 0" );
   WALBERLA_ASSERT            ( isinf(I_), "Moment of inertia of plane is not infinity\n" << I_ );
   WALBERLA_ASSERT_FLOAT_EQUAL( Iinv_    , Mat3(real_c(0)), "Inverse moment of inertia of plane is not 0" );

   // Checking the plane normal
   // Since the plane constructor is never directly called but only used in a small number
   // of functions that already check the plane arguments, only asserts are used here to
   // double check the arguments.
   WALBERLA_ASSERT_FLOAT_EQUAL( normal.sqrLength(), real_c(1) , "Invalid plane normal" );

   // Setting the global position (anchor point) of the plane
   setPosition(gpos);

   // Calculating the orientation and rotation
   // The default normal of a plane is <0,0,1>. The rotation of the plane is calculated
   // as the rotation of this default normal to the specified normal.
   if( normal[0]*normal[0] + normal[1]*normal[1] < math::Limits<real_t>::accuracy() )
      setOrientation( Quat(Vec3( 1, 0, 0 ), std::acos( normal[2] ) ) );
   else
      setOrientation( Quat( Vec3( -normal[1], normal[0], real_c(0) ), std::acos( normal[2] ) ) );

   // Setting the axis-aligned bounding box
   Plane::calcBoundingBox();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Destructor for the PlaneBase class.
 */
Plane::~Plane()
{
   // Logging the destruction of the plane
   WALBERLA_LOG_DETAIL( "Destroyed plane " << sid_ );
}
//*************************************************************************************************




//=================================================================================================
//
//  UTILITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Calculation of the bounding box of the plane.
 *
 * \return void
 *
 * This function updates the axis-aligned bounding box of the plane primitive according to the
 * current position and orientation of the plane. Note that the bounding box is increased in
 * all dimensions by pe::contactThreshold to guarantee that rigid bodies in close proximity of
 * the plane are also considered during the collision detection process.
 */
void Plane::calcBoundingBox()
{
   aabb_ = AABB( ( normal_[0] <  real_c(0) && floatIsEqual(normal_[1], real_c(0)) && floatIsEqual(normal_[2], real_c(0)) )
                 ? ( -d_ - contactThreshold ) : ( -math::Limits<real_t>::inf() ),
                 ( floatIsEqual(normal_[0], real_c(0)) && normal_[1] <  real_c(0) && floatIsEqual(normal_[2], real_c(0)) )
                 ? ( -d_ - contactThreshold ) : ( -math::Limits<real_t>::inf() ),
                 ( floatIsEqual(normal_[0], real_c(0)) && floatIsEqual(normal_[1], real_c(0)) && normal_[2] <  real_c(0) )
                 ? ( -d_ - contactThreshold ) : ( -math::Limits<real_t>::inf() ),
                 ( normal_[0] >  real_c(0) && floatIsEqual(normal_[1], real_c(0)) && floatIsEqual(normal_[2], real_c(0)) )
                 ? (  d_ + contactThreshold ) : (  math::Limits<real_t>::inf() ),
                 ( floatIsEqual(normal_[0], real_c(0)) && normal_[1] >  real_c(0) && floatIsEqual(normal_[2], real_c(0)) )
                 ? (  d_ + contactThreshold ) : (  math::Limits<real_t>::inf() ),
                 ( floatIsEqual(normal_[0], real_c(0)) && floatIsEqual(normal_[1], real_c(0)) && normal_[2] >  real_c(0) )
                 ? (  d_ + contactThreshold ) : (  math::Limits<real_t>::inf() ) );
}
//*************************************************************************************************

//=================================================================================================
//
//  SET FUNCTIONS
//
//=================================================================================================


//*************************************************************************************************
/*!\brief Setting the global position of the plane.
 *
 * \param px The x-component of the global position.
 * \param py The y-component of the global position.
 * \param pz The z-component of the global position.
 * \return void
 * \exception std::logic_error Invalid translation of a plane inside an exclusive section.
 *
 * This function sets the global position (anchor point) of the plane.
 *
 * \b Note:
 * - Setting the position of a plane contained in a union changes the geometry of the union.
 *   Therefore this may cause an invalidation of links contained in the union.
 */
void Plane::setPositionImpl( real_t px, real_t py, real_t pz )
{
   RigidBody::setPositionImpl(px,py,pz);
   d_ = normal_ * getPosition();
   Plane::calcBoundingBox();    // Updating the axis-aligned bounding box of the plane
#if MOBILE_INFINITE
   wake();               // Waking the box from sleep mode
#else
   // As an immobile, infinite body a plane doesn't have to be awakened from sleep mode!
#endif
   signalTranslation();  // Signaling the position change to the superordinate body
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Setting the global orientation of the plane.
 *
 * \param r The quaternion scalar part.
 * \param i,j,k The quaternion vector part.
 * \return void
 * \exception std::logic_error Invalid rotation of a plane inside an exclusive section.
 *
 * This function sets the normal of the plane corresponding to the given global orientation,
 * where the initial orientation of the plane's normal is (0,0,1). This change of orientation
 * corresponds to a rotation around the anchor point of the plane, which is consequently not
 * changed. However, this rotation changes the distance (displacement) to the origin of the
 * global world frame.
 *
 * \b Note:
 * - Setting the orientation of a plane contained in a union changes the geometry of the union.
 *   Therefore this may cause an invalidation of links contained in the union.
 */
void Plane::setOrientationImpl( real_t r, real_t i, real_t j, real_t k )
{
   RigidBody::setOrientationImpl(r,i,j,k);
   Mat3 R = getRotation();
   // Updating the normal of the plane ( R * <0,0,1> )
   normal_[0] = R[2];
   normal_[1] = R[5];
   normal_[2] = R[8];

   // Updating the displacement from the origin
   d_ = normal_ * getPosition();

   Plane::calcBoundingBox();  // Updating the axis-aligned bounding box of the plane
#if MOBILE_INFINITE
   wake();               // Waking the box from sleep mode
#else
   // As an immobile, infinite body a plane doesn't have to be awakened from sleep mode!
#endif
   signalRotation();   // Signaling the change of orientation to the superordinate body
}
//*************************************************************************************************




//=================================================================================================
//
//  TRANSLATION FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Translation of the global position of the plane by the displacement vector
 * \brief (\a dx,\a dy,\a dz).
 *
 * \param dx The x-component of the translation/displacement.
 * \param dy The y-component of the translation/displacement.
 * \param dz The z-component of the translation/displacement.
 * \return void
 * \exception std::logic_error Invalid translation of a plane inside an exclusive section.
 *
 * \b Note:
 * - Translating a plane contained in a union changes the geometry of the union. Therefore this
 *   may cause an invalidation of links contained in the union.
 */
void Plane::translateImpl( real_t dx, real_t dy, real_t dz )
{

   setPosition(getPosition() + Vec3(dx,dy,dz));
   d_ = normal_ * getPosition();

   Plane::calcBoundingBox();    // Updating the axis-aligned bounding box
#if MOBILE_INFINITE
   wake();               // Waking the box from sleep mode
#else
   // As an immobile, infinite body a plane doesn't have to be awakened from sleep mode!
#endif
   signalTranslation();  // Signaling the position change to the superordinate body
}
//*************************************************************************************************




//=================================================================================================
//
//  ROTATION FUNCTIONS
//
//=================================================================================================


//*************************************************************************************************
/*!\brief Rotation of the plane by the quaternion \a dq.
 *
 * \param dq The quaternion for the rotation.
 * \return void
 * \exception std::logic_error Invalid rotation of a plane inside an exclusive section.
 *
 * Changing the orientation/rotation of the plane. The plane is rotated around its anchor point
 * (its current global position) by the quaternion \a dq. This rotation changes the normal of
 * the plane and its distance (displacement) to the origin of the global world frame.\n
 *
 * \b Note:
 * - Rotating a plane contained in a union changes the geometry of the union. Therefore this may
 *   cause an invalidation of links contained in the union.
 */
void Plane::rotateImpl( const Quat& dq )
{
//   if( ExclusiveSection::isActive() )
//      throw std::logic_error( "Invalid rotation of a plane inside an exclusive section" );

   setOrientation(dq * getQuaternion());

   // Updating the normal of the plane ( R * <0,0,1> )
   Mat3 R = getRotation();
   normal_[0] = R[2];
   normal_[1] = R[5];
   normal_[2] = R[8];

   // Updating the displacement from the origin
   d_ = normal_ * getPosition();

   Plane::calcBoundingBox();  // Updating the axis-aligned bounding box of the plane
#if MOBILE_INFINITE
   wake();               // Waking the box from sleep mode
#else
   // As an immobile, infinite body a plane doesn't have to be awakened from sleep mode!
#endif
   signalRotation();   // Signaling the change of orientation to the superordinate body
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Rotation of the plane around the origin of the global world frame.
 *
 * \param dq The quaternion for the rotation.
 * \return void
 * \exception std::logic_error Invalid rotation of a plane inside an exclusive section.
 *
 * Changing the orientation/rotation of the plane. The plane is rotated around the origin of
 * the global frame around by the quaternion \a dq.
 * Therefore the anchor point (global position) and the normal of the plane are changed, not
 * its distance (displacement) to the origin.\n
 *
 * \b Note:
 * - Rotating a plane contained in a union changes the geometry of the union. Therefore this may
 *   cause an invalidation of links contained in the union.
 */
void Plane::rotateAroundOriginImpl( const Quat& dq )
{

   setPosition(dq.rotate( getPosition() ));     // Updating the global position of the plane
   setOrientation(dq * getQuaternion());                // Updating the orientation of the plane


   // Updating the normal of the plane ( R * <0,0,1> )
   Mat3 R = getRotation();
   normal_[0] = R[2];
   normal_[1] = R[5];
   normal_[2] = R[8];

   Plane::calcBoundingBox();  // Updating the axis-aligned bounding box of the plane
   signalRotation();   // Signaling the change of orientation to the superordinate body

#if MOBILE_INFINITE
   wake();               // Waking the box from sleep mode
#else
   // As an immobile, infinite body a plane doesn't have to be awakened from sleep mode!
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Rotation of the plane around a specific global coordinate.
 *
 * \param point The global center of the rotation.
 * \param dq The global rotation.
 * \return void
 * \exception std::logic_error Invalid rotation of a plane inside an exclusive section.
 *
 * This function rotates the plane around the given global coordinate \a point and changes
 * the global position, the displacement from the origin, and the orientation/rotation of the
 * plane. The plane is rotated around the given axis \a axis by \a angle degrees (radian
 * measure).\n
 *
 * \b Note:
 * - Rotating a plane contained in a union changes the geometry of the union. Therefore this may
 *   cause an invalidation of links contained in the union.
 */
void Plane::rotateAroundPointImpl( const Vec3& point, const Quat& dq )
{
   const Vec3 dp( getPosition() - point );

   setPosition(point + dq.rotate( dp ) );  // Updating the global position of the box
   setOrientation(dq * getQuaternion());   // Updating the orientation of the box

   // Updating the normal of the plane ( R * <0,0,1> )
   Mat3 R = getRotation();
   normal_[0] = R[2];
   normal_[1] = R[5];
   normal_[2] = R[8];

   // Updating the displacement from the origin
   d_ = normal_ * getPosition();

   Plane::calcBoundingBox();  // Updating the axis-aligned bounding box of the plane
#if MOBILE_INFINITE
   wake();               // Waking the box from sleep mode
#else
   // As an immobile, infinite body a plane doesn't have to be awakened from sleep mode!
#endif
   signalRotation();   // Signaling the change of orientation to the superordinate body
}
//*************************************************************************************************




//=================================================================================================
//
//  UTILITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Checks, whether a point in body relative coordinates lies inside the plane.
 *
 * \param px The x-component of the relative coordinate.
 * \param py The y-component of the relative coordinate.
 * \param pz The z-component of the relative coordinate.
 * \return \a true if the point lies inside the plane, \a false if not.
 */
bool Plane::containsRelPointImpl( real_t /*px*/, real_t /*py*/, real_t pz ) const
{
   return pz <= real_c(0);
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checks, whether a point in body relative coordinates lies on the plane.
 *
 * \param px The x-component of the relative coordinate.
 * \param py The y-component of the relative coordinate.
 * \param pz The z-component of the relative coordinate.
 * \return \a true if the point lies on the plane, \a false if not.
 *
 * The tolerance level of the check is pe::surfaceThreshold.
 */
bool Plane::isSurfaceRelPointImpl( real_t /*px*/, real_t /*py*/, real_t pz ) const
{
   return std::fabs( pz ) <= surfaceThreshold;
}
//*************************************************************************************************


//=================================================================================================
//
//  OUTPUT FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Output of the current state of a plane.
 *
 * \param os Reference to the output stream.
 * \param tab Indentation in front of every line of the plane output.
 * \return void
 */
void Plane::print( std::ostream& os, const char* tab ) const
{
   using std::setw;

   os << tab << " Plane " << uid_ << " with normal " << normal_ << " and displacement " << d_ << "\n"
      << tab << "   System ID         = " << sid_ << "\n"
      << tab << "   Material          = " << Material::getName( material_ ) << "\n"
      << tab << "   Global position   = " << getPosition() << "\n";

   os << tab << "   Relative position = " << getRelPosition() << "\n"
      << tab << "   Bounding box      = " << aabb_ << "\n"
      << tab << "   Quaternion        = " << getQuaternion()<< "\n";
   Mat3 R = getRotation();
   os << tab << "   Rotation matrix   = ( " << setw(9) << R[0] << " , " << setw(9) << R[1] << " , " << setw(9) << R[2] << " )\n"
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
/*!\brief Global output operator for planes.
 *
 * \param os Reference to the output stream.
 * \param p Reference to a constant plane object.
 * \return Reference to the output stream.
 */
std::ostream& operator<<( std::ostream& os, const Plane& p )
{
   os << "--" << "PLANE PARAMETERS"
      << "--------------------------------------------------------------\n";
   p.print( os, "" );
   os << "--------------------------------------------------------------------------------\n"
      << std::endl;
   return os;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Global output operator for plane handles.
 *
 * \param os Reference to the output stream.
 * \param p Constant plane handle.
 * \return Reference to the output stream.
 */
std::ostream& operator<<( std::ostream& os, ConstPlaneID p )
{
   os << "--" << "PLANE PARAMETERS"
      << "--------------------------------------------------------------\n";
   p->print( os, "" );
   os << "--------------------------------------------------------------------------------\n"
      << std::endl;
   return os;
}
//*************************************************************************************************

id_t Plane::staticTypeID_ = std::numeric_limits<id_t>::max();

} // namespace pe
} // namespace walberla
