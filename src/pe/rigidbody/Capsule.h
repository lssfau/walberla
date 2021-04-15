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
//! \file Capsule.h
//! \author Klaus Iglberger
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
#include <core/math/Constants.h>
#include <core/math/Limits.h>
#include <core/math/Shims.h>
#include <core/math/Utility.h>
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
 * \brief Capsule geometry.
 *
 * The Capsule class represents the geometric primitive capsule, which is one of the basic
 * geometric primitives of the \b pe physics module. The class is derived from the GeomPrimitive
 * base class, which makes the capsule both a geometric primitive and a rigid body.\n
 * A capsule is the combination of a cylinder and two hemisphere caps at both ends of the cylinder.
 * This combination allows to calculate a unique normal on each point of the capsule's surface.
 * In order to setup a capsule, the two values radius and length are required: radius specifies
 * the radius of both the cylinder part and the two hemispheres, length is the length of the
 * cylinder part. A capsule is created axis-aligned with the x-axis (the length of the cylinder
 * part is parallel to the x-axis).
 */
class Capsule : public GeomPrimitive
{
public:
   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit Capsule( id_t sid, id_t uid, const Vec3& gpos,  const Quat& q,
                     real_t  radius, real_t  length, MaterialID material,
                     const bool global, const bool communicating, const bool infiniteMass );
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   /*!\name Destructor */
   //@{
   ~Capsule() override;
   //@}
   //**********************************************************************************************

public:
   //**Get functions*******************************************************************************
   /*!\name Get functions */
   //@{
   inline real_t  getRadius() const;
   inline real_t  getLength() const;
   inline real_t  getVolume() const override;
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   inline Vec3 support( const Vec3& d ) const override;
   //@}
   //**********************************************************************************************

   //**Volume, mass and density functions**********************************************************
   /*!\name Volume, mass and density functions */
   //@{
   static inline real_t  calcVolume( real_t  radius, real_t  length );
   static inline real_t  calcMass( real_t  radius, real_t  length, real_t  density );
   static inline real_t  calcDensity( real_t  radius, real_t  length, real_t  mass );
   static        Mat3    calcInertia( const real_t radius, const real_t length, const real_t density);
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   inline real_t  getRelDepth   ( real_t  px, real_t  py, real_t  pz ) const;
   inline real_t  getRelDepth   ( const Vec3& rpos )          const;
   inline real_t  getDepth      ( real_t  px, real_t  py, real_t  pz ) const;
   inline real_t  getDepth      ( const Vec3& gpos )          const;
   inline real_t  getRelDistance( real_t  px, real_t  py, real_t  pz ) const;
   inline real_t  getRelDistance( const Vec3& rpos )          const;
   inline real_t  getDistance   ( real_t  px, real_t  py, real_t  pz ) const;
   inline real_t  getDistance   ( const Vec3& gpos )          const;
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
   bool containsRelPointImpl ( real_t  px, real_t  py, real_t  pz ) const override;
   bool isSurfaceRelPointImpl( real_t  px, real_t  py, real_t  pz ) const override;

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   void calcBoundingBox() override;  // Calculation of the axis-aligned bounding box
   //@}
   //**********************************************************************************************

   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   real_t  radius_;  //!< The radius of the cylinder part and the caps on both ends of the cylinder.
   real_t  length_;  //!< The length of the cylinder part.
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
/*!\brief Returns the radius of the capsule.
 *
 * \return The radius of the capsule.
 */
inline real_t  Capsule::getRadius() const
{
   return radius_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the length of the cylinder part.
 *
 * \return The length of the cylinder part.
 */
inline real_t  Capsule::getLength() const
{
   return length_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the volume of the capsule.
 *
 * \return The volume of the capsule.
 */
inline real_t  Capsule::getVolume() const
{
   return Capsule::calcVolume( getRadius(), getLength() );
}
//*************************************************************************************************




//=================================================================================================
//
//  VOLUME, MASS AND DENSITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Calculates the volume of a capsule for a given radius and length.
 *
 * \param radius The radius of the cylinder part and the caps on both ends of the cylinder.
 * \param length The length of the cylinder part.
 * \return The volume of the capsule.
 */
inline real_t  Capsule::calcVolume( real_t  radius, real_t  length )
{
   return math::pi*radius*radius*( ( static_cast<real_t >( 4 ) / static_cast<real_t >( 3 ) )*radius + length );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Calculates the mass of a capsule for a given radius, length and density.
 *
 * \param radius The radius of the cylinder part and the caps on both ends of the cylinder.
 * \param length The length of the cylinder part.
 * \param density The density of the capsule.
 * \return The total mass of the capsule.
 */
inline real_t  Capsule::calcMass( real_t  radius, real_t  length, real_t  density )
{
   return math::pi*radius*radius*( ( static_cast<real_t >( 4 ) / static_cast<real_t >( 3 ) )*radius + length ) * density;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Calculates the density of a capsule for a given radius, length and mass.
 *
 * \param radius The radius of the cylinder part and the caps on both ends of the cylinder.
 * \param length The length of the cylinder part.
 * \param mass The total mass of the capsule.
 * \return The density of the capsule.
 */
inline real_t  Capsule::calcDensity( real_t  radius, real_t  length, real_t  mass )
{
   return mass / ( math::pi*radius*radius*( ( static_cast<real_t >( 4 ) / static_cast<real_t >( 3 ) )*radius + length ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Estimates the point which is farthest in direction \a d.
 *
 * \param d The normalized search direction in world-frame coordinates.
 * \return The support point in world-frame coordinates in direction \a d.
 */
inline Vec3 Capsule::support( const Vec3& d ) const
{
   auto len = d.sqrLength();
   if (math::equal(len, real_t(0)))
      return Vec3(0,0,0);

   Vec3 dnorm = d / sqrt(len);

   const Vec3 bfD = vectorFromWFtoBF(dnorm); //d in body frame coordinates

   const Vec3 supportSegment = Vec3( math::sign(bfD[0])*length_*real_t (0.5), real_t (0.0), real_t (0.0));
   const Vec3 supportSphere = radius_ * dnorm;

   return getPosition() + vectorFromBFtoWF(supportSegment) + supportSphere;
}
//*************************************************************************************************



//=================================================================================================
//
//  UTILITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Calculates the depth of a point relative to the capsule's frame of reference.
 *
 * \param px The x-component of the relative coordinate.
 * \param py The y-component of the relative coordinate.
 * \param pz The z-component of the relative coordinate.
 * \return Depth of the relative point.
 *
 * Returns a positive value, if the point lies inside the capsule and a negative value,
 * if the point lies outside the capsule.
 */
inline real_t  Capsule::getRelDepth( real_t  px, real_t  py, real_t  pz ) const
{
   const real_t  xabs( std::fabs( px ) );         // Absolute x-distance
   const real_t  hlength( real_t (0.5) * length_ );  // Capsule half length

   if( xabs > hlength ) {
      return ( radius_ - std::sqrt( math::sq(xabs-hlength) + math::sq(py) + math::sq(pz) ) );
   }
   else {
      return ( radius_ - std::sqrt( math::sq(py) + math::sq(pz) ) );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Calculates the depth of a point relative to the capsule's frame of reference.
 *
 * \param rpos The relative coordinate.
 * \return Depth of the relative point.
 *
 * Returns a positive value, if the point lies inside the capsule and a negative value,
 * if the point lies outside the capsule.
 */
inline real_t  Capsule::getRelDepth( const Vec3& rpos ) const
{
   return getRelDepth( rpos[0], rpos[1], rpos[2] );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Calculates the depth of a point in global coordinates.
 *
 * \param px The x-component of the global coordinate.
 * \param py The y-component of the global coordinate.
 * \param pz The z-component of the global coordinate.
 * \return Depth of the global point.
 *
 * Returns a positive value, if the point lies inside the capsule and a negative value,
 * if the point lies outside the capsule.
 */
inline real_t  Capsule::getDepth( real_t  px, real_t  py, real_t  pz ) const
{
   return getRelDepth( pointFromWFtoBF( px, py, pz ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Calculates the depth of a point in global coordinates.
 *
 * \param gpos The global coordinate.
 * \return Depth of the global point.
 *
 * Returns a positive value, if the point lies inside the capsule and a negative value,
 * if the point lies outside the capsule.
 */
inline real_t  Capsule::getDepth( const Vec3& gpos ) const
{
   return getRelDepth( pointFromWFtoBF( gpos ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Calculates the distance of a point relative to the capsule's frame of reference.
 *
 * \param px The x-component of the relative coordinate.
 * \param py The y-component of the relative coordinate.
 * \param pz The z-component of the relative coordinate.
 * \return Distance of the relative point.
 *
 * Returns a positive value, if the point lies outside the capsule and a negative value,
 * if the point lies inside the capsule.
 */
inline real_t  Capsule::getRelDistance( real_t  px, real_t  py, real_t  pz ) const
{
   return -getRelDepth( px, py, pz );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Calculates the distance of a point relative to the capsule's frame of reference.
 *
 * \param rpos The relative coordinate.
 * \return Distance of the relative point.
 *
 * Returns a positive value, if the point lies outside the capsule and a negative value,
 * if the point lies inside the capsule.
 */
inline real_t  Capsule::getRelDistance( const Vec3& rpos ) const
{
   return -getRelDepth( rpos );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Calculates the distance of a point in global coordinates.
 *
 * \param px The x-component of the global coordinate.
 * \param py The y-component of the global coordinate.
 * \param pz The z-component of the global coordinate.
 * \return Distance of the global point.
 *
 * Returns a positive value, if the point lies outside the capsule and a negative value,
 * if the point lies inside the capsule.
 */
inline real_t  Capsule::getDistance( real_t  px, real_t  py, real_t  pz ) const
{
   return -getRelDepth( pointFromWFtoBF( px, py, pz ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Calculates the distance of a point in global coordinates.
 *
 * \param gpos The global coordinate.
 * \return Distance of the global point.
 *
 * Returns a positive value, if the point lies outside the capsule and a negative value,
 * if the point lies inside the capsule.
 */
inline real_t  Capsule::getDistance( const Vec3& gpos ) const
{
   return -getRelDepth( pointFromWFtoBF( gpos ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns unique type id of this type
 *
 * \return geometry specific type id
 */
inline id_t Capsule::getStaticTypeID()
{
   return staticTypeID_;
}




//=================================================================================================
//
//  GLOBAL OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\name Capsule operators */
//@{
std::ostream& operator<<( std::ostream& os, const Capsule& cc );
std::ostream& operator<<( std::ostream& os, ConstCapsuleID cc );
//@}
//*************************************************************************************************

} // namespace pe
} // namespace walberla
