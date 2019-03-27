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
//! \file Sphere.h
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
 * \brief Base class for the sphere geometry.
 *
 * The Sphere class represents the base class for the sphere geometry. It provides
 * the basic functionality of a sphere. For a full description of the sphere geometry,
 * see the Sphere class description.
 */
class Sphere : public GeomPrimitive
{
public:
   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit Sphere( id_t sid, id_t uid, const Vec3& gpos, const Vec3& rpos, const Quat& q,
                    real_t radius, MaterialID material,
                    const bool global, const bool communicating, const bool infiniteMass );
   explicit Sphere( id_t const typeID, id_t sid, id_t uid, const Vec3& gpos, const Vec3& rpos, const Quat& q,
                    real_t radius, MaterialID material,
                    const bool global, const bool communicating, const bool infiniteMass );
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   /*!\name Destructor */
   //@{
   virtual ~Sphere();
   //@}
   //**********************************************************************************************
   //**********************************************************************************************

public:
   //**Get functions*******************************************************************************
   /*!\name Get functions */
   //@{
   inline real_t getRadius() const;
   virtual inline real_t getVolume()         const;
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   inline real_t getRelDepth   ( real_t px, real_t py, real_t pz ) const;
   inline real_t getRelDepth   ( const Vec3& rpos )          const;
   inline real_t getDepth      ( real_t px, real_t py, real_t pz ) const;
   inline real_t getDepth      ( const Vec3& gpos )          const;
   inline real_t getRelDistance( real_t px, real_t py, real_t pz ) const;
   inline real_t getRelDistance( const Vec3& rpos )          const;
   inline real_t getDistance   ( real_t px, real_t py, real_t pz ) const;
   inline real_t getDistance   ( const Vec3& gpos )          const;
   static inline id_t getStaticTypeID();
   //@}
   //**********************************************************************************************

   //**Output functions****************************************************************************
   /*!\name Output functions */
   //@{
   virtual void print( std::ostream& os, const char* tab ) const;
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   inline virtual Vec3 support( const Vec3& d ) const;
   inline virtual Vec3 supportContactThreshold( const Vec3& d ) const;
   //@}
   //**********************************************************************************************

   //**Volume, mass and density functions**********************************************************
   /*!\name Volume, mass and density functions */
   //@{
   static inline real_t calcVolume( real_t radius );
   static inline real_t calcMass( real_t radius, real_t density );
   static inline real_t calcDensity( real_t radius, real_t mass );
   static inline Mat3   calcInertia( const real_t mass, const real_t radius );
   //@}
   //**********************************************************************************************

protected:
   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   virtual bool containsRelPointImpl ( real_t px, real_t py, real_t pz ) const;
   virtual bool isSurfaceRelPointImpl( real_t px, real_t py, real_t pz ) const;
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   inline virtual void calcBoundingBox();  // Calculation of the axis-aligned bounding box
   //@}
   //**********************************************************************************************

   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   real_t radius_;  //!< Radius of the sphere.
                  /*!< The radius is constrained to values larger than 0.0. */
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
/*!\brief Returns the radius of the sphere.
 *
 * \return The radius of the sphere.
 */
inline real_t Sphere::getRadius() const
{
   return radius_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the volume of the sphere.
 *
 * \return The volume of the sphere.
 */
inline real_t Sphere::getVolume() const
{
   return Sphere::calcVolume(getRadius());
}
//*************************************************************************************************




//=================================================================================================
//
//  VOLUME, MASS AND DENSITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Calculates the volume of a sphere for a given radius.
 *
 * \param radius The radius of the sphere.
 * \return The volume of the sphere.
 */
inline real_t Sphere::calcVolume( real_t radius )
{
   return real_c(4.0)/real_c(3.0) * math::M_PI * radius * radius * radius;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Calculates the mass of a sphere for a given radius and density.
 *
 * \param radius The radius of the sphere.
 * \param density The density of the sphere.
 * \return The total mass of the sphere.
 */
inline real_t Sphere::calcMass( real_t radius, real_t density )
{
   return real_c(4.0)/real_c(3.0) * math::M_PI * radius * radius * radius * density;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Calculates the density of a sphere for a given radius and mass.
 *
 * \param radius The radius of the sphere.
 * \param mass The total mass of the sphere.
 * \return The density of the sphere.
 */
inline real_t Sphere::calcDensity( real_t radius, real_t mass )
{
   return real_c(0.75) * mass / ( math::M_PI * radius * radius * radius );
}
//*************************************************************************************************




//=================================================================================================
//
//  UTILITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Calculation of the bounding box of the sphere.
 *
 * \return void
 *
 * This function updates the axis-aligned bounding box of the sphere primitive according to the
 * current position and orientation of the sphere. Note that the bounding box is increased in
 * all dimensions by pe::contactThreshold to guarantee that rigid bodies in close proximity of
 * the sphere are also considered during the collision detection process.
 */
inline void Sphere::calcBoundingBox()
{
   const real_t length( radius_ + contactThreshold );

   aabb_    = AABB(
              gpos_[0] - length,
              gpos_[1] - length,
              gpos_[2] - length,
              gpos_[0] + length,
              gpos_[1] + length,
              gpos_[2] + length
   );

//   WALBERLA_ASSERT( aabb_.isValid()        , "Invalid bounding box detected" );
   WALBERLA_ASSERT( aabb_.contains( gpos_ ), "Invalid bounding box detected("<< getSystemID() <<")\n" << "pos: " << gpos_ << "\nlength: " << length << "\nvel: " << getLinearVel() << "\nbox: " << aabb_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Calculation of the moment of inertia in reference to the body frame of the sphere.
 *
 * \return void
 */
inline Mat3 Sphere::calcInertia( const real_t mass, const real_t radius )
{
   return Mat3::makeDiagonalMatrix( real_c(0.4) * mass * radius * radius );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Estimates the point which is farthest in direction \a d.
 *
 * \param d The normalized search direction in world-frame coordinates.
 * \return The support point in world-frame coordinates in direction a\ d.
 */
inline Vec3 Sphere::support( const Vec3& d ) const
{
   auto len = d.sqrLength();
   if (!math::equal(len, real_t(0)))
   {
      //WALBERLA_ASSERT_FLOAT_EQUAL( len, real_t(1), "search direction not normalized!");
      const Vec3 s = getPosition() + (getRadius() / sqrt(len)) * d;
      //std::cout << "Support in direction " << d << " with center " << getPosition() << " (r=" << getRadius() << ") is " << s << std::endl;
      return s;
   } else
   {
      return Vec3(0,0,0);
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Estimates the point which is farthest in direction \a d.
 *
 * \param d The normalized search direction in world-frame coordinates
 * \return The support point in world-frame coordinates in direction a\ d extended by a vector in
 *         direction \a d of length \a pe::contactThreshold.
 */
inline Vec3 Sphere::supportContactThreshold( const Vec3& d ) const
{
   auto len = d.sqrLength();
   if (!math::equal(len, real_t(0)))
   {
      //WALBERLA_ASSERT_FLOAT_EQUAL( len, real_t(1), "search direction not normalized!");
      const Vec3 s = getPosition() + (getRadius() / sqrt(len) + contactThreshold) * d;
      //std::cout << "Support in direction " << d << " with center " << getPosition() << " (r=" << getRadius() << ") is " << s << std::endl;
      return s;
   } else
   {
      return Vec3(0,0,0);
   }
}
//*************************************************************************************************

//=================================================================================================
//
//  UTILITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Calculates the depth of a point relative to the sphere's geometric center.
 *
 * \param px The x-component of the relative coordinate.
 * \param py The y-component of the relative coordinate.
 * \param pz The z-component of the relative coordinate.
 * \return Depth of the relative point.
 *
 * Returns a positive value, if the point lies inside the sphere and a negative value,
 * if the point lies outside the sphere.
 */
inline real_t Sphere::getRelDepth( real_t px, real_t py, real_t pz ) const
{
   return getRelDepth( Vec3( px, py, pz ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Calculates the depth of a point relative to the sphere's geometric center.
 *
 * \param rpos The relative coordinate.
 * \return Depth of the relative point.
 *
 * Returns a positive value, if the point lies inside the sphere and a negative value,
 * if the point lies outside the sphere.
 */
inline real_t Sphere::getRelDepth( const Vec3& rpos ) const
{
   return ( radius_ - rpos.length() );
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
 * Returns a positive value, if the point lies inside the sphere and a negative value,
 * if the point lies outside the sphere.
 */
inline real_t Sphere::getDepth( real_t px, real_t py, real_t pz ) const
{
   return getDepth( Vec3( px, py, pz ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Calculates the depth of a point in global coordinates.
 *
 * \param gpos The global coordinate.
 * \return Depth of the global point.
 *
 * Returns a positive value, if the point lies inside the sphere and a negative value,
 * if the point lies outside the sphere.
 */
inline real_t Sphere::getDepth( const Vec3& gpos ) const
{
   return ( radius_ - ( gpos - gpos_ ).length() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Calculates the distance of a point relative to the sphere's frame of reference.
 *
 * \param px The x-component of the relative coordinate.
 * \param py The y-component of the relative coordinate.
 * \param pz The z-component of the relative coordinate.
 * \return Distance of the relative point.
 *
 * Returns a positive value, if the point lies outside the sphere and a negative value, if the
 * point lies inside the sphere.
 */
inline real_t Sphere::getRelDistance( real_t px, real_t py, real_t pz ) const
{
   return getRelDistance( Vec3( px, py, pz ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Calculates the distance of a point relative to the sphere's frame of reference.
 *
 * \param rpos The relative coordinate.
 * \return Distance of the relative point.
 *
 * Returns a positive value, if the point lies outside the sphere and a negative value, if the
 * point lies inside the sphere.
 */
inline real_t Sphere::getRelDistance( const Vec3& rpos ) const
{
   return ( rpos.length() - radius_ );
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
 * Returns a positive value, if the point lies outside the sphere and a negative value, if the
 * point lies inside the sphere.
 */
inline real_t Sphere::getDistance( real_t px, real_t py, real_t pz ) const
{
   return getDistance( Vec3( px, py, pz ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Calculates the distance of a point in global coordinates.
 *
 * \param gpos The global coordinate.
 * \return Distance of the global point.
 *
 * Returns a positive value, if the point lies outside the sphere and a negative value, if the
 * point lies inside the sphere.
 */
inline real_t Sphere::getDistance( const Vec3& gpos ) const
{
   return ( ( gpos - gpos_ ).length() - radius_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns unique type id of this type
 *
 * \return geometry specific type id
 */
inline id_t Sphere::getStaticTypeID()
{
   return staticTypeID_;
}




//=================================================================================================
//
//  GLOBAL OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\name Sphere operators */
//@{
std::ostream& operator<<( std::ostream& os, const Sphere& s );
std::ostream& operator<<( std::ostream& os, ConstSphereID s );
//@}
//*************************************************************************************************


} // namespace pe
}  // namespace walberla
