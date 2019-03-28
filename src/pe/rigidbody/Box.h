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
//! \file Box.h
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

/**
 * \ingroup pe
 * \brief Box geometry.
 *
 * The Box class represents the geometric primitive box, which is one of the basic geometric
 * primitives of the \b pe physics module. The class is derived from the GeomPrimitive base class,
 * which makes the box both a geometric primitive and a rigid body.\n
 * A box is created axis-aligned with the global coordinate system, where its geometric center is
 * exactly in the middle of the box, the x-side is aligned with the x-axis, the y-side is aligned
 * with the y-axis and the z-side is aligned with the z-axis.
 */
class Box : public GeomPrimitive
{
public:
   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit Box( id_t sid, id_t uid, const Vec3& gpos, const Vec3& rpos, const Quat& q,
                 const Vec3& lengths, MaterialID material,
                 const bool global, const bool communicating, const bool infiniteMass );
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   /*!\name Destructor */
   //@{
   virtual ~Box();
   //@}
   //**********************************************************************************************

public:
   //**Get functions*******************************************************************************
   /*!\name Get functions */
   //@{
   inline const Vec3& getLengths() const;
   virtual inline real_t getVolume() const;
   //@}
   //**********************************************************************************************

   //**Volume, mass and density functions**********************************************************
   /*!\name Volume, mass and density functions */
   //@{
   static inline real_t calcVolume( real_t x, real_t y, real_t z );
   static inline real_t calcVolume( const Vec3& l );
   static inline real_t calcMass( real_t x, real_t y, real_t z, real_t density );
   static inline real_t calcMass( const Vec3& l, real_t density );
   static inline real_t calcDensity( real_t x, real_t y, real_t z, real_t mass );
   static inline real_t calcDensity( const Vec3& l, real_t mass );
   static inline Mat3   calcInertia( const Vec3& length, const real_t mass);      // Calculation of the moment of inertia
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
          real_t getRelDepth   ( real_t px, real_t py, real_t pz ) const;
          real_t getRelDepth   ( const Vec3& rpos )          const;
   inline real_t getDepth      ( real_t px, real_t py, real_t pz ) const;
   inline real_t getDepth      ( const Vec3& gpos )          const;
   inline real_t getRelDistance( real_t px, real_t py, real_t pz ) const;
   inline real_t getRelDistance( const Vec3& rpos )          const;
   inline real_t getDistance   ( real_t px, real_t py, real_t pz ) const;
   inline real_t getDistance   ( const Vec3& gpos )          const;

   static inline id_t getStaticTypeID();
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   inline virtual Vec3 support( const Vec3& d ) const;
   inline virtual Vec3 supportContactThreshold( const Vec3& d ) const;
   //@}
   //**********************************************************************************************

   //**Output functions****************************************************************************
   /*!\name Output functions */
   //@{
   virtual void print( std::ostream& os, const char* tab ) const;
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
   virtual void calcBoundingBox();  // Calculation of the axis-aligned bounding box
   //@}
   //**********************************************************************************************

   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   Vec3 lengths_;  //!< Lengths of the x-, y- and z-sides of the box.
                   /*!< All side lengths are constrained to values larger than 0. */
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
//  UTILITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Calculates the depth of a point in global coordinates.
 *
 * \param px The x-component of the global coordinate.
 * \param py The y-component of the global coordinate.
 * \param pz The z-component of the global coordinate.
 * \return Depth of the global point.
 *
 * Returns a positive value, if the point lies inside the box and a negative value, if the point
 * lies outside the box. The returned depth is calculated relative to the closest side of the box.
 */
inline real_t Box::getDepth( real_t px, real_t py, real_t pz ) const
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
 * Returns a positive value, if the point lies inside the box and a negative value, if the point
 * lies outside the box. The returned depth is calculated relative to the closest side of the box.
 */
inline real_t Box::getDepth( const Vec3& gpos ) const
{
   return getDepth( gpos[0], gpos[1], gpos[2] );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Calculates the distance of a point relative to the box's frame of reference.
 *
 * \param px The x-component of the relative coordinate.
 * \param py The y-component of the relative coordinate.
 * \param pz The z-component of the relative coordinate.
 * \return Distance of the relative point.
 *
 * Returns a positive value, if the point lies outside the box and a negative value, if the
 * point lies inside the box. The returned distance is calculated relative to the closest
 * side of the box.
 */
inline real_t Box::getRelDistance( real_t px, real_t py, real_t pz ) const
{
   return -getRelDepth( px, py, pz );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Calculates the distance of a point relative to the box's frame of reference.
 *
 * \param rpos The relative coordinate.
 * \return Distance of the relative point.
 *
 * Returns a positive value, if the point lies outside the box and a negative value, if the
 * point lies inside the box. The returned distance is calculated relative to the closest
 * side of the box.
 */
inline real_t Box::getRelDistance( const Vec3& rpos ) const
{
   return getRelDistance( rpos[0], rpos[1], rpos[2] );
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
 * Returns a positive value, if the point lies outside the box and a negative value, if the
 * point lies inside the box. The returned distance is calculated relative to the closest
 * side of the box.
 */
inline real_t Box::getDistance( real_t px, real_t py, real_t pz ) const
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
 * Returns a positive value, if the point lies outside the box and a negative value, if the
 * point lies inside the box. The returned distance is calculated relative to the closest
 * side of the box.
 */
inline real_t Box::getDistance( const Vec3& gpos ) const
{
   return getDistance( gpos[0], gpos[1], gpos[2] );
}
//*************************************************************************************************

//*************************************************************************************************
/*!\brief Returns unique type id of this type
 *
 * \return geometry specific type id
 */
inline id_t Box::getStaticTypeID()
{
   return staticTypeID_;
}

//=================================================================================================
//
//  GET FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Returns the side lengths of the box.
 *
 * \return The side lengths of the box.
 */
inline const Vec3& Box::getLengths() const
{
   return lengths_;
}
//*************************************************************************************************

//*************************************************************************************************
/*!\brief Returns the volume of the box.
 *
 * \return The volume of the box.
 */
inline real_t Box::getVolume() const
{
   return Box::calcVolume(getLengths());
}
//*************************************************************************************************


//=================================================================================================
//
//  VOLUME, MASS AND DENSITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Calculates the volume of a box for given side lengths.
 *
 * \param x The x-length of the box.
 * \param y The y-length of the box.
 * \param z The z-length of the box.
 * \return The volume of the box.
 */
inline real_t Box::calcVolume( real_t x, real_t y, real_t z )
{
   return x * y * z;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Calculates the volume of a box for given side lengths.
 *
 * \param l The side lengths of the box.
 * \return The volume of the box.
 */
inline real_t Box::calcVolume( const Vec3& l )
{
   return calcVolume(l[0], l[1], l[2]);
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Calculates the mass of a box for given side lengths and density.
 *
 * \param x The x-length of the box.
 * \param y The y-length of the box.
 * \param z The z-length of the box.
 * \param density The density of the box.
 * \return The total mass of the box.
 */
inline real_t Box::calcMass( real_t x, real_t y, real_t z, real_t density )
{
   return x * y * z * density;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Calculates the mass of a box for given side lengths and density.
 *
 * \param l The side lengths of the box.
 * \param density The density of the box.
 * \return The total mass of the box.
 */
inline real_t Box::calcMass( const Vec3& l, real_t density )
{
   return calcMass(l[0], l[1], l[2], density);
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Calculates the density of a box for given side lengths and mass.
 *
 * \param x The x-length of the box.
 * \param y The y-length of the box.
 * \param z The z-length of the box.
 * \param mass The total mass of the box.
 * \return The density of the box.
 */
inline real_t Box::calcDensity( real_t x, real_t y, real_t z, real_t mass )
{
   return mass / ( x * y * z );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Calculates the density of a box for given side lengths and mass.
 *
 * \param l The side lengths of the box.
 * \param mass The total mass of the box.
 * \return The density of the box.
 */
inline real_t Box::calcDensity( const Vec3& l, real_t mass )
{
   return mass / ( l[0] * l[1] * l[2] );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Estimates the point which is farthest in direction \a d.
 *
 * \param d The normalized search direction in world-frame coordinates.
 * \return The support point in world-frame coordinates in direction a\ d.
 */
inline Vec3 Box::support( const Vec3& d ) const
{
   auto len = d.sqrLength();
   if (math::equal(len, real_t(0)))
      return Vec3(0,0,0);

   const Vec3 bfD = vectorFromWFtoBF(d / sqrt(len)); //d in body frame coordinates

   /*as there is no signum funktion in the std-lib I found this
   http://stackoverflow.com/questions/1903954/is-there-a-standard-sign-function-signum-sgn-in-c-c
   template <typename T>
   int sgn(T val) {
      return (T(0) < val) - (val < T(0));
   }
   */

   //As it is save to say we have atleast one component of the d-vector != 0 we can use
   Vec3 relativSupport = Vec3( math::sign(bfD[0])*lengths_[0]*real_t(0.5),
                               math::sign(bfD[1])*lengths_[1]*real_t(0.5),
                               math::sign(bfD[2])*lengths_[2]*real_t(0.5) );

   return gpos_ + vectorFromBFtoWF(relativSupport);
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Estimates the point which is farthest in direction \a d.
 *
 * \param d The normalized search direction in world-frame coordinates.
 * \return The support point in world-frame coordinates in direction a\ d extended by a vector in
 *         direction \a d of length \a pe::contactThreshold.
 */
inline Vec3 Box::supportContactThreshold( const Vec3& d ) const
{
   auto len = d.sqrLength();
   if (math::equal(len, real_t(0)))
      return Vec3(0,0,0);

   return support(d) + d*contactThreshold / sqrt(len);
}
//*************************************************************************************************


//=================================================================================================
//
//  GLOBAL OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\name Box operators */
//@{
std::ostream& operator<<( std::ostream& os, const Box& b );
std::ostream& operator<<( std::ostream& os, ConstBoxID b );
//@}
//*************************************************************************************************

} // namespace pe
} // namespace walberla
