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
//! \file Ellipsoid.h
//! \author Klaus Iglberger
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//! \author Tobias Leemann <tobias.leemann@fau.de>
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
 * \brief Base class for the Ellipsoid geometry.
 *
 * The Ellipsoid class represents the base class for the Ellipsoid geometry. It provides
 * the basic functionality of a Ellipsoid. An Ellipsoid is obtained from a sphere by deforming it by means
 * of directional scalings. Its three semi-axes corresponding to the x, y, z direction are labeled
 * a, b, c.
 */
class Ellipsoid : public GeomPrimitive
{
public:
   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit Ellipsoid( id_t sid, id_t uid, const Vec3& gpos, const Quat& q,
                       const Vec3& semiAxes, MaterialID material,
                       const bool global, const bool communicating, const bool infiniteMass );
   explicit Ellipsoid( id_t const typeID, id_t sid, id_t uid, const Vec3& rpos, const Quat& q,
                       const Vec3& semiAxes, MaterialID material,
                       const bool global, const bool communicating, const bool infiniteMass );
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   /*!\name Destructor */
   //@{
   ~Ellipsoid() override;
   //@}
   //**********************************************************************************************
   //**********************************************************************************************

public:
   //**Get functions*******************************************************************************
   /*!\name Get functions */
   //@{
   inline const Vec3& getSemiAxes() const;
   inline real_t getVolume()         const override;
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   static inline id_t getStaticTypeID();
   //@}
   //**********************************************************************************************

   //**Output functions****************************************************************************
   /*!\name Output functions */
   //@{
   void print( std::ostream& os, const char* tab ) const override;
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
   static inline real_t calcVolume( const Vec3& semiAxes );
   static inline real_t calcMass( const Vec3& semiAxes, real_t density );
   static inline real_t calcDensity( const Vec3& semiAxes, real_t mass );
   static inline Mat3   calcInertia( const real_t mass, const Vec3& semiAxes );
   //@}
   //**********************************************************************************************

protected:
   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   bool containsRelPointImpl ( real_t px, real_t py, real_t pz ) const override;
   bool isSurfaceRelPointImpl( real_t px, real_t py, real_t pz ) const override;
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   inline void calcBoundingBox() override;  // Calculation of the axis-aligned bounding box
   //@}
   //**********************************************************************************************

   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   Vec3 semiAxes_;  //!< Semi-axes of the Ellipsoid.
   /*!< The radius is constrained to values larger than 0.0. */
   //@}
   //**********************************************************************************************
private:
   static id_t staticTypeID_;  //< type id of Ellipsoid, will be set by SetBodyTypeIDs
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
/*!\brief Returns the semi-axes of the Ellipsoid.
 *
 * \return The semi-axes of the Ellipsoid.
 */
inline const Vec3& Ellipsoid::getSemiAxes() const
{
   return semiAxes_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the volume of the Ellipsoid.
 *
 * \return The volume of the Ellipsoid.
 */
inline real_t Ellipsoid::getVolume() const
{
   return Ellipsoid::calcVolume(getSemiAxes());
}
//*************************************************************************************************




//=================================================================================================
//
//  VOLUME, MASS AND DENSITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Calculates the volume of a Ellipsoid for a given vector of semi-axes.
 *
 * \param semiAxes The vector of semi-axes of the Ellipsoid.
 * \return The volume of the Ellipsoid.
 */
inline real_t Ellipsoid::calcVolume(const Vec3& semiAxes ) 
{
   return real_c(4.0)/real_c(3.0) * math::pi * semiAxes[0] * semiAxes[1] * semiAxes[2];
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Calculates the mass of a Ellipsoid for a given a given vector of semi-axes and density.
 *
 * \param semiAxes The vector of semi-axes of the Ellipsoid.
 * \param density The density of the Ellipsoid.
 * \return The total mass of the Ellipsoid.
 */
inline real_t Ellipsoid::calcMass(const Vec3& semiAxes, real_t density )
{
   return real_c(4.0)/real_c(3.0) * math::pi * semiAxes[0] * semiAxes[1] * semiAxes[2] * density;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Calculates the density of a Ellipsoid for a given vector of semi-axes and mass.
 *
 * \param semiAxes The vector of semi-axes of the Ellipsoid.
 * \param mass The total mass of the Ellipsoid.
 * \return The density of the Ellipsoid.
 */
inline real_t Ellipsoid::calcDensity(const Vec3& semiAxes, real_t mass )
{
   return real_c(0.75) * mass / ( math::pi * semiAxes[0] * semiAxes[1] * semiAxes[2] );
}
//*************************************************************************************************




//=================================================================================================
//
//  UTILITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Calculation of the bounding box of the Ellipsoid.
 *
 * \return void
 *
 * This function updates the axis-aligned bounding box of the Ellipsoid primitive according to the
 * current position and orientation of the Ellipsoid. Note that the bounding box is increased in
 * all dimensions by pe::contactThreshold to guarantee that rigid bodies in close proximity of
 * the Ellipsoid are also considered during the collision detection process.
 * Algorithm: Use a non-axes-aligned bounding box of the ellipse (box
 * with sides twice the semi-axes long) and calc its AABB.
 */
inline void Ellipsoid::calcBoundingBox()
{	
   using std::fabs;
   Mat3 R(getRotation());
   Vec3 gpos = getPosition();
   const real_t xlength( fabs(R[0]*semiAxes_[0]) + fabs(R[1]*semiAxes_[1]) + fabs(R[2]*semiAxes_[2])  + contactThreshold );
   const real_t ylength( fabs(R[3]*semiAxes_[0]) + fabs(R[4]*semiAxes_[1]) + fabs(R[5]*semiAxes_[2])  + contactThreshold );
   const real_t zlength( fabs(R[6]*semiAxes_[0]) + fabs(R[7]*semiAxes_[1]) + fabs(R[8]*semiAxes_[2])  + contactThreshold );
   aabb_ = math::AABB(
           gpos[0] - xlength,
           gpos[1] - ylength,
           gpos[2] - zlength,
           gpos[0] + xlength,
           gpos[1] + ylength,
           gpos[2] + zlength
         );
   //   WALBERLA_ASSERT( aabb_.isValid()        , "Invalid bounding box detected" );
   WALBERLA_ASSERT( aabb_.contains( gpos ), "Invalid bounding box detected("<< getSystemID() <<")\n" << "pos: " << getPosition() << "\nlengths: " << xlength << "," << ylength << ", " << zlength<< "\nvel: " << getLinearVel() << "\nbox: " << aabb_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Calculation of the moment of inertia in reference to the body frame of the Ellipsoid.
 *
 * \return void
 */
inline Mat3 Ellipsoid::calcInertia( const real_t mass, const Vec3& semiAxes )
{
   return Mat3::makeDiagonalMatrix(
            real_c(0.2) * mass * (semiAxes[1] * semiAxes[1] + semiAxes[2] * semiAxes[2]),
            real_c(0.2) * mass * (semiAxes[2] * semiAxes[2] + semiAxes[0] * semiAxes[0]),
            real_c(0.2) * mass * (semiAxes[0] * semiAxes[0] + semiAxes[1] * semiAxes[1]));
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Estimates the point which is farthest in direction \a d.
 *
 * \param d The normalized search direction in world-frame coordinates.
 * \return The support point in world-frame coordinates in direction \a d.
 */
inline Vec3 Ellipsoid::support( const Vec3& d ) const
{
   auto len = d.sqrLength();
   if (!math::equal(len, real_t(0)))
   {
      Vec3 d_loc = vectorFromWFtoBF(d);
      Vec3 norm_vec(d_loc[0] * semiAxes_[0], d_loc[1] * semiAxes_[1], d_loc[2] * semiAxes_[2]);
      real_t norm_length = norm_vec.length();
      Vec3 local_support = (real_t(1.0) / norm_length) * Vec3(semiAxes_[0] * semiAxes_[0] * d_loc[0],
            semiAxes_[1] * semiAxes_[1] * d_loc[1], semiAxes_[2] * semiAxes_[2] * d_loc[2]);
      return pointFromBFtoWF(local_support);
   } else
   {
      return Vec3(0,0,0);
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns unique type id of this type
 *
 * \return geometry specific type id
 */
inline id_t Ellipsoid::getStaticTypeID()
{
   return staticTypeID_;
}
//*************************************************************************************************


//=================================================================================================
//
//  GLOBAL OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\name Ellipsoid operators */
//@{
std::ostream& operator<<( std::ostream& os, const Ellipsoid& s );
std::ostream& operator<<( std::ostream& os, ConstEllipsoidID s );
//@}
//*************************************************************************************************


} // namespace pe
}  // namespace walberla
