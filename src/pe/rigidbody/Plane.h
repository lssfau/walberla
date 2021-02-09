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
//! \file Plane.h
//! \author Klaus Iglberger
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#pragma once


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include "pe/rigidbody/GeomPrimitive.h"
#include "pe/Types.h"

#include <limits>

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
 * \brief Plane geometry.
 *
 * The Plane class represents the geometric primitive plane, which is one of the basic geometric
 * primitives of the pe module. The class is derived from the GeomPrimitive base
 * class, which makes the plane both a geometric primitive and a rigid body.\n
 * The plane geometry is an infinite rigid body dividing the global space in two half spaces.
 * One of these half spaces is considered to be inside the plane. Bodies entering this half
 * space are therefore colliding with the plane. The other half space is considered to be
 * outside the plane. The plane is represented by the following equation:
 *
 *                                \f[ ax + by + cz = d , \f]
 *
 * where \a a, \a b and \a c are the x, y and z component of the normal vector. The normal
 * \a n of the plane is a normalized vector that points towards the half space outside the
 * plane. \a d is the distance/displacement from the origin of the global world frame to the
 * plane. A positive value of \a d indicates that the origin of the global world frame is
 * inside the plane, whereas a negative value of \a d indicates that the origin is outside
 * the plane. A value of 0 therefore indicates that the origin is on the surface of the plane.
 * The global position \f$ (x,y,z) \f$ of the plane can be considered the anchor point of the
 * plane. Rotations that are performed via the setOrientation() or the rotate() functions
 * rotate the plane around this anchor point.
 */
class Plane : public GeomPrimitive
{
public:
   //**Constructor*********************************************************************************
   /*!\name Constructor */
   //@{
   explicit Plane( id_t sid, id_t uid, const Vec3& gpos, const Vec3& normal,
                   real_t d, MaterialID material );
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   /*!\name Destructor */
   //@{
   ~Plane() override;
   //@}
   //**********************************************************************************************

public:
   //**Get functions*******************************************************************************
   /*!\name Get functions */
   //@{
   inline real_t getVolume()       const override;
   inline const Vec3&    getNormal()       const;
   inline real_t         getDisplacement() const;
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
   void print( std::ostream& os, const char* tab ) const override;
   //@}
   //**********************************************************************************************

protected:
   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   void setPositionImpl       ( real_t px, real_t py, real_t pz ) override;
   void setOrientationImpl    ( real_t r, real_t i, real_t j, real_t k ) override;
   void translateImpl         ( real_t dx, real_t dy, real_t dz ) override;
   void rotateImpl            ( const Quat& dq ) override;
   void rotateAroundOriginImpl( const Quat& dq ) override;
   void rotateAroundPointImpl ( const Vec3& point, const Quat& dq ) override;
   bool containsRelPointImpl   ( real_t px, real_t py, real_t pz ) const override;
   bool isSurfaceRelPointImpl  ( real_t px, real_t py, real_t pz ) const override;
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   void calcBoundingBox() override;  // Calculation of the axis-aligned bounding box
   //@}
   //**********************************************************************************************

   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   Vec3 normal_;  //!< Normal of the plane in reference to the global world frame.
                  /*!< The normal of the plane is always pointing towards the halfspace
                       outside the plane. */
   real_t d_;       //!< Plane displacement from the origin.
                  /*!< The displacement can be categorized in the following way:\n
                        - > 0: The global origin is inside the plane\n
                        - < 0: The global origin is outside the plane\n
                        - = 0: The global origin is on the surface of the plane */
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
/*!\brief Returns the volume of the plane.
 *
 * \return The volume is always infinity.
 */
inline real_t Plane::getVolume() const
{
   return std::numeric_limits<real_t>::infinity();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the normal of the plane in reference to the global world frame.
 *
 * \return The normal of the plane.
 */
inline const Vec3& Plane::getNormal() const
{
   return normal_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the displacement/distance from the origin of the global world frame.
 *
 * \return The displacement of the plane.
 *
 * A positive displacement value indicates that the origin of the global world frame is contained
 * in the plane, whereas a negative value indicates that the origin is not contained in the plane.
 */
inline real_t Plane::getDisplacement() const
{
   return d_;
}
//*************************************************************************************************

//=================================================================================================
//
//  UTILITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Calculates the depth of a point relative to the plane's frame of reference.
 *
 * \param px The x-component of the relative coordinate.
 * \param py The y-component of the relative coordinate.
 * \param pz The z-component of the relative coordinate.
 * \return Depth of the relative point.
 *
 * Returns a positive value, if the point lies inside the plane and a negative value,
 * if the point lies outside the plane.
 */
inline real_t Plane::getRelDepth( real_t /*px*/, real_t /*py*/, real_t pz ) const
{
   return -pz;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Calculates the depth of a point relative to the plane's frame of reference.
 *
 * \param rpos The relative coordinate.
 * \return Depth of the relative point.
 *
 * Returns a positive value, if the point lies inside the plane and a negative value,
 * if the point lies outside the plane.
 */
inline real_t Plane::getRelDepth( const Vec3& rpos ) const
{
   return -rpos[2];
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
 * Returns a positive value, if the point lies inside the plane and a negative value,
 * if the point lies outside the plane.
 */
inline real_t Plane::getDepth( real_t px, real_t py, real_t pz ) const
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
 * Returns a positive value, if the point lies inside the plane and a negative value,
 * if the point lies outside the plane.
 */
inline real_t Plane::getDepth( const Vec3& gpos ) const
{
   return d_ - ( gpos * normal_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Calculates the distance of a point relative to the plane's frame of reference.
 *
 * \param px The x-component of the relative coordinate.
 * \param py The y-component of the relative coordinate.
 * \param pz The z-component of the relative coordinate.
 * \return Distance of the relative point.
 *
 * Returns a positive value, if the point lies outside the plane and a negative value, if the
 * point lies inside the plane.
 */
inline real_t Plane::getRelDistance( real_t px, real_t py, real_t pz ) const
{
   return -getRelDepth( px, py, pz );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Calculates the distance of a point relative to the plane's frame of reference.
 *
 * \param rpos The relative coordinate.
 * \return Distance of the relative point.
 *
 * Returns a positive value, if the point lies outside the plane and a negative value, if the
 * point lies inside the plane.
 */
inline real_t Plane::getRelDistance( const Vec3& rpos ) const
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
 * Returns a positive value, if the point lies outside the plane and a negative value, if the
 * point lies inside the plane.
 */
inline real_t Plane::getDistance( real_t px, real_t py, real_t pz ) const
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
 * Returns a positive value, if the point lies outside the plane and a negative value, if the
 * point lies inside the plane.
 */
inline real_t Plane::getDistance( const Vec3& gpos ) const
{
   return -getRelDepth( pointFromWFtoBF( gpos ) );
}
//*************************************************************************************************

//*************************************************************************************************
/*!\brief Returns unique type id of this type
 *
 * \return geometry specific type id
 */
inline id_t Plane::getStaticTypeID()
{
   return staticTypeID_;
}

//=================================================================================================
//
//  GLOBAL OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\name Plane operators */
//@{
std::ostream& operator<<( std::ostream& os, const Plane& p );
std::ostream& operator<<( std::ostream& os, ConstPlaneID p );
//@}
//*************************************************************************************************


} // namespace pe
} // namespace walberla
