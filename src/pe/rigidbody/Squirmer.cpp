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
//! \file Squirmer.cpp
//! \author Christian Burkard <buch@icp.uni-stuttgart.de>
//! \author Michael Kuron <mkuron@icp.uni-stuttgart.de>
//
//======================================================================================================================

#include "Squirmer.h"

namespace walberla {
namespace pe {

//=================================================================================================
//
//  CONSTRUCTOR
//
//=================================================================================================


//*************************************************************************************************
//*************************************************************************************************
/*!\brief Constructor for the Squirmer class.
 *
 * \param sid Unique system-specific ID for the sphere.
 * \param uid User-specific ID for the sphere.
 * \param gpos Global geometric center of the sphere.
 * \param q The orientation of the sphere's body frame in the global world frame.
 * \param radius The radius of the sphere \f$ (0..\infty) \f$.
 * \param squirmerVelocity The velocity of the squirmer.
 * \param squirmerBeta The dipolar characteristic of the squirmer.
 * \param material The material of the sphere.
 * \param global specifies if the sphere should be created in the global storage
 * \param communicating specifies if the sphere should take part in synchronization (syncNextNeighbour, syncShadowOwner)
 * \param infiniteMass specifies if the sphere has infinite mass and will be treated as an obstacle
 */
Squirmer::Squirmer( id_t sid, id_t uid, const Vec3& gpos, const Quat& q,
                    real_t radius, real_t squirmerVelocity, real_t squirmerBeta, MaterialID material,
                    const bool global, const bool communicating, const bool infiniteMass )
   : Sphere( getStaticTypeID(), sid, uid, gpos, q, radius, material, global, communicating, infiniteMass )  // Initialization of the parent class
   , squirmerVelocity_(squirmerVelocity), squirmerBeta_(squirmerBeta)
{
}

//=================================================================================================
//
//  DESTRUCTOR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Destructor for the Squirmer class.
 */
Squirmer::~Squirmer()
{
   // Logging the destruction of the squirmer
   WALBERLA_LOG_DETAIL( "Destroyed squirmer " << sid_ );
}
//*************************************************************************************************


//=================================================================================================
//
//  GET FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Calculation of the global velocity of a relative point.
 *
 * \param px The x-component of the relative coordinate.
 * \param py The y-component of the relative coordinate.
 * \param pz The z-component of the relative coordinate.
 * \return The global velocity.
 *
 * The function calculates the global velocity of a point relative to the body's center of mass.
 */
Vec3 Squirmer::velFromBF( real_t px, real_t py, real_t pz ) const
{
   return velFromBF( Vec3( px, py, pz ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Calculation of the global velocity of a relative point.
 *
 * \param rpos The relative coordinate.
 * \return The global velocity.
 *
 * The function calculates the global velocity of a point relative to the body's center of mass.
 */
Vec3 Squirmer::velFromBF( const Vec3& rpos ) const
{
   return Sphere::velFromBF( rpos ) + getSquirmerVelocity( getRotation() * rpos );
}
//*************************************************************************************************

//*************************************************************************************************
/*!\brief Calculation of the global velocity of a point in global coordinates.
 *
 * \param px The x-component of the global coordinate.
 * \param py The y-component of the global coordinate.
 * \param pz The z-component of the global coordinate.
 * \return The global velocity.
 *
 * The function calculates the global velocity of a point in global coordinates.
 */
Vec3 Squirmer::velFromWF( real_t px, real_t py, real_t pz ) const
{
   return velFromWF( Vec3( px, py, pz ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the squirmer's surface velocity at a given relative position.
 *
 * \param rpos The relative global coordinate.
 * \return The local surface velocity of the squirmer.
 */
Vec3 Squirmer::getSquirmerVelocity( const Vec3& rpos ) const
{
   const auto rs = rpos.getNormalized();
   const auto e = getQuaternion().rotate(Vec3(0.0,0.0,1.0)).getNormalized();
   const auto v0 = getSquirmerVelocity();
   const auto beta = getSquirmerBeta();
   return ((e * rs) * rs - e) * real_c(1.5) * v0 * (1 + beta * (e * rs));
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Calculation of the global velocity of a point in global coordinates.
 *
 * \param gpos The global coordinate.
 * \return The global velocity.
 *
 * The function calculates the global velocity of a point in global coordinates.
 */
Vec3 Squirmer::velFromWF( const Vec3& gpos ) const
{
   return Sphere::velFromWF( gpos ) + getSquirmerVelocity( gpos - getPosition() );
}
//*************************************************************************************************

id_t Squirmer::staticTypeID_ = std::numeric_limits<id_t>::max();

} // namespace pe
} // namespace walberla
