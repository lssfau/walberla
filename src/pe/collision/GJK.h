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
//! \file GJK.h
//! \author Tobias Scharpff
//! \author Tobias Leemann
//
//======================================================================================================================

#pragma once

//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <vector>

#include <pe/rigidbody/GeomPrimitive.h>
#include <pe/Thresholds.h>

#include <core/Abort.h>
#include <core/math/Limits.h>
#include <core/math/Vector3.h>

namespace walberla {
namespace pe {
namespace fcd {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Implementation of the Gilbert-Johnson-Keerthi Algorithm.
 */
class GJK
{
public:

   //**Constructor*********************************************************************************
   /*! \name Constructor */
   //@{
   explicit inline GJK();
   //@}
   //**********************************************************************************************

   //**Query functions*****************************************************************************
   /*! \name Query functions */
   //@{
   real_t doGJK( GeomPrimitive &geom1, GeomPrimitive &geom2, Vec3& normal, Vec3& contactPoint );

   bool doGJKmargin( GeomPrimitive &geom1, GeomPrimitive &geom2, const real_t margin = contactThreshold);
   //@}
   //**********************************************************************************************

   //**Get functions*******************************************************************************
   /*! \name Get functions */
   //@{
   inline const std::vector<Vec3>& getSimplex()     const;
   inline size_t                   getSimplexSize() const;
   inline const std::vector<Vec3>& getSupportA()    const;
   inline const std::vector<Vec3>& getSupportB()    const;
   //@}
   //**********************************************************************************************

private:
   //**Utility functions***************************************************************************
   /*! \name Utility functions */
   //@{
   bool simplex2(Vec3& d);
   bool simplex3(Vec3& d);
   bool simplex4(Vec3& d);

   inline bool sameDirection   ( const Vec3& vec1, const Vec3& vec2 ) const;
   inline bool zeroLengthVector( const Vec3& vec )                     const;
   real_t calcDistance    ( Vec3& normal, Vec3& contactPoint );
   inline const Vec3 putSupport(const GeomPrimitive &geom1, const GeomPrimitive &geom2, const Vec3& dir, const real_t margin,
                                std::vector<Vec3> &simplex, std::vector<Vec3> &supportA, std::vector<Vec3> &supportB, size_t index);
   //@}
   //**********************************************************************************************

   //**Member variables****************************************************************************
   /*! \name Member variables */
   //@{
   std::vector<Vec3> simplex_;   //<! Container to hold the simplex.
   std::vector<Vec3> supportA_;  //<! Container to hold the support points generated in triangle mesh mA
   std::vector<Vec3> supportB_;  //<! Container to hold the support points generated in triangle mesh mB
   unsigned char     numPoints_; //<! Current number of points in the simplex.
   Vec3              d_;         //<! The next search direction.
   //@}
   //**********************************************************************************************
};
//*************************************************************************************************

//=================================================================================================
//
//  CONSTRUCTOR
//
//=================================================================================================

//*************************************************************************************************
inline GJK::GJK() : simplex_(4), supportA_(4), supportB_(4), numPoints_(0)
{
   d_ = Vec3(real_t(0.0),real_t(0.6),real_t(0.8)); // just start with any vector of length 1
}
//*************************************************************************************************


//=================================================================================================
//
//  GET FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
inline const std::vector<Vec3>& GJK::getSimplex() const
{
   return simplex_;
}
//*************************************************************************************************


//*************************************************************************************************
inline size_t GJK::getSimplexSize() const
{
   return numPoints_;
}
//*************************************************************************************************


//*************************************************************************************************
inline const std::vector<Vec3>& GJK::getSupportA() const
{
   return supportA_;
}
//*************************************************************************************************


//*************************************************************************************************
inline const std::vector<Vec3>& GJK::getSupportB() const
{
   return supportB_;
}
//*************************************************************************************************


//=================================================================================================
//
//  UTILITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Checks if two vectors roughly point in the same directionTODO
 */
inline bool GJK::sameDirection(const Vec3& vec1, const Vec3& vec2) const
{
   return vec1 * vec2 > real_t(0.0);
}
//*************************************************************************************************


//*************************************************************************************************
/* Checks if the length of a vector is zero or as close to zero that it can not be distinguished form zero
 */
inline bool GJK::zeroLengthVector(const Vec3& vec) const
{
   return vec.sqrLength() < math::Limits<real_t>::fpuAccuracy();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Calculate a support point of a body extended by a threshold.
 * \param geom1 The body A.
 * \param geom2 The body B.
 * \param dir The support point direction.
 * \param margin Extension of the Body.
 */
inline const Vec3 GJK::putSupport(const GeomPrimitive &geom1, const GeomPrimitive &geom2, const Vec3& dir, const real_t margin, 
                                  std::vector<Vec3> &simplex, std::vector<Vec3> &supportA, std::vector<Vec3> &supportB, size_t index){
   supportA[index] = geom1.support(dir);
   supportB[index] = geom2.support(-dir);
   Vec3 supp = supportA[index]- supportB[index] + (real_t(2.0) * dir * margin);
   simplex[index] = supp;
   return supp;
}
//*************************************************************************************************


} // namespace fcd

} // namespace pe

} // namespace walberla
