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
//! \file
//! \author Tobias Scharpff
//! \author Tobias Leemann
//
//======================================================================================================================

#pragma once

#include <mesa_pd/collision_detection/Support.h>
#include <mesa_pd/data/DataTypes.h>

#include <vector>

#include <core/Abort.h>
#include <core/math/Limits.h>
#include <core/math/Vector3.h>

namespace walberla {
namespace mesa_pd {
namespace collision_detection {

//*************************************************************************************************
/*!\brief Impelementation of the Gilbert-Johnson-Keerthi Algorithm.
 */
class GJK
{
public:

   //**Constructor*********************************************************************************
   /*! \name Constructor */
   //@{
   explicit GJK();
   //@}
   //**********************************************************************************************

   //**Query functions*****************************************************************************
   /*! \name Query functions */
   //@{
   real_t doGJK( const Support &geom1, const Support &geom2, Vec3& normal, Vec3& contactPoint );

   bool doGJKmargin( const Support &geom1, const Support &geom2, const real_t margin = real_t(0));
   //@}
   //**********************************************************************************************

   //**Get functions*******************************************************************************
   /*! \name Get functions */
   //@{
   inline const std::vector<Vec3>& getSimplex()     const { return simplex_; }
   inline size_t                   getSimplexSize() const { return numPoints_; }
   inline const std::vector<Vec3>& getSupportA()    const { return supportA_; }
   inline const std::vector<Vec3>& getSupportB()    const { return supportB_; }
   //@}
   //**********************************************************************************************

private:
   //**Utility functions***************************************************************************
   /*! \name Utility functions */
   //@{
   bool simplex2(Vec3& d);
   bool simplex3(Vec3& d);
   bool simplex4(Vec3& d);

   /// Checks if two vectors roughly point in the same direction
   inline bool sameDirection   ( const Vec3& vec1, const Vec3& vec2 ) const { return vec1 * vec2 > real_t(0.0); }
   /// Checks if the length of a vector is zero or as close to zero that it can not be distinguished form zero
   inline bool zeroLengthVector( const Vec3& vec )                    const { return vec.sqrLength() < math::Limits<real_t>::fpuAccuracy(); }
   real_t calcDistance    ( Vec3& normal, Vec3& contactPoint );

   inline Vec3 putSupport(const Support &geom1,
                                const Support &geom2,
                                const Vec3& dir,
                                const real_t margin,
                                std::vector<Vec3> &simplex,
                                std::vector<Vec3> &supportA,
                                std::vector<Vec3> &supportB,
                                size_t index);
   //@}
   //**********************************************************************************************

   //**Member variables****************************************************************************
   /*! \name Member variables */
   //@{
   std::vector<Vec3> simplex_;   ///< Container to hold the simplex.
   std::vector<Vec3> supportA_;  ///< Container to hold the support points generated in triangle mesh mA
   std::vector<Vec3> supportB_;  ///< Container to hold the support points generated in triangle mesh mB
   unsigned char     numPoints_; ///< Current number of points in the simplex.
   Vec3              d_;         ///< The next search direction.
   //@}
   //**********************************************************************************************
};
//*************************************************************************************************

} // namespace collision_detection
} // namespace mesa_pd
} // namespace walberla
