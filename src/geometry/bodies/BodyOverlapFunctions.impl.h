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
//! \file BodyOverlapFunctions.impl.h
//! \ingroup geometry
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================


namespace walberla {
namespace geometry {

   template <typename Body>
   FastOverlapResult fastOverlapCheck ( const Body & /*body*/, const AABB & /*box*/ )
   {
      // Default implementation has to fastOverlapCheck
      return DONT_KNOW;
   }


   template <typename Body>
   FastOverlapResult fastOverlapCheck ( const Body & /*body*/,
                                        const Vector3<real_t> & /*cellMidpoint*/,
                                        const Vector3<real_t> & /*dx*/ )
   {
      // Default implementation has to fastOverlapCheck
      return DONT_KNOW;
   }


   template< typename Body>
   real_t cellSupersampling( const Vector3<real_t> & cellMidpoint, const Vector3<real_t> & dx, const Body & body, uint_t maxDepth=4, uint_t depth = uint_t(0u) )
   {
      FastOverlapResult r = fastOverlapCheck( body, cellMidpoint, dx );
      if ( r == CONTAINED_INSIDE_BODY )
         return real_t(1);
      else if ( r == COMPLETELY_OUTSIDE )
         return real_t(0);

      uint_t nrCornerPointsInBody(0u);
      for( int signX = -1; signX <= 1; signX += 2 )
         for( int signY = -1; signY <= 1; signY += 2 )
            for( int signZ = -1; signZ <= 1; signZ += 2 )
            {
               // epsilon is subtracted due to symmetry reasons ( i.e. a sphere on a cell boundary should be symmetric)
               const Vector3<real_t> corner( cellMidpoint[0] + real_c(signX) * dx[0] * (real_t(0.5) - real_comparison::Epsilon<real_t>::value ),
                                             cellMidpoint[1] + real_c(signY) * dx[1] * (real_t(0.5) - real_comparison::Epsilon<real_t>::value ),
                                             cellMidpoint[2] + real_c(signZ) * dx[2] * (real_t(0.5) - real_comparison::Epsilon<real_t>::value ) );
               if ( contains( body, corner ) )
                  ++nrCornerPointsInBody;
            }

      if ( nrCornerPointsInBody == uint_t(8u) )
         return real_t(1);
      else if ( nrCornerPointsInBody == uint_t(0u) && !contains( body, cellMidpoint ) )
         return real_t(0);
      else if ( depth == maxDepth )
          return real_c(nrCornerPointsInBody) * real_t(0.125);

      // Recursive calls for 8 sub-cubes
      real_t fraction(0);
      for( int signX = -1; signX <= 1; signX += 2 )
         for( int signY = -1; signY <= 1; signY += 2 )
            for( int signZ = -1; signZ <= 1; signZ += 2 )
            {
               const Vector3<real_t> offsetVec ( real_c(signX) * real_t(0.25) * dx[0], real_c(signY) * real_t(0.25) * dx[1], real_c(signZ) * real_t(0.25) * dx[2] );
               fraction += cellSupersampling( cellMidpoint + offsetVec, dx*real_t(0.5), body, maxDepth, depth+uint_t(1u) );
            }
      fraction *= real_t(0.125);

      return fraction;
   }



   template < typename Body >
   real_t overlapFraction ( const Body & body, const Vector3<real_t> & cellMidpoint, real_t dx, uint_t maxDepth )
   {
      return overlapFraction<Body>(body, cellMidpoint, Vector3<real_t>(dx), maxDepth);
   }

   template < typename Body >
   real_t overlapFraction ( const Body & body, const Vector3<real_t> & cellMidpoint, real_t dx, int maxDepth )
   {
      if ( maxDepth >= 0 )
         return overlapFraction( body, cellMidpoint, dx, uint_c(maxDepth));
      if( contains( body, cellMidpoint ) )
         return real_t(1);
      return real_t(0);
   }

   template < typename Body >
   real_t overlapFraction ( const Body & body, const Vector3<real_t> & cellMidpoint, const Vector3<real_t> & dx, uint_t maxDepth )
   {
      FastOverlapResult r = fastOverlapCheck( body, cellMidpoint, dx );
      if ( r == CONTAINED_INSIDE_BODY )
         return real_t(1);
      else if ( r == COMPLETELY_OUTSIDE )
         return real_t(0);

      // default: fall-back to super-sampling
      real_t overlapFractionBySuperSampling = cellSupersampling( cellMidpoint, dx, body, maxDepth );
      WALBERLA_ASSERT_GREATER_EQUAL(overlapFractionBySuperSampling, real_t(0));
      WALBERLA_ASSERT_LESS_EQUAL(overlapFractionBySuperSampling, real_t(1));
      return overlapFractionBySuperSampling;
   }


} // namespace geometry
} // namespace walberla


