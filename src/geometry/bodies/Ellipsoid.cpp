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
//! \file Ellipsoid.cpp
//! \ingroup geometry
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#include "Ellipsoid.h"


namespace walberla {
namespace geometry {

   Ellipsoid::Ellipsoid( const Vector3<real_t> & midp,
            Vector3<real_t> axis1,
            Vector3<real_t> axis2,
            const Vector3<real_t>& radii )
      : midpoint_( midp ), radii_( radii )
   {
      normalize( axis1 );
      axis2 = axis2 - ( axis1 * axis2) * axis1;
      normalize( axis2 );
      WALBERLA_ASSERT_FLOAT_EQUAL( axis1*axis2, 0.0);

      Vector3<real_t> axis3 = axis1 % axis2;

      for( uint_t i =0; i< 3; ++i )
      {
         rotationMatrix_(i,0) = axis1[i];
         rotationMatrix_(i,1) = axis2[i];
         rotationMatrix_(i,2) = axis3[i];
      }

      Matrix3<real_t> diagonalMatrix ( 0.0 );
      diagonalMatrix(0,0) = real_t(1) / ( radii_[0] * radii_[0] );
      diagonalMatrix(1,1) = real_t(1) / ( radii_[1] * radii_[1] );
      diagonalMatrix(2,2) = real_t(1) / ( radii_[2] * radii_[2] );


      mat_ = rotationMatrix_ *  diagonalMatrix * rotationMatrix_.getTranspose() ;

      minRadius_ = std::min( radii_[0], std::min( radii_[1], radii_[2] ) );
      maxRadius_ = std::max( radii_[0], std::max( radii_[1], radii_[2] ) );

      updateBoundingBox( );
   }




   bool Ellipsoid::contains( const Vector3<real_t> & point ) const
   {
      const Vector3<real_t> p = point - midpoint_;
      return (  p *  (mat_ * p) ) < 1;
   }


   void Ellipsoid::updateBoundingBox()
   {
      Matrix3<real_t> diagonalMatrix ( 0.0 );
      diagonalMatrix(0,0) = radii_[0];
      diagonalMatrix(1,1) = radii_[1];
      diagonalMatrix(2,2) = radii_[2];
      Matrix3<real_t> M = rotationMatrix_.rotate( diagonalMatrix );
      real_t xMax = (M * ( M.multTranspose( Vector3<real_t>(1, 0, 0 ) ).getNormalized() ) )[0];
      real_t yMax = (M * ( M.multTranspose( Vector3<real_t>(0, 1, 0 ) ).getNormalized() ) )[1];
      real_t zMax = (M * ( M.multTranspose( Vector3<real_t>(0, 0, 1 ) ).getNormalized() ) )[2];
      boundingBox_ = AABB::createFromMinMaxCorner( midpoint_[0] - xMax, midpoint_[1] - yMax, midpoint_[2] - zMax,
                                                   midpoint_[0] + xMax, midpoint_[1] + yMax, midpoint_[2] + zMax  );
   }


   //===================================================================================================================
   //
   //  Body concept implementation
   //
   //===================================================================================================================


   template<>
   FastOverlapResult fastOverlapCheck ( const Ellipsoid & ellipsoid, const AABB & box )
   {
      if ( ! ellipsoid.boundingBox().intersects( box ) ) return COMPLETELY_OUTSIDE;
      else                                               return DONT_KNOW;
   }

   template<>
   FastOverlapResult fastOverlapCheck ( const Ellipsoid & ellipsoid, const Vector3<real_t> & cellMidpoint, const Vector3<real_t> & dx )
   {
      AABB box = AABB::createFromMinMaxCorner( cellMidpoint[0] - real_t(0.5)*dx[0], cellMidpoint[1] - real_t(0.5)*dx[1], cellMidpoint[2] - real_t(0.5)*dx[2],
                                               cellMidpoint[0] + real_t(0.5)*dx[0], cellMidpoint[1] + real_t(0.5)*dx[1], cellMidpoint[2] + real_t(0.5)*dx[2]);

      if ( ! ellipsoid.boundingBox().intersects( box ) )
         return COMPLETELY_OUTSIDE;


      // Check against inner circle
      static const real_t sqrt3half = std::sqrt( real_t(3) ) / real_t(2);

      const real_t midPointDistSq = (ellipsoid.midpoint() - cellMidpoint).sqrLength();

      // Check against inner circle of box
      const real_t dxMax = dx.max();
      const real_t dist2 = ellipsoid.minRadius() - sqrt3half * dxMax;
      if ( midPointDistSq < dist2 * dist2 )
         return CONTAINED_INSIDE_BODY;

      return DONT_KNOW;
   }

   template<>
   bool contains ( const Ellipsoid & ellipsoid, const Vector3<real_t> & point )
   {
      return ellipsoid.contains( point );
   }


} // namespace geometry
} // namespace walberla


