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
//! \file Torus.cpp
//! \ingroup geometry
//! \author Matthias Markl <matthias.markl@fau.de>
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#include "Torus.h"

namespace walberla {
namespace geometry {

   Torus::Torus( const Vector3<real_t> & _midPnt, const Vector3<real_t> & _normal, real_t _radius, real_t _distance )
      : midPnt_( _midPnt ), normal_( _normal ), radius_( _radius ), distnc_( _distance )
   {
      update();
   }
   
   void Torus::update()
   {
      WALBERLA_ASSERT_GREATER( distnc_, real_t(0) );
      WALBERLA_ASSERT_GREATER( radius_, distnc_   );
      WALBERLA_ASSERT_GREATER( normal_.sqrLength(), real_t(1e-16) );

      math::normalize( normal_ );
      math::normals( normal_, defTan_, comTan_ );

      const real_t sum = radius_ + distnc_;

      boundingBox_.initMinMaxCorner( midPnt_[0] - sum, midPnt_[1] - sum, midPnt_[2] - sum,
                                     midPnt_[0] + sum, midPnt_[1] + sum, midPnt_[2] + sum );
   }


   //===================================================================================================================
   //
   //  Body concept implementation
   //
   //===================================================================================================================
   
   template<>
   FastOverlapResult fastOverlapCheck ( const Torus & torus, const AABB & box )
   {
      if ( ! torus.boundingBox().intersects( box ) )  return COMPLETELY_OUTSIDE;
      else                                            return DONT_KNOW;
   }

   template<>
   bool contains ( const Torus & torus, const Vector3<real_t> & point )
   {
      const auto dir = point - torus.midPnt();
      const auto tan = (dir*torus.defTan())*torus.defTan() + (dir*torus.comTan())*torus.comTan();
      if( tan.sqrLength() < real_t(1e-16) )
         return false;
      const auto pnt = torus.midPnt() + tan.getNormalized()*torus.radius();
      const auto len = ( pnt - point ).sqrLength();
      return len < torus.distance()*torus.distance();
   }
   
} // namespace geometry
} // namespace walberla
