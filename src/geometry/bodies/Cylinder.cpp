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
//! \file Cylinder.cpp
//! \ingroup geometry
//! \author Matthias Markl <matthias.markl@fau.de>
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#include "Cylinder.h"

namespace walberla {
namespace geometry {

   Cylinder::Cylinder( const Vector3<real_t> & _start, const Vector3<real_t> & _end, real_t _rad )
      : start_( _start ), end_( _end ), rad_ ( _rad )
   {
      update();
   }
   
   void Cylinder::update()
   {
      WALBERLA_ASSERT_GREATER( rad_, real_t(0) );

      dir_ = end_ - start_;
      len_ = dir_.length();
      WALBERLA_ASSERT_GREATER( len_, real_t(1e-16) );
      math::normalize(dir_);

      boundingBox_.initMinMaxCorner( std::min( start_[0], end_[0] ) - rad_,
                                     std::min( start_[1], end_[1] ) - rad_,
                                     std::min( start_[2], end_[2] ) - rad_,
                                     std::max( start_[0], end_[0] ) + rad_,
                                     std::max( start_[1], end_[1] ) + rad_,
                                     std::max( start_[2], end_[2] ) + rad_ );
   }


   //===================================================================================================================
   //
   //  Body concept implementation
   //
   //===================================================================================================================
   
   template<>
   FastOverlapResult fastOverlapCheck ( const Cylinder & cyl, const AABB & box )
   {
      if ( ! cyl.boundingBox().intersects( box ) )  return COMPLETELY_OUTSIDE;
      else                                          return DONT_KNOW;
   }

   template<>
   bool contains ( const Cylinder & cyl, const Vector3<real_t> & point )
   {
      const real_t ps = ( point - cyl.start() ) * cyl.dir();

      // before start point or after end point
      if( ps < real_t(0) || ps > cyl.length() )
         return false;

      // outside radius
      if( ( cyl.start() + ps*cyl.dir() - point ).sqrLength() > cyl.radius()*cyl.radius() )
         return false;

      return true;
   }
   
} // namespace geometry
} // namespace walberla
