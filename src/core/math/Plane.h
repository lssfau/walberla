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
//! \ingroup core
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/DataTypes.h"
#include "core/math/Vector3.h"

#include <iostream>


namespace walberla {
namespace math {

class Plane
{
public:
   using Vec3Real = Vector3<real_t>;

   inline Plane();
   inline Plane( const Vec3Real & origin, const Vec3Real & _normal );
   inline Plane( const Vec3Real & _normal, const real_t _d );

   inline bool operator==( const Plane & other ) const;
   inline bool operator!=( const Plane & other ) const;

   inline void assign( const Vec3Real & origin, const Vec3Real & _normal );
   inline void assign( const Vec3Real & _normal, const real_t _d );

   inline real_t signedDistance( const Vec3Real & point ) const;
   inline real_t distance( const Vec3Real & point ) const;

   inline bool isInHalfSpace( const Vec3Real & point ) const;

   inline void invert();
   inline void shift( const real_t distance );
   inline void scale( const real_t factor );

   real_t                d() const { return d_; }
   const Vec3Real & normal() const { return normal_; }

   inline friend std::ostream& operator<<( std::ostream& os, const Plane & p );
   inline friend std::istream& operator>>( std::istream& is, Plane & plane );

private:
   Vec3Real normal_;
   real_t   d_;
};

Plane::Plane() :
   normal_( Vec3Real( real_t(1), real_t(0), real_t(0) ) ), d_( real_t(0) )
{
}


Plane::Plane( const Vec3Real & origin, const Vec3Real & _normal ) :
   normal_( _normal.getNormalized() ), d_( normal_ * origin )
{
}


Plane::Plane( const Vec3Real & _normal, const real_t _d ) :
   normal_( _normal.getNormalized() ), d_( _d )
{
}


bool Plane::operator==( const Plane & other ) const
{
   return isIdentical(normal_, other.normal_ ) && isIdentical( d_, other.d_ );
}


bool Plane::operator!=( const Plane & other ) const
{
   return !( *this == other );
}


void Plane::assign( const Vec3Real & origin, const Vec3Real & _normal )
{
   normal_ = _normal.getNormalized();
   d_ = normal_ * origin;
}

void Plane::assign( const Vec3Real & _normal, const real_t _d )
{
   normal_ = _normal.getNormalized();
   d_ = _d;
}


real_t Plane::signedDistance( const Vec3Real & point ) const
{
   return normal_ * point - d_;
}


real_t Plane::distance( const Vec3Real & point ) const
{
   return std::fabs( normal_ * point - d_);
}


bool Plane::isInHalfSpace( const Vec3Real & point ) const
{
   return normal_ * point <= d_;
}


void Plane::invert()
{
   normal_ = -normal_;
   d_      = -d_;
}


void Plane::shift( const real_t dist )
{
   d_ += dist;
}


void Plane::scale( const real_t factor )
{
   d_ *= factor;
}


std::ostream& operator<<( std::ostream& os, const Plane & plane )
{
   os << "(" << plane.normal_ << ", " << plane.d_ << ")";

   return os;
}


std::istream& operator>>( std::istream& is, Plane & plane )
{
   if( !is ) return is;

   char bracket1, bracket2, comma;
   Plane::Vec3Real normal;
   real_t   d;

   const std::istream::pos_type pos( is.tellg() );
   const std::istream::fmtflags oldFlags( is.flags() );

   // Setting the 'skip whitespaces' flag
   is >> std::skipws;

   bool success = false;

   // Extracting the plane in form (<x, y, z>, d)
   if( !(is >> bracket1 >> normal >> comma >> d >> bracket2) ||
      bracket1 != '(' || comma != ',' || bracket2 != ')' )
   {
      is.clear();
      is.seekg( pos );
      is.flags( oldFlags );
   }
   else
   {
      plane.assign( normal, d );
      success = true;
   }

   if( !success )
   {
      Plane::Vec3Real origin;
      if( !(is >> bracket1 >> origin >> comma >> normal >> bracket2) ||
         bracket1 != '(' || comma != ',' || bracket2 != ')' )
      {
         is.clear();
         is.seekg( pos );
         is.flags( oldFlags );
      }
      else
      {
         plane.assign( origin, normal );
         success = true;
      }
   }

   // Resetting the flags
   is.flags( oldFlags );

   if( !success )
      is.setstate( std::istream::failbit );

   return is;
}


} // namespace math
} // namespace walberla
