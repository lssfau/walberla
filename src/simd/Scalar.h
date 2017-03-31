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
//! \file Scalar.h
//! \ingroup simd
//! \author Martin Bauer <martin.bauer@fau.de>
//! \brief Fallback to scalar operations
//
//======================================================================================================================

#pragma once

#include "waLBerlaDefinitions.h"

#include "core/math/FastInvSqrt.h"
#include <cmath>
#include <limits>

#include "core/logging/Logging.h"

namespace walberla {
namespace simd {
namespace scalar {

const double   TRUE_MASK  = -std::numeric_limits<double>::signaling_NaN();
const double   FALSE_MASK = static_cast<double>( 0x0ul );
const uint64_t MSB        = 0x8000000000000000ul;


class double4_t {
public:
   double4_t() {}

   double4_t( double v3, double v2, double v1, double v0 ) {
      data[3] = v3;
      data[2] = v2;
      data[1] = v1;
      data[0] = v0;
   }

   double  operator[]( unsigned int i ) const { return data[i]; }
   double& operator[]( unsigned int i )       { return data[i]; }


   uint64_t  asUInt( unsigned int i ) const { return uintData[i]; }
   uint64_t& asUInt( unsigned int i )       { return uintData[i]; }

   double4_t operator+( const double4_t & o ) const { double4_t r;  for (unsigned int i=0; i<4; ++i) r[i] = (*this)[i] + o[i]; return r; }
   double4_t operator-( const double4_t & o ) const { double4_t r;  for (unsigned int i=0; i<4; ++i) r[i] = (*this)[i] - o[i]; return r; }
   double4_t operator*( const double4_t & o ) const { double4_t r;  for (unsigned int i=0; i<4; ++i) r[i] = (*this)[i] * o[i]; return r; }
   double4_t operator/( const double4_t & o ) const { double4_t r;  for (unsigned int i=0; i<4; ++i) r[i] = (*this)[i] / o[i]; return r; }

private:
   union {
      double   data[4];
      uint64_t uintData[4];
   };
};

inline const char * usedInstructionSet() { return "scalar emulation"; }

inline double4_t make_double4   ( double d, double c, double b, double a ) { return double4_t ( d,c,b,a ); }
inline double4_t make_double4_r ( double a, double b, double c, double d ) { return double4_t ( d,c,b,a ); }

inline double4_t make_double4  ( double a                               )  { return double4_t ( a,a,a,a ); }

inline double4_t make_zero()   { return make_double4(0.0); }

inline double4_t load_aligned    ( double const * m )         { return make_double4_r( m[0], m[1],m[2],m[3] ); }
inline double4_t load_unaligned  ( double const * m )         { return make_double4_r( m[0], m[1],m[2],m[3] ); }
inline void      store_aligned ( double * m, double4_t a )  { m[0]=a[0]; m[1]=a[1]; m[2]=a[2]; m[3]=a[3];    }

inline void loadNeighbors( const double * p, double4_t & r_left, double4_t & r_center, double4_t & r_right )
{
   r_left   = load_unaligned( p-1 );
   r_center = load_aligned  ( p   );
   r_right  = load_unaligned( p+1 );
}


inline double getComponent ( const double4_t & v, int i )           { return v[(unsigned int)(i)]; }
inline double getComponent ( const double4_t & v, unsigned long i ) { return v[(unsigned int)(i)]; }

inline bool   getBoolComponent ( const double4_t & v, int i           ) { return (v.asUInt((unsigned int)(i))) != 0; }
inline bool   getBoolComponent ( const double4_t & v, unsigned long i ) { return (v.asUInt((unsigned int)(i))) != 0; }



inline double4_t hadd( double4_t a,  double4_t b )
{
   double4_t r;
   r[0] = a[1] + a[0];
   r[1] = b[1] + b[0];
   r[2] = a[3] + a[2];
   r[3] = b[3] + b[2];
   return r;
}


inline double4_t horizontalSum ( double4_t a ){
   return make_double4( a[0]+a[1]+a[2]+a[3] );
}

inline double4_t exchangeLowerUpperHalf ( double4_t a )
{
   return make_double4( a[1], a[0], a[3], a[2] );
}

inline void extract( double4_t in, double4_t &d, double4_t &c, double4_t &b, double4_t & a )
{
   a = make_double4( in[0] );
   b = make_double4( in[1] );
   c = make_double4( in[2] );
   d = make_double4( in[3] );
}

inline double4_t rotateRight( double4_t a )
{
   return make_double4( a[0], a[3], a[2], a[1] );
}

inline double4_t rotateLeft( double4_t a )
{
   return make_double4( a[2], a[1], a[0], a[3] );
}

#ifdef WALBERLA_CXX_COMPILER_IS_GNU
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-equal"
#endif


inline double4_t compareEQ( double4_t a, double4_t b )
{
   return make_double4( a[3] == b[3] ? TRUE_MASK : FALSE_MASK,
                        a[2] == b[2] ? TRUE_MASK : FALSE_MASK,
                        a[1] == b[1] ? TRUE_MASK : FALSE_MASK,
                        a[0] == b[0] ? TRUE_MASK : FALSE_MASK
                      );
}
inline double4_t compareNEQ( double4_t a, double4_t b )
{
   return make_double4( (a[3] != b[3]) ? TRUE_MASK : FALSE_MASK,
                        (a[2] != b[2]) ? TRUE_MASK : FALSE_MASK,
                        (a[1] != b[1]) ? TRUE_MASK : FALSE_MASK,
                        (a[0] != b[0]) ? TRUE_MASK : FALSE_MASK
                      );
}

#ifdef WALBERLA_CXX_COMPILER_IS_GNU
#pragma GCC diagnostic pop
#endif


inline double4_t compareGE( double4_t a, double4_t b )
{
   return make_double4( (a[3] >= b[3]) ? TRUE_MASK : FALSE_MASK,
                        (a[2] >= b[2]) ? TRUE_MASK : FALSE_MASK,
                        (a[1] >= b[1]) ? TRUE_MASK : FALSE_MASK,
                        (a[0] >= b[0]) ? TRUE_MASK : FALSE_MASK
                      );
}

inline double4_t compareLE( double4_t a, double4_t b )
{
   return make_double4( (a[3] <= b[3]) ? TRUE_MASK : FALSE_MASK,
                        (a[2] <= b[2]) ? TRUE_MASK : FALSE_MASK,
                        (a[1] <= b[1]) ? TRUE_MASK : FALSE_MASK,
                        (a[0] <= b[0]) ? TRUE_MASK : FALSE_MASK
                      );
}

inline double4_t logicalAND( double4_t a, double4_t b )
{
   double4_t result;
   result.asUInt(3) = a.asUInt(3) & b.asUInt(3);
   result.asUInt(2) = a.asUInt(2) & b.asUInt(2);
   result.asUInt(1) = a.asUInt(1) & b.asUInt(1);
   result.asUInt(0) = a.asUInt(0) & b.asUInt(0);
   return result;
}
inline double4_t logicalOR( double4_t a, double4_t b ) {
   double4_t result;
   result.asUInt(3) = a.asUInt(3) | b.asUInt(3);
   result.asUInt(2) = a.asUInt(2) | b.asUInt(2);
   result.asUInt(1) = a.asUInt(1) | b.asUInt(1);
   result.asUInt(0) = a.asUInt(0) | b.asUInt(0);
   return result;
}



inline int movemask( double4_t a )
{
   int result = 0;
   if ( a.asUInt(3) & MSB ) result |= ( 1 << 3);
   if ( a.asUInt(2) & MSB ) result |= ( 1 << 2);
   if ( a.asUInt(1) & MSB ) result |= ( 1 << 1);
   if ( a.asUInt(0) & MSB ) result |= ( 1 << 0);
   return result;
}

inline double4_t blendv( double4_t a, double4_t b, double4_t mask)
{
   return make_double4( (mask.asUInt(3) & MSB ) ? b[3] : a[3],
                        (mask.asUInt(2) & MSB ) ? b[2] : a[2],
                        (mask.asUInt(1) & MSB ) ? b[1] : a[1],
                        (mask.asUInt(0) & MSB ) ? b[0] : a[0]
                      );
}

template<int mask>
inline double4_t blend( double4_t a, double4_t b )
{
   return make_double4( ( mask & (1<<3) ) ? b[3] : a[3],
                        ( mask & (1<<2) ) ? b[2] : a[2],
                        ( mask & (1<<1) ) ? b[1] : a[1],
                        ( mask & (1<<0) ) ? b[0] : a[0]
                      );
}


inline double4_t sqrt( double4_t a)
{
   return make_double4( std::sqrt( a[3] ), std::sqrt( a[2] ),  std::sqrt( a[1] ),  std::sqrt( a[0] ) );
}

template< unsigned int numIter >
inline double4_t invSqrt( double4_t a )
{
   return make_double4( math::fastInvSqrt<numIter>( a[3] ),
                        math::fastInvSqrt<numIter>( a[2] ),
                        math::fastInvSqrt<numIter>( a[1] ),
                        math::fastInvSqrt<numIter>( a[0] ) );
}



} // namespace scalar
} // namespace simd
} // namespace walberla


