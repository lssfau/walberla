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
//! \file QPX.h
//! \ingroup phasefield
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#pragma once


namespace walberla {
namespace simd {
namespace qpx {

struct double4_t
{
   double4_t() {}
   double4_t( vector4double init ) : v(init) {}
   
   operator vector4double() const { return v; }  
   double   operator[]( int i ) const { return v[i]; }
   double & operator[]( int i )       { return v[i]; }
   
   vector4double v;
};
//typedef vector4double double4_t;

inline const char * usedInstructionSet() { return "QPX"; }

inline double4_t make_double4   ( double d, double c, double b, double a ) { return double4_t( (vector4double) {a,b,c,d} ); }
inline double4_t make_double4_r ( double a, double b, double c, double d ) { return double4_t( (vector4double) {a,b,c,d} ); }

inline double4_t make_double4   ( double a                               ) { return vec_splats(a); }

inline double4_t make_zero()   { return vec_splats(0.0); }

inline double4_t load_aligned  ( const double * mem_addr )          { return vec_ld(0ul, const_cast<double*>(mem_addr) ); }
inline void      store_aligned ( double * mem_addr, double4_t a )   { vec_st( a, 0ul, mem_addr);                          }

inline double4_t load_unaligned  ( const double * mem_addr )
{
   double * p =  const_cast<double*>( mem_addr );
   vector4double v1,v2,xv,pctl;
   v1 = vec_ld ( 0ul, p   );
   v2 = vec_ld ( 0ul, p+4 );

   pctl = vec_lvsl( 0, p );
   xv   = vec_perm( v1, v2, pctl );
   return xv;
}

inline void loadNeighbors( const double * mem_addr, double4_t & r_left, double4_t & r_center, double4_t & r_right )
{
    double * p =  const_cast<double*>( mem_addr );

    vector4double left   = vec_ld(0, const_cast<double*>( p-4 ) );
    vector4double center = vec_ld(0, const_cast<double*>( p   ) );
    vector4double right  = vec_ld(0, const_cast<double*>( p+4 ) );

    vector4double ctrlLeft  = vec_lvsl( 0, p-1 );
    vector4double ctrlRight = vec_lvsl( 0, p+1 );
    r_left   = vec_perm( left, center, ctrlLeft  );
    r_center = center;
    r_right  = vec_perm( center,right, ctrlRight );
}

inline double getComponent ( const double4_t & v, int i )           { return v[i]; }
inline double getComponent ( const double4_t & v, unsigned long i ) { return v[i]; }

inline double4_t hadd( double4_t a,  double4_t b )
{
   //TODO is there an instruction for this?
   return make_double4_r( a[0] + a[1],
                          b[0] + b[1],
                          a[2] + a[3],
                          b[2] + b[3] );
}

inline double4_t horizontalSum ( double4_t a )
{
   //TODO is there an instruction for this?
   return vec_splats( a[0] + a[1] + a[2] + a[3] );
}

inline double4_t exchangeLowerUpperHalf ( double4_t a )
{
   static const vector4double permCode = vec_gpci( 0x4c1 ); // hex for octal 2301
   return vec_perm( a, a, permCode );
}


inline void extract( double4_t in, double4_t &d, double4_t &c, double4_t &b, double4_t & a )
{
   a = vec_splat( in, 0) ;
   b = vec_splat( in, 1 );
   c = vec_splat( in, 2 );
   d = vec_splat( in, 3 );
}

inline double4_t rotateRight( double4_t a )
{
   static const vector4double permCode = vec_gpci( 0x298 ); // hex for octal 1230
   return vec_perm( a, a, permCode );
}

inline double4_t rotateLeft( double4_t a )
{
   static const vector4double permCode = vec_gpci( 0x60a ); // hex for octal 3012
   return vec_perm( a, a, permCode );
}


inline double4_t compareEQ( double4_t a, double4_t b ) {
   return vec_cmpeq( a,b );
}
inline double4_t compareNEQ( double4_t a, double4_t b ) {
   return vec_not( vec_cmpeq(a,b) );
}
inline double4_t compareGE( double4_t a, double4_t b ) {
   return vec_not( vec_cmplt(a,b) );
}
inline double4_t compareLE( double4_t a, double4_t b ) {
   return vec_not( vec_cmpgt(a,b) );
}

inline double4_t logicalAND( double4_t a, double4_t b ) {
   return vec_and( a,b );
}
inline double4_t logicalOR( double4_t a, double4_t b ) {
   return vec_or( a,b );
}



inline int movemask( double4_t a )
{
   int result = 0;
   if ( a[3] == 1.0 ) result |= ( 1 << 3);
   if ( a[2] == 1.0 ) result |= ( 1 << 2);
   if ( a[1] == 1.0 ) result |= ( 1 << 1);
   if ( a[0] == 1.0 ) result |= ( 1 << 0);
   return result;
}

inline double4_t blendv( double4_t a, double4_t b, double4_t mask)
{
   return vec_sel( a,b, mask );
}

template<int mask>
inline double4_t blend( double4_t a, double4_t b )
{
   const double m0 = mask & ( 1<<0 ) ? 1.0 : -1.0;
   const double m1 = mask & ( 1<<1 ) ? 1.0 : -1.0;
   const double m2 = mask & ( 1<<2 ) ? 1.0 : -1.0;
   const double m3 = mask & ( 1<<3 ) ? 1.0 : -1.0;

   const vector4double vec_mask = { m0,m1,m2,m3};
   return vec_sel( a,b,vec_mask);
}


inline double4_t sqrt( double4_t a) {
   return vec_swsqrt_nochk (a );
}

template< unsigned int numIter>
inline double4_t invSqrt( double4_t y )
{
   return vec_rsqrte( y);
}


inline double4_t operator+( const double4_t & a, const double4_t & b ) { return vec_add ( a, b); }
inline double4_t operator-( const double4_t & a, const double4_t & b ) { return vec_sub ( a, b); }
inline double4_t operator*( const double4_t & a, const double4_t & b ) { return vec_mul  ( a, b); }
inline double4_t operator/( const double4_t & a, const double4_t & b ) { return vec_swdivs_nochk ( a, b); }



} // namespace qpx
} // namespace simd
} // namespace walberla



