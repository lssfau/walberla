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
//! \file AVX.h
//! \ingroup avx
//! \author Martin Bauer <martin.bauer@fau.de>
//! \brief Wrapper functions around AVX intrinsics
//
//======================================================================================================================

#pragma once

#include "emmintrin.h" //SSE2
#include "waLBerlaDefinitions.h"

#include "core/DataTypes.h"

#ifdef __SSE2__

namespace walberla {
namespace simd {
namespace sse2 {

   struct double4_t
   {
      __m128d low;
      __m128d high;

      inline double4_t operator+( const double4_t & o ) const {
         double4_t res;
         res.low  = _mm_add_pd( this->low,  o.low );
         res.high = _mm_add_pd( this->high, o.high );
         return res;
      }
      inline double4_t operator-( const double4_t & o ) const {
         double4_t res;
         res.low  = _mm_sub_pd( this->low,  o.low );
         res.high = _mm_sub_pd( this->high, o.high );
         return res;
      }
      inline double4_t operator*( const double4_t & o ) const {
         double4_t res;
         res.low  = _mm_mul_pd( this->low,  o.low );
         res.high = _mm_mul_pd( this->high, o.high );
         return res;
      }
      inline double4_t operator/( const double4_t & o ) const {
         double4_t res;
         res.low  = _mm_div_pd( this->low,  o.low );
         res.high = _mm_div_pd( this->high, o.high );
         return res;
      }
   };

   inline const char * usedInstructionSet() { return "SSE2"; }


   inline double4_t make_double4   ( double d, double c, double b, double a ) {
      double4_t res;
      res.low  = _mm_set_pd( b,a );
      res.high = _mm_set_pd( d,c );
      return res;
   }
   inline double4_t make_double4_r ( double a, double b, double c, double d )  {
      double4_t res;
      res.low  = _mm_setr_pd( a,b );
      res.high = _mm_setr_pd( c,d );
      return res;
   }

   inline double4_t make_double4  ( double a  ) {
      return make_double4 ( a,a,a,a );
   }

   inline double4_t make_zero()   {
      double4_t res;
      res.low  = _mm_setzero_pd();
      res.high = _mm_setzero_pd();
      return res;
   }

   inline double4_t load_aligned  ( double const * mem_addr ){
      double4_t res;
      res.low  = _mm_load_pd( mem_addr   );
      res.high = _mm_load_pd( mem_addr+2 );
      return res;
   }

   inline double4_t load_unaligned  ( double const * mem_addr ){
      double4_t res;
      res.low  = _mm_loadu_pd( mem_addr   );
      res.high = _mm_loadu_pd( mem_addr+2 );
      return res;
   }

   inline void loadNeighbors( const double * p, double4_t & r_left, double4_t & r_center, double4_t & r_right )
   {
      r_left   = load_unaligned( p-1 );
      r_center = load_aligned  ( p   );
      r_right  = load_unaligned( p+1 );
   }

   inline void  store_aligned ( double * mem_addr, double4_t a ) {
      _mm_store_pd( mem_addr  , a.low  );
      _mm_store_pd( mem_addr+2, a.high );
   }

   inline double getComponent ( const double4_t & v, int i           ) { return reinterpret_cast<const double*>(&v)[i]; }
   inline double getComponent ( const double4_t & v, unsigned long i ) { return reinterpret_cast<const double*>(&v)[i]; }

   inline bool   getBoolComponent ( const double4_t & v, int i           ) { return (reinterpret_cast<const uint64_t*>(&v)[i]) != 0; }
   inline bool   getBoolComponent ( const double4_t & v, unsigned long i ) { return (reinterpret_cast<const uint64_t*>(&v)[i]) != 0; }

   inline double4_t hadd( double4_t a,  double4_t b ) {
      double4_t res;
      res.low  = _mm_set_pd( getComponent(b,0) + getComponent(b,1),
                             getComponent(a,0) + getComponent(a,1) );
      res.high = _mm_set_pd( getComponent(b,2) + getComponent(b,3),
                             getComponent(a,2) + getComponent(a,3) );
      return res;
   }

   inline double4_t horizontalSum ( double4_t a ) {
     const double sum = getComponent(a,0) + getComponent(a,1) + getComponent(a,2) + getComponent(a,3);
     return make_double4( sum );
   }

   inline double4_t exchangeLowerUpperHalf ( double4_t a )
   {
      double4_t res;
      res.low  = a.high;
      res.high = a.low;
      return res;
   }


   inline void extract( double4_t in, double4_t &d, double4_t &c, double4_t &b, double4_t & a )
   {
      double const * vals = reinterpret_cast<double const * >(&in);
      a.low  = _mm_set1_pd( vals[0] );
      a.high = _mm_set1_pd( vals[0] );

      b.low  = _mm_set1_pd( vals[1] );
      b.high = _mm_set1_pd( vals[1] );

      c.low  = _mm_set1_pd( vals[2] );
      c.high = _mm_set1_pd( vals[2] );

      d.low  = _mm_set1_pd( vals[3] );
      d.high = _mm_set1_pd( vals[3] );
   }

   inline double4_t rotateRight( const double4_t & a )
   {
      double const * v = reinterpret_cast<double const * >(&a);
      return make_double4( v[0], v[3], v[2], v[1]  );
   }

   inline double4_t rotateLeft( const double4_t & a )
   {
      double const * v = reinterpret_cast<double const * >(&a);
      return make_double4( v[2], v[1], v[0], v[3]  );
   }


   inline double4_t compareEQ( double4_t a, double4_t b ) {
      double4_t result;
      result.low  = _mm_cmpeq_pd( a.low ,b.low  );
      result.high = _mm_cmpeq_pd( a.high,b.high );
      return result;
   }

   inline double4_t compareNEQ( double4_t a, double4_t b ) {
      double4_t result;
      result.low  = _mm_cmpneq_pd( a.low ,b.low  );
      result.high = _mm_cmpneq_pd( a.high,b.high );
      return result;
   }
   inline double4_t compareGE( double4_t a, double4_t b ) {
      double4_t result;
      result.low  = _mm_cmpge_pd( a.low ,b.low  );
      result.high = _mm_cmpge_pd( a.high,b.high );
      return result;
   }
   inline double4_t compareLE( double4_t a, double4_t b ) {
      double4_t result;
      result.low  = _mm_cmple_pd( a.low ,b.low  );
      result.high = _mm_cmple_pd( a.high,b.high );
      return result;
   }

   inline double4_t logicalAND( double4_t a, double4_t b ) {
      double4_t result;
      result.low  = _mm_and_pd( a.low ,b.low  );
      result.high = _mm_and_pd( a.high,b.high );
      return result;
   }
   inline double4_t logicalOR( double4_t a, double4_t b ) {
      double4_t result;
      result.low  = _mm_or_pd( a.low ,b.low  );
      result.high = _mm_or_pd( a.high,b.high );
      return result;
   }



   inline int movemask( double4_t m ) {
      int res;
      res = _mm_movemask_pd( m.high );
      res <<= 2;
      res |= _mm_movemask_pd( m.low );
      return res;
   }

   inline double4_t blendv( double4_t a, double4_t b, double4_t mask) {
      return make_double4 (  getBoolComponent(mask,3) ? getComponent(b,3) : getComponent(a,3),
                             getBoolComponent(mask,2) ? getComponent(b,2) : getComponent(a,2),
                             getBoolComponent(mask,1) ? getComponent(b,1) : getComponent(a,1),
                             getBoolComponent(mask,0) ? getComponent(b,0) : getComponent(a,0)
                           );
   }

   template<int mask>
   inline double4_t blend( double4_t a, double4_t b ) {
      return make_double4 (  (uint64_t)( mask & ( 1 << 3)) ? getComponent(b,3) : getComponent(a,3),
                             (uint64_t)( mask & ( 1 << 2)) ? getComponent(b,2) : getComponent(a,2),
                             (uint64_t)( mask & ( 1 << 1)) ? getComponent(b,1) : getComponent(a,1),
                             (uint64_t)( mask & ( 1 << 0)) ? getComponent(b,0) : getComponent(a,0)
                           );
   }


   inline double4_t sqrt( double4_t a) {
     double4_t result;
     result.high = _mm_sqrt_pd(a.high);
     result.low  = _mm_sqrt_pd(a.low);
      return result;
   }

   template< unsigned int numIter = 3 >
   inline double4_t invSqrt( double4_t y )
   {
      double4_t yHalf = make_double4( 0.5 ) * y;

      __m128i magicConstant = _mm_set_epi64x ( 0x5fe6ec85e7de30daLL, 0x5fe6ec85e7de30daLL );

      // lower part
      {
         __m128i shifted = _mm_srli_epi64( ( __m128i ) y.low, 1 );
         __m128i i = _mm_sub_epi64( magicConstant, shifted );
         y.low = ( __m128d ) i;
      }
      // upper part
      {
         __m128i shifted = _mm_srli_epi64( ( __m128i ) y.high, 1 );
         __m128i i = _mm_sub_epi64( magicConstant, shifted );
         y.high = ( __m128d ) i;
      }

      double4_t onePointFive = make_double4( 1.5 );
      for( unsigned int k=0; k < numIter; ++k )
      y = y * ( onePointFive - yHalf * y * y );

      return y;
   }

} // namespace sse2
} // namespace simd
} // namespace walberla


#endif // __SSE2__
