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

#include "immintrin.h"

#include "waLBerlaDefinitions.h"
#include "core/DataTypes.h"
#include "IntelVecTypesCppOperators.h"


namespace walberla {
namespace simd {
namespace avx2 {

   typedef __m256d double4_t;

   inline const char * usedInstructionSet() { return "AVX2"; }

   inline double4_t make_double4   ( double d, double c, double b, double a ) { return _mm256_set_pd  ( d,c,b,a ); }
   inline double4_t make_double4_r ( double a, double b, double c, double d ) { return _mm256_setr_pd ( a,b,c,d ); }

   inline double4_t make_double4  ( double a                               ) { return _mm256_set_pd ( a,a,a,a ); }

   inline double4_t make_zero()   { return _mm256_setzero_pd(); }

   inline double4_t load_aligned  ( double const * mem_addr )                { return _mm256_load_pd (mem_addr); }
   inline double4_t load_unaligned ( double const * mem_addr )               { return _mm256_loadu_pd (mem_addr); }
   inline void      store_aligned ( double * mem_addr, double4_t a )         { _mm256_store_pd ( mem_addr, a) ;  }

   inline void loadNeighbors( const double * p, double4_t & r_left, double4_t & r_center, double4_t & r_right )
   {
      r_left   = load_unaligned( p-1 );
      r_center = load_aligned  ( p   );
      r_right  = load_unaligned( p+1 );
   }

   inline double getComponent ( const double4_t & v, int           i ) { return reinterpret_cast<const double*>(&v)[i]; }
   inline double getComponent ( const double4_t & v, unsigned long i ) { return reinterpret_cast<const double*>(&v)[i]; }

   inline bool   getBoolComponent ( const double4_t & v, int i           ) { return (reinterpret_cast<const uint64_t*>(&v)[i]) != 0; }
   inline bool   getBoolComponent ( const double4_t & v, unsigned long i ) { return (reinterpret_cast<const uint64_t*>(&v)[i]) != 0; }

   inline double4_t hadd( double4_t a,  double4_t b ) { return _mm256_hadd_pd ( a,b); }

   inline double4_t horizontalSum ( double4_t a )
   {
     double4_t t1 = avx2::hadd(a,a);
     double4_t t2 = _mm256_permute2f128_pd(t1, t1, 0x01);  // exchange lower and upper half
     return  _mm256_add_pd(t1,t2);
   }

   inline double4_t exchangeLowerUpperHalf ( double4_t a )
   {
      return _mm256_permute2f128_pd(a, a, 0x01);
   }


   inline void extract( double4_t in, double4_t &d, double4_t &c, double4_t &b, double4_t & a )
   {
      double const * vals = reinterpret_cast<double const * >(&in);
      a = _mm256_broadcast_sd ( &vals[0] );
      b = _mm256_broadcast_sd ( &vals[1] );
      c = _mm256_broadcast_sd ( &vals[2] );
      d = _mm256_broadcast_sd ( &vals[3] );
   }

   inline double4_t rotateRight( double4_t a )
   {
//      //Quelle: https://software.intel.com/en-us/forums/topic/343530
//                                                              // [3  2  1  0]  ->a
//      __m256d t0 = _mm256_permute2f128_pd(a, a,  0x01 );       // [1  0  3  2]  ->b
//      return _mm256_castsi256_pd(_mm256_alignr_epi8(  _mm256_castpd_si256(t0), _mm256_castpd_si256(a), 0x08 ));
//                                                              // ->a   [ x x 1 0 ]
//                                                              // ->a   [ 1 0 x x ]
//                                                              // ->b   [ x x 3 2 ]
//                                                              // a|b   [ 1 0 3 2 ]
//                                                              // a|b>> [ 0 3 2 1 ]
//                                                              //       [ x x 2 1 ] -> temp
//                                                              // ->a   [ x x 3 2 ]
//                                                              // ->a   [ 3 2 x x ]
//                                                              // ->b   [ x x 1 0 ]
//                                                              // a|b   [ 3 2 1 0 ]
//                                                              // a|b>> [ 2 1 0 3 ]
//                                                              //       [ x x 0 3 ] -> temp1
//                                                              //  ==>  [ 0 3 2 1 ]
      return _mm256_permute4x64_pd(a, 0x39);
   }

   inline double4_t rotateLeft( double4_t a )
   {
//                                                              // [3  2  1  0]  ->a
//      __m256d t0 = _mm256_permute2f128_pd(a, a,  0x01 );       // [1  0  3  2]  ->b
//      return _mm256_castsi256_pd(_mm256_alignr_epi8(_mm256_castpd_si256(a), _mm256_castpd_si256(t0), 0x08 ));
//                                                              // ->a   [ x x 1 0 ]
//                                                              // ->a   [ 1 0 x x ]
//                                                              // ->b   [ x x 3 2 ]
//                                                              // a|b   [ 1 0 3 2 ]
//                                                              // a|b>> [ 2 1 0 3 ]
//                                                              //       [ x x 0 3 ] -> temp
//                                                              // ->a   [ x x 3 2 ]
//                                                              // ->a   [ 3 2 x x ]
//                                                              // ->b   [ x x 1 0 ]
//                                                              // a|b   [ 3 2 1 0 ]
//                                                              // a|b>> [ 0 3 2 1 ]
//                                                              //       [ x x 2 1 ] -> temp1
//                                                              //  ==>  [ 2 1 0 3 ]
      return _mm256_permute4x64_pd(a, 0x93);
   }



   inline double4_t compareEQ( double4_t a, double4_t b ) {
      return _mm256_cmp_pd ( a, b, _CMP_EQ_UQ );
   }
   inline double4_t compareNEQ( double4_t a, double4_t b ) {
      return _mm256_cmp_pd ( a, b, _CMP_NEQ_UQ );
   }
   inline double4_t compareGE( double4_t a, double4_t b ) {
      return _mm256_cmp_pd ( a, b, _CMP_GE_OQ );
   }
   inline double4_t compareLE( double4_t a, double4_t b ) {
      return _mm256_cmp_pd ( a, b, _CMP_LE_OQ );
   }

   inline double4_t logicalAND( double4_t a, double4_t b ) {
      return _mm256_and_pd ( a, b );
   }
   inline double4_t logicalOR( double4_t a, double4_t b ) {
      return _mm256_or_pd ( a, b );
   }


   inline int movemask( double4_t m ) {
      return _mm256_movemask_pd ( m );
   }
   inline double4_t blendv( double4_t a, double4_t b, double4_t mask) {
      return _mm256_blendv_pd (a,b,mask);
   }
   template<int mask>
   inline double4_t blend( double4_t a, double4_t b ) {
      return _mm256_blend_pd (a,b,mask);
   }

   inline double4_t sqrt( double4_t a) {
      return  _mm256_sqrt_pd (a );
   }

   template< unsigned int numIter = 3 >
   inline double4_t invSqrt( double4_t y )
   {
      //  (Add, Mul, Div ): (numIter, numIter*3+1,0)

      const long long MAGIC = 0x5fe6ec85e7de30daLL;

      double4_t yHalf = make_double4( 0.5 ) * y;

      auto yAsInt = reinterpret_cast<const long long*>(&y);

      __m256i i = _mm256_setr_epi64x( MAGIC - (yAsInt[0] >> 1),
                                      MAGIC - (yAsInt[1] >> 1),
                                      MAGIC - (yAsInt[2] >> 1),
                                      MAGIC - (yAsInt[3] >> 1) );

      y = (double4_t) i;
      double4_t onePointFive = make_double4( 1.5 );
      for( unsigned int k=0; k < numIter; ++k )
         y = y * ( onePointFive - yHalf * y * y );

      return y;
   }


} // namespace avx2
} // namespace simd
} // namespace walberla



