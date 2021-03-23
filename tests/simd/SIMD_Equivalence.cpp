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
//! \file SSE_AVX_Equivalence.h
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#include "core/debug/TestSubsystem.h"
#include "core/Environment.h"
#include "core/logging/Logging.h"

#include "simd/SIMD.h"


#define QUIT_ON_ERROR 1


//===================================================================================================================
//
//  Conditional include of relevant SIMD Header
//
//===================================================================================================================


#ifdef WALBERLA_SIMD_AVX2_AVAILABLE
#include "simd/AVX2.h"
#endif

#ifdef WALBERLA_SIMD_AVX_AVAILABLE
#include "simd/AVX.h"
#endif


#ifdef WALBERLA_SIMD_SSE4_AVAILABLE
#include "simd/SSE4.h"
#endif

#ifdef WALBERLA_SIMD_SSE2_AVAILABLE
#include "simd/SSE2.h"
#endif

#include "simd/Scalar.h"

#include "simd/StreamOutput.h"


//===================================================================================================================
//
//  Fallback to scalar emulation if instruction set not available
//
//===================================================================================================================

// ---------------- AVX2 ------------
#ifdef IS0_AVX2
#ifdef WALBERLA_SIMD_AVX2_AVAILABLE
#define is0 avx2
#else
#define is0 scalar
#endif
#endif

#ifdef IS1_AVX2
#ifdef WALBERLA_SIMD_AVX2_AVAILABLE
#define is1 avx2
#else
#define is1 scalar
#endif
#endif

// ---------------- AVX  ------------
#ifdef IS0_AVX
#ifdef WALBERLA_SIMD_AVX_AVAILABLE
#define is0 avx
#else
#define is0 scalar
#endif
#endif

#ifdef IS1_AVX
#ifdef WALBERLA_SIMD_AVX_AVAILABLE
#define is1 avx
#else
#define is1 scalar
#endif
#endif

// ---------------- SSE4  ------------
#ifdef IS0_SSE4
#ifdef WALBERLA_SIMD_SSE4_AVAILABLE
#define is0 sse4
#else
#define is0 scalar
#endif
#endif

#ifdef IS1_SSE4
#ifdef WALBERLA_SIMD_SSE4_AVAILABLE
#define is1 sse4
#else
#define is1 scalar
#endif
#endif


// ---------------- SSE2  ------------
#ifdef IS0_SSE2
#ifdef WALBERLA_SIMD_SSE2_AVAILABLE
#define is0 sse2
#else
#define is0 scalar
#endif
#endif

#ifdef IS1_SSE2
#ifdef WALBERLA_SIMD_SSE2_AVAILABLE
#define is1 sse2
#else
#define is1 scalar
#endif
#endif


// ---------------- Scalar -----------
#ifdef IS0_SCALAR
#define is0 scalar
#endif

#ifdef IS1_SCALAR
#define is1 scalar
#endif



using namespace walberla;
using namespace simd;


void print0( const is0::double4_t& vec)
{
   WALBERLA_LOG_DEVEL( is0::getComponent(vec, 0 ) );
   WALBERLA_LOG_DEVEL( is0::getComponent(vec, 1 ) );
   WALBERLA_LOG_DEVEL( is0::getComponent(vec, 2 ) );
   WALBERLA_LOG_DEVEL( is0::getComponent(vec, 3 ) );
}

void print1( const is1::double4_t& vec)
{
   WALBERLA_LOG_DEVEL( is1::getComponent(vec, 0 ) );
   WALBERLA_LOG_DEVEL( is1::getComponent(vec, 1 ) );
   WALBERLA_LOG_DEVEL( is1::getComponent(vec, 2 ) );
   WALBERLA_LOG_DEVEL( is1::getComponent(vec, 3 ) );
}

void checkVecEqual ( is0::double4_t a, is1::double4_t b, const std::string & description = "" )
{
   if ( description.size() > 0) 
   {
      double meanError = std::abs( is0::getComponent(a, 0 )- is1::getComponent(b,0 ) ) +
			 std::abs( is0::getComponent(a, 1 )- is1::getComponent(b,1 ) ) +
			 std::abs( is0::getComponent(a, 2 )- is1::getComponent(b,2 ) ) +
			 std::abs( is0::getComponent(a, 3 )- is1::getComponent(b,3 ) );
      WALBERLA_LOG_RESULT("Checking " << description << " mean error " << meanError );                          
   }
#if QUIT_ON_ERROR
   WALBERLA_CHECK_FLOAT_EQUAL( is0::getComponent(a, 0 ), is1::getComponent(b,0 )  );
   WALBERLA_CHECK_FLOAT_EQUAL( is0::getComponent(a, 1 ), is1::getComponent(b,1 )  );
   WALBERLA_CHECK_FLOAT_EQUAL( is0::getComponent(a, 2 ), is1::getComponent(b,2 )  );
   WALBERLA_CHECK_FLOAT_EQUAL( is0::getComponent(a, 3 ), is1::getComponent(b,3 )  );
#endif
  
}


void basic()
{
   auto a1 = is0::make_double4( 1.0, 2.0, 3.0, 4.0 );
   auto a2 = is0::make_double4( 6.0, 7.0, 8.0, 9.0 );
   auto b1 = is1::make_double4( 1.0, 2.0, 3.0, 4.0 );
   auto b2 = is1::make_double4( 6.0, 7.0, 8.0, 9.0 );
   checkVecEqual( a1, b1 );

   checkVecEqual( a1+a2, b1+b2, "Addition" );
   checkVecEqual( a1-a2, b1-b2, "Subtraction" );
   checkVecEqual( a1*a2, b1*b2, "Multiplication" );
   checkVecEqual( a1/a2, b1/b2, "Division" );

   checkVecEqual( is0::rotateLeft (a1), is1::rotateLeft (b1), "Rotate reft" );
   checkVecEqual( is0::rotateRight(a1), is1::rotateRight(b1), "Rotate right" );

   checkVecEqual( is0::horizontalSum(a1), is1::horizontalSum(b1), "HorizontalSum" );
   checkVecEqual( is0::exchangeLowerUpperHalf(a1), is1::exchangeLowerUpperHalf(b1), "exchangeLowerUpperHalf" );

}

void extract()
{
  auto inA = is0::make_double4( 1.0, 2.0, 3.0, 4.0 );
  auto inB = is1::make_double4( 1.0, 2.0, 3.0, 4.0 );

  is0::double4_t a0, a1, a2, a3;
  is1::double4_t b0, b1, b2, b3;
  is0::extract( inA, a3, a2, a1, a0 );
  is1::extract( inB, b3, b2, b1, b0 );
  checkVecEqual( a0, b0, "Extract 0" );
  checkVecEqual( a1, b1, "Extract 1" );
  checkVecEqual( a2, b2, "Extract 2" );
  checkVecEqual( a3, b3, "Extract 3" );
}


void comparisonAndBlend()
{
  auto inA = is0::make_double4_r( 1.0, 2.0, 3.0, 4.0 );
  auto inB = is1::make_double4_r( 1.0, 2.0, 3.0, 4.0 );

  is0::double4_t maskvA = is0::compareGE( inA, is0::make_double4( 3.0)  );
  is1::double4_t maskvB = is1::compareGE( inB, is1::make_double4( 3.0)  );

  WALBERLA_LOG_DEVEL("-------------------");
  print0(maskvA);
  print1(maskvB);
  WALBERLA_LOG_DEVEL("-------------------");
  print0(is0::blendv(inA, is0::make_zero(), maskvA ));
  print1(is1::blendv(inB, is1::make_zero(), maskvB ));
  WALBERLA_LOG_DEVEL("-------------------");

  checkVecEqual( is0::blendv(inA, is0::make_zero(), maskvA ),
                 is1::blendv(inB, is1::make_zero(), maskvB ), "comparisonAndBlend");

  WALBERLA_CHECK_EQUAL( is0::movemask( maskvA ), is1::movemask( maskvB ) );
}

void blendInteger()
{
   auto inA1  = is0::make_double4( 4.0, 3.0, 2.0, 1.0 );
   auto inA2  = is0::make_double4( 8.0, 7.0, 6.0, 5.0 );

   auto inB1  = is1::make_double4( 4.0, 3.0, 2.0, 1.0 );
   auto inB2  = is1::make_double4( 8.0, 7.0, 6.0, 5.0 );

   const int mask = 10;
   auto refA = is0::make_double4( 8.0, 3.0, 6.0, 1.0 );
   auto refB = is1::make_double4( 8.0, 3.0, 6.0, 1.0 );

   auto resA = is0::blend<mask>( inA1, inA2 );
   auto resB = is1::blend<mask>( inB1, inB2 );

   checkVecEqual( resA, resB );
   for( int i=0; i< 4; ++i )
   {
      WALBERLA_CHECK_FLOAT_EQUAL( is0::getComponent(resA,i), is0::getComponent(refA,i) );
      WALBERLA_CHECK_FLOAT_EQUAL( is1::getComponent(resB,i), is1::getComponent(refB,i) );
   }
}


void sqrtTest()
{
   auto inA = is0::make_double4_r( 1.0, 2.0, 3.0, 4.0 );
   auto inB = is1::make_double4_r( 1.0, 2.0, 3.0, 4.0 );

   checkVecEqual( is0::invSqrt<1>(inA), is1::invSqrt<1>(inB), "invSqrt1" );
   checkVecEqual( is0::invSqrt<3>(inA), is1::invSqrt<3>(inB), "invSqrt3" );
}

void compareAndMaskTest()
{
   auto a_v1 = is0::make_double4_r( 1.0, 2.0, 3.0, 4.0 );
   auto a_v2 = is0::make_double4_r( 1.0, 2.0, 3.0, 4.0 );

   auto b_v1 = is1::make_double4_r( 1.0, 2.0, 3.0, 4.0 );
   auto b_v2 = is1::make_double4_r( 1.0, 2.0, 3.0, 4.0 );

   auto a_compareRes = is0::compareEQ( a_v1, a_v2 );
   auto a_mask       = is0::movemask( a_compareRes );

   auto b_compareRes = is1::compareEQ( b_v1, b_v2 );
   auto b_mask       = is1::movemask( b_compareRes );

   WALBERLA_CHECK_EQUAL( a_mask, b_mask );
}



int main( int argc, char ** argv )
{
   debug::enterTestMode();
   walberla::Environment walberlaEnv( argc, argv );

   WALBERLA_LOG_RESULT("Comparing " << is0::usedInstructionSet() << " and " << is1::usedInstructionSet() );

   if ( std::string( is0::usedInstructionSet() ) ==  std::string( is1::usedInstructionSet() ) )
   {
      WALBERLA_LOG_RESULT("Nothing to compare since instruction set not available on this machine" );
      return 0;
   }

   basic();
   extract();
   comparisonAndBlend();
   sqrtTest();
   blendInteger();
   compareAndMaskTest();
   return 0;
}


