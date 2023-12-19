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
//! \file Float16SupportTest.cpp
//! \ingroup core
//! \author Michael Zikeli <michael.zikeli@fau.de>
//
//======================================================================================================================

#include <memory>
#include <numeric>

#include "core/DataTypes.h"
#include "core/Environment.h"
#include "core/logging/Logging.h"

namespace walberla::simple_Float16_test {
using walberla::floatIsEqual;
using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;

// === Choosing Accuracy ===
//+++ Precision : fp16 +++
using walberla::float16;
using walberla::float32;
using walberla::float64;
using dst_t                         = float16;
using src_t                         = real_t;
constexpr real_t     precisionLimit = walberla::float16( 1e-3 );
const std::string    precisionType  = "float16";
constexpr const auto maxLevel       = uint_t( 3 );

void simple_array_test()
{
   auto fpSrc = std::make_shared< src_t[] >( 10 );
   auto fpDst = std::make_shared< dst_t[] >( 10 );

   std::fill_n( fpSrc.get(), 10, 17. );
   std::fill_n( fpDst.get(), 10, (dst_t) 17. );

   fpSrc[5] = 8.;
   fpDst[5] = (dst_t) 8.;

   // Test equality with custom compare
   WALBERLA_CHECK_LESS( std::fabs( fpSrc[9] - (src_t) fpDst[9] ), precisionLimit );
   WALBERLA_CHECK_LESS( std::fabs( fpSrc[5] - (src_t) fpDst[5] ), precisionLimit );
   // Test specialized floatIsEqual
   WALBERLA_CHECK( floatIsEqual( fpSrc[9], (src_t) fpDst[9], (src_t) precisionLimit ) );
   WALBERLA_CHECK( floatIsEqual( (dst_t) fpSrc[9], fpDst[9], (dst_t) precisionLimit ) );
   WALBERLA_CHECK_FLOAT_EQUAL( (dst_t) fpSrc[9], fpDst[9] );

   // Test std::fill_n
   auto other_fpDst = std::make_shared< dst_t[] >( 10 );
   std::fill_n( other_fpDst.get(), 10, (dst_t) 2. );
   WALBERLA_CHECK_FLOAT_EQUAL( (dst_t) 2., other_fpDst[9] );
   WALBERLA_CHECK_FLOAT_EQUAL( (dst_t) 2., other_fpDst[5] );

   // Test std::swap
   std::swap( fpDst, other_fpDst );
   fpDst[5] = (dst_t) 9.;

   WALBERLA_CHECK_FLOAT_EQUAL( (dst_t) fpSrc[9], other_fpDst[9] );
   WALBERLA_CHECK_FLOAT_EQUAL( (dst_t) fpSrc[5], other_fpDst[5] );
   WALBERLA_CHECK_FLOAT_EQUAL( (dst_t) 2., fpDst[9] );
   WALBERLA_CHECK_FLOAT_EQUAL( (dst_t) 9., fpDst[5] );

} // simple_Float16_test::simple_array_test()

void vector_test()
{
   auto fpSrc      = std::vector< src_t >( 10 );
   auto fpDst_cast = std::vector< dst_t >( 10 );
   auto fp32       = std::vector< walberla::float32 >( 10 );
   auto fpDst      = std::vector< dst_t >( 10 );

   fpSrc.assign( 10, 1.5 );
   fpDst_cast.assign( 10, (dst_t) 1.5 );
   fp32.assign( 10, 1.5f );
   std::copy( fpSrc.begin(), fpSrc.end(), fpDst.begin() );
   WALBERLA_LOG_WARNING_ON_ROOT(
       " std::vector.assign is not able to assign "
       << typeid( src_t ).name() << " values to container of type " << precisionType << ".\n"
       << " Therefore, the floating point value for assign must be cast beforehand or std::copy must be used, since it uses a static_cast internally." );

   fpSrc[5]      = 2.3;
   fpDst_cast[5] = (dst_t) 2.3;
   fp32[5]       = 2.3f;
   fpDst[5]      = (dst_t) 2.3;

   WALBERLA_CHECK_FLOAT_EQUAL( (walberla::float32) fpSrc[0], fp32[0] );
   WALBERLA_CHECK_FLOAT_EQUAL( (walberla::float32) fpSrc[9], fp32[9] );
   WALBERLA_CHECK_FLOAT_EQUAL( (walberla::float32) fpSrc[5], fp32[5] );
   WALBERLA_CHECK_FLOAT_EQUAL( (dst_t) fpSrc[0], fpDst_cast[0] );
   WALBERLA_CHECK_FLOAT_EQUAL( (dst_t) fpSrc[9], fpDst_cast[9] );
   WALBERLA_CHECK_FLOAT_EQUAL( (dst_t) fpSrc[5], fpDst_cast[5] );
   WALBERLA_CHECK_FLOAT_EQUAL( (dst_t) fpSrc[0], fpDst[0] );
   WALBERLA_CHECK_FLOAT_EQUAL( (dst_t) fpSrc[9], fpDst[9] );
   WALBERLA_CHECK_FLOAT_EQUAL( (dst_t) fpSrc[5], fpDst[5] );
   WALBERLA_CHECK_EQUAL( typeid( fpDst ), typeid( fpDst_cast ) );

   // Add up all elements of the vector to check whether the result is sufficiently correct.
   {
      const auto sumSrc = std::reduce(fpSrc.begin(), fpSrc.end());
      const auto sumDst = std::reduce(fpDst.begin(), fpDst.end());
      WALBERLA_CHECK_FLOAT_EQUAL( (dst_t)sumSrc, sumDst );
   }
   {
      fpSrc.assign( 13, 1.3 );
      std::copy( fpSrc.begin(), fpSrc.end(), fpDst.begin() );
      const auto sumSrc = std::reduce(fpSrc.begin(), fpSrc.end());
      const auto sumDst = std::reduce(fpDst.begin(), fpDst.end());
      WALBERLA_CHECK_FLOAT_UNEQUAL( (dst_t)sumSrc, sumDst );
   }
} // simple_Float16_test::vector_test()

int main( int argc, char** argv )
{
   // This check only works since C++23 and is used in many implementations, so it's important, that it works.
   WALBERLA_CHECK( std::is_arithmetic< dst_t >::value );

   walberla::Environment env( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::INFO );
   walberla::MPIManager::instance()->useWorldComm();

   WALBERLA_LOG_INFO_ON_ROOT( " This run is executed with " << precisionType );
   WALBERLA_LOG_INFO_ON_ROOT( " machine precision limit is " << precisionLimit );
   const std::string stringLine( 125, '=' );
   WALBERLA_LOG_INFO_ON_ROOT( stringLine );

   WALBERLA_LOG_INFO_ON_ROOT( " Start a test with shared_pointer<float16[]>." );
   simple_array_test();

   WALBERLA_LOG_INFO_ON_ROOT( " Start a test with std::vector<float16>." );
   vector_test();

   WALBERLA_LOG_INFO_ON_ROOT( " Start a where float32 is sufficient but float16 is not." );
   WALBERLA_CHECK_FLOAT_UNEQUAL( dst_t(1.0)-dst_t(0.3), 1.0-0.3 );
   WALBERLA_CHECK_FLOAT_EQUAL( 1.0f-0.3f, 1.0-0.3 );

   return 0;
} // simple_Float16_test::main()

} // namespace walberla::simple_Float16_test

int main( int argc, char** argv )
{
   walberla::simple_Float16_test::main( argc, argv );

   return EXIT_SUCCESS;
} // main()
