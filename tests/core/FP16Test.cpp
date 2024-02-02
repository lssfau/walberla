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
//! \file FP16Test.cpp
//! \ingroup core
//! \author Nils Kohl <nils.kohl@fau.de>
//
//======================================================================================================================

#include "core/DataTypes.h"
#include "core/debug/Debug.h"
#include "core/debug/TestSubsystem.h"
#include "core/logging/Logging.h"
#include "core/Environment.h"

#include <cstdlib>
#include <iostream>

namespace walberla {

void fp16Test( int argc, char ** argv )
{
   Environment const env( argc, argv );

   WALBERLA_LOG_INFO_ON_ROOT("-------------")
   WALBERLA_LOG_INFO_ON_ROOT(" FP16 checks ")
   WALBERLA_LOG_INFO_ON_ROOT("-------------")

#ifndef WALBERLA_BUILD_WITH_HALF_PRECISION_SUPPORT
   WALBERLA_LOG_INFO_ON_ROOT(" - Test does nothing as it was not built with fp16 support.")
   WALBERLA_LOG_INFO_ON_ROOT(" - Apparently you have not enabled half precision support.")
   WALBERLA_LOG_INFO_ON_ROOT(" - Reconfigure by setting the respective CMake variable "
                             "(at the time of writing this it's called WALBERLA_BUILD_WITH_HALF_PRECISION_SUPPORT) "
                             "to ON.")
#else

   WALBERLA_LOG_INFO_ON_ROOT(" - Half precision support enabled via CMake!")

   WALBERLA_LOG_INFO_ON_ROOT(" - Sizeof checks: ")
   auto sfloat64 = sizeof(float64);
   auto sfloat32 = sizeof(float32);
   auto sfloat16 = sizeof(float16);
   WALBERLA_CHECK_EQUAL( sfloat64, 8, "Your types don't seem to have the expected sizes." );
   WALBERLA_CHECK_EQUAL( sfloat32, 4, "Your types don't seem to have the expected sizes." );
   WALBERLA_CHECK_EQUAL( sfloat16, 2, "Your types don't seem to have the expected sizes." );

   WALBERLA_LOG_INFO_ON_ROOT(" - Casting checks (promotion is required to format strings): ")
   const float64 a64 = 42;
   const float32 a32 = 42;
   const float16 a16 = 42;
   WALBERLA_LOG_INFO_ON_ROOT("   + float64: " << a64)
   WALBERLA_LOG_INFO_ON_ROOT("   + float32: " << a32)
   WALBERLA_LOG_INFO_ON_ROOT("   + float16: " << (double) a16)
   WALBERLA_LOG_INFO_ON_ROOT("   Casting and output compiles.")

   WALBERLA_LOG_INFO_ON_ROOT(" - Basic arithmetic check: ")
   const float16 x = 1.2f16;
   const float16 y = -1.8f16;
   const float64 z = -0.6;
   WALBERLA_LOG_INFO_ON_ROOT("   + " << (double) x << " + " << (double) y << " == " << (float64) (x + y) << " ? ")
   WALBERLA_CHECK_FLOAT_EQUAL( (x + y), (float16) z, "float16 addition does not work correctly.");
#endif
}

}


int main( int argc, char** argv )
{
   walberla::debug::enterTestMode();
   walberla::fp16Test( argc, argv );
   return EXIT_SUCCESS;
}
