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
//! \file CheckFP16.cpp
//! \brief Checks the availability of float16 (half precision) and verifies some properties.
//! \author Nils Kohl <nils.kohl@fau.de>
//
//======================================================================================================================

#include <core/DataTypes.h>
#include <core/Environment.h>
#include <core/logging/Logging.h>
#include <core/perf_analysis/extern/likwid.h>

namespace walberla
{

template< typename T >
void kernel(T* v, T* vv, T* r, size_t vsize)
{
   for (size_t i = 0; i < vsize; i++)
   {
      r[i] = v[i] + vv[i];
   }
}

int main(int argc, char** argv)
{
   Environment const env(argc, argv);

   WALBERLA_LOG_INFO_ON_ROOT("-------------")
   WALBERLA_LOG_INFO_ON_ROOT(" FP16 checks ")
   WALBERLA_LOG_INFO_ON_ROOT("-------------")

#ifndef WALBERLA_BUILD_WITH_HALF_PRECISION_SUPPORT
   WALBERLA_LOG_INFO_ON_ROOT(" - Apparently you have not enabled half precision support.")
   WALBERLA_LOG_INFO_ON_ROOT("   Reconfigure by setting the respective CMake variable to ON.")
   WALBERLA_LOG_INFO_ON_ROOT("   At the time of writing this it's called WALBERLA_BUILD_WITH_HALF_PRECISION_SUPPORT.")

   return EXIT_FAILURE;
#else
   WALBERLA_LOG_INFO_ON_ROOT(" - Half precision support enabled via CMake!")

   WALBERLA_LOG_INFO_ON_ROOT(" - Sizeof checks: ")
   const auto sfloat64 = sizeof(float64);
   const auto sfloat32 = sizeof(float32);
   const auto sfloat16 = sizeof(float16);
   WALBERLA_LOG_INFO_ON_ROOT("   + sizeof( float64 ) == " << sfloat64)
   WALBERLA_LOG_INFO_ON_ROOT("   + sizeof( float32 ) == " << sfloat32)
   WALBERLA_LOG_INFO_ON_ROOT("   + sizeof( float16 ) == " << sfloat16)
   if (sfloat64 != 8 || sfloat32 != 4 || sfloat16 != 2)
   {
      WALBERLA_LOG_INFO_ON_ROOT("   Your types don't seem to have the expected sizes.")
      return EXIT_FAILURE;
   }
   WALBERLA_LOG_INFO_ON_ROOT("   -> works out!")

   WALBERLA_LOG_INFO_ON_ROOT(" - Casting checks (promotion is required to format strings): ")
   const float64 a64 = 42;
   const float32 a32 = 42;
   const float16 a16 = 42;
   WALBERLA_LOG_INFO_ON_ROOT("   + float64: " << a64)
   WALBERLA_LOG_INFO_ON_ROOT("   + float32: " << a32)
   WALBERLA_LOG_INFO_ON_ROOT("   + float16: " << (double) a16)
   WALBERLA_LOG_INFO_ON_ROOT("   Casting and output compiles.")

   WALBERLA_LOG_INFO_ON_ROOT(" - Basic arithmetic check: ")
   const auto x   = float16(1.2);
   const auto y   = float16(-1.8);
   const float64 z   = -0.6;
   const float16 sum = x + y;
   WALBERLA_LOG_INFO_ON_ROOT("     " << (double) x << " + " << (double) y << " == " << (float64) (x + y) << "")
   WALBERLA_CHECK(std::abs((float64) sum - z) < 1e-3, "Float16 arithmetic is broken.");
   WALBERLA_LOG_INFO_ON_ROOT("")

#   ifdef WALBERLA_BUILD_WITH_LIKWID_MARKERS
   WALBERLA_LOG_INFO_ON_ROOT(" - Memory traffic test. You have built with likwid enabled. Make sure to run ")
   WALBERLA_LOG_INFO_ON_ROOT("     $ likwid-perfctr -g MEM_DP    -m ./CheckFP16")
   WALBERLA_LOG_INFO_ON_ROOT("   to compare the memory traffic, and")
   WALBERLA_LOG_INFO_ON_ROOT("     $ likwid-perfctr -g FLOPS_AVX -m ./CheckFP16")
   WALBERLA_LOG_INFO_ON_ROOT(
      "   for the stream-triad-like benchmark to check whether automatic float32 vectorization works.")
   WALBERLA_LOG_INFO_ON_ROOT("")
   WALBERLA_LOG_INFO_ON_ROOT("   The only real benefit of using float16 is reduced memory traffic since internally,\n"
                             "all arithmetic operations are preceded by promotions to float32 (likely - depends on "
                             "the machine).")
   WALBERLA_LOG_INFO_ON_ROOT("   + Stream test ... ")

   LIKWID_MARKER_INIT;
   LIKWID_MARKER_THREADINIT;

   LIKWID_MARKER_REGISTER("float64-mem");
   LIKWID_MARKER_REGISTER("float32-mem");
   LIKWID_MARKER_REGISTER("float16-mem");

   LIKWID_MARKER_REGISTER("float64-vec");
   LIKWID_MARKER_REGISTER("float32-vec");
   LIKWID_MARKER_REGISTER("float16-vec");

   size_t vsize = 100000000;

   std::vector< float64 > v64(vsize, 0.01);
   std::vector< float32 > v32(vsize, 0.01f);
   std::vector< float16 > v16(vsize, float16(0.01));

   std::vector< float64 > vv64(vsize, 0.02);
   std::vector< float32 > vv32(vsize, 0.02f);
   std::vector< float16 > vv16(vsize, float16(0.02));

   std::vector< float64 > r64(vsize);
   std::vector< float32 > r32(vsize);
   std::vector< float16 > r16(vsize);

   LIKWID_MARKER_START("float64-mem");
   float64 sum64 = 0;
   for (size_t j = 0; j < vsize; j++)
   {
      if (0 == j % 2) { sum64 += v64[j]; }
      else { sum64 -= v64[j]; }
   }
   WALBERLA_LOG_INFO_ON_ROOT(
      "   + Printing sum of float64 vector entries. Should be zero up to rounding errors: " << sum64);
   LIKWID_MARKER_STOP("float64-mem");

   // Start measurements
   LIKWID_MARKER_START("float32-mem");
   float32 sum32 = 0;
   for (size_t j = 0; j < vsize; j++)
   {
      if (0 == j % 2) { sum32 += v32[j]; }
      else { sum32 -= v32[j]; }
   }
   WALBERLA_LOG_INFO_ON_ROOT(
      "   + Printing sum of float32 vector entries. Should be zero up to rounding errors: " << sum32);
   LIKWID_MARKER_STOP("float32-mem");

   // Start measurements
   LIKWID_MARKER_START("float16-mem");
   float16 sum16 = 0;
   for (size_t j = 0; j < vsize; j++)
   {
      if (0 == j % 2) { sum16 += v16[j]; }
      else { sum16 -= v16[j]; }
   }
   WALBERLA_LOG_INFO_ON_ROOT(
      "   + Printing sum of float16 vector entries. Should be zero up to rounding errors: " << (double) sum16);
   LIKWID_MARKER_STOP("float16-mem");

   WALBERLA_LOG_INFO_ON_ROOT("   + Vectorization test ... ")

   float64* v64_ptr  = v64.data();
   float64* vv64_ptr = vv64.data();
   float64* r64_ptr  = r64.data();
   LIKWID_MARKER_START("float64-vec");
   kernel(v64_ptr, vv64_ptr, r64_ptr, vsize);
   WALBERLA_LOG_INFO_ON_ROOT("   + Printing entry of float64 vector sum: " << r64[vsize / 2]);
   LIKWID_MARKER_STOP("float64-vec");

   float32* v32_ptr  = v32.data();
   float32* vv32_ptr = vv32.data();
   float32* r32_ptr  = r32.data();
   LIKWID_MARKER_START("float32-vec");
   kernel(v32_ptr, vv32_ptr, r32_ptr, vsize);
   WALBERLA_LOG_INFO_ON_ROOT("   + Printing entry of float32 vector sum: " << r32[vsize / 2]);
   LIKWID_MARKER_STOP("float32-vec");

   float16* v16_ptr  = v16.data();
   float16* vv16_ptr = vv16.data();
   float16* r16_ptr  = r16.data();
   LIKWID_MARKER_START("float16-vec");
   kernel(v16_ptr, vv16_ptr, r16_ptr, vsize);
   WALBERLA_LOG_INFO_ON_ROOT("   + Printing entry of float16 vector sum: " << (double) r16[vsize / 2]);
   LIKWID_MARKER_STOP("float16-vec");

   LIKWID_MARKER_CLOSE;

#   else
   WALBERLA_LOG_INFO_ON_ROOT(" - Build and run with likwid to run memory traffic test.")
#   endif
#endif
   return EXIT_SUCCESS;
}
} // namespace walberla

int main(int argc, char** argv) { return walberla::main(argc, argv); }
