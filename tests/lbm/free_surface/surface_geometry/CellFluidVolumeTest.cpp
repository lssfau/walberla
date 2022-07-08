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
//! \file CellFluidVolumeTest.cpp
//! \ingroup lbm/free_surface/surface_geometry
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//! \brief Calculate fluid volume within a cell and compare with results obtained with ParaView.
//
//======================================================================================================================

#include "core/Environment.h"
#include "core/debug/Debug.h"
#include "core/debug/TestSubsystem.h"

#include "lbm/free_surface/surface_geometry/CurvatureSweep.h"

namespace walberla
{
namespace free_surface
{
namespace CellFluidVolumeTest
{
inline void test(const Vector3< real_t >& normal, real_t offset, real_t expectedVolume, real_t tolerance)
{
   real_t volume = computeCellFluidVolume(normal, offset);
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(expectedVolume, volume, tolerance);
}

int main(int argc, char** argv)
{
   debug::enterTestMode();
   Environment env(argc, argv);

   // allowed deviation from the expected result
   real_t tolerance = real_c(1e-3);

   // test case where simple formula is applied (Figure 2.12, p. 24 in dissertation of Thomas Pohl)
   const Vector3< real_t > normalSimpleCase(real_c(0), real_c(0.5), real_c(0.866));
   const real_t offsetSimpleCase         = real_c(-0.1);
   const real_t expectedVolumeSimpleCase = real_c(0.3845); // obtained with ParaView
   test(normalSimpleCase, offsetSimpleCase, expectedVolumeSimpleCase, tolerance);

   // test cases as in dissertation of Thomas Pohl, page 25
   const Vector3< real_t > normal(real_c(0.37), real_c(0.61), real_c(0.7));
   const std::vector< real_t > offsetList{ real_c(-0.52), real_c(-0.3), real_c(-0.17), real_c(0) };
   const std::vector< real_t > volumeList{ real_c(0.0347), real_c(0.1612), real_c(0.2887),
                                           real_c(0.5) }; // obtained with ParaView

   WALBERLA_ASSERT_EQUAL(offsetList.size(), volumeList.size());

   for (size_t i = 0; i != offsetList.size(); ++i)
   {
      test(normal, offsetList[i], volumeList[i], tolerance);
   }

   return EXIT_SUCCESS;
}
} // namespace CellFluidVolumeTest
} // namespace free_surface
} // namespace walberla

int main(int argc, char** argv) { return walberla::free_surface::CellFluidVolumeTest::main(argc, argv); }
