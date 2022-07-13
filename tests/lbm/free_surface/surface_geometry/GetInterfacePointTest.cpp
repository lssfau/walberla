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
//! \file GetInterfacePointTest.cpp
//! \ingroup lbm/free_surface/surface_geometry
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//! \brief Compute the interface point from normal and fill level, and compare with results obtained with ParaView.
//
//======================================================================================================================

#include "core/Environment.h"
#include "core/debug/TestSubsystem.h"

#include "lbm/free_surface/surface_geometry/CurvatureSweep.h"

namespace walberla
{
namespace free_surface
{
namespace GetInterfacePointTest
{
inline void test(const Vector3< real_t >& normal, real_t fillLevel, Vector3< real_t > expectedInterfacePoint,
                 real_t tolerance)
{
   const Vector3< real_t > interfacePoint = getInterfacePoint(normal, fillLevel);

   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(interfacePoint[0], expectedInterfacePoint[0], tolerance);
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(interfacePoint[1], expectedInterfacePoint[1], tolerance);
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(interfacePoint[2], expectedInterfacePoint[2], tolerance);
}

int main(int argc, char** argv)
{
   debug::enterTestMode();
   Environment env(argc, argv);

   // allowed deviation from the expected result
   real_t tolerance = real_c(1e-2);

   // tests identical to those in CellFluidVolumeTest.cpp, results obtained with ParaView
   Vector3< real_t > normal(real_c(0), real_c(0.5), real_c(0.866));
   const real_t offset                      = real_c(-0.1);
   const real_t fillLevel                   = real_c(0.3845);
   Vector3< real_t > expectedInterfacePoint = Vector3< real_t >(real_c(0.5)) + offset * normal;
   test(normal, fillLevel, expectedInterfacePoint, tolerance);

   normal = Vector3< real_t >(real_c(0.37), real_c(0.61), real_c(0.7));
   const std::vector< real_t > offsetList{ real_c(-0.52), real_c(-0.3), real_c(-0.17), real_c(0) };
   const std::vector< real_t > fillLevelList{ real_c(0.0347), real_c(0.1612), real_c(0.2887), real_c(0.5) };

   WALBERLA_ASSERT_EQUAL(offsetList.size(), fillLevelList.size());
   for (size_t i = 0; i != offsetList.size(); ++i)
   {
      expectedInterfacePoint = Vector3< real_t >(real_c(0.5)) + offset * normal;
      test(normal, fillLevel, expectedInterfacePoint, tolerance);
   }

   return EXIT_SUCCESS;
}
} // namespace GetInterfacePointTest
} // namespace free_surface
} // namespace walberla

int main(int argc, char** argv) { return walberla::free_surface::GetInterfacePointTest::main(argc, argv); }
