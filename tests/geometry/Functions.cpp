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
//! \file Functions.cpp
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include <core/DataTypes.h>
#include <core/debug/TestSubsystem.h>
#include <core/logging/Logging.h>
#include <core/mpi/MPIManager.h>
#include <geometry/GeometricalFunctions.h>

using namespace walberla;
using namespace walberla::geometry;

using Vec3 = math::Vector3<real_t>;

int main( int argc, char** argv )
{
    debug::enterTestMode();

    MPIManager::instance()->initializeMPI( &argc, &argv );

    Vec3 a,b;
    geometry::getClosestLineSegmentPoints(Vec3(real_c(7.41657),real_c( 9.82461),real_c(17.28197)),
                                          Vec3(real_c(7.21657),real_c( 9.82461),real_c(17.28197)),
                                          Vec3(real_c(7.36474),real_c(10.20490),real_c(17.45072)),
                                          Vec3(real_c(7.16474),real_c(10.20490),real_c(17.45072)), a, b);
    WALBERLA_LOG_DEVEL(a << "\n" << b);
    // might fail due to numerically instable algorithm
    geometry::getClosestLineSegmentPoints(Vec3(real_c(7.416573195643714),real_c(9.82461330720222),real_c(17.28197302209214)),
                                          Vec3(real_c(7.216573195643715),real_c(9.82461330720222),real_c(17.28197302209214)),
                                          Vec3(real_c(7.364742740420613),real_c(10.2049088375177),real_c(17.45072353608904)),
                                          Vec3(real_c(7.164742740420614),real_c(10.2049088375177),real_c(17.45072353608904)), a, b);
    WALBERLA_LOG_DEVEL(a << "\n" << b);

    return EXIT_SUCCESS;
}
