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
//! \file   DistanceCalculation.cpp
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include <mesa_pd/domain/BlockForestDomain.h>

#include <core/Environment.h>
#include <core/logging/Logging.h>

#include <iostream>

namespace walberla {

using namespace walberla::mesa_pd;

void main( int argc, char ** argv )
{
   Environment env(argc, argv);
   WALBERLA_UNUSED(env);
   mpi::MPIManager::instance()->useWorldComm();

   using namespace walberla::mesa_pd;

   WALBERLA_CHECK_FLOAT_EQUAL( sqDistanceLineToPoint(real_t(0.7), real_t(0.9), real_t(1.2)) + real_t(1),
                               real_t(0.04) + real_t(1) );
   WALBERLA_CHECK_FLOAT_EQUAL( sqDistanceLineToPoint(real_t(1.0), real_t(0.9), real_t(1.2)) + real_t(1),
                               real_t(0) + real_t(1) );
   WALBERLA_CHECK_FLOAT_EQUAL( sqDistanceLineToPoint(real_t(1.5), real_t(0.9), real_t(1.2)) + real_t(1),
                               real_t(0.09) + real_t(1) );

   WALBERLA_CHECK_FLOAT_EQUAL( sqDistancePointToAABB( Vec3(1,1,1),
                                                      math::AABB(real_t(2),real_t(2),real_t(2),real_t(3),real_t(3),real_t(3)) ) + real_t(1),
                               real_t(3) + real_t(1) );

   WALBERLA_CHECK_FLOAT_EQUAL( sqDistancePointToAABB( Vec3(real_t(0.5),real_t(0.5),real_t(0.5)),
                                                      math::AABB(real_t(2),real_t(2),real_t(2),real_t(3),real_t(3),real_t(3)) ) + real_t(1),
                               real_t(6.75) + real_t(1) );

   WALBERLA_CHECK_FLOAT_EQUAL( sqDistancePointToAABBPeriodic( Vec3(real_t(0.5),real_t(0.5),real_t(0.5)),
                                                              math::AABB(real_t(2),real_t(2),real_t(2),real_t(3),real_t(3),real_t(3)),
                                                              math::AABB(real_t(1),real_t(1),real_t(1),real_t(3),real_t(3),real_t(3)),
                                                              std::array< bool, 3 >{{false, false, false}}) + real_t(1),
                               real_t(6.75) + real_t(1) );

   WALBERLA_CHECK_FLOAT_EQUAL( sqDistancePointToAABBPeriodic( Vec3(real_t(0.5),real_t(0.5),real_t(0.5)),
                                                              math::AABB(real_t(2),real_t(2),real_t(2),real_t(3),real_t(3),real_t(3)),
                                                              math::AABB(real_t(1),real_t(1),real_t(1),real_t(3),real_t(3),real_t(3)),
                                                              std::array< bool, 3 >{{true, true, true}}) + real_t(1),
                               real_t(0) + real_t(1) );

   WALBERLA_CHECK_FLOAT_EQUAL( sqDistancePointToAABBPeriodic( Vec3(real_t(1.1),real_t(1.1),real_t(1.1)),
                                                              math::AABB(real_t(2),real_t(2),real_t(2),real_t(3),real_t(3),real_t(3)),
                                                              math::AABB(real_t(1),real_t(1),real_t(1),real_t(3),real_t(3),real_t(3)),
                                                              std::array< bool, 3 >{{true, true, true}}) + real_t(1),
                               real_t(0.03) + real_t(1) );
}

}

int main( int argc, char ** argv )
{
   walberla::main(argc, argv);

   return EXIT_SUCCESS;
}
