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
//! \file ParallelTimerTest.cpp
//! \author Sebastian Eibl <sebastian.eibl@web.de>
//
//======================================================================================================================


#include "core/debug/TestSubsystem.h"
#include "core/Environment.h"
#include "core/math/Constants.h"
#include "core/timing/StaticPolicy.h"
#include "core/timing/Timer.h"

#include <cmath>

namespace walberla {

void reductionTest()
{
   using namespace timing;
   StaticPolicy::setTime(0.0);
   Timer<StaticPolicy> t0;

   t0.start();
   StaticPolicy::setTime(2.0);
   t0.end();

   t0.start();
   StaticPolicy::setTime(6.0);
   t0.end();

   WALBERLA_CHECK_FLOAT_EQUAL( t0.min(), 2.0 );
   WALBERLA_CHECK_FLOAT_EQUAL( t0.max(), 4.0 );
   WALBERLA_CHECK_FLOAT_EQUAL( t0.average(), 3.0 );
   WALBERLA_CHECK_EQUAL( t0.getCounter(), 2 );

   auto timer_reduced = getReduced(t0, REDUCE_TOTAL, -1); //allreduce
   WALBERLA_CHECK_FLOAT_EQUAL( timer_reduced->min(), 6.0 );
   WALBERLA_CHECK_FLOAT_EQUAL( timer_reduced->max(), 6.0 );
   WALBERLA_CHECK_FLOAT_EQUAL( timer_reduced->average(), 6.0 );
   WALBERLA_CHECK_EQUAL( timer_reduced->getCounter(), 2 );

   timer_reduced = getReduced(t0, REDUCE_TOTAL, 0);
   if (mpi::MPIManager::instance()->worldRank() == 0)
   {
      WALBERLA_CHECK_FLOAT_EQUAL( timer_reduced->min(), 6.0 );
      WALBERLA_CHECK_FLOAT_EQUAL( timer_reduced->max(), 6.0 );
      WALBERLA_CHECK_FLOAT_EQUAL( timer_reduced->average(), 6.0 );
      WALBERLA_CHECK_EQUAL( timer_reduced->getCounter(), 2 );
   } else
   {
      WALBERLA_CHECK_EQUAL( timer_reduced, nullptr );
   }
}

}

int main( int argc, char ** argv )
{
   walberla::debug::enterTestMode();
   walberla::Environment env(argc, argv);
   WALBERLA_UNUSED(env);

   walberla::reductionTest();
}
