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
//! \file TimingPoolTest.cpp
//! \ingroup core
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================


#include "core/debug/TestSubsystem.h"
#include "core/timing/StaticPolicy.h"
#include "core/timing/Timer.h"

namespace walberla {

void mergeTest() 
{
   timing::Timer<timing::StaticPolicy> t0, t1, t2, t_all;

   t0.start();
   t_all.start();
   timing::StaticPolicy::addTime(1.23);
   t0.end();
   t_all.end();

   t1.start();
   t_all.start();
   timing::StaticPolicy::addTime(1.23);
   t1.end();
   t_all.end();

   t2.start();
   t_all.start();
   timing::StaticPolicy::addTime(1.23);
   t2.end();
   t_all.end();

   t0.merge( t1 );
   t0.merge( t2 );

   WALBERLA_CHECK_FLOAT_EQUAL( t0.average(), t_all.average() );
   WALBERLA_CHECK_FLOAT_EQUAL( t1.average(), t_all.average() );
   WALBERLA_CHECK_FLOAT_EQUAL( t2.average(), t_all.average() );

   WALBERLA_CHECK_FLOAT_EQUAL( t0.sumOfSquares(), t_all.sumOfSquares() );
   WALBERLA_CHECK_FLOAT_EQUAL( t1.sumOfSquares(), t2.sumOfSquares() );

   WALBERLA_CHECK_FLOAT_EQUAL( t1.total(), t2.total() );
   WALBERLA_CHECK_FLOAT_EQUAL( t0.total(), t_all.total() );
   WALBERLA_CHECK_EQUAL( t0.getCounter(), t_all.getCounter() );
   WALBERLA_CHECK_EQUAL( t1.getCounter(), t2.getCounter() );

   WALBERLA_CHECK_FLOAT_EQUAL( t0.max(), t_all.max() );
   WALBERLA_CHECK_FLOAT_EQUAL( t0.min(), t_all.min() );
   WALBERLA_CHECK_FLOAT_EQUAL( t1.max(), t_all.max() );
   WALBERLA_CHECK_FLOAT_EQUAL( t1.min(), t_all.min() );
   WALBERLA_CHECK_FLOAT_EQUAL( t2.max(), t_all.max() );
   WALBERLA_CHECK_FLOAT_EQUAL( t2.min(), t_all.min() );

   WALBERLA_CHECK_FLOAT_EQUAL( t0.variance(), t_all.variance() );
   WALBERLA_CHECK_FLOAT_EQUAL( t1.variance(), t_all.variance() );
   WALBERLA_CHECK_FLOAT_EQUAL( t2.variance(), t_all.variance() );
}

}

int main()
{
  walberla::debug::enterTestMode();
  walberla::mergeTest();
}
