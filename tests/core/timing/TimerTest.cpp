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
//
//======================================================================================================================


#include "core/debug/TestSubsystem.h"
#include "core/math/Constants.h"
#include "core/timing/Timer.h"

#include <cmath>

namespace walberla {

static double burnTime()
{
   double sum = 0.0;
   for( double d = 0.0; d < math::M_PI; d += 0.000001 )
   {
      sum += std::atan( std::tan( d ) );
      sum += std::asin( std::sin( d ) );
      sum += std::acos( std::cos( d ) );
   }

   return sum;
}

void mergeTest() 
{
   WcTimer t0, t1, t2, t_all;

   t0.start();
   t_all.start();
   burnTime();
   t0.end();
   t_all.end();

   t1.start();
   t_all.start();
   burnTime();
   t1.end();
   t_all.end();

   t2.start();
   t_all.start();
   burnTime();
   t2.end();
   t_all.end();

   t0.merge( t1 );
   t0.merge( t2 );

   WALBERLA_CHECK_LESS( t0.average() - t_all.average(), 1e-4 );
   WALBERLA_CHECK_LESS( t0.sumOfSquares() - t_all.sumOfSquares(), 3e-2 );
   WALBERLA_CHECK_LESS( t0.total() - t_all.total(), 3e-4 );
   WALBERLA_CHECK_EQUAL( t0.getCounter(), t_all.getCounter() );
   WALBERLA_CHECK_LESS( t0.max() - t_all.max(), 1e-4 );
   WALBERLA_CHECK_LESS( t0.min() - t_all.min(), 1e-4 );
   WALBERLA_CHECK_LESS( t0.variance() - t_all.variance(), 1e-4 );
}

}

int main()
{
  walberla::debug::enterTestMode();
  walberla::mergeTest();
}
