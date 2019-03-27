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
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#include "core/debug/TestSubsystem.h"
#include "core/logging/Logging.h"
#include "core/math/Constants.h"
#include "core/mpi/MPIManager.h"
#include "core/timing/Timer.h"
#include "core/timing/TimingPool.h"

#include <iostream>


namespace walberla {

using std::cout;
using std::endl;

void simpleTiming()
{
   WcTimingPool pool;

   pool["All Together"];

   double sum=0;
   for(int i=0 ;i<30; i++)
   {
      pool["In the loop"].start();
      for( int j=0; j< 10000; ++j)
         sum += 1;
      pool["In the loop"].end();
   }

   pool["All Together"].end();

   cout << pool << endl;
}

void scopedTimer()
{
   WcTimingPool pool;
   {
      pool["normal timer"].start();
      auto scopedTimer = pool.getScopeTimer( "scope timer" );

      double sum = 0.0;
      for( double d = 0.0; d < math::M_PI; d += 0.00001 )
      {
         sum += std::atan( std::tan( d ) );
         sum += std::asin( std::sin( d ) );
         sum += std::acos( std::cos( d ) );
      }

      pool["normal timer"].end();
   }

   WALBERLA_CHECK_LESS( pool["normal timer"].last() - pool["scope timer"].last(), 1e-4 ); // less than 0.1 ms
}


void reduction()
{
   WALBERLA_CHECK( MPIManager::instance()->numProcesses() == 3 );

   int myRank = MPIManager::instance()->worldRank();

   WcTimingPool pool;

   if ( myRank == 0 )
   {
      // measured 2,4,6
      pool["test"] = WcTimer( 3, 2, 6, 12, 56 );
   }
   else if ( myRank == 1 )
   {
      // measured 1,2,4
      pool["test"] = WcTimer( 3, 1, 4, 7, 21 );
   }
   else if ( myRank == 2 )
   {
      // measured 1,2,3,4
      pool["test"] = WcTimer( 4, 1, 4, 10, 30 );
   }

   WcTimer timerBackup = pool["test"];

   // Test minimum reduction
   auto red = pool.getReduced( timing::REDUCE_MIN, 0 );
   WALBERLA_ROOT_SECTION() {
      WcTimer & t = (*red)["test"];
      WALBERLA_CHECK_FLOAT_EQUAL( t.min(), 1.0 );
      WALBERLA_CHECK_FLOAT_EQUAL( t.max(), 2.0 );
      WALBERLA_CHECK_FLOAT_EQUAL( t.average(), (1 + 1 + 2)/3.0 );
      WALBERLA_CHECK_FLOAT_EQUAL( t.variance(), (2*2+1+1)/3.0 - t.average()*t.average() );
      WALBERLA_CHECK_EQUAL( t.getCounter(), 3 );
      pool["test"] = timerBackup;
   }

   // Test maximum reduction
   red = pool.getReduced( timing::REDUCE_MAX, 0 );
   WALBERLA_ROOT_SECTION() {
      WcTimer & t = (*red)["test"];
      WALBERLA_CHECK_FLOAT_EQUAL( t.min(), 4.0 );
      WALBERLA_CHECK_FLOAT_EQUAL( t.max(), 6.0 );
      WALBERLA_CHECK_FLOAT_EQUAL( t.average(), (4 + 4 + 6)/3.0 );
      WALBERLA_CHECK_FLOAT_EQUAL( t.variance(), (16+16+36)/3.0 - t.average()*t.average() );
      WALBERLA_CHECK_EQUAL( t.getCounter(), 3 );
      pool["test"] = timerBackup;
   }

   // Test complete reduction
   red = pool.getReduced( timing::REDUCE_TOTAL, 0 );
   WALBERLA_ROOT_SECTION() {
      WcTimer & t = (*red)["test"];
      WALBERLA_CHECK_FLOAT_EQUAL( t.min(), 1.0 );
      WALBERLA_CHECK_FLOAT_EQUAL( t.max(), 6.0 );
      WALBERLA_CHECK_EQUAL( t.getCounter(), 10 );
      pool["test"] = timerBackup;
   }


   red.reset();
   red = pool.getReduced(  timing::REDUCE_TOTAL, -1 );
   WALBERLA_CHECK_NOT_NULLPTR( red );
   WALBERLA_CRITICAL_SECTION_START
   cout << *red << endl;
   WALBERLA_CRITICAL_SECTION_END


   bool enablePrinting = false;

   if(enablePrinting)
   {
      WALBERLA_CRITICAL_SECTION_START
      cout << "Pool on process " << myRank << endl;
      cout << pool << endl;
      WALBERLA_CRITICAL_SECTION_END

      red = pool.getReduced( timing::REDUCE_TOTAL, 0 );
      WALBERLA_ROOT_SECTION() {
         cout << "Reduced System" << endl << *red << endl;
      }
   }
}


void merging()
{
   WcTimingPool pool1;
   WcTimingPool pool2;

   WcTimer timer1 = WcTimer( 3, 2, 6, 12, 56 ); // measured 2,4,6
   WcTimer timer2 = WcTimer( 3, 1, 4,  7, 21 ); // measured 1,2,4

   pool1["test"] = timer1;
   pool2["test"] = timer2;

   pool1.merge( pool2, false );
   WALBERLA_CHECK_EQUAL ( pool1["test"].getCounter(), 3 );

   pool1.merge( pool2, true );
   WALBERLA_CHECK_EQUAL ( pool1["test"].getCounter(), 6 );
}
}// namespace walberla


int main( int argc, char ** argv )
{
   walberla::debug::enterTestMode();
   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );

   walberla::simpleTiming();
   walberla::scopedTimer();
   walberla::reduction();
   walberla::merging();

   return 0;
}




