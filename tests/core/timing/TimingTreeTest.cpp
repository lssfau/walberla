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
//! \file TimingTreeTest.cpp
//! \ingroup core
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#include "core/debug/TestSubsystem.h"
#include "core/logging/Logging.h"
#include "core/mpi/Environment.h"
#include "core/timing/StaticPolicy.h"
#include "core/timing/TimingTree.h"

#include <iostream>

void mssleep(unsigned int ms)
{
   walberla::timing::StaticPolicy::addTime( ms * 1e-3);
}

int main( int argc, char ** argv )
{
   walberla::debug::enterTestMode();

   walberla::mpi::Environment mpiEnv(argc, argv);
   WALBERLA_UNUSED( mpiEnv );

   const unsigned int rank = static_cast<unsigned int> ( walberla::MPIManager::instance()->worldRank() );

   walberla::timing::TimingTree<walberla::timing::StaticPolicy> tt;
   walberla::timing::StaticPolicy::setTime(0);

   tt.start("A");
   mssleep(100 * rank);
   tt.start("AA");
   mssleep(100 * rank);
   tt.stop("AA");
   tt.start("AB");
   tt.start("ABA");
   mssleep(100 * rank);
   tt.stop("ABA");
   tt.start("ABB");
   mssleep(100 * rank);
   tt.stop("ABB");
   tt.stop("AB");
   tt.start("AC");
   mssleep(100 * rank);
   tt.start("ACA");
   mssleep(100 * rank);
   tt.stop("ACA");
   tt.stop("AC");
   tt.stop("A");

   WALBERLA_ASSERT(tt.timerExists("A"));
   WALBERLA_ASSERT(tt.timerExists("A.AA"));
   WALBERLA_ASSERT(tt.timerExists("A.AB.ABA"));
   WALBERLA_ASSERT(tt.timerExists("A.AB.ABB"));
   WALBERLA_ASSERT(tt.timerExists("A.AC.ACA"));
   WALBERLA_ASSERT(!tt.timerExists("AAC"));
   WALBERLA_ASSERT(!tt.timerExists("A.AA.C"));

   // check copy constructor
   walberla::timing::TimingTree<walberla::timing::StaticPolicy> tt2(tt);
   // check assignment operator
   walberla::timing::TimingTree<walberla::timing::StaticPolicy> tt3;
   tt3 = tt;

   WALBERLA_ASSERT(tt2.timerExists("A"));
   WALBERLA_ASSERT(tt2.timerExists("A.AA"));
   WALBERLA_ASSERT(tt2.timerExists("A.AB.ABA"));
   WALBERLA_ASSERT(tt2.timerExists("A.AB.ABB"));
   WALBERLA_ASSERT(tt2.timerExists("A.AC.ACA"));
   WALBERLA_ASSERT(!tt2.timerExists("AAC"));
   WALBERLA_ASSERT(!tt2.timerExists("A.AA.C"));

   WALBERLA_ASSERT(tt3.timerExists("A"));
   WALBERLA_ASSERT(tt3.timerExists("A.AA"));
   WALBERLA_ASSERT(tt3.timerExists("A.AB.ABA"));
   WALBERLA_ASSERT(tt3.timerExists("A.AB.ABB"));
   WALBERLA_ASSERT(tt3.timerExists("A.AC.ACA"));
   WALBERLA_ASSERT(!tt3.timerExists("AAC"));
   WALBERLA_ASSERT(!tt3.timerExists("A.AA.C"));

   tt2 = tt.getReduced( walberla::timing::REDUCE_TOTAL, 0 );
   tt2 = tt.getReduced( walberla::timing::REDUCE_TOTAL, 1 );

   tt2 = tt.getReduced( walberla::timing::REDUCE_MIN, -1 );
   {
   const auto& data = tt2.getRawData();
   WALBERLA_CHECK_FLOAT_EQUAL( data.tree_.at("A").timer_.total(), (1.8) );
   WALBERLA_CHECK_FLOAT_EQUAL( data.tree_.at("A").tree_.at("AA").timer_.total(), (0.3) );

   WALBERLA_CHECK_FLOAT_EQUAL( tt["A.AB.ABB"].total(), (0.100 * rank), "total time: " << tt["A.AB.ABB"].total() );
   }

   tt2 = tt.getReduced( walberla::timing::REDUCE_MAX, -1 );
   {
   const auto& data = tt2.getRawData();
   WALBERLA_CHECK_FLOAT_EQUAL( data.tree_.at("A").timer_.total(), (1.8) );
   WALBERLA_CHECK_FLOAT_EQUAL( data.tree_.at("A").tree_.at("AA").timer_.total(), (0.3) );

   WALBERLA_CHECK_FLOAT_EQUAL( tt["A.AB.ABB"].total(), (0.100 * rank), "total time: " << tt["A.AB.ABB"].total() );
   }

   tt2 = tt.getReduced( walberla::timing::REDUCE_AVG, -1 );
   {
   const auto& data = tt2.getRawData();
   WALBERLA_CHECK_FLOAT_EQUAL( data.tree_.at("A").timer_.total(), (1.8) );
   WALBERLA_CHECK_FLOAT_EQUAL( data.tree_.at("A").tree_.at("AA").timer_.total(), (0.3) );

   WALBERLA_CHECK_FLOAT_EQUAL( tt["A.AB.ABB"].total(), (0.100 * rank), "total time: " << tt["A.AB.ABB"].total() );
   }

   tt2 = tt.getReduced( walberla::timing::REDUCE_TOTAL, -1 );
   {
   const auto& data = tt2.getRawData();
   WALBERLA_CHECK_FLOAT_EQUAL( data.tree_.at("A").timer_.total(), (1.8) );
   WALBERLA_CHECK_FLOAT_EQUAL( data.tree_.at("A").tree_.at("AA").timer_.total(), (0.3) );

   WALBERLA_CHECK_FLOAT_EQUAL( tt["A.AB.ABB"].total(), (0.100 * rank), "total time: " << tt["A.AB.ABB"].total() );
   }

   WALBERLA_ROOT_SECTION()
   {
//      std::cout << tt;
//      std::cout << tt2;
//      std::cout << tt3;
   }

   return 0;
}




