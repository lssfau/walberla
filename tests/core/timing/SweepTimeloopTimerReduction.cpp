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
//! \file SweepTimeloopTimerReduction.cpp
//! \ingroup core
//! \author Martin Bauer <martin.bauer@fau.de>
//!
//! This test checks if a timing pool measured by a SweepTimeloop can be reduced
//! even if on process does not have any block. See ticket #289
//
//======================================================================================================================

#include "core/Environment.h"
#include "core/debug/TestSubsystem.h"
#include "core/logging/Logging.h"
#include "core/timing/TimingPool.h"

#include "blockforest/Initialization.h"
#include "blockforest/SetupBlockForest.h"
#include "blockforest/StructuredBlockForest.h"

#include "domain_decomposition/IBlock.h"

#include "timeloop/SweepTimeloop.h"


void dummySweep ( walberla::IBlock * block )
{
   WALBERLA_LOG_DEVEL("DummySweep on block " << block );
}

int main( int argc, char ** argv )
{
   walberla::debug::enterTestMode();
   walberla::Environment walberlaEnv( argc, argv );

   // 9 processes and 8 blocks:   one process does not have a block
   WALBERLA_CHECK_EQUAL( walberla::MPIManager::instance()->numProcesses(), 9 );

   using namespace walberla::blockforest;

   SetupBlockForest sforest;
   sforest.addWorkloadMemorySUIDAssignmentFunction( uniformWorkloadAndMemoryAssignment );
   sforest.init( walberla::AABB(0,0,0, 2,2,2), 2, 2, 2, false, false, false );

   // calculate process distribution
   const memory_t memoryLimit = 1;

   GlobalLoadBalancing::MetisConfiguration< SetupBlock > metisConfig( false, false, uniformFacesDominantCommunication );
   sforest.calculateProcessDistribution_Default( walberla::uint_c( walberla::MPIManager::instance()->numProcesses() ), memoryLimit, "hilbert", 10, false, metisConfig );
   walberla::MPIManager::instance()->useWorldComm();

   // create StructuredBlockForest (encapsulates a newly created BlockForest)
   auto bf = std::make_shared< BlockForest >( walberla::uint_c( walberla::MPIManager::instance()->rank() ), sforest, true );
   auto sbf = std::make_shared< StructuredBlockForest >( bf, 10, 10, 10 );
   sbf->createCellBoundingBoxes();


   walberla::SweepTimeloop tl( sbf, 1 );
   tl.add() << walberla::Sweep( dummySweep, "DummySweep" );

   walberla::WcTimingPool timingPool;
   tl.run( timingPool );

   timingPool.logResultOnRoot();

   return 0;
}
