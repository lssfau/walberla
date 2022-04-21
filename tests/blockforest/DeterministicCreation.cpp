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
//! \file DeterministicCreation.cpp
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include "blockforest/Initialization.h"
#include "core/DataTypes.h"
#include "core/debug/TestSubsystem.h"
#include "core/logging/Logging.h"
#include "core/mpi/Environment.h"

namespace walberla{

void TwoBlockForestsTest()
{
   Vector3<uint_t> cells   (10,10,10);
   Vector3<uint_t> blocks  (2,2,2);
   Vector3<uint_t> procs   (2,2,2);
   Vector3<bool>   periodic(true, true, true);

   uint_t peDomainFactor = uint_t(20);

   AABB lbmDomainAABB( real_c(0), real_c(0), real_c(0), real_c(cells[0]), real_c(cells[1]), real_c(cells[2]) );
   AABB peDomainAABB ( real_c(0), real_c(0), real_c(0), real_c(cells[0]), real_c(cells[1]), real_c(cells[2]*peDomainFactor) );

   mpi::MPIManager::instance()->resetMPI();

   auto lbmBlocks = blockforest::createUniformBlockGrid(
      lbmDomainAABB,
      blocks[0],          blocks[1],          blocks[2], // number of blocks in x/y/z direction
      cells[0]/blocks[0], cells[1]/blocks[1], cells[2]/blocks[2], // number of cells per block in x/y/z direction
      procs[0],           procs[1],           procs[2], // number of processes in x/y/z direction
      periodic[0],        periodic[1],        periodic[2], true );       // periodic in every direction

   mpi::MPIManager::instance()->resetMPI();

   auto peBlocks = blockforest::createUniformBlockGrid(
       peDomainAABB,
       blocks[0], blocks[1], blocks[2],           // number of blocks in x/y/z direction
       cells[0] / blocks[0], cells[1] / blocks[1], peDomainFactor * cells[2] / blocks[2],  // number of cells per block in x/y/z direction
       procs[0], procs[1], procs[2],            // number of processes in x/y/z direction
       periodic[0], periodic[1], periodic[2], true);       // periodic in every direction

   WALBERLA_LOG_DEVEL(lbmBlocks->begin()->getAABB())
   WALBERLA_LOG_DEVEL(peBlocks->begin()->getAABB())

   WALBERLA_CHECK_EQUAL( lbmBlocks->begin()->getId(), peBlocks->begin()->getId() )
}

int main( int argc, char** argv )
{
   debug::enterTestMode();

   mpi::Environment env( argc, argv );
   //MPIManager::instance()->useWorldComm();

   TwoBlockForestsTest();

   return EXIT_SUCCESS;
}
}

int main( int argc, char** argv )
{
   return walberla::main(argc,argv);
}