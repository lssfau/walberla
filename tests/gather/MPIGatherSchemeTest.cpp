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
//! \file MPIGatherSchemeTest.cpp
//! \ingroup gather
//! \author Martin Bauer <martin.bauer@fau.de>
//! \brief Basic test for MPIGatherScheme
//!
//! In this test half of the processes ( that with ranks >= nrProcs / 2 ) take place in a
//! gather operation. Every process sends its own global rank r, r times
//! such that the messages have different size.
//
//======================================================================================================================

#include "gather/GatherPackInfo.h"
#include "gather/MPIGatherScheme.h"

#include "blockforest/Initialization.h"

#include "core/Environment.h"
#include "core/debug/TestSubsystem.h"
#include "core/mpi/MPIManager.h"


namespace walberla {

class TestGatherPackInfo : public gather::GatherPackInfo
{
public:
   TestGatherPackInfo() = default;
   ~TestGatherPackInfo() override = default;


   void packData  ( const IBlock *, mpi::GenericSendBuffer<unsigned char> & outBuffer ) override
   {
      auto mpi = MPIManager::instance();

      if ( mpi->rank() >= mpi->numProcesses() / 2 )
      {
         outBuffer << mpi->rank();
         for( int i=0; i <  mpi->rank(); ++i )
            outBuffer << i;
      }
   }


   void unpackData( mpi::GenericRecvBuffer<unsigned char> & buffer ) override
   {
      int nrOfItems;
      buffer >> nrOfItems;
      for( int i=0 ;i < nrOfItems; ++i ) {
         int curItem= -1;
         buffer >> curItem;
         WALBERLA_CHECK_EQUAL( curItem, i );
      }
      receivedRanks_.insert( nrOfItems );
   }

   void gatherFinished() override
   {
      auto mpi = MPIManager::instance();

      int startRank = mpi->numProcesses() /2;

      for ( int i=startRank; i < mpi->numProcesses(); ++i )
      {
         auto found = receivedRanks_.find( i );
         WALBERLA_CHECK( found != receivedRanks_.end() );
      }

      receivedRanks_.clear();
   }

private:
   std::set<int> receivedRanks_;
};
}// namespace walberla



int main( int argc, char ** argv )
{
   walberla::Environment env( argc, argv );

   walberla::debug::enterTestMode();

   walberla::uint_t numProcesses = walberla::uint_c( walberla::MPIManager::instance()->numProcesses() );
   WALBERLA_CHECK( numProcesses >= 2 );

   using walberla::blockforest::createUniformBlockGrid;
   auto blocks = createUniformBlockGrid( numProcesses, 1, 1,
                                         1, 1, 1,
                                         1,
                                         true,                // one block per process
                                         false, false, false, // periodicity
                                         false );             // do NOT keep global information

   walberla::gather::MPIGatherScheme gatherScheme( blocks->getBlockStorage(), 0 );
   gatherScheme.addPackInfo( walberla::make_shared<walberla::TestGatherPackInfo>() );

   for(int i=0; i<3; ++i )
      gatherScheme();

   return 0;
}
