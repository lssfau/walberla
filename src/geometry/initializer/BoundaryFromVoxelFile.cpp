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
//! \file BoundaryFromVoxelFile.cpp
//! \ingroup geometry
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//
//======================================================================================================================

#include "BoundaryFromVoxelFile.h"
#include "core/mpi/BufferDataTypeExtensions.h"
#include "core/mpi/BufferSystem.h"


namespace walberla {
namespace geometry {
namespace initializer {


bool IBlockIDPtrCompare::operator()( const IBlockID * lhs, const IBlockID * rhs ) const
{
   return *lhs < *rhs;
}

CellIntervalDataMap readCellIntervalsOnRoot( const std::string & geometryFile,
                                             const Cell & offset, const CellIntervalMap & cellIntervals )
{
   CellIntervalDataMap result;

   mpi::BufferSystem allToRootBs( MPI_COMM_WORLD );
   mpi::BufferSystem rootToAllBs( MPI_COMM_WORLD );

   WALBERLA_ROOT_SECTION() {
      allToRootBs.setReceiverInfo( mpi::BufferSystem::allRanksButRoot(), true );
      rootToAllBs.setReceiverInfo( std::set<int>(), true );
   }
   else {
      allToRootBs.setReceiverInfo( std::set<int>(), true );
      rootToAllBs.setReceiverInfo( mpi::BufferSystem::onlyRoot(), true );
   }


   // send regions of geometry file to read to root
   WALBERLA_NON_ROOT_SECTION()
   {
      auto & sendBuffer = allToRootBs.sendBuffer(0);

      sendBuffer << uint8_t(0); // dummy byte to make sure a message is sent even if cellIntervals is empty

      for( auto ciIt = cellIntervals.begin(); ciIt != cellIntervals.end(); ++ciIt )
         sendBuffer << ciIt->second;

   }


   allToRootBs.sendAll();

   VoxelFileReader<uint8_t> reader( geometryFile );
   std::vector<uint8_t> data;

   WALBERLA_ROOT_SECTION()
   {
      // Read root cell intervals from geometry file
      for( auto ciIt = cellIntervals.begin(); ciIt != cellIntervals.end(); ++ciIt )
      {
         CellInterval shifted = ciIt->second;
         shifted.shift( -offset );
         reader.read(shifted, data);
         result[ ciIt->first ] = std::make_pair( ciIt->second, data );
      }
   }


   // Receive and read cell intervals from other processes, send data back
   for( auto it = allToRootBs.begin(); it != allToRootBs.end(); ++it )
   {
      // only root should receive messages
      WALBERLA_ASSERT_EQUAL( MPIManager::instance()->worldRank(), 0 );

      int rank = it.rank();
      auto & recvBuffer = it.buffer();
      auto & sendBuffer = rootToAllBs.sendBuffer( rank );

      uint8_t dummy;
      recvBuffer >> dummy;
      sendBuffer << uint8_t(0); // dummy byte to make sure a message is sent even if recvBuffer is empty

      while( !recvBuffer.isEmpty() )
      {
         CellInterval ci;
         recvBuffer >> ci;
         ci.shift(-offset);
         reader.read(ci, data);
         sendBuffer << data;
      }
   }


   rootToAllBs.sendAll();


   // receive read data from root
   for( auto it = rootToAllBs.begin(); it != rootToAllBs.end(); ++it )
   {
      WALBERLA_ASSERT_UNEQUAL( MPIManager::instance()->worldRank(), 0 );

      auto & recvBuffer = it.buffer();

      uint8_t dummy;
      recvBuffer >> dummy;

      for( auto ciIt = cellIntervals.begin(); ciIt != cellIntervals.end(); ++ciIt )
      {
         WALBERLA_ASSERT( !recvBuffer.isEmpty() );
         recvBuffer >> data;
         result[ ciIt->first ] = std::make_pair( ciIt->second, data );
      }

      WALBERLA_ASSERT( recvBuffer.isEmpty() );
   }

   return result;
}

CellVector findCellsWithFlag( const CellInterval & cellInterval, const std::vector<uint8_t> & data, uint8_t flag )
{
   WALBERLA_ASSERT_EQUAL( cellInterval.numCells(), data.size() );

   CellVector result;

   auto cellIt = cellInterval.begin();
   auto dataIt = data.begin();

   for( ; cellIt != cellInterval.end(); ++cellIt, ++dataIt )
   {
      if( *dataIt == flag )
         result.push_back( *cellIt );
   }

   return result;
}


} // namespace initializer
} // namespace geometry
} // namespace walberla
