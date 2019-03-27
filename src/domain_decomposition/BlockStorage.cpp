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
//! \file BlockStorage.cpp
//! \ingroup domain_decomposition
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#include "BlockStorage.h"
#include "MapPointToPeriodicDomain.h"
#include "PeriodicIntersect.h"

#include "core/Abort.h"
#include "core/debug/CheckFunctions.h"
#include "core/debug/Debug.h"
#include "core/logging/Logging.h"
#include "core/mpi/MPIIO.h"
#include "core/uid/GlobalState.h"
#include "core/mpi/Gather.h"


namespace walberla {
namespace domain_decomposition {



//**********************************************************************************************************************
/*!
*   For documentation, see documentation of free function
*   'void mapPointToPeriodicDomain( const std::array< bool, 3 > & periodic, const AABB & domain,
*                                   real_t & x, real_t & y, real_t & z )'
*/
//**********************************************************************************************************************
void BlockStorage::mapToPeriodicDomain( real_t & x, real_t & y, real_t & z ) const
{
   std::array< bool, 3 > periodic;
   periodic[0] = periodic_[0];
   periodic[1] = periodic_[1];
   periodic[2] = periodic_[2];

   domain_decomposition::mapPointToPeriodicDomain( periodic, domain_, x, y, z );
}

//**********************************************************************************************************************
/*!
*   For documentation, see documentation of free function
*   'bool periodicIntersect( const std::array< bool, 3 > & periodic, const math::AABB & domain, const math::AABB & box1, const math::AABB & box2 )'
*/
//**********************************************************************************************************************
bool BlockStorage::periodicIntersect( const math::AABB & box1, const math::AABB & box2 ) const
{
   std::array< bool, 3 > periodic;
   periodic[0] = periodic_[0];
   periodic[1] = periodic_[1];
   periodic[2] = periodic_[2];

   return domain_decomposition::periodicIntersect( periodic, domain_, box1, box2);
}

//**********************************************************************************************************************
/*!
*   For documentation, see documentation of free function
*   'bool periodicIntersect( const std::array< bool, 3 > & periodic, const math::AABB & domain, const math::AABB & box1, const math::AABB & box2 )'
*/
//**********************************************************************************************************************
bool BlockStorage::periodicIntersect( const math::AABB & box1, const math::AABB & box2, const real_t dx ) const
{
   std::array< bool, 3 > periodic;
   periodic[0] = periodic_[0];
   periodic[1] = periodic_[1];
   periodic[2] = periodic_[2];

   return domain_decomposition::periodicIntersect( periodic, domain_, box1, box2, dx);
}



//**********************************************************************************************************************
/*!
*   After registering one/multiple block data handling objects, this function is called for actually adding
*   data to the blocks.
*/
//**********************************************************************************************************************
BlockDataID BlockStorage::addBlockData( const internal::SelectableBlockDataHandlingWrapper & dataHandling,
                                        const std::string& identifier )
{
   WALBERLA_LOG_PROGRESS( "Adding block data (\"" << identifier << "\")" );

   BlockDataID id( blockDataItem_.size() );
   internal::BlockDataItem item( id, identifier, dataHandling );
   blockDataItem_.push_back( item );

   for( auto block = begin(); block != end(); ++block )
   {
      auto dh = item.getDataHandling( block.get() );
      if( dh )
         block->addData( id, dh->initialize( block.get() ) );
      else
         block->addData( id, nullptr );
   }

   return id;
}



namespace internal {
   bool sortBlocksByID( IBlock * lhs, IBlock * rhs ) { return lhs->getId() < rhs->getId(); }
}

//**********************************************************************************************************************
/*!
*   After registering one/multiple block data handling objects, this function is called for actually adding data to
*   the blocks by loading the data from file.
*   ATTENTION: If the data was previously saved to file during a parallel run, the current run must also be in parallel
*              with the exact same block structure (same blocks, same block states, etc.) and the same number of
*              processes. Analogously, data saved during a serial run can only be loaded during another serial run.
*/
//**********************************************************************************************************************
BlockDataID BlockStorage::loadBlockData( const std::string & file, const internal::SelectableBlockDataHandlingWrapper & dataHandling,
                                         const std::string & identifier )
{
   WALBERLA_LOG_PROGRESS( "Adding block data (\"" << identifier << "\"), loading data from file \"" << file << "\" ..." );
   
   BlockDataID id( blockDataItem_.size() );
   internal::BlockDataItem item( id, identifier, dataHandling );
   blockDataItem_.push_back( item );
   
   mpi::RecvBuffer buffer;
   mpi::readMPIIO(file, buffer);

   std::vector< IBlock * > blocks;
   for( auto block = begin(); block != end(); ++block )
      blocks.push_back( block.get() );
   std::sort( blocks.begin(), blocks.end(), internal::sortBlocksByID );
   
   for( auto it = blocks.begin(); it != blocks.end(); ++it )
   {
      IBlock * block = *it;
      auto dh = item.getDataHandling( block );
      if( dh )
      {
         block->addData( id, dh->deserialize( block ) );
         dh->deserialize( block, id, buffer );
      }
      else
         block->addData( id, nullptr );
   }

   return id;   
}



//**********************************************************************************************************************
/*!
*   This function can be used to store the data that corresponds to 'id' to file, so that later (probably when
*   restarting a simulation) the data can be loaded by using the member function 'loadBlockData'.
*/
//**********************************************************************************************************************
void BlockStorage::saveBlockData( const std::string & file, const BlockDataID & id )
{
   WALBERLA_CHECK_LESS( uint_t(id), blockDataItem_.size() );

   mpi::SendBuffer buffer;
   serializeBlockData(id, buffer);
   mpi::writeMPIIO(file, buffer);
}


//**********************************************************************************************************************
/*!
*   This function serializes the block data into a sendbuffer
*/
//**********************************************************************************************************************
void BlockStorage::serializeBlockData( const BlockDataID & id, mpi::SendBuffer & buffer )
{
   WALBERLA_CHECK_LESS( uint_t(id), blockDataItem_.size() );

   std::vector< IBlock * > blocks;
   for( auto block = begin(); block != end(); ++block )
      blocks.push_back( block.get() );
   std::sort( blocks.begin(), blocks.end(), internal::sortBlocksByID );

   auto & item = blockDataItem_[ uint_t(id) ];

   for( auto it = blocks.begin(); it != blocks.end(); ++it )
   {
      IBlock * block = *it;
      auto dh = item.getDataHandling( block );
      if( dh )
         dh->serialize( block, id, buffer );
   }
}


//**********************************************************************************************************************
/*!
*   Deserializes data form a recv buffer into existing block data
 *  Only deserialize data, that has be serialized with serializeBlockData
*/
//**********************************************************************************************************************
void BlockStorage::deserializeBlockData( const BlockDataID & id, mpi::RecvBuffer & buffer )
{
   WALBERLA_CHECK_LESS( uint_t(id), blockDataItem_.size() );

   std::vector< IBlock * > blocks;
   for( auto block = begin(); block != end(); ++block )
      blocks.push_back( block.get() );
   std::sort( blocks.begin(), blocks.end(), internal::sortBlocksByID );

   auto & item = blockDataItem_[ uint_t(id) ];

   for( auto it = blocks.begin(); it != blocks.end(); ++it )
   {
      IBlock * block = *it;
      auto dh = item.getDataHandling( block );
      if( dh )
         dh->deserialize( block, id, buffer );
   }
}


//**********************************************************************************************************************
/*!
*   Returns a MPI communicator that only contains processes that possess blocks
*
*   ATTENTION: This function must always be called by ALL processes, otherwise the mpi::allGather call will stall!
*              Also, the MPI communicator returned by this function must not be stored outside of this BlockStorage
*              class. Every time a MPI communicator that only contains processes that possess blocks is needed, this
*              member function of class BlockStorage must be called - only by doing so you are guaranteed to receive a
*              valid MPI communicator. If blocks have moved between processes, an updated MPI communicator will be
*              returned, otherwise a still valid, internally stored communicator is returned.
*/
//**********************************************************************************************************************
MPI_Comm BlockStorage::processesWithBlocksCommunicator()
{
   WALBERLA_MPI_SECTION()
   {
      if( rebuildProcessesWithBlocksCommunicator_ )
      {
         int8_t hasBlocks = ( getNumberOfBlocks() > uint_t(0) ) ? int8_t(1) : int8_t(0);

         std::vector<int8_t> processHasBlocks = mpi::allGather( hasBlocks );

         std::vector<int> excludedProcesses;
         for( uint_t i=0; i < processHasBlocks.size(); ++i )
            if( processHasBlocks[i] == int8_t(0) )
               excludedProcesses.push_back( int_c(i) );

         if( processesWithBlocksCommunicator_ != MPI_COMM_NULL )
            MPI_Comm_free( &processesWithBlocksCommunicator_ );

         if( excludedProcesses.empty() )
         {
            MPI_Comm_dup( MPI_COMM_WORLD, &processesWithBlocksCommunicator_ );
         }
         else
         {
            MPI_Group worldGroup;
            MPI_Comm_group( MPI_COMM_WORLD, &worldGroup );

            MPI_Group newGroup;
            MPI_Group_excl( worldGroup, int_c( excludedProcesses.size() ), &excludedProcesses[0], &newGroup );
            MPI_Comm_create( MPI_COMM_WORLD, newGroup, &processesWithBlocksCommunicator_ );
            MPI_Group_free( &newGroup );

            MPI_Group_free( &worldGroup );
         }

         rebuildProcessesWithBlocksCommunicator_ = false;
      }
   }

   return processesWithBlocksCommunicator_;
}



} // namespace domain_decomposition
} // namespace walberla
