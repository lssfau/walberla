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
#include "core/mpi/Gather.h"
#include "core/mpi/Reduce.h"
#include "core/uid/GlobalState.h"


namespace walberla {
namespace domain_decomposition {



//**********************************************************************************************************************
/*!
*   For documentation, see documentation of free function
*   'void mapPointToPeriodicDomain( const boost::array< bool, 3 > & periodic, const AABB & domain,
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
*   'bool periodicIntersect( const boost::array< bool, 3 > & periodic, const math::AABB & domain, const math::AABB & box1, const math::AABB & box2 )'
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

   WALBERLA_NON_MPI_SECTION()
   {
      std::ifstream ifile( file.c_str(), std::ifstream::binary );
      
      ifile.seekg( 0, std::ios::end );
      const uint_t length = uint_c( static_cast< std::streamoff >( ifile.tellg() ) );
      ifile.seekg( 0, std::ios::beg );
      
      buffer.resize( length / sizeof( mpi::RecvBuffer::ElementType ) );
      
      ifile.read( reinterpret_cast< char* >( buffer.ptr() ), numeric_cast< std::streamsize >( length ) );
      ifile.close();
   }
   
   WALBERLA_MPI_SECTION()
   {
      MPI_File mpiFile = MPI_FILE_NULL;
      int result = MPI_SUCCESS;
      result = MPI_File_open( MPIManager::instance()->comm(), const_cast<char*>( file.c_str() ), MPI_MODE_RDONLY, MPI_INFO_NULL, &mpiFile );

      if( result != MPI_SUCCESS )
         WALBERLA_ABORT( "Error while opening file \"" << file << "\" for writing. MPI Error is \"" << MPIManager::instance()->getMPIErrorString( result ) << "\"" );
      
      int uintSize;
      MPI_Type_size( MPITrait< uint_t >::type(), &uintSize );
      
      MPI_Datatype offsettype;
      MPI_Type_contiguous( 2, MPITrait< uint_t >::type(), &offsettype );
      MPI_Type_commit( &offsettype );
      
      // offsets to processes data (+ buffer size)
      
      result = MPI_File_set_view( mpiFile, numeric_cast<MPI_Offset>( mpi::MPIManager::instance()->rank() * 2 * uintSize ),
                                  MPITrait< uint_t >::type(), offsettype, const_cast<char*>( "native" ), MPI_INFO_NULL );
      
      if( result != MPI_SUCCESS )
         WALBERLA_ABORT( "Internal MPI-IO error! MPI Error is \"" << MPIManager::instance()->getMPIErrorString( result ) << "\"" );

      uint_t offsetData[2];
      
      result = MPI_File_read_all( mpiFile, reinterpret_cast<char*>( offsetData ), 2, MPITrait< uint_t >::type(), MPI_STATUS_IGNORE );
      
      if( result != MPI_SUCCESS )
         WALBERLA_ABORT( "Error while writing to file \"" << file << "\". MPI Error is \"" << MPIManager::instance()->getMPIErrorString( result ) << "\"" );

      // the data
      
      MPI_Datatype arraytype;
      MPI_Type_contiguous( int_c( offsetData[1] ), MPITrait< mpi::RecvBuffer::ElementType >::type(), &arraytype );
      MPI_Type_commit( &arraytype );

      result = MPI_File_set_view( mpiFile, numeric_cast<MPI_Offset>( offsetData[0] ), MPITrait< mpi::RecvBuffer::ElementType >::type(), arraytype,
                                  const_cast<char*>( "native" ), MPI_INFO_NULL );
      
      if( result != MPI_SUCCESS )
         WALBERLA_ABORT( "Internal MPI-IO error! MPI Error is \"" << MPIManager::instance()->getMPIErrorString( result ) << "\"" );
      
      buffer.resize( offsetData[1] );
      
      result = MPI_File_read_all( mpiFile, reinterpret_cast<char*>( buffer.ptr() ), int_c( buffer.size() ), MPITrait< mpi::SendBuffer::ElementType >::type(), MPI_STATUS_IGNORE );

      if( result != MPI_SUCCESS )
         WALBERLA_ABORT( "Error while writing to file \"" << file << "\". MPI Error is \"" << MPIManager::instance()->getMPIErrorString( result ) << "\"" );
      
      result = MPI_File_close( &mpiFile );

      if( result != MPI_SUCCESS )
         WALBERLA_ABORT( "Error while closing file \"" << file << "\". MPI Error is \"" << MPIManager::instance()->getMPIErrorString( result ) << "\"" );

      MPI_Type_free( &arraytype );
      MPI_Type_free( &offsettype );
   }
   
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
   
   std::vector< IBlock * > blocks;
   for( auto block = begin(); block != end(); ++block )
      blocks.push_back( block.get() );
   std::sort( blocks.begin(), blocks.end(), internal::sortBlocksByID );
   
   auto & item = blockDataItem_[ uint_t(id) ];
   mpi::SendBuffer buffer;
   
   for( auto it = blocks.begin(); it != blocks.end(); ++it )
   {
      IBlock * block = *it;
      auto dh = item.getDataHandling( block );
      if( dh )
         dh->serialize( block, id, buffer );
   }
   
   uint_t dataSize = sizeof( mpi::SendBuffer::ElementType ) * buffer.size();
   
   WALBERLA_NON_MPI_SECTION()
   {
      std::ofstream ofile( file.c_str(), std::ofstream::binary );
      ofile.write( reinterpret_cast< const char* >( buffer.ptr() ), numeric_cast< std::streamsize >( dataSize ) );
      ofile.close();
   }
   
   WALBERLA_MPI_SECTION()
   {
      MPI_File mpiFile = MPI_FILE_NULL;
      int result = MPI_SUCCESS;
      result = MPI_File_open( MPIManager::instance()->comm(), const_cast<char*>( file.c_str() ), MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &mpiFile );

      if( result != MPI_SUCCESS )
         WALBERLA_ABORT( "Error while opening file \"" << file << "\" for writing. MPI Error is \"" << MPIManager::instance()->getMPIErrorString( result ) << "\"" );
      
      int uintSize;
      MPI_Type_size( MPITrait< uint_t >::type(), &uintSize );
      
      const MPI_Offset filesize = numeric_cast<MPI_Offset>( mpi::MPIManager::instance()->numProcesses() * 2 * uintSize ) +
                                  numeric_cast<MPI_Offset>( mpi::allReduce( dataSize, mpi::SUM, MPIManager::instance()->comm() ) );
      MPI_File_set_size( mpiFile, filesize );
      
      uint_t exscanResult;
      MPI_Exscan( &dataSize, &exscanResult, 1, MPITrait<uint_t>::type(), MPI_SUM, MPIManager::instance()->comm() );
      if( MPIManager::instance()->rank() == 0 )
         exscanResult = uint_t( 0 );
      const uint_t offset = uint_c( mpi::MPIManager::instance()->numProcesses() * 2 * uintSize ) + exscanResult;

      MPI_Datatype offsettype;
      MPI_Type_contiguous( 2, MPITrait< uint_t >::type(), &offsettype );
      MPI_Type_commit( &offsettype );
      
      MPI_Datatype arraytype;
      MPI_Type_contiguous( int_c( buffer.size() ), MPITrait< mpi::SendBuffer::ElementType >::type(), &arraytype );
      MPI_Type_commit( &arraytype );
      
      // offsets to processes data (+ buffer size)
      
      result = MPI_File_set_view( mpiFile, numeric_cast<MPI_Offset>( mpi::MPIManager::instance()->rank() * 2 * uintSize ),
                                  MPITrait< uint_t >::type(), offsettype, const_cast<char*>( "native" ), MPI_INFO_NULL );
      
      if( result != MPI_SUCCESS )
         WALBERLA_ABORT( "Internal MPI-IO error! MPI Error is \"" << MPIManager::instance()->getMPIErrorString( result ) << "\"" );

      uint_t offsetData[2];
      offsetData[0] = offset;
      offsetData[1] = uint_c( buffer.size() );
      
      result = MPI_File_write_all( mpiFile, reinterpret_cast<char*>( offsetData ), 2, MPITrait< uint_t >::type(), MPI_STATUS_IGNORE );
      
      if( result != MPI_SUCCESS )
         WALBERLA_ABORT( "Error while writing to file \"" << file << "\". MPI Error is \"" << MPIManager::instance()->getMPIErrorString( result ) << "\"" );

      // the data

      result = MPI_File_set_view( mpiFile, numeric_cast<MPI_Offset>( offset ), MPITrait< mpi::SendBuffer::ElementType >::type(), arraytype,
                                  const_cast<char*>( "native" ), MPI_INFO_NULL );
      
      if( result != MPI_SUCCESS )
         WALBERLA_ABORT( "Internal MPI-IO error! MPI Error is \"" << MPIManager::instance()->getMPIErrorString( result ) << "\"" );
      
      result = MPI_File_write_all( mpiFile, reinterpret_cast<char*>( buffer.ptr() ), int_c( buffer.size() ), MPITrait< mpi::SendBuffer::ElementType >::type(), MPI_STATUS_IGNORE );

      if( result != MPI_SUCCESS )
         WALBERLA_ABORT( "Error while writing to file \"" << file << "\". MPI Error is \"" << MPIManager::instance()->getMPIErrorString( result ) << "\"" );

      result = MPI_File_close( &mpiFile );

      if( result != MPI_SUCCESS )
         WALBERLA_ABORT( "Error while closing file \"" << file << "\". MPI Error is \"" << MPIManager::instance()->getMPIErrorString( result ) << "\"" );

      MPI_Type_free( &arraytype );
      MPI_Type_free( &offsettype );
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
