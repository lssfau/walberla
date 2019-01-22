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
//! \file MPIIO.cppp
//! \ingroup core
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#include "core/DataTypes.h"
#include "core/mpi/RecvBuffer.h"
#include "core/mpi/SendBuffer.h"
#include "core/mpi/MPIManager.h"
#include "core/mpi/Reduce.h"
#include <fstream>


namespace walberla {
namespace mpi  {

   void writeMPIIO(const std::string & file, SendBuffer & buffer)
   {
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


   void readMPIIO(const std::string & file, RecvBuffer & buffer)
   {
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
   }


} // namespace mpi
} // namespace walberla


