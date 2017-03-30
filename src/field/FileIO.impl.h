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
//! \file FileIO.impl.h
//! \ingroup field
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//
//======================================================================================================================

#pragma once

namespace walberla {
namespace field {

namespace internal {


template< typename FieldT >
class FieldWriter
{
public:
   FieldWriter( const std::string & filename, const BlockDataID & fieldID,
                const Set<SUID> & requiredSelectors = Set<SUID>::emptySet(), const Set<SUID> & incompatibleSelectors = Set<SUID>::emptySet() )
      : filename_( filename ), fieldID_( fieldID ), requiredSelectors_( requiredSelectors ), incompatibleSelectors_( incompatibleSelectors )
   {
   }

   void writeToFile( const BlockStorage & blockStorage ) const;
   void readFromFile( BlockStorage & blockStorage ) const;
private:

   std::vector< IBlock * > getBlocks( BlockStorage & blockStorage ) const;
   std::vector< const IBlock * > getBlocks( const BlockStorage & blockStorage ) const;

   void writeToFileNonMPI( const std::vector< const IBlock * > & blocks ) const;
   void readFromFileNonMPI( const std::vector< IBlock * > & blocks ) const;

   std::vector< uint_t > computeBlockOffsets( const std::vector< const IBlock * > & blocks ) const;
   std::vector< uint_t > computeBlockOffsets( const std::vector< IBlock * > & blocks ) const;

   MPI_Offset computeProcessByteOffset( uint_t processNumElements ) const;


   std::string filename_;
   BlockDataID fieldID_;
   
   Set<SUID> requiredSelectors_;
   Set<SUID> incompatibleSelectors_;
};



template< typename FieldT >
void FieldWriter<FieldT>::writeToFile( const BlockStorage & blockStorage ) const
{
   std::vector< const IBlock * > blocks = getBlocks( blockStorage );
   
   WALBERLA_NON_MPI_SECTION()
   {
      writeToFileNonMPI( blocks );
      return;
   }

   std::vector< uint_t > blockOffsets = computeBlockOffsets( blocks );

   if( blockOffsets.back() > uint_c( std::numeric_limits<int>::max() ) )
      WALBERLA_ABORT( "writeToFile does not support writing more than " << std::numeric_limits<int>::max() << " field elements per process!" );

   const MPI_Offset filesize = numeric_cast<MPI_Offset>( mpi::allReduce( blockOffsets.back(), mpi::SUM, MPIManager::instance()->comm() ) * sizeof( typename FieldT::value_type ) );

   MPI_Datatype arraytype;
   MPI_Type_contiguous( int_c( blockOffsets.back() ), MPITrait< typename FieldT::value_type >::type(), &arraytype );
   MPI_Type_commit( &arraytype );

   const MPI_Offset offset = computeProcessByteOffset( blockOffsets.back() );

   MPI_File mpiFile = MPI_FILE_NULL;
   int result = MPI_SUCCESS;
   result = MPI_File_open( MPIManager::instance()->comm(), const_cast<char*>( filename_.c_str() ), MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &mpiFile );

   if( result != MPI_SUCCESS )
      WALBERLA_ABORT( "Error while opening file \"" << filename_ << "\" for writing. MPI Error is \"" << MPIManager::instance()->getMPIErrorString( result ) << "\"" );

   MPI_File_set_size( mpiFile, filesize );

   result = MPI_File_set_view( mpiFile, offset, MPITrait< typename FieldT::value_type >::type(), arraytype, const_cast<char*>( "native" ), MPI_INFO_NULL );

   if( result != MPI_SUCCESS )
      WALBERLA_ABORT( "Internal MPI-IO error! MPI Error is \"" << MPIManager::instance()->getMPIErrorString( result ) << "\"" );

   size_t blockIdx = 0;
   for( auto block = blocks.begin(); block != blocks.end(); ++block, ++blockIdx )
   {
      const FieldT * field = (*block)->template getData<FieldT>( fieldID_ );

      std::vector< typename FieldT::value_type > data( field->xSize() * field->ySize() * field->zSize() * field->fSize() );

      auto dataIt = data.begin();
      for( auto fieldIt = field->begin(); fieldIt != field->end(); ++fieldIt, ++dataIt )
      {
         WALBERLA_ASSERT( dataIt != data.end() );
         *dataIt = *fieldIt;
      }
      WALBERLA_ASSERT_EQUAL( dataIt, data.end() );

      result = MPI_File_write_at( mpiFile, numeric_cast<MPI_Offset>( blockOffsets[blockIdx] ), &data[0], int_c( data.size() ),
                                  MPITrait<typename FieldT::value_type>::type(), MPI_STATUS_IGNORE );

      if( result != MPI_SUCCESS )
         WALBERLA_ABORT( "Error while writing to file \"" << filename_ << "\". MPI Error is \"" << MPIManager::instance()->getMPIErrorString( result ) << "\"" );
   }

   result = MPI_File_close( &mpiFile );

   if( result != MPI_SUCCESS )
      WALBERLA_ABORT( "Error while closing file \"" << filename_ << "\". MPI Error is \"" << MPIManager::instance()->getMPIErrorString( result ) << "\"" );

   MPI_Type_free( &arraytype );
}



template< typename FieldT >
void FieldWriter<FieldT>::readFromFile( BlockStorage & blockStorage ) const
{
   std::vector< IBlock * > blocks = getBlocks( blockStorage );
   
   WALBERLA_NON_MPI_SECTION()
   {
      readFromFileNonMPI( blocks );
      return;
   }

   std::vector< uint_t > blockOffsets = computeBlockOffsets( blocks );

   if( blockOffsets.back() > uint_c( std::numeric_limits<int>::max() ) )
      WALBERLA_ABORT( "readFromFile does not support reading more than " << std::numeric_limits<int>::max() << " field elements per process!" );

   MPI_Datatype arraytype;
   MPI_Type_contiguous( int_c( blockOffsets.back() ), MPITrait< typename FieldT::value_type >::type(), &arraytype );
   MPI_Type_commit( &arraytype );

   const MPI_Offset offset = computeProcessByteOffset( blockOffsets.back() );

   MPI_File mpiFile;
   int result = MPI_SUCCESS;
   result = MPI_File_open( MPIManager::instance()->comm(), const_cast<char*>( filename_.c_str() ), MPI_MODE_RDONLY, MPI_INFO_NULL, &mpiFile );

   if( result != MPI_SUCCESS )
      WALBERLA_ABORT( "Error while opening file \"" << filename_ << "\" for reading. MPI Error is \"" << MPIManager::instance()->getMPIErrorString( result ) << "\"" );

   result = MPI_File_set_view( mpiFile, offset, MPITrait< typename FieldT::value_type >::type(), arraytype, const_cast<char*>( "native" ), MPI_INFO_NULL );

   if( result != MPI_SUCCESS )
      WALBERLA_ABORT( "Internal MPI-IO error! MPI Error is \"" << MPIManager::instance()->getMPIErrorString( result ) << "\"" );

   size_t blockIdx = 0;
   for( auto block = blocks.begin(); block != blocks.end(); ++block, ++blockIdx )
   {
      FieldT * field = (*block)->template getData<FieldT>( fieldID_ );

      std::vector< typename FieldT::value_type > data( field->xSize() * field->ySize() * field->zSize() * field->fSize() );

      result = MPI_File_read_at( mpiFile, numeric_cast<MPI_Offset>( blockOffsets[blockIdx] ), &data[0], int_c( data.size() ),
                                 MPITrait<typename FieldT::value_type>::type(), MPI_STATUS_IGNORE );

      auto dataIt = data.begin();
      for( auto fieldIt = field->begin(); fieldIt != field->end(); ++fieldIt, ++dataIt )
      {
         WALBERLA_ASSERT_UNEQUAL( dataIt, data.end() );
         *fieldIt = *dataIt;
      }
      WALBERLA_ASSERT_EQUAL( dataIt, data.end() );

      if( result != MPI_SUCCESS )
         WALBERLA_ABORT( "Error while reading from file \"" << filename_ << "\". MPI Error is \"" << MPIManager::instance()->getMPIErrorString( result ) << "\"" );
   }

   result = MPI_File_close( &mpiFile );

   if( result != MPI_SUCCESS )
      WALBERLA_ABORT( "Error while closing file \"" << filename_ << "\". MPI Error is \"" << MPIManager::instance()->getMPIErrorString( result ) << "\"" );

   MPI_Type_free( &arraytype );
}



inline bool sortBlocksByID( IBlock * lhs, IBlock * rhs ) { return lhs->getId() < rhs->getId(); }

template< typename FieldT >
std::vector< IBlock * > FieldWriter<FieldT>::getBlocks( BlockStorage & blockStorage ) const
{
   std::vector< IBlock * > blocks;
   for( auto it = blockStorage.begin( requiredSelectors_, incompatibleSelectors_ ); it != blockStorage.end(); ++it )
      blocks.push_back( it.get() );
   std::sort( blocks.begin(), blocks.end(), sortBlocksByID );
   return blocks;
}



inline bool sortConstBlocksByID( const IBlock * lhs, const IBlock * rhs ) { return lhs->getId() < rhs->getId(); }

template< typename FieldT >
std::vector< const IBlock * > FieldWriter<FieldT>::getBlocks( const BlockStorage & blockStorage ) const
{
   std::vector< const IBlock * > blocks;
   for( auto it = blockStorage.begin( requiredSelectors_, incompatibleSelectors_ ); it != blockStorage.end(); ++it )
      blocks.push_back( it.get() );
   std::sort( blocks.begin(), blocks.end(), sortConstBlocksByID );
   return blocks;
}



template< typename FieldT >
void FieldWriter<FieldT>::writeToFileNonMPI( const std::vector< const IBlock * > & blocks ) const
{
   std::ofstream ofs( filename_.c_str(), std::ofstream::out | std::ofstream::binary );

   for( auto block = blocks.begin(); block != blocks.end(); ++block )
   {
      const FieldT * field = (*block)->template getData<FieldT>( fieldID_ );

      for( auto fieldIt = field->begin(); fieldIt != field->end(); ++fieldIt )
      {
         ofs.write( reinterpret_cast<const char*>( &( *fieldIt ) ), sizeof( typename FieldT::value_type ) );
      }
   }

   ofs.close();
}



template< typename FieldT >
void FieldWriter<FieldT>::readFromFileNonMPI( const std::vector< IBlock * > & blocks ) const
{
   std::ifstream ifs( filename_.c_str(), std::ifstream::in | std::ifstream::binary );

   for( auto block = blocks.begin(); block != blocks.end(); ++block )
   {
      FieldT * field = (*block)->template getData<FieldT>( fieldID_ );

      for( auto fieldIt = field->begin(); fieldIt != field->end(); ++fieldIt )
      {
         ifs.read( reinterpret_cast<char*>( &( *fieldIt ) ), sizeof( typename FieldT::value_type ) );
      }
   }

   ifs.close();
}



template< typename FieldT >
std::vector< uint_t > FieldWriter<FieldT>::computeBlockOffsets( const std::vector< const IBlock * > & blocks ) const
{
   std::vector< uint_t > blockOffsets;
   blockOffsets.push_back( 0 );

   for( auto block = blocks.begin(); block != blocks.end(); ++block )
   {
      const FieldT * field = (*block)->template getData<FieldT>( fieldID_ );
      const uint_t offset = blockOffsets.back() + field->xSize() * field->ySize() * field->zSize() * field->fSize();
      blockOffsets.push_back( offset );
   }

   return blockOffsets;
}



template< typename FieldT >
std::vector< uint_t > FieldWriter<FieldT>::computeBlockOffsets( const std::vector< IBlock * > & blocks ) const
{
   std::vector< uint_t > blockOffsets;
   blockOffsets.push_back( 0 );

   for( auto block = blocks.begin(); block != blocks.end(); ++block )
   {
      const FieldT * field = (*block)->template getData<FieldT>( fieldID_ );
      const uint_t offset = blockOffsets.back() + field->xSize() * field->ySize() * field->zSize() * field->fSize();
      blockOffsets.push_back( offset );
   }

   return blockOffsets;
}



template< typename FieldT >
MPI_Offset FieldWriter<FieldT>::computeProcessByteOffset( uint_t processNumElements ) const
{
   uint_t exscanResult;
   MPI_Exscan( &processNumElements, &exscanResult, 1, MPITrait<uint_t>::type(), MPI_SUM, MPIManager::instance()->comm() );
   if( MPIManager::instance()->rank() == 0 )
      exscanResult = uint_t( 0 );

   return numeric_cast<MPI_Offset>( exscanResult * sizeof( typename FieldT::value_type ) );
}


} // namespace internal



template< typename FieldT >
void writeToFile( const std::string & filename, const BlockStorage & blockStorage, const BlockDataID & fieldID,
                  const Set<SUID> & requiredSelectors, const Set<SUID> & incompatibleSelectors )
{
   internal::FieldWriter<FieldT> writer( filename, fieldID, requiredSelectors, incompatibleSelectors );
   writer.writeToFile( blockStorage );
}



template< typename FieldT >
void readFromFile( const std::string & filename, BlockStorage & blockStorage, const BlockDataID & fieldID,
                   const Set<SUID> & requiredSelectors, const Set<SUID> & incompatibleSelectors )
{
   internal::FieldWriter<FieldT> writer( filename, fieldID, requiredSelectors, incompatibleSelectors );
   writer.readFromFile( blockStorage );
}

} // namespace walberla
} // namespace field
