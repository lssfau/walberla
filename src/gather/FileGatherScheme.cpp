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
//! \file FileGatherScheme.cpp
//! \ingroup gather
//! \author Martin Bauer <martin.bauer@fau.de>
//! \brief Gathering Data using Files
//
//======================================================================================================================

#include "CommonSchemeFunctions.h"
#include "FileGatherScheme.h"
#include "core/debug/Debug.h"
#include "core/mpi/MPIManager.h"
#include "core/mpi/RecvBuffer.h"
#include "core/mpi/SendBuffer.h"

#include "core/Filesystem.h"


namespace walberla {
namespace gather {


size_t FileGatherScheme::nextId_ = 0;


FileGatherScheme::FileGatherScheme( BlockStorage & blockStorage, uint_t everyNTimestep )
   : blocks_( blockStorage ), myId_(nextId_++), everyNTimestep_(everyNTimestep)
{

}

FileGatherScheme::~FileGatherScheme()
{
   collectFromFiles();

   if ( fileStream_.is_open() )
      fileStream_.close();
}



void FileGatherScheme::deleteTemporaryFiles()
{
   using namespace filesystem;

   for(int rank=0; rank < MPIManager::instance()->numProcesses(); ++rank)
   {
      std::string fileName = "tmp_collection_" + std::to_string(myId_) + "_"
                                               + std::to_string(rank);

      if( exists(fileName) )
         remove(fileName);
   }
}



void FileGatherScheme::writeToFile()
{
   static mpi::SendBuffer buffer;

   buffer.clear();

   internal::packData( blocks_, packInfos_, buffer );

   const size_t & bufSize = buffer.size();

   if ( bufSize > 0 && !fileStream_.is_open() )
   {
      // not important if rank() or worldRank() we need only a distinct number
      // for each process
      int rank = MPIManager::instance()->worldRank();
      std::string fileName = "tmp_collection_" + std::to_string(myId_) + "_"
                                               + std::to_string(rank);

      fileStream_.open( fileName.c_str(), std::ios::binary );
   }

   fileStream_.write( reinterpret_cast<const char*>( &bufSize     ), numeric_cast< std::streamsize >( sizeof( bufSize ) ) );
   fileStream_.write( reinterpret_cast<const char*>( buffer.ptr() ), numeric_cast< std::streamsize >( buffer.size()     ) );

   // All contents have to be physically in file, since collectFromFiles might be called
   // in next timestep
   fileStream_.flush();
}





void FileGatherScheme::collectFromFiles()
{
   using namespace filesystem;
   using namespace std;

   if ( fileStream_.is_open() )
      fileStream_.close();

   // Wait until all processes have closed their file
   WALBERLA_MPI_WORLD_BARRIER();

   // Only the root process collects all files together
   if(MPIManager::instance()->worldRank() != 0)
      return;

   vector< std::ifstream* > inputStreams;

   //look for files from all processes
   for(int rank=0; rank < MPIManager::instance()->numProcesses(); ++rank)
   {
      string fileName = "tmp_collection_" + to_string(myId_) + "_"
                                          + to_string(rank);

      if(exists(fileName)) {
		  std::ifstream * str = new std::ifstream(fileName.c_str(),ios::in|ios::binary);
         inputStreams.push_back(str);
      }
   }

   bool finished=false;
   while(!finished) // loop over all timesteps
   {
      // loop over all process-files ( only read a chunk corresponding to a timestep )
      for(size_t i=0; i<inputStreams.size(); ++i)
      {
		  std::ifstream & iStream = *(inputStreams[i]);

         size_t bufferSize=0;
         iStream.read( reinterpret_cast<char*>(&bufferSize),sizeof( size_t ) );


         mpi::GenericRecvBuffer<unsigned char> recvBuffer;
         recvBuffer.resize(bufferSize);
         iStream.read(reinterpret_cast<char*>( recvBuffer.ptr()), numeric_cast< std::streamsize >( bufferSize ) );

         internal::unpackData( packInfos_, recvBuffer );

         if( iStream.peek() == std::ifstream::traits_type::eof() )
         {
            // this is hit when the first file was read completely
            // then the finished flag is set to true, so that the outer timestep loop terminates
            // however the inner file loop goes over all remaining files, reads them and checks
            // here that they are also at the end
            if( i==0 )
               finished=true;
            else
               WALBERLA_ASSERT( finished );
         }
      }

      // Notify all packInfos that this timestep is finished
      for( size_t  s =0; s < packInfos_.size() ; ++s )
         packInfos_[s]->gatherFinished( );

   }



   for(size_t i=0; i<inputStreams.size(); ++i)
      delete inputStreams[i];

   deleteTemporaryFiles();
}



} // namespace gather
} // namespace walberla



