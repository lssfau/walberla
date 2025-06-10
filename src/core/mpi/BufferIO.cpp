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
//! \file BufferIO.cpp
//! \ingroup core
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//! \author Martin Bauer <martin.bauer@fau.de>
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//
//======================================================================================================================

#include "BufferIO.h"

#include "core/mpi/RecvBuffer.h"
#include "core/mpi/SendBuffer.h"

#include <fstream>

namespace walberla::mpi
{

constexpr uint_t arraySize = 2;
using uintArray = std::array< uint_t, arraySize >;

void writeBufferNoMPI(const std::string& file, SendBuffer& buffer)
{
   WALBERLA_MPI_SECTION(){ WALBERLA_ABORT("This function can only be used if waLBerla is build without MPI"); }
   
   uintArray metaData;
   metaData[0] = uint_c(0);
   metaData[1] = uint_c(0);
   const std::streamsize metaDataSize = numeric_cast< std::streamsize >(arraySize) * numeric_cast< std::streamsize >( sizeof(uint_t) );

   std::ofstream ofile(file.c_str(), std::ofstream::binary);
   if (ofile.fail()) { WALBERLA_ABORT("Error while opening file \"" << file << "\" for writing."); }

   // write metaData to file
   ofile.write(reinterpret_cast< const char* >(metaData.data()), metaDataSize);
   if (ofile.fail()) { WALBERLA_ABORT("Error while writing to file \"" << file << "\"."); }

   ofile.write(reinterpret_cast< const char* >(buffer.ptr()), numeric_cast< std::streamsize >(buffer.dataSize()));
   if (ofile.fail()) { WALBERLA_ABORT("Error while writing to file \"" << file << "\"."); }

   ofile.close();
   if (ofile.fail()) { WALBERLA_ABORT("Error while closing file \"" << file << "\" for reading."); }
}

void writeBufferSerialIO(const std::string& file, SendBuffer& buffer)
{
   WALBERLA_NON_MPI_SECTION(){ WALBERLA_ABORT("This function can only be used if waLBerla is build with MPI"); }

   std::ofstream ofile;

   uintArray metaData;
   metaData[0] = uint_c(1);
   metaData[1] = uint_c(MPIManager::instance()->numProcesses());
   const std::streamsize metaDataSize = numeric_cast< std::streamsize >(arraySize) * numeric_cast< std::streamsize >(MPITrait< uintArray::value_type >::size());

   uintArray offsetData;
   offsetData[0] = internal::getProcessOffset(buffer.dataSize());
   offsetData[1] = uint_c(buffer.size());
   const std::streamsize offsetDataSize = numeric_cast< std::streamsize >(arraySize) * numeric_cast< std::streamsize >(MPITrait< uintArray::value_type >::size());

   // write each process' offset and buffer size to file
   for (int i = 0; i != MPIManager::instance()->numProcesses(); ++i)
   {
      if (i == MPIManager::instance()->rank())
      {
         if (MPIManager::instance()->rank() == 0)
         {
            ofile.open(file.c_str(), std::ofstream::binary);
            if (ofile.fail()) { WALBERLA_ABORT("Error while opening file \"" << file << "\" for writing."); }
            // write metaData to file
            ofile.write(reinterpret_cast< const char* >(metaData.data()), metaDataSize);
            if (ofile.fail()) { WALBERLA_ABORT("Error while writing to file \"" << file << "\"."); }
         }
         else
         {
            ofile.open(file.c_str(), std::ofstream::binary | std::ofstream::app);
         }
         if (ofile.fail()) { WALBERLA_ABORT("Error while opening file \"" << file << "\" for writing."); }

         ofile.write(reinterpret_cast< const char* >(offsetData.data()), offsetDataSize);
         if (ofile.fail()) { WALBERLA_ABORT("Error while writing to file \"" << file << "\"."); }

         ofile.close();
         if (ofile.fail()) { WALBERLA_ABORT("Error while closing file \"" << file << "\" for writing."); }
      }
      WALBERLA_MPI_BARRIER();
   }

   // write each process' data to file (starting at offset)
   for (int i = 0; i != MPIManager::instance()->numProcesses(); ++i)
   {
      if (i == MPIManager::instance()->rank())
      {
         ofile.open(file.c_str(), std::ofstream::binary | std::ofstream::app);
         if (ofile.fail()) { WALBERLA_ABORT("Error while opening file \"" << file << "\" for writing."); }

         ofile.write(reinterpret_cast< const char* >(buffer.ptr()), numeric_cast< std::streamsize >(buffer.dataSize()));
         if (ofile.fail()) { WALBERLA_ABORT("Error while writing to file \"" << file << "\"."); }

         ofile.close();
         if (ofile.fail()) { WALBERLA_ABORT("Error while closing file \"" << file << "\" for writing."); }
      }
      WALBERLA_MPI_BARRIER();
   }
}

void writeBufferMPIIO(const std::string& file, SendBuffer& buffer)
{
   WALBERLA_NON_MPI_SECTION(){ WALBERLA_ABORT("This function can only be used if waLBerla is build with MPI"); }

   uintArray metaData;
   metaData[0] = uint_c(2);
   metaData[1] = uint_c(MPIManager::instance()->numProcesses());
   const MPI_Offset metaDataSize = numeric_cast< MPI_Offset >(arraySize) * numeric_cast< MPI_Offset >(MPITrait< uintArray::value_type >::size());

   uintArray offsetData;
   offsetData[0] = internal::getProcessOffset(buffer.dataSize());
   offsetData[1] = uint_c(buffer.size());
   const MPI_Offset offsetDataSize = numeric_cast< MPI_Offset >(arraySize) * numeric_cast< MPI_Offset >(MPITrait< uintArray::value_type >::size());

   MPI_File mpiFile = MPI_FILE_NULL;
   int result = MPI_File_open(MPIManager::instance()->comm(), const_cast< char* >(file.c_str()),
                          MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &mpiFile);

   if (result != MPI_SUCCESS)
      WALBERLA_ABORT("Error while opening file \"" << file << "\" for writing. MPI Error is \""
                                                   << MPIManager::instance()->getMPIErrorString(result) << "\"");

   if (MPIManager::instance()->rank() == 0)
   {
      result = MPI_File_write_at(mpiFile, 0, reinterpret_cast< char* >(metaData.data()), metaData.size(), MPITrait< uint_t >::type(), MPI_STATUS_IGNORE);
      if (result != MPI_SUCCESS)
         WALBERLA_ABORT("Error while writing metaData to file \"" << file << "\". MPI Error is \""
                                                         << MPIManager::instance()->getMPIErrorString(result)
                                                         << "\"");
   }

   WALBERLA_MPI_BARRIER()

   MPI_Datatype offsettype;
   MPI_Type_contiguous(offsetData.size(), MPITrait< uint_t >::type(), &offsettype);
   MPI_Type_commit(&offsettype);

   MPI_Datatype arraytype;
   MPI_Type_contiguous(int_c(buffer.size()), MPITrait< mpi::SendBuffer::ElementType >::type(), &arraytype);
   MPI_Type_commit(&arraytype);

   // write each process' offset and buffer size to file
   const auto displacement = numeric_cast< MPI_Offset >(mpi::MPIManager::instance()->rank()) * offsetDataSize + metaDataSize;
   result = MPI_File_set_view(mpiFile, displacement, MPITrait< uint_t >::type(), offsettype, const_cast< char* >("native"), MPI_INFO_NULL);

   if (result != MPI_SUCCESS)
      WALBERLA_ABORT("Internal MPI-IO error! MPI Error is \"" << MPIManager::instance()->getMPIErrorString(result)
                                                              << "\"");

   result = MPI_File_write_all(mpiFile, reinterpret_cast< char* >(offsetData.data()), offsetData.size(), MPITrait< uint_t >::type(),
                               MPI_STATUS_IGNORE);

   if (result != MPI_SUCCESS)
      WALBERLA_ABORT("Error while writing to file \"" << file << "\". MPI Error is \""
                                                      << MPIManager::instance()->getMPIErrorString(result)
                                                      << "\"");

   // write each process' data to file (starting at offset)
   result = MPI_File_set_view(mpiFile, numeric_cast< MPI_Offset >( offsetData[0] ) + metaDataSize,
                              MPITrait< mpi::SendBuffer::ElementType >::type(), arraytype,
                              const_cast< char* >("native"), MPI_INFO_NULL);


   if (result != MPI_SUCCESS)
      WALBERLA_ABORT("Internal MPI-IO error! MPI Error is \"" << MPIManager::instance()->getMPIErrorString(result)
                                                              << "\"");

   result = MPI_File_write_all(mpiFile, reinterpret_cast< char* >(buffer.ptr()), int_c(buffer.size()),
                               MPITrait< mpi::SendBuffer::ElementType >::type(), MPI_STATUS_IGNORE);

   if (result != MPI_SUCCESS)
      WALBERLA_ABORT("Error while writing to file \"" << file << "\". MPI Error is \""
                                                      << MPIManager::instance()->getMPIErrorString(result)
                                                      << "\"");

   result = MPI_File_close(&mpiFile);

   if (result != MPI_SUCCESS)
      WALBERLA_ABORT("Error while closing file \"" << file << "\". MPI Error is \""
                                                   << MPIManager::instance()->getMPIErrorString(result) << "\"");

   MPI_Type_free(&arraytype);
   MPI_Type_free(&offsettype);


}

void writeBuffer(const std::string& file, SendBuffer& buffer, const bool forceSerialIO)
{

   WALBERLA_NON_MPI_SECTION()
   {
      writeBufferNoMPI(file, buffer);
   }

   WALBERLA_MPI_SECTION()
   {
      // use serial I/O for versions of OpenMPI that produce segmentation faults when using MPI-IO with a 3D
      // Cartesian MPI communicator (see waLBerla issue #73)
      if (forceSerialIO || !(MPIManager::instance()->isCommMPIIOValid()) )
      {
         writeBufferSerialIO(file, buffer);
      }
      else // use MPI-IO
      {
         writeBufferMPIIO(file, buffer);
      }
   }
}

void readBufferNoMPI(const std::string& file, RecvBuffer& buffer)
{
   WALBERLA_MPI_SECTION(){ WALBERLA_ABORT("This function can only be used if waLBerla is build without MPI"); }

   std::ifstream ifile(file.c_str(), std::ifstream::binary);
   if (ifile.fail()) { WALBERLA_ABORT("Error while opening file \"" << file << "\" for reading."); }

   // check that file was written using writeBufferSerialIO (metaData[0] must be zero)
   uintArray metaData;
   const std::streamoff metaDataSize = numeric_cast< std::streamoff >( arraySize ) * numeric_cast< std::streamoff >( sizeof(uint_t) );

   ifile.seekg(numeric_cast< std::streamoff >(0));
   ifile.read(reinterpret_cast< char* >(metaData.data()), numeric_cast< std::streamsize >(metaDataSize));

   if (ifile.fail()) { WALBERLA_ABORT("Error while reading from file \"" << file << "\"."); }

   if (metaData[0] != uint_c(0)){
      WALBERLA_ABORT("Error file was not written without MPI. MetaData[0] = " << metaData[0] << " != 0");
   }

   ifile.seekg(0, std::ios::end);
   const uint_t length = uint_c(static_cast< std::streamoff >(ifile.tellg())) - uint_c(metaDataSize);
   ifile.seekg(numeric_cast< std::streamsize >(metaDataSize), std::ios::beg);

   buffer.resize(length / sizeof(mpi::RecvBuffer::ElementType));

   ifile.read(reinterpret_cast< char* >(buffer.ptr()), numeric_cast< std::streamsize >(length));
   if (ifile.fail()) { WALBERLA_ABORT("Error while reading from file \"" << file << "\"."); }

   ifile.close();
   if (ifile.fail()) { WALBERLA_ABORT("Error while closing file \"" << file << "\" for reading."); }

}


void readBufferSerialIO(const std::string& file, RecvBuffer& buffer)
{
   WALBERLA_NON_MPI_SECTION(){ WALBERLA_ABORT("This function can only be used if waLBerla is build with MPI"); }

   std::ifstream ifile;
   uintArray offsetData;

   const std::streamoff metaDataSize   = numeric_cast< std::streamoff >(arraySize ) * numeric_cast< std::streamoff >(MPITrait< uintArray::value_type >::size());
   const std::streamoff offsetDataSize = numeric_cast< std::streamoff >(arraySize ) * numeric_cast< std::streamoff >(MPITrait< uintArray::value_type >::size());

   // read each process' offset and buffer size from file
   for (int i = 0; i != MPIManager::instance()->numProcesses(); ++i)
   {
      if (i == MPIManager::instance()->rank())
      {
         ifile.open(file.c_str(), std::ofstream::binary);
         if (ifile.fail()) { WALBERLA_ABORT("Error while opening file \"" << file << "\" for reading."); }

         if (i == 0)
         {
            // check that file was written using writeBufferSerialIO (metaData[0] must be one)
            // check that file was witten using same process distribution (metaData[1] must be equal to number of MPI processes)
            uintArray metaData;
            ifile.seekg(numeric_cast< std::streamoff >(0));
            ifile.read(reinterpret_cast< char* >(metaData.data()), numeric_cast< std::streamsize >(metaDataSize));
            if (ifile.fail()) { WALBERLA_ABORT("Error while reading from file \"" << file << "\"."); }

            if (metaData[0] != uint_c(1)){
               WALBERLA_ABORT("Error file was not written using serial I/O. MetaData[0] = " << metaData[0] << " != 1");
            }

            if (metaData[1] != uint_c(MPIManager::instance()->numProcesses())){
               WALBERLA_ABORT("Error file was written with different process distribution. Number of processes written file = " << metaData[1] << " != " << MPIManager::instance()->numProcesses())
            }
         }

         const auto displacement = numeric_cast< std::streamoff >(mpi::MPIManager::instance()->rank()) * offsetDataSize + metaDataSize;
         ifile.seekg(displacement);
         ifile.read(reinterpret_cast< char* >(offsetData.data()), offsetDataSize);
         if (ifile.fail()) { WALBERLA_ABORT("Error while reading from file \"" << file << "\"."); }

         ifile.close();
         if (ifile.fail()) { WALBERLA_ABORT("Error while closing file \"" << file << "\" for reading."); }
      }
      WALBERLA_MPI_BARRIER();
   }

   buffer.resize(offsetData[1] / sizeof(mpi::RecvBuffer::ElementType));

   // read each process' data from file (starting at offset)
   for (int i = 0; i != MPIManager::instance()->numProcesses(); ++i)
   {
      if (i == MPIManager::instance()->rank())
      {
         ifile.open(file.c_str(), std::ofstream::binary);
         if (ifile.fail()) { WALBERLA_ABORT("Error while opening file \"" << file << "\" for reading."); }

         ifile.seekg(numeric_cast< std::streamoff >( offsetData[0] ) + metaDataSize);
         ifile.read(reinterpret_cast< char* >(buffer.ptr()), numeric_cast< std::streamsize >(offsetData[1]));
         if (ifile.fail()) { WALBERLA_ABORT("Error while reading from file \"" << file << "\"."); }

         ifile.close();
         if (ifile.fail()) { WALBERLA_ABORT("Error while closing file \"" << file << "\" for reading."); }
      }
      WALBERLA_MPI_BARRIER();
   }

}


void readBufferMPIIO(const std::string& file, RecvBuffer& buffer)
{
   WALBERLA_NON_MPI_SECTION(){ WALBERLA_ABORT("This function can only be used if waLBerla is build with MPI"); }
   
   uintArray offsetData;

   const MPI_Offset metaDataSize   = numeric_cast< MPI_Offset >( arraySize ) * numeric_cast< MPI_Offset >( MPITrait< uintArray::value_type >::size() );
   const MPI_Offset offsetDataSize = numeric_cast< MPI_Offset >( arraySize ) * numeric_cast< MPI_Offset >( MPITrait< uintArray::value_type >::size() );


   MPI_File mpiFile = MPI_FILE_NULL;
   int result = MPI_File_open(MPIManager::instance()->comm(), const_cast< char* >(file.c_str()), MPI_MODE_RDONLY,
                              MPI_INFO_NULL, &mpiFile);

   if (result != MPI_SUCCESS)
      WALBERLA_ABORT("Error while opening file \"" << file << "\" for reading. MPI Error is \""
                                                   << MPIManager::instance()->getMPIErrorString(result) << "\"");

   MPI_Datatype offsettype;
   MPI_Type_contiguous(arraySize, MPITrait< uint_t >::type(), &offsettype);
   MPI_Type_commit(&offsettype);

   if (MPIManager::instance()->rank() == 0)
   {
      uintArray metaData;
      result = MPI_File_read_at(mpiFile, 0, reinterpret_cast< char* >(metaData.data()), metaData.size(), MPITrait< uint_t >::type(), MPI_STATUS_IGNORE);
      if (result != MPI_SUCCESS)
         WALBERLA_ABORT("Error while reading metaData from file \"" << file << "\". MPI Error is \""
                                                         << MPIManager::instance()->getMPIErrorString(result)
                                                         << "\"");
      if (metaData[0] != uint_c(2)){
         WALBERLA_ABORT("Error file was not written using MPI I/O. MetaData[0] = " << metaData[0] << " != 2");
      }

      if (metaData[1] != uint_c(MPIManager::instance()->numProcesses())){
         WALBERLA_ABORT("Error file was written with different process distribution. Number of processes written file = " << metaData[1] << " != " << MPIManager::instance()->numProcesses())
      }
   }

   // read each process' offset and buffer size from file
   const auto displacement = numeric_cast< MPI_Offset >(mpi::MPIManager::instance()->rank()) * offsetDataSize + metaDataSize;
   result = MPI_File_set_view(mpiFile, displacement, MPITrait< uint_t >::type(), offsettype, const_cast< char* >("native"), MPI_INFO_NULL);

   if (result != MPI_SUCCESS)
      WALBERLA_ABORT("Internal MPI-IO error! MPI Error is \"" << MPIManager::instance()->getMPIErrorString(result)
                                                              << "\"");

   result = MPI_File_read_all(mpiFile, reinterpret_cast< char* >(offsetData.data()), offsetData.size(), MPITrait< uint_t >::type(),
                              MPI_STATUS_IGNORE);

   if (result != MPI_SUCCESS)
      WALBERLA_ABORT("Error while reading from file \"" << file << "\". MPI Error is \""
                                                        << MPIManager::instance()->getMPIErrorString(result)
                                                        << "\"");

   // read each process' data from file (starting at offset)
   MPI_Datatype arraytype;
   MPI_Type_contiguous(int_c(offsetData[1]), MPITrait< mpi::RecvBuffer::ElementType >::type(), &arraytype);
   MPI_Type_commit(&arraytype);

   result = MPI_File_set_view(mpiFile, numeric_cast< MPI_Offset >(offsetData[0] ) + metaDataSize,
                              MPITrait< mpi::RecvBuffer::ElementType >::type(), arraytype,
                              const_cast< char* >("native"), MPI_INFO_NULL);

   if (result != MPI_SUCCESS)
      WALBERLA_ABORT("Internal MPI-IO error! MPI Error is \"" << MPIManager::instance()->getMPIErrorString(result)
                                                              << "\"");

   buffer.resize(offsetData[1]);

   result = MPI_File_read_all(mpiFile, reinterpret_cast< char* >(buffer.ptr()), int_c(buffer.size()),
                              MPITrait< mpi::SendBuffer::ElementType >::type(), MPI_STATUS_IGNORE);

   if (result != MPI_SUCCESS)
      WALBERLA_ABORT("Error while reading from file \"" << file << "\". MPI Error is \""
                                                        << MPIManager::instance()->getMPIErrorString(result)
                                                        << "\"");

   result = MPI_File_close(&mpiFile);

   if (result != MPI_SUCCESS)
      WALBERLA_ABORT("Error while closing file \"" << file << "\". MPI Error is \""
                                                   << MPIManager::instance()->getMPIErrorString(result) << "\"");

   MPI_Type_free(&arraytype);
   MPI_Type_free(&offsettype);

}


void readBuffer(const std::string& file, RecvBuffer& buffer, const bool forceSerialIO)
{
   WALBERLA_NON_MPI_SECTION()
   {
      readBufferNoMPI(file, buffer);
   }

   WALBERLA_MPI_SECTION()
   {

      // use serial I/O for versions of OpenMPI that produce segmentation faults when using MPI-IO with a 3D
      // Cartesian MPI communicator (see waLBerla issue #73)
      if (forceSerialIO || !(MPIManager::instance()->isCommMPIIOValid()) )
      {
         readBufferSerialIO(file, buffer);
      }
      else // use MPI-IO
      {
         readBufferMPIIO(file, buffer);
      }
   }
}

} // namespace walberla::mpi
