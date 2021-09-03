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
//! \file MPITextFile.cpp
//! \ingroup core
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//
//======================================================================================================================

#include "MPITextFile.h"

#include "core/logging/Logging.h"
#include "core/mpi/Reduce.h"

#include <fstream>

namespace walberla
{
namespace mpi
{
//======================================================================================================================
/*!
 *  \brief Writes file using MPI IO with each process providing a part of it
 *
 *  This method has the be called collectively by all the processes in comm. The file will be assembled in the order of
 *  the ranks of the calling processes.
 *
 *  \param filename          The name of the file to be written
 *  \param processLocalPart  The part of the file belonging to the calling process (size may differ among processes)
 *  \param comm              The MPI communicator used for communication
 */
//======================================================================================================================
void writeMPITextFile(const std::string& filename, const std::string& processLocalPart,
                      const MPI_Comm comm /*= MPI_COMM_WORLD*/)
{
   WALBERLA_NON_MPI_SECTION()
   {
      std::ofstream ofs(filename.c_str());
      ofs << processLocalPart;
      ofs.close();
   }

   WALBERLA_MPI_SECTION()
   {
      int rank;
      int numProcesses;
      MPI_Comm_rank(comm, &rank);
      MPI_Comm_size(comm, &numProcesses);

      // use serial I/O for versions of OpenMPI that produce segmentation faults when using MPI-IO with a 3D Cartesian
      // MPI communicator (see waLBerla issue #73)
      if (!MPIManager::instance()->isCommMPIIOValid())
      {
         std::ofstream ofs;

         for (int i = 0; i != numProcesses; ++i)
         {
            if (i == rank)
            {
               if (rank == 0) { ofs.open(filename.c_str()); }
               else
               {
                  ofs.open(filename.c_str(), std::ofstream::app);
               }
               ofs << processLocalPart;
               ofs.close();
            }
            WALBERLA_MPI_BARRIER();
         }
      }
      else
      {
         if (processLocalPart.size() > numeric_cast< std::string::size_type >(std::numeric_limits< int >::max()))
            WALBERLA_ABORT("writeMPITextFile does not support more than " << std::numeric_limits< int >::max()
                                                                          << " characters per process!");

         MPI_File mpiFile;
         int result = MPI_SUCCESS;
         result     = MPI_File_open(comm, const_cast< char* >(filename.c_str()), MPI_MODE_WRONLY | MPI_MODE_CREATE,
                                MPI_INFO_NULL, &mpiFile);
         if (result != MPI_SUCCESS)
            WALBERLA_ABORT("Error while opening file \"" << filename << "\" for writing. MPI Error is \""
                                                         << MPIManager::instance()->getMPIErrorString(result) << "\"");

         const MPI_Offset filesize = numeric_cast< MPI_Offset >(allReduce(processLocalPart.size(), mpi::SUM, comm));

         size_t exscanResult;
         size_t input = processLocalPart.size();
         MPI_Exscan(&input, &exscanResult, 1, MPITrait< size_t >::type(), MPI_SUM, comm);
         if (rank == 0) exscanResult = size_t(0);

         const MPI_Offset offset = numeric_cast< MPI_Offset >(exscanResult);

         MPI_File_set_size(mpiFile, filesize);

         MPI_Datatype arraytype;
         MPI_Type_contiguous(int_c(processLocalPart.size()), MPITrait< char >::type(), &arraytype);
         MPI_Type_commit(&arraytype);

         result = MPI_File_set_view(mpiFile, offset, MPITrait< char >::type(), arraytype, const_cast< char* >("native"),
                                    MPI_INFO_NULL);

         if (result != MPI_SUCCESS)
            WALBERLA_ABORT("Internal MPI-IO error! MPI Error is \"" << MPIManager::instance()->getMPIErrorString(result)
                                                                    << "\"");

         result = MPI_File_write_all(mpiFile, const_cast< char* >(processLocalPart.c_str()),
                                     int_c(processLocalPart.size()), MPITrait< char >::type(), MPI_STATUS_IGNORE);

         if (result != MPI_SUCCESS)
            WALBERLA_ABORT("Error while writing to file \"" << filename << "\". MPI Error is \""
                                                            << MPIManager::instance()->getMPIErrorString(result)
                                                            << "\"");

         result = MPI_File_close(&mpiFile);

         if (result != MPI_SUCCESS)
            WALBERLA_ABORT("Error while closing file \"" << filename << "\". MPI Error is \""
                                                         << MPIManager::instance()->getMPIErrorString(result) << "\"");

         MPI_Type_free(&arraytype);
      }
   }
}

} // namespace mpi
} // namespace walberla