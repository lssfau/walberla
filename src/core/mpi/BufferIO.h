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
//! \file BufferIO.h
//! \ingroup core
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/DataTypes.h"
#include "core/logging/Logging.h"

#include "core/mpi/MPIManager.h"
#include "core/mpi/MPIWrapper.h"
#include "core/mpi/RecvBuffer.h"
#include "core/mpi/SendBuffer.h"

namespace walberla::mpi::internal {
  inline uint_t getProcessOffset(uint_t bufferDataSize){
    // get the position at which a process writes its data (offset from the beginning of the file)
    uint_t exscanResult;
    MPI_Exscan(&bufferDataSize, &exscanResult, 1, MPITrait< uint_t >::type(), MPI_SUM, MPIManager::instance()->comm());
    if (MPIManager::instance()->rank() == 0) exscanResult = uint_c(0);
    return uint_c(mpi::MPIManager::instance()->numProcesses() * 2 * MPITrait< uint_t >::size()) + exscanResult;
  }
}


namespace walberla::mpi {
 
//======================================================================================================================
/*!
 *  \brief Writes contents of local buffer to a single binary file without MPI. Has to be used when waLBerla is compiled without MPI.
 *
 *  \param file     Name of the output file
 *  \param buffer   A SendBuffer which contains the data which should be written
 */
//======================================================================================================================
  void writeBufferNoMPI(const std::string& file, SendBuffer& buffer);

//======================================================================================================================
/*!
 *  \brief Writes contents of local buffer to a single binary file. The data is written process by process.
 *
 *  \param file     Name of the output file
 *  \param buffer   A SendBuffer which contains the data which should be written
 */
//======================================================================================================================
  void writeBufferSerialIO(const std::string& file, SendBuffer& buffer);

//======================================================================================================================
/*!
 *  \brief Writes contents of local buffer to a single binary file using MPIIO.
 *
 *  \param file     Name of the output file
 *  \param buffer   A SendBuffer which contains the data which should be written
 */
//======================================================================================================================
  void writeBufferMPIIO(const std::string& file, SendBuffer& buffer);

//======================================================================================================================
/*!
 *  \brief Writes contents of local buffer to a single binary file. It is internally decided weather to use MPIIO or not
 *
 *  \param file           Name of the output file
 *  \param buffer         A SendBuffer which contains the data which should be written
 *  \param forceSerialIO  This parameter can be set to prevent the usage of MPIIO
 */
//======================================================================================================================
  void writeBuffer(const std::string& file, SendBuffer& buffer, const bool forceSerialIO = true);


//======================================================================================================================
/*!
 *  \brief Counterpart to writeBufferNoMPI
 *
 *  \param file           Name of the input file
 *  \param buffer         A RecvBuffer to which the content of the file will be written to
 */
//======================================================================================================================
  void readBufferNoMPI(const std::string& file, RecvBuffer& buffer);

//======================================================================================================================
/*!
 *  \brief Counterpart to writeBufferSerialIO. Both functions need to be called with the same process distribution
 *
 *  \param file           Name of the input file
 *  \param buffer         A RecvBuffer to which the content of the file will be written to
 */
//======================================================================================================================
  void readBufferSerialIO(const std::string& file, RecvBuffer& buffer);

//======================================================================================================================
/*!
 *  \brief Counterpart to writeBufferMPIIO. Both functions need to be called with the same process distribution
 *
 *  \param file           Name of the input file
 *  \param buffer         A RecvBuffer to which the content of the file will be written to
 */
//======================================================================================================================
  void readBufferMPIIO(const std::string& file, RecvBuffer& buffer);

//======================================================================================================================
/*!
 *  \brief Counterpart to writeBuffer. Both functions need to be called with the same process distribution
 *
 *  \param file           Name of the input file
 *  \param buffer         A RecvBuffer to which the content of the file will be written to
 *  \param forceSerialIO  This parameter can be set to prevent the usage of MPIIO
 */
//======================================================================================================================
  void readBuffer(const std::string& file, RecvBuffer& buffer, const bool forceSerialIO = true);


} // namespace walberla::mpi


