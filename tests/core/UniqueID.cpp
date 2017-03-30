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
//! \file UniqueID.cpp
//! \ingroup field
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include "core/debug/TestSubsystem.h"
#include "core/mpi/MPIManager.h"
#include "core/mpi/BufferDataTypeExtensions.h"
#include "core/mpi/SendBuffer.h"
#include "core/mpi/RecvBuffer.h"
#include "core/mpi/MPIWrapper.h"
#include "core/UniqueID.h"

#include <algorithm>

using namespace walberla;

int main( int argc, char** argv )
{
    walberla::debug::enterTestMode();

    WALBERLA_MPI_SECTION(){

    walberla::MPIManager::instance()->initializeMPI( &argc, &argv );

    {
       ///check for unequal ids

        const int size = 100;
        std::vector<size_t> ids;
        ids.reserve(size);
        for (int i = 0; i < size; ++i){
           ids.push_back(UniqueID<int>::create());
        }

        using namespace walberla::mpi;
        SendBuffer buffer;
        buffer << ids;

        WALBERLA_NON_ROOT_SECTION(){
           MPI_Send(buffer.ptr(), static_cast<int>(buffer.size()), MPI_BYTE, 0, 0, MPI_COMM_WORLD);
        }

        WALBERLA_ROOT_SECTION(){
           for (int i = 1; i < walberla::MPIManager::instance()->numProcesses(); ++i){
              RecvBuffer rbuf;
              rbuf.resize(buffer.size());
              MPI_Status status;
              MPI_Recv(rbuf.ptr(), static_cast<int>(buffer.size()), MPI_BYTE, i, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
              std::vector<size_t> vec;
              rbuf >> vec;
              ids.insert(ids.end(), vec.begin(), vec.end());
           }
        }

        std::sort(ids.begin(), ids.end());
        auto it = std::adjacent_find(ids.begin(), ids.end());
        WALBERLA_CHECK_EQUAL(it, ids.end());

    }

    {
       ///check for equal global ids

        const int size = 100;
        std::vector<size_t> ids;
        ids.reserve(size);
        for (int i = 0; i < size; ++i){
           ids.push_back(UniqueID<int>::createGlobal());
        }

        using namespace walberla::mpi;
        SendBuffer buffer;
        buffer << ids;

        WALBERLA_NON_ROOT_SECTION(){
           MPI_Send(buffer.ptr(), static_cast<int>(buffer.size()), MPI_BYTE, 0, 0, MPI_COMM_WORLD);
        }

        WALBERLA_ROOT_SECTION(){
           for (int i = 1; i < walberla::MPIManager::instance()->numProcesses(); ++i){
              RecvBuffer rbuf;
              rbuf.resize(buffer.size());
              MPI_Status status;
              MPI_Recv(rbuf.ptr(), static_cast<int>(buffer.size()), MPI_BYTE, i, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
              std::vector<size_t> vec;
              rbuf >> vec;
              WALBERLA_CHECK_EQUAL(ids, vec);
           }
        }

    }
    }

    return EXIT_SUCCESS;
}
