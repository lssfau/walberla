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
//! \file DeviceSelectMPI.cpp
//! \ingroup cuda
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#include "DeviceSelectMPI.h"
#include "core/mpi/MPIWrapper.h"
#include "cuda/ErrorChecking.h"
#include "core/logging/Logging.h"

namespace walberla {
namespace cuda {

#if MPI_VERSION == 2 || MPI_VERSION == 1

void selectDeviceBasedOnMpiRank() {
   WALBERLA_ABORT("Your MPI implementation is tool old - it does not support CUDA device selection based on MPI rank");
}

#else

void selectDeviceBasedOnMpiRank()
{
#ifdef WALBERLA_BUILD_WITH_MPI
   int deviceCount;
   WALBERLA_CUDA_CHECK( cudaGetDeviceCount( &deviceCount ));
   WALBERLA_LOG_INFO_ON_ROOT( "Selecting CUDA device depending on MPI Rank" );

   MPI_Info info;
   MPI_Info_create( &info );
   MPI_Comm newCommunicator;
   MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, info, &newCommunicator );

   int processesOnNode;
   int rankOnNode;
   MPI_Comm_size( newCommunicator, &processesOnNode );
   MPI_Comm_rank( newCommunicator, &rankOnNode );

   if ( deviceCount == processesOnNode )
   {
      WALBERLA_CUDA_CHECK( cudaSetDevice( rankOnNode ));
   }
   else if ( deviceCount > processesOnNode )
   {
      WALBERLA_LOG_WARNING( "Not using all available GPUs on node. Processes on node "
                               << processesOnNode << " available GPUs on node " << deviceCount );
      WALBERLA_CUDA_CHECK( cudaSetDevice( rankOnNode ));
   }
   else
   {
      WALBERLA_LOG_WARNING( "Too many processes started per node - should be one per GPU. Number of processes per node "
                               << processesOnNode << ", available GPUs on node " << deviceCount );
      WALBERLA_CUDA_CHECK( cudaSetDevice( rankOnNode % deviceCount ));
   }
#endif
}

#endif


} // namespace cuda
} // namespace walberla