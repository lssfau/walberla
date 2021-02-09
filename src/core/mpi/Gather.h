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
//! \file Gather.h
//! \ingroup core
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/DataTypes.h"
#include "core/debug/Debug.h"
#include "core/mpi/MPIManager.h"
#include "core/mpi/MPIWrapper.h"

#include <vector>


namespace walberla {
namespace mpi {


//======================================================================================================================
/*!
 *  \brief Gathers values from MPI processes and stores them into a std::vector
 *
 *  T has to be a native MPI_Datatype
 *
 *  \param value      The value gathered from the process
 *  \param recvRank   The rank of the process receiving the gathered information
 *  \param comm       The MPI communicator used for communication
 *
 *  \returns          A std::vector in with the result on recvRank, else an empty vector
 */
//======================================================================================================================
template< typename T >
std::vector<T> gather( T value, int recvRank = 0, MPI_Comm comm = MPI_COMM_WORLD )
{
   WALBERLA_NON_MPI_SECTION()
   {
      WALBERLA_ASSERT_EQUAL( recvRank, 0 );
      return std::vector<T>( 1u, value );
   }

   int myRank;
   MPI_Comm_rank( comm, &myRank );
   int numProcesses;
   MPI_Comm_size( comm, &numProcesses );
   WALBERLA_ASSERT_GREATER( numProcesses, 0 );

   if( numProcesses == 1 )
   {
      WALBERLA_ASSERT_EQUAL( recvRank, 0 );
      return std::vector<T>( 1u, value );
   }

   if( recvRank == myRank )
   {
      std::vector<T> result( static_cast<size_t>( numProcesses ) );
      MPI_Gather( &value, 1, MPITrait<T>::type(), &(result[0]), 1, MPITrait<T>::type(),
                  recvRank, comm );

      WALBERLA_ASSERT( isIdentical( result[ uint_c(recvRank) ], value ) );

      return result;
   }

   MPI_Gather( &value, 1, MPITrait<T>::type(), nullptr, 1, MPITrait<T>::type(),
               recvRank, comm );

   return std::vector<T>();

}


//======================================================================================================================
/*!
 *  \brief Gathers values from MPI processes and stores them into a std::vector on all Processes
 *
 *  T has to be a native MPI_Datatype
 *
 *  \param value      The value gathered from the process
 *  \param comm       The MPI communicator used for communication
 *
 *  \returns          A std::vector in with the result
 */
//======================================================================================================================
template< typename T >
std::vector<T> allGather( T value, MPI_Comm comm = MPI_COMM_WORLD )
{
   WALBERLA_NON_MPI_SECTION()
   {
      return std::vector<T>( 1u, value );
   }

   int myRank;
   MPI_Comm_rank( comm, &myRank );
   int numProcesses;
   MPI_Comm_size( comm, &numProcesses );
   WALBERLA_ASSERT_GREATER( numProcesses, 0 );

   if( numProcesses == 1 )
   {
      return std::vector<T>( 1u, value );
   }

   std::vector<T> result( static_cast<size_t>( numProcesses ) );
   MPI_Allgather( &value, 1, MPITrait<T>::type(), &(result[0]), 1, MPITrait<T>::type(),
                  comm );

   WALBERLA_ASSERT( isIdentical( result[ uint_c(myRank) ], value ) );

   return result;
}

} // namespace mpi
} // namespace walberla
