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
//! \file Gatherv.h
//! \ingroup core
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "Gather.h"
#include "MPIManager.h"
#include "MPIWrapper.h"
#include "SendBuffer.h"
#include "RecvBuffer.h"
#include "core/DataTypes.h"
#include "core/debug/Debug.h"

#include <string>
#include <vector>



namespace walberla {
namespace mpi {



//======================================================================================================================
/*!
 *  \brief Gathers values from MPI processes and stores them into a std::vector
 *
 *  T has to be a native MPI_Datatype
 *
 *  \param values     The values gathered from the process
 *  \param recvRank   The rank of the process receiving the gathered information
 *  \param comm       The MPI communicator used for communication
 *
 *  \returns          A std::vector in with the result on recvRank, else an empty vector
 */
//======================================================================================================================
template< typename T >
std::vector<T> gatherv( const std::vector<T> & values, int recvRank = 0, MPI_Comm comm = MPI_COMM_WORLD )
{
   WALBERLA_NON_MPI_SECTION()
   {
      WALBERLA_ASSERT_EQUAL( recvRank, 0 );
      return std::vector<T>( values );
   }

   int myRank;
   MPI_Comm_rank( comm, &myRank );
   int numProcesses;
   MPI_Comm_size( comm, &numProcesses );
   WALBERLA_ASSERT_GREATER( numProcesses, 0 );

   if( numProcesses == 1 )
   {
      WALBERLA_ASSERT_EQUAL( recvRank, 0 );
      return std::vector<T>( values );
   }

   if( recvRank == myRank )
   {
      std::vector<int> recvCounts = gather( int_c( values.size() ), recvRank, comm );
      std::vector<int> displacements( recvCounts.size(), 0 );

      size_t resultSize =  numeric_cast<size_t>( recvCounts[0] );
      for( uint_t i = 1; i < recvCounts.size(); ++i )
      {
         displacements[i] = displacements[i-1] + recvCounts[i-1];
         resultSize += numeric_cast<size_t>( recvCounts[i] );
      }

      std::vector<T> result( resultSize );
      if( values.empty() )
      {
         if( !result.empty() )
            MPI_Gatherv( nullptr, 0, MPITrait<T>::type(), &(result[0]), &(recvCounts[0]), &(displacements[0]), MPITrait<T>::type(), recvRank, comm );
         else
            MPI_Gatherv( nullptr, 0, MPITrait<T>::type(), nullptr, &(recvCounts[0]), &(displacements[0]), MPITrait<T>::type(), recvRank, comm );
      }
      else
      {
         MPI_Gatherv( const_cast<T*>( &(values[0]) ), int_c( values.size() ), MPITrait<T>::type(), &(result[0]), &(recvCounts[0]),
                      &(displacements[0]), MPITrait<T>::type(), recvRank, comm );
      }

      WALBERLA_ASSERT( std::equal( values.begin(), values.end(), result.begin() + displacements[ uint_c(myRank) ] ) );

      return result;
   }

   gather( int_c( values.size() ), recvRank, comm );

   if( values.empty() )
   {
      MPI_Gatherv( nullptr, 0, MPITrait<T>::type(), nullptr, nullptr, nullptr, MPITrait<T>::type(), recvRank, comm );
   }
   else
   {
      MPI_Gatherv( const_cast<T*>( &(values[0]) ), int_c( values.size() ), MPITrait<T>::type(), nullptr, nullptr, nullptr, MPITrait<T>::type(),
                   recvRank, comm );
   }

   return std::vector<T>();
}

template<>
std::vector< std::string > gatherv( const std::vector< std::string > & values, int recvRank, MPI_Comm comm );



//======================================================================================================================
/*!
 *  \brief Gathers values from MPI processes and stores them into a std::vector on all Processes
 *
 *  T has to be a native MPI_Datatype
 *
 *  \param values     The values gathered from the process
 *  \param comm       The MPI communicator used for communication
 *
 *  \returns          A std::vector in with the result
 */
//======================================================================================================================
template< typename T >
std::vector<T> allGatherv( const std::vector<T> & values, MPI_Comm comm = MPI_COMM_WORLD )
{
   WALBERLA_NON_MPI_SECTION()
   {
      return std::vector<T>( values );
   }

   int myRank;
   MPI_Comm_rank( comm, &myRank );
   int numProcesses;
   MPI_Comm_size( comm, &numProcesses );
   WALBERLA_ASSERT_GREATER( numProcesses, 0 );

   if( numProcesses == 1 )
   {
      return std::vector<T>( values );
   }

   std::vector<int> recvCounts = allGather( int_c( values.size() ), comm );
   std::vector<int> displacements( recvCounts.size(), 0 );

   size_t resultSize =  numeric_cast<size_t>( recvCounts[0] );
   for( uint_t i = 1; i < recvCounts.size(); ++i )
   {
      displacements[i] = displacements[i-1] + recvCounts[i-1];
      resultSize += numeric_cast<size_t>( recvCounts[i] );
   }

   std::vector<T> result( resultSize );
   if( !result.empty() )
   {
      if( values.empty() )
      {
         MPI_Allgatherv( nullptr, 0, MPITrait<T>::type(), &(result[0]), &(recvCounts[0]), &(displacements[0]),
                         MPITrait<T>::type(), comm );
      }
      else
      {
         MPI_Allgatherv( const_cast<T*>( &(values[0]) ), int_c( values.size() ), MPITrait<T>::type(), &(result[0]), &(recvCounts[0]),
                         &(displacements[0]), MPITrait<T>::type(), comm );
      }

      WALBERLA_ASSERT( std::equal( values.begin(), values.end(), result.begin() + displacements[ uint_c(myRank) ] ) );
   }

   return result;
}

template<>
std::vector< std::string > allGatherv( const std::vector< std::string > & values, MPI_Comm comm );



//*******************************************************************************************************************
/*! Gathers the buffer content on a single target process
*
* - every process holds one mpi::SendBuffer ( can have different size on each process )
* - the buffer contents are gathered on process with targetRank
* - buffer contents are sorted by rank and stored consecutively in a mpi::RecvBuffer
*
* \param sendBuffer [in]  sendBuffer with (possibly) different size on each process
* \param recvBuffer [out] recvBuffer which is left unchanged on all processes but targetRank
*                         on targetRank  recvBuffer holds the gathered result
* \param targetRank [in]  rank of the process where data is gathered
* \param comm       [in]  mpi communicator to use
*/
//*******************************************************************************************************************
void gathervBuffer( const mpi::SendBuffer & sendBuffer, mpi::RecvBuffer & recvBuffer, int targetRank = 0, MPI_Comm comm = MPI_COMM_WORLD );



//*******************************************************************************************************************
/*! Almost identical to mpi::gathervBuffer, the only difference: The result is stored on every process.
*/
//*******************************************************************************************************************
void allGathervBuffer( const mpi::SendBuffer & sendBuffer, mpi::RecvBuffer & recvBuffer, MPI_Comm comm = MPI_COMM_WORLD );



} // namespace mpi
} // namespace walberla
