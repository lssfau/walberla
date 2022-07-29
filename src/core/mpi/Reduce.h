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
//! \file Reduce.h
//! \ingroup core
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//
//======================================================================================================================

#pragma once

#include "BufferDataTypeExtensions.h"

#include "core/Abort.h"
#include "core/DataTypes.h"
#include "core/debug/Debug.h"
#include "core/mpi/MPIManager.h"
#include "core/mpi/MPIWrapper.h"

#include "core/math/Vector3.h"

#include <type_traits>
#include <vector>


namespace walberla {
namespace mpi {



enum Operation { MIN, MAX, SUM, PRODUCT, LOGICAL_AND, BITWISE_AND, LOGICAL_OR, BITWISE_OR, LOGICAL_XOR, BITWISE_XOR };

inline MPI_Op toMPI_Op( Operation operation )
{
   switch( operation )
   {
   case MIN:         return MPI_MIN;
   case MAX:         return MPI_MAX;
   case SUM:         return MPI_SUM;
   case PRODUCT:     return MPI_PROD;
   case LOGICAL_AND: return MPI_LAND;
   case BITWISE_AND: return MPI_BAND;
   case LOGICAL_OR:  return MPI_LOR;
   case BITWISE_OR:  return MPI_BOR;
   case LOGICAL_XOR: return MPI_LXOR;
   case BITWISE_XOR: return MPI_BXOR;
   default:          WALBERLA_ABORT( "Unknown operation!" );
   }
#ifdef __IBMCPP__
   return MPI_SUM; // never reached, helps to suppress a warning from the IBM compiler
#endif
}

//======================================================================================================================
/*!
 *  \brief Reduces a value over all processes in-place
 *
 *  T has to be an integer or floating point value
 *
 *  \param value      The value to be reduced
 *  \param operation  The operation to be performed
 *  \param recvRank   The rank of the process receiving the reduced value
 *  \param comm       The MPI communicator used for communication
 */
//======================================================================================================================
template< typename T >
void reduceInplace( T & value, Operation operation, int recvRank = 0, MPI_Comm comm = MPI_COMM_WORLD )
{
   static_assert( std::is_arithmetic<T>::value, "reduceInplace(...) may only by called with integer or floating point types!" );

   WALBERLA_NON_MPI_SECTION()
   {
      WALBERLA_ASSERT_EQUAL( recvRank, 0 );
      return;
   }

   int myRank;
   MPI_Comm_rank( comm, &myRank );

   if( myRank == recvRank )
   {
      MPI_Reduce( MPI_IN_PLACE, &value, 1, MPITrait<T>::type(), toMPI_Op(operation), recvRank, comm );
   }
   else
   {
      MPI_Reduce( &value, nullptr, 1, MPITrait<T>::type(), toMPI_Op(operation), recvRank, comm );
   }
}



//======================================================================================================================
/*!
 *  \brief Reduces a boolean value over all processes in-place
 *
 *  \param value      The boolean value to be reduced
 *  \param operation  The operation to be performed (one of LOGICAL_AND, LOGICAL_OR or LOGICAL_XOR)
 *  \param recvRank   The rank of the process receiving the reduced value
 *  \param comm       The MPI communicator used for communication
 */
//======================================================================================================================
inline void reduceInplace( bool & value, Operation operation, int recvRank = 0, MPI_Comm comm = MPI_COMM_WORLD )
{
   WALBERLA_ASSERT( operation == LOGICAL_AND || operation == LOGICAL_OR || operation == LOGICAL_XOR );

   WALBERLA_NON_MPI_SECTION()
   {
      WALBERLA_ASSERT_EQUAL( recvRank, 0 );
      return;
   }

   int myRank;
   MPI_Comm_rank( comm, &myRank );

   int intValue  = value ? 1 : 0;

   if( myRank == recvRank )
   {
      MPI_Reduce( MPI_IN_PLACE, &intValue, 1, MPITrait<int>::type(), toMPI_Op(operation), recvRank, comm );
   }
   else
   {
      MPI_Reduce( &intValue, nullptr, 1, MPITrait<int>::type(), toMPI_Op(operation), recvRank, comm );
   }

   value = intValue != 0;
}



//======================================================================================================================
/*!
 *  \brief Reduces a value over all processes
 *
 *  T has to be an integer or floating point value
 *
 *  \param value      The value to be reduced
 *  \param operation  The operation to be performed
 *  \param recvRank   The rank of the process receiving the reduced value
 *  \param comm       The MPI communicator used for communication
 *
 *  \returns          The reduced value on recvRank, 0 on all other ranks.
 */
//======================================================================================================================
template< typename T >
T reduce( const T value, Operation operation, int recvRank = 0, MPI_Comm comm = MPI_COMM_WORLD )
{
   static_assert( std::is_arithmetic<T>::value, "reduce(...) may only by called with integer or floating point types!" );

   WALBERLA_NON_MPI_SECTION()
   {
      WALBERLA_ASSERT_EQUAL( recvRank, 0 );
      return value;
   }

   int myRank;
   MPI_Comm_rank( comm, &myRank );

   T result = T(0);

   if( myRank == recvRank )
   {
      MPI_Reduce( const_cast<T*>( &value ), &result, 1, MPITrait<T>::type(), toMPI_Op(operation), recvRank, comm );
   }
   else
   {
      MPI_Reduce( const_cast<T*>( &value ), nullptr, 1, MPITrait<T>::type(), toMPI_Op(operation), recvRank, comm );
   }

   return result;
}

//======================================================================================================================
/*!
 *  \brief Reduces a boolean value over all processes
 *
 *  \param value      The boolean value to be reduced
 *  \param operation  The operation to be performed (one of LOGICAL_AND, LOGICAL_OR or LOGICAL_XOR)
 *  \param recvRank   The rank of the process receiving the reduced value
 *  \param comm       The MPI communicator used for communication
 *
 *  \returns          The reduced boolean value on recvRank, false on all other ranks.
 */
//======================================================================================================================
inline bool reduce( const bool value, Operation operation, int recvRank = 0, MPI_Comm comm = MPI_COMM_WORLD )
{
   WALBERLA_ASSERT( operation == LOGICAL_AND || operation == LOGICAL_OR || operation == LOGICAL_XOR );

   WALBERLA_NON_MPI_SECTION()
   {
      WALBERLA_ASSERT_EQUAL( recvRank, 0 );
      return value;
   }

   int myRank;
   MPI_Comm_rank( comm, &myRank );

   int intValue  = value ? 1 : 0;

   int intResult = 0;

   if( myRank == recvRank )
   {
      MPI_Reduce( &intValue, &intResult, 1, MPITrait<int>::type(), toMPI_Op(operation), recvRank, comm );
   }
   else
   {
      MPI_Reduce( &intValue, nullptr, 1, MPITrait<int>::type(), toMPI_Op(operation), recvRank, comm );
   }

   return intResult != 0;
}



/*!
 *  \brief Reduces values in a std::vector<T> over all processes in-place
 *
 *  T has to be an integer or floating point value
 *
 *  \param values     The values to be reduced
 *  \param operation  The operation to be performed
 *  \param recvRank   The rank of the process receiving the reduced values
 *  \param comm       The MPI communicator used for communication
 */
//======================================================================================================================
template< typename T >
void reduceInplace( std::vector<T> & values, Operation operation, int recvRank = 0, MPI_Comm comm = MPI_COMM_WORLD )
{
   static_assert( std::is_arithmetic<T>::value, "reduceInplace(...) may only by called with integer or floating point types!" );

   WALBERLA_NON_MPI_SECTION()
   {
      WALBERLA_ASSERT_EQUAL( recvRank, 0 );
      return;
   }

   int myRank;
   MPI_Comm_rank( comm, &myRank );

   if( myRank == recvRank )
   {
      MPI_Reduce( MPI_IN_PLACE, values.empty() ? nullptr : &values[0], int_c( values.size() ), MPITrait<T>::type(), toMPI_Op(operation), recvRank, comm );
   }
   else
   {
      MPI_Reduce( values.empty() ? nullptr : &values[0], nullptr, int_c( values.size() ), MPITrait<T>::type(), toMPI_Op(operation), recvRank, comm );
   }
}


//======================================================================================================================
/*!
 *  \brief Reduces boolean values in a std::vector<bool> over all processes in-place
 *
 *  Specialization of reduceInplace<T>
 *
 *  \param values     The boolean values to be reduced
 *  \param operation  The operation to be performed (one of  BITWISE_AND, BITWISE_OR or BITWISE_XOR)
 *  \param recvRank   The rank of the process receiving the reduced values
 *  \param comm       The MPI communicator used for communication
 */
//======================================================================================================================
inline void reduceInplace( std::vector<bool> & values, Operation operation, int recvRank = 0, MPI_Comm comm = MPI_COMM_WORLD )
{
   WALBERLA_ASSERT( operation == BITWISE_AND || operation == BITWISE_OR || operation == BITWISE_XOR );

   WALBERLA_NON_MPI_SECTION()
   {
      WALBERLA_ASSERT_EQUAL( recvRank, 0 );
      return;
   }

   int myRank;
   MPI_Comm_rank( comm, &myRank );

   std::vector<uint8_t> sendBuffer;

   convert( values, sendBuffer );

   if( myRank == recvRank )
   {
      MPI_Reduce( MPI_IN_PLACE, sendBuffer.empty() ? nullptr : &sendBuffer[0], int_c( sendBuffer.size() ), MPITrait<uint8_t>::type(), toMPI_Op(operation), recvRank, comm );
      size_t size = values.size();
      convert( sendBuffer, values );
      values.resize(size);
   }
   else
   {
      MPI_Reduce( sendBuffer.empty() ? nullptr : &sendBuffer[0], nullptr, int_c( sendBuffer.size() ), MPITrait<uint8_t>::type(), toMPI_Op(operation), recvRank, comm );
   }
}


/*!
 *  \brief Reduces values in a math::Vector3<T> over all processes in-place
 *
 *  T has to be an integer or floating point value
 *
 *  \param values     The values to be reduced
 *  \param operation  The operation to be performed
 *  \param recvRank   The rank of the process receiving the reduced values
 *  \param comm       The MPI communicator used for communication
 */
//======================================================================================================================
template< typename T >
void reduceInplace( math::Vector3<T> & values, Operation operation, int recvRank = 0, MPI_Comm comm = MPI_COMM_WORLD )
{
   static_assert( std::is_arithmetic<T>::value, "reduceInplace(...) may only by called with integer or floating point types!" );

   WALBERLA_NON_MPI_SECTION()
   {
      WALBERLA_ASSERT_EQUAL( recvRank, 0 );
      return;
   }

   int myRank;
   MPI_Comm_rank( comm, &myRank );

   if( myRank == recvRank )
   {
      MPI_Reduce( MPI_IN_PLACE, values.data(), 3, MPITrait<T>::type(), toMPI_Op(operation), recvRank, comm );
   }
   else
   {
      MPI_Reduce( values.data(), nullptr, 3, MPITrait<T>::type(), toMPI_Op(operation), recvRank, comm );
   }
}

/*!
 *  \brief Reduces booleans in a math::Vector3 over all processes in-place
 *
 *  T has to be an integer or floating point value
 *
 *  \param values     The values to be reduced
 *  \param operation  The operation to be performed
 *  \param recvRank   The rank of the process receiving the reduced values
 *  \param comm       The MPI communicator used for communication
 */
//======================================================================================================================
inline void reduceInplace( math::Vector3<bool> & values, Operation operation, int recvRank = 0, MPI_Comm comm = MPI_COMM_WORLD )
{
   WALBERLA_ASSERT( operation == LOGICAL_AND || operation == LOGICAL_OR || operation == LOGICAL_XOR );

   WALBERLA_NON_MPI_SECTION()
   {
      WALBERLA_ASSERT_EQUAL( recvRank, 0 );
      return;
   }

   int myRank;
   MPI_Comm_rank( comm, &myRank );

   math::Vector3<int> intValues{values[0] ? 1 : 0, values[1] ? 1 : 0, values[2] ? 1 : 0};

   if( myRank == recvRank )
   {
      MPI_Reduce( MPI_IN_PLACE, intValues.data(), 3, MPITrait<int>::type(), toMPI_Op(operation), recvRank, comm );
   }
   else
   {
      MPI_Reduce( intValues.data(), nullptr, 3, MPITrait<int>::type(), toMPI_Op(operation), recvRank, comm );
   }

   for(uint_t i = 0; i < 3; ++i)
      values[i] = intValues[i] != 0;
}


//======================================================================================================================
/*!
 *  \brief Reduces values in a math::Vector3<T> over all processes
 *
 *  T has to be an integer or floating point value
 *
 *  \param value      The value to be reduced
 *  \param operation  The operation to be performed
 *  \param recvRank   The rank of the process receiving the reduced value
 *  \param comm       The MPI communicator used for communication
 *
 *  \returns          The reduced value on recvRank, 0 on all other ranks.
 */
//======================================================================================================================
template< typename T >
math::Vector3<T> reduce( const math::Vector3<T> & values, Operation operation, int recvRank = 0, MPI_Comm comm = MPI_COMM_WORLD )
{
   static_assert( std::is_arithmetic<T>::value, "reduce(...) may only by called with integer or floating point types!" );

   WALBERLA_NON_MPI_SECTION()
   {
      WALBERLA_ASSERT_EQUAL( recvRank, 0 );
      return values;
   }

   int myRank;
   MPI_Comm_rank( comm, &myRank );

   math::Vector3<T> result{T(0)};

   if( myRank == recvRank )
   {
      MPI_Reduce( const_cast<T*>( values.data() ), result.data(), 3, MPITrait<T>::type(), toMPI_Op(operation), recvRank, comm );
   }
   else
   {
      MPI_Reduce( const_cast<T*>( values.data() ), nullptr, 3, MPITrait<T>::type(), toMPI_Op(operation), recvRank, comm );
   }

   return result;
}

//======================================================================================================================
/*!
 *  \brief Reduces boolean values in a Vector3 over all processes
 *
 *  \param value      The boolean value to be reduced
 *  \param operation  The operation to be performed (one of LOGICAL_AND, LOGICAL_OR or LOGICAL_XOR)
 *  \param recvRank   The rank of the process receiving the reduced value
 *  \param comm       The MPI communicator used for communication
 *
 *  \returns          The reduced boolean value on recvRank, false on all other ranks.
 */
//======================================================================================================================
inline math::Vector3<bool> reduce( const math::Vector3<bool> & values, Operation operation, int recvRank = 0, MPI_Comm comm = MPI_COMM_WORLD )
{
   WALBERLA_ASSERT( operation == LOGICAL_AND || operation == LOGICAL_OR || operation == LOGICAL_XOR );

   WALBERLA_NON_MPI_SECTION()
   {
      WALBERLA_ASSERT_EQUAL( recvRank, 0 );
      return values;
   }

   int myRank;
   MPI_Comm_rank( comm, &myRank );

   math::Vector3<int> intValues{values[0] ? 1 : 0, values[1] ? 1 : 0, values[2] ? 1 : 0};

   math::Vector3<bool> results(false);

   if( myRank == recvRank )
   {
      MPI_Reduce( MPI_IN_PLACE, intValues.data(), 3, MPITrait<int>::type(), toMPI_Op(operation), recvRank, comm );

      for(uint_t i = 0; i < 3; ++i)
         results[i] = intValues[i] != 0;
   }
   else
   {
      MPI_Reduce( intValues.data(), nullptr, 3, MPITrait<int>::type(), toMPI_Op(operation), recvRank, comm );
   }

   return results;
}


//======================================================================================================================
/*!
 *  \brief Reduces a value over all processes
 *
 *  T has to be an integer or floating point value
 *
 *  \param value      The value to be reduced
 *  \param operation  The operation to be performed
 *  \param comm       The MPI communicator used for communication
 *
 *  \returns          The reduced value on recvRank, 0 on all other ranks.
 */
//======================================================================================================================
template< typename T >
T allReduce( const T & value, Operation operation, MPI_Comm comm = MPI_COMM_WORLD )
{
   static_assert( std::is_arithmetic<T>::value, "allReduce(...) may only by called with integer or floating point types!" );

   WALBERLA_NON_MPI_SECTION() { return value; }

   T result;
   MPI_Allreduce( const_cast<T*>( &value ), &result, 1, MPITrait<T>::type(), toMPI_Op(operation), comm );
   return result;
}



//======================================================================================================================
/*!
 *  \brief Reduces a boolean value over all processes
 *
 *  T has to be a boolean value
 *
 *  \param value      The boolean value to be reduced
 *  \param operation  The operation to be performed
 *  \param comm       The MPI communicator used for communication
 *
 *  \returns          The reduced value on recvRank, 0 on all other ranks.
 */
//======================================================================================================================
inline bool allReduce( const bool value, Operation operation, MPI_Comm comm = MPI_COMM_WORLD )
{
   WALBERLA_ASSERT( operation == LOGICAL_AND || operation == LOGICAL_OR || operation == LOGICAL_XOR );

   WALBERLA_NON_MPI_SECTION() { return value; }

   int intValue = value ? 1 : 0;

   MPI_Allreduce( MPI_IN_PLACE, &intValue, 1, MPITrait<int>::type(), toMPI_Op(operation), comm );

   return intValue != 0;
}



//======================================================================================================================
/*!
 *  \brief Reduces a value over all processes in-place
 *
 *  T has to be an integer or floating point value
 *
 *  \param value      The value to be reduced
 *  \param operation  The operation to be performed
 *  \param comm       The MPI communicator used for communication
 */
//======================================================================================================================
template< typename T >
void allReduceInplace( T & value, Operation operation, MPI_Comm comm = MPI_COMM_WORLD )
{
   static_assert( std::is_arithmetic<T>::value, "allReduceInplace(...) may only by called with integer or floating point types!" );

   WALBERLA_NON_MPI_SECTION() { return; }

   MPI_Allreduce( MPI_IN_PLACE, &value, 1, MPITrait<T>::type(), toMPI_Op(operation), comm );
}


//======================================================================================================================
/*!
 *  \brief Reduces a boolean value over all processes in-place
 *
 *  T has to be a boolean value
 *
 *  \param value      The boolean value to be reduced
 *  \param operation  The operation to be performed (one of LOGICAL_AND, LOGICAL_OR or LOGICAL_XOR)
 *  \param comm       The MPI communicator used for communication
 */
//======================================================================================================================
inline void allReduceInplace( bool & value, Operation operation, MPI_Comm comm = MPI_COMM_WORLD )
{
   WALBERLA_ASSERT( operation == LOGICAL_AND || operation == LOGICAL_OR || operation == LOGICAL_XOR );

   WALBERLA_NON_MPI_SECTION() { return; }

   int intValue = value ? 1 : 0;

   MPI_Allreduce( MPI_IN_PLACE, &intValue, 1, MPITrait<int>::type(), toMPI_Op(operation), comm );

   value = intValue != 0;
}



//======================================================================================================================
/*!
 *  \brief Reduces values in a std::vector<T> over all processes in-place
 *
 *  T has to be an integer or floating point value
 *
 *  \param values     The values to be reduced
 *  \param operation  The operation to be performed
 *  \param comm       The MPI communicator used for communication
 */
//======================================================================================================================
template< typename T >
void allReduceInplace( std::vector<T> & values, Operation operation, MPI_Comm comm = MPI_COMM_WORLD )
{
   static_assert( std::is_arithmetic<T>::value, "allReduceInplace(...) may only by called with integer or floating point types!" );

   WALBERLA_NON_MPI_SECTION() { return; }

   MPI_Allreduce( MPI_IN_PLACE, values.empty() ? nullptr : &values[0], int_c( values.size() ), MPITrait<T>::type(), toMPI_Op(operation), comm );
}



//======================================================================================================================
/*!
 *  \brief Reduces values in a std::vector<bool> over all processes in-place
 *
 *  Specialization of allReduceInplace<T>
 *
 *  \param bools      The boolean values to be reduced
 *  \param operation  The operation to be performed (one of BITWISE_AND, BITWISE_OR or BITWISE_XOR)
 *  \param comm       The MPI communicator used for communication
 */
//======================================================================================================================
inline void allReduceInplace( std::vector<bool> & bools, Operation operation, MPI_Comm comm = MPI_COMM_WORLD  )
{
   WALBERLA_ASSERT( operation == BITWISE_AND || operation == BITWISE_OR || operation == BITWISE_XOR );

   WALBERLA_NON_MPI_SECTION() { return; }

   std::vector<uint8_t> sendBuffer;

   convert( bools, sendBuffer );
   MPI_Allreduce( MPI_IN_PLACE, sendBuffer.empty() ? nullptr : &sendBuffer[0], int_c( sendBuffer.size() ), MPITrait<uint8_t>::type(), toMPI_Op(operation), comm );
   auto size = bools.size();
   convert(sendBuffer, bools);
   bools.resize(size);
}


//======================================================================================================================
/*!
 *  \brief Reduces values in math::Vector3<T> over all processes in-place
 *
 *  T has to be an integer or floating point value
 *
 *  \param values     The values to be reduced
 *  \param operation  The operation to be performed
 *  \param comm       The MPI communicator used for communication
 */
//======================================================================================================================
template< typename T >
void allReduceInplace( math::Vector3<T> & values, Operation operation, MPI_Comm comm = MPI_COMM_WORLD )
{
   static_assert( std::is_arithmetic<T>::value, "allReduceInplace(...) may only by called with integer or floating point types!" );

   WALBERLA_NON_MPI_SECTION() { return; }

   MPI_Allreduce( MPI_IN_PLACE, values.data(), 3, MPITrait<T>::type(), toMPI_Op(operation), comm );
}


//======================================================================================================================
/*!
 *  \brief Reduces boolean values in math::Vector3 over all processes in-place
 *
 *  T has to be an integer or floating point value
 *
 *  \param values     The values to be reduced
 *  \param operation  The operation to be performed
 *  \param comm       The MPI communicator used for communication
 */
//======================================================================================================================
inline void allReduceInplace( math::Vector3<bool> & bools, Operation operation, MPI_Comm comm = MPI_COMM_WORLD )
{
   WALBERLA_ASSERT( operation == LOGICAL_AND || operation == LOGICAL_OR || operation == LOGICAL_XOR );

   WALBERLA_NON_MPI_SECTION()
   {
      return;
   }

   math::Vector3<int> intValues{bools[0] ? 1 : 0, bools[1] ? 1 : 0, bools[2] ? 1 : 0};

   MPI_Allreduce( MPI_IN_PLACE, intValues.data(), 3, MPITrait<int>::type(), toMPI_Op(operation), comm );

   for(uint_t i = 0; i < 3; ++i)
   {
      bools[i] = intValues[i] != 0;
   }
}

} // namespace mpi
} // namespace walberla
