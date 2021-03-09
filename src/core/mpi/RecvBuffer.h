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
//! \file RecvBuffer.h
//! \ingroup core
//! \author Klaus Iglberger
//! \brief Implementation of a MPI receive buffer.
//!
//! Copyright (C) 2009 Klaus Iglberger
//! Taken from "pe Physics Engine" with small changes
//
//======================================================================================================================

#pragma once

#include "waLBerlaDefinitions.h"
#include "BufferSizeTrait.h"
#include "SendBuffer.h"

#include "core/debug/Debug.h"
#include "core/Sanitizer.h"

#include <algorithm>
#include <cstring>
#include <type_traits>

namespace walberla {
namespace mpi {

//======================================================================================================================
//
//  CLASS DEFINITION
//
//======================================================================================================================

//**********************************************************************************************************************
/*!\brief Implementation of a MPI receive buffer.
// \ingroup mpi
//
// The RecvBuffer class is a special purpose implementation for the MPI communication
// functionality. It offers a convenient and safe access to the received data values
// even for mixed-type communication. The following example gives an impression of the
// usage of the RecvBuffer class:

   \code
   using namespace mpi;

   // Preparing a receive buffer for an incoming message of 32 bytes
   RecvBuffer buffer;
   buffer.resize( 32 );

   // Receiving a MPI message from process 0 in a blocking MPI_Recv() function
   MPI_Status status;
   MPI_Recv( buffer.ptr(), 32, MPI_BYTE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, status );

   // Extracting a double and an integer from the buffer
   double d;
   int i;
   buffer >> d >> i;
   \endcode
// For another example see also the Buffer Unit Test in File BufferTest.cpp
//
// Note that the order of data values in the receive buffer is depending on the message
// sent via MPI. See also the SendBuffer class description for the sender side of the MPI
// communication.
*/
template< typename T = unsigned char >  // Element type
class GenericRecvBuffer
{
public:
   //**Type definitions*************************************************************************************************
   using ElementType = T;  //!< Type of the receive buffer elements.
   //*******************************************************************************************************************

   //**Constructors*****************************************************************************************************
   /*!\name Constructors */
   //@{
   explicit inline GenericRecvBuffer();
            inline GenericRecvBuffer( const GenericRecvBuffer & rb );
   template< typename G >
   explicit inline GenericRecvBuffer( GenericSendBuffer<T,G> & sb );
   //@}
   //*******************************************************************************************************************

   //**Destructor*******************************************************************************************************
   /*!\name Destructor */
   //@{
   inline ~GenericRecvBuffer();
   //@}
   //*******************************************************************************************************************

   //**Assignment operator**********************************************************************************************
   /*!\name Assignment operator */
   //@{
   GenericRecvBuffer& operator=( const GenericRecvBuffer& sb );
   template< typename G >
   GenericRecvBuffer& operator=( const GenericSendBuffer<T,G> & sb );
   //@}
   //*******************************************************************************************************************

   //**Get functions****************************************************************************************************
   /*!\name Get functions */
   //@{
   inline size_t maxSize () const;
   inline size_t size    () const;
   inline size_t capacity() const;
   inline bool   isEmpty () const;
   //@}
   //*******************************************************************************************************************

   //**Operators********************************************************************************************************
   /*!\name Operators */
   //@{
   template< typename V >
   typename std::enable_if< std::is_arithmetic<V>::value || std::is_enum<V>::value,
                              GenericRecvBuffer& >::type
   operator>>( V& value );
   //@}
   //*******************************************************************************************************************

   //**Utility functions************************************************************************************************
   /*!\name Utility functions */
   //@{
                          inline T*   ptr    () const;
                          inline void reserve( size_t newCapacity );
                          inline void resize ( size_t newSize     );
   template< typename V >        void peek   ( V& value ) const;
                          inline T *  skip   ( size_t elements    );
                          inline void clear  ();
                          inline void reset  ();
                          inline void readDebugMarker( const char * marker );
   //@}
   //*******************************************************************************************************************

private:
   //**Member variables*************************************************************************************************
   /*!\name Member variables */
   //@{
   size_t capacity_;  //!< The current size of the receive buffer.
   T* begin_;         //!< Pointer to the first element of the receive buffer.
   T* cur_;           //!< Pointer to the current element of the receive buffer.
   T* end_;           //!< Pointer to the last element of the receive buffer.
   //@}
   //*******************************************************************************************************************

   //**Utility functions************************************************************************************************
   /*!\name Utility functions */
   //@{
   template< typename V >
   typename std::enable_if< std::is_arithmetic<V>::value || std::is_enum<V>::value,
                              GenericRecvBuffer& >::type
   get( V& value );
   //@}
   //*******************************************************************************************************************



   //**Compile time checks**********************************************************************************************
   static_assert( std::is_arithmetic<T>::value, "SendBuffer<T>: T has to be native datatype" ) ;
   //*******************************************************************************************************************
};
//**********************************************************************************************************************

using RecvBuffer = GenericRecvBuffer<>;



//======================================================================================================================
//
//  CONSTRUCTORS
//
//======================================================================================================================

//**********************************************************************************************************************
/*!\brief Standard constructor for RecvBuffer.
*/
template< typename T >  // Element type
inline GenericRecvBuffer<T>::GenericRecvBuffer()
   : capacity_( 0    )  // Capacity of the receive buffer
   , begin_   ( nullptr )  // Pointer to the first element
   , cur_     ( nullptr )  // Pointer to the current element
   , end_     ( nullptr )  // Pointer to the last element
{}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\brief Copy constructor for RecvBuffer.
//
// \param rb The receive buffer to be copied.
*/
template< typename T >  // Element type
inline GenericRecvBuffer<T>::GenericRecvBuffer( const GenericRecvBuffer& rb )
   : capacity_( rb.size()        )  // Capacity of the receive buffer
   , begin_   ( new T[capacity_] )  // Pointer to the first element
   , cur_     ( begin_           )  // Pointer to the current element
   , end_     ( begin_+capacity_ )  // Pointer to the last element
{
   for( size_t i=0; i<capacity_; ++i )
      cur_[i] = rb.cur_[i];
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\brief Constructor for RecvBuffer.
//
// \param sb The send buffer whose content is transfered to this receive buffer.
*/
template< typename T >  // Element type
template< typename G >
inline GenericRecvBuffer<T>::GenericRecvBuffer( GenericSendBuffer<T,G> & sb )
   : capacity_( sb.capacity() )
   , begin_( sb.begin_ )
   , cur_( sb.begin_ )
   , end_( sb.end_ )
{
   sb.begin_ = new T[0];
   sb.cur_   = sb.begin_;
   sb.end_   = sb.begin_;
}
//**********************************************************************************************************************



//======================================================================================================================
//
//  DESTRUCTOR
//
//======================================================================================================================

//**********************************************************************************************************************
/*!\brief Destructor for RecvBuffer.
*/
template< typename T >  // Element type
inline GenericRecvBuffer<T>::~GenericRecvBuffer()
{
   delete [] begin_;
}
//**********************************************************************************************************************




//======================================================================================================================
//
//  ASSIGNMENT OPERATOR
//
//======================================================================================================================

//**********************************************************************************************************************
/*!\brief Copy assignment operator for RecvBuffer.
//
// \param rb The receive buffer to be copied.
// \return Reference to the assigned receive buffer.
*/
template< typename T >  // Element type
GenericRecvBuffer<T>& GenericRecvBuffer<T>::operator=( const GenericRecvBuffer& rb )
{
   if( &rb == this ) return *this;

   if( rb.size() > capacity_ ) {
      T* newBegin( new T[rb.size()] );
      end_ = std::copy( rb.cur_, rb.end_, newBegin );
      std::swap( begin_, newBegin );
      delete [] newBegin;

      capacity_ = rb.size();
      cur_ = begin_;
   }
   else {
      end_ = std::copy( rb.cur_, rb.end_, begin_ );
      cur_ = begin_;
   }

   return *this;
}

template< typename T >  // Element type
template< typename G >
GenericRecvBuffer<T>& GenericRecvBuffer<T>::operator=( const GenericSendBuffer<T,G> & sb )
{
   if( sb.size() > capacity_ ) {
      T* newBegin( new T[sb.size()] );
      end_ = std::copy( sb.begin_, sb.cur_, newBegin );
      std::swap( begin_, newBegin );
      delete [] newBegin;

      capacity_ = sb.size();
      cur_ = begin_;
   }
   else {
      end_ = std::copy( sb.begin_, sb.cur_, begin_ );
      cur_ = begin_;
   }

   return *this;
}
//**********************************************************************************************************************




//======================================================================================================================
//
//  GET FUNCTIONS
//
//======================================================================================================================

//**********************************************************************************************************************
/*!\brief Returns the maximum possible size of the receive buffer.
//
// \return The maximum possible size.
*/
template< typename T >  // Element type
inline size_t GenericRecvBuffer<T>::maxSize() const
{
   return size_t(-1) / sizeof(T);
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\brief Returns the current size of the receive buffer.
//
// \return The current size.
*/
template< typename T >  // Element type
inline size_t GenericRecvBuffer<T>::size() const
{
   return numeric_cast< size_t >( end_ - cur_ );
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\brief Returns the capacity of the receive buffer.
//
// \return The capacity.
*/
template< typename T >  // Element type
inline size_t GenericRecvBuffer<T>::capacity() const
{
   return numeric_cast< size_t >( capacity_ );
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\brief Returns \a true if the receive buffer is empty.
//
// \return \a true if the receive buffer is empty, \a false if it is not.
*/
template< typename T >  // Element type
inline bool GenericRecvBuffer<T>::isEmpty() const
{
   return cur_ >= end_;
}
//**********************************************************************************************************************




//======================================================================================================================
//
//  OPERATORS
//
//======================================================================================================================
//**********************************************************************************************************************
/*!\brief Implements operator>> without debugging system.
//
// \return Reference to the receive buffer.
*/
template< typename T >  // Element type
template< typename V >  // Type of the built-in data value
typename std::enable_if< std::is_arithmetic<V>::value || std::is_enum<V>::value,
                           GenericRecvBuffer<T> & >::type
GenericRecvBuffer<T>::get( V& value )
{
   // Compile time check that V is built-in data type
   static_assert( std::is_arithmetic<V>::value || std::is_enum<V>::value,
                            "RecvBuffer accepts only built-in data types");


   static_assert( sizeof(V) >= sizeof(T), "V has to be bigger than T" );
   static_assert( sizeof(V) %  sizeof(T) == 0, "V has to be divisible by T" );



   // Checking the validity of the read operation
   WALBERLA_ASSERT_LESS_EQUAL( cur_ + (sizeof(V) / sizeof(T)), end_ );

   // Extracting the data value
   std::memcpy( &value, cur_, sizeof(V) );
   cur_ += sizeof(V) / sizeof(T);

   // Invariants check
   WALBERLA_ASSERT_LESS_EQUAL( cur_, end_);

   return *this;
}
//**********************************************************************************************************************



//**********************************************************************************************************************
/*!\brief Reads a built-in data value from the receive buffer.
//
// \return Reference to the receive buffer.
//
// This function extracts one data value of built-in data type from the receive buffer.
//
// \b Note: This operator may only be used for built-in data types. The attempt to use
// a user-defined data type results in a compile time error!
*/
template< typename T >  // Element type
template< typename V >  // Type of the built-in data value
typename std::enable_if< std::is_arithmetic<V>::value || std::is_enum<V>::value,
                           GenericRecvBuffer<T> & >::type
GenericRecvBuffer<T>::operator>>( V& value )
{
   readDebugMarker( typeid(V).name() );
   return get( value );

}
//**********************************************************************************************************************




//======================================================================================================================
//
//  UTILITY FUNCTIONS
//
//======================================================================================================================

//**********************************************************************************************************************
/*!\brief Returns a pointer to the first element of the receive buffer.
//
// \return Pointer to the first element of the receive buffer.
//
// This utility function enables the RecvBuffer to be used directly as receive buffer in MPI
// receive functions. Note however, that this operation is only allowed for a reinitialized
// receive buffer (for instance via the resize() function).

   \code
   using namespace pe;

   // Preparing a receive buffer for a message of 50 bytes
   RecvBuffer<byte> buffer;
   buffer.resize( 50 );

   // Receiving an MPI message from process 0 in a blocking MPI_Recv() function
   MPI_Status status;
   MPI_Recv( buffer.ptr(), buffer.size(), MPI_BYTE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, status );
   \endcode
*/
template< typename T >  // Element type
inline T* GenericRecvBuffer<T>::ptr() const
{
   WALBERLA_ASSERT_EQUAL( begin_, cur_);
   return begin_;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\brief Setting the minimum capacity of the receive buffer.
//
// \param newCapacity The new minimum capacity of the receive buffer.
// \return void
//
// This function reserves at least \a newCapacity elements of data type \a T for the receive
// buffer. Note that this operation involves a complete reset of the receive buffer!
*/
template< typename T >  // Element type
void GenericRecvBuffer<T>::reserve( size_t newCapacity )
{
   if( newCapacity > capacity_ )
   {
      // Allocating a new array
      T* tmp = new T[newCapacity];

      // Replacing the old array
      std::swap( tmp, begin_ );
      delete [] tmp;

      // Adjusting the members
      capacity_ = newCapacity;
   }

   // Clearing the receive buffer
   clear();
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\brief Changing the size of the receive buffer.
//
// \param newSize The new size of the receive buffer.
// \return void
//
// This function resizes the receive buffer to the given size \a newSize. Note that this
// operation does not preserve the current contents of the receive buffer!
*/
template< typename T >  // Element type
void GenericRecvBuffer<T>::resize( size_t newSize )
{
   if( newSize > capacity_ )
   {
      // Allocating a new array
      T* tmp = new T[newSize];

      // Replacing the old array
      std::swap( tmp, begin_ );
      delete [] tmp;

      // Adjusting the capacity
      capacity_ = newSize;
   }

   cur_ = begin_;
   end_ = begin_+newSize;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\brief Reads a built-in data value from the receive buffer without extracting it.
//
// \return void
//
// This function reads the next built-in data value of type \a V from the receive buffer
// without extracting/removing it from the buffer.
*/
template< typename T >  // Element type
template< typename V >  // Type of the built-in data value
inline void GenericRecvBuffer<T>::peek( V& value ) const
{
   WALBERLA_STATIC_ASSERT( std::is_arithmetic<V>::value );

   WALBERLA_STATIC_ASSERT( sizeof(V) > sizeof(T) );
   WALBERLA_STATIC_ASSERT( sizeof(V) % sizeof(T) == 0);

   // Checking the validity of the read operation
   WALBERLA_ASSERT_LESS_EQUAL( cur_ + BUFFER_DEBUG_OVERHEAD + sizeof(V)/sizeof(T), end_);

   // Extracting the data value
   V* tmp( reinterpret_cast<V*>( cur_ + BUFFER_DEBUG_OVERHEAD ) );
   value = *tmp;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\brief Skipping the given number of elements.
//
// \param elements The number of elements to be skipped.
// \return void
//
// This function skips \a element receive buffer elements of type \a T.
*/
template< typename T >  // Element type
T * GenericRecvBuffer<T>::skip( size_t elements )
{
   auto previous = cur_;
   cur_ += elements;

   // Invariants check
   WALBERLA_ASSERT_LESS_EQUAL( cur_, end_ );

   return previous;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\brief Clearing the receive buffer.
//
// \return void
//
// This function performs a complete reset of the receive buffer.
*/
template< typename T >  // Element type
inline void GenericRecvBuffer<T>::clear()
{
   cur_ = begin_;
   end_ = begin_;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\brief Clearing the receive buffer.
//
// \return void
//
// This function performs a complete reset of the receive buffer - including the deletion of allocated memory!
*/
template< typename T >  // Element type
inline void GenericRecvBuffer<T>::reset()
{
   delete [] begin_;
   capacity_ = 0;
   begin_    = nullptr;
   cur_      = nullptr;
   end_      = nullptr;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\brief Reads debug marker and raises assert when marker does not match
//
//
*/
#ifdef WALBERLA_BUFFER_DEBUG

template< typename T >  // Element type
inline void GenericRecvBuffer<T>::readDebugMarker( const char * expectedMarker )
{
   bool valid = true;
   uint_t len = std::strlen( expectedMarker );
   // read the first BUFFER_DEBUG_OVERHEAD chars of the marker
   for( uint_t i = 0; i < len && i < BUFFER_DEBUG_OVERHEAD; ++i )
   {
      T readMarker;
      get( readMarker );
      if ( readMarker != T( expectedMarker[i] ))
         valid = false;
   }
   // If marker was shorter then BUFFER_DEBUG_OVERHEAD the rest should be '-'
   // to always have the same length
   for( uint_t i = len; i < BUFFER_DEBUG_OVERHEAD; ++i ) {
      T readMarker;
      get( readMarker );
      if ( readMarker != T('-') )
         valid = false;
   }

   WALBERLA_ASSERT( valid , "Packed and unpacked (" << expectedMarker << ") datatype do not match." );

}
#else

template< typename T >  // Element type
inline void GenericRecvBuffer<T>::readDebugMarker( const char *  ){}

#endif

//**********************************************************************************************************************




} // namespace mpi
} // namespace walberla
