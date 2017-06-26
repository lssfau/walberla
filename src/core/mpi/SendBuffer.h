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
//! \file SendBuffer.h
//! \ingroup core
//! \author Klaus Iglberger
//! \brief Implementation of a MPI send buffer.
//!
//! Copyright (C) 2009 Klaus Iglberger
//! Taken from "pe Physics Engine" with small changes
//
//======================================================================================================================

#pragma once

#include "waLBerlaDefinitions.h"
#include "BufferSizeTrait.h"
#include "growPolicies/ConstantGrowth.h"
#include "growPolicies/LinearGrowth.h"
#include "growPolicies/OptimalGrowth.h"

#include "core/debug/Debug.h"

#include <boost/mpl/logical.hpp>
#include <boost/type_traits/is_enum.hpp>
#include <boost/type_traits/is_arithmetic.hpp>
#include <boost/utility/enable_if.hpp>

#include <algorithm>
#include <typeinfo>


namespace walberla {
namespace mpi {

//======================================================================================================================
//
//  CLASS DEFINITION
//
//======================================================================================================================

//**********************************************************************************************************************
/*!\brief Implementation of a MPI send buffer.
// \ingroup mpi
//
// The SendBuffer class is a special purpose implementation for the MPI communication
// functionality. It offers a convenient and safe setup of MPI messages even for mixed-type
// communication. The following example gives an impression of the usage of the SendBuffer
// class:

   \code
   using namespace mpi;

   // Adding a double and an integer to a send buffer
   SendBuffer buffer;
   double d( 2.1 );
   Vec3 v( 3.5, -2.1, 0.7 );
   buffer << d << v;

   // Sending a MPI message to process 0 with the contents of the send buffer
   MPI_Send( buffer.ptr(), buffer.size(), MPI_BYTE, 0, 0, MPI_COMM_WORLD );
   \endcode
// For another example see also the Buffer Unit Test in File BufferTest.cpp
//
// Note that the order of data values in the send buffer is depending on the order the
// elements are added to the buffer. This order is not changed during the MPI communication
// and the same order has to be used during the extraction of the sent data elements. See
// also the RecvBuffer class description for the receiver side of the MPI communication.
*/
template< typename T = unsigned char    // Element type
        , typename G = OptimalGrowth >  // Growth policy
class GenericSendBuffer
{
public:
   //**Type definitions*************************************************************************************************
   typedef T  ElementType;  //!< Type of the receive buffer elements.
   //*******************************************************************************************************************

   //**Constructors*****************************************************************************************************
   /*!\name Constructors */
   //@{
   explicit inline GenericSendBuffer( size_t initCapacity = 0 );
            inline GenericSendBuffer( const GenericSendBuffer& sb );
   //@}
   //*******************************************************************************************************************

   //**Destructor*******************************************************************************************************
   /*!\name Destructor */
   //@{
   inline ~GenericSendBuffer();
   //@}
   //*******************************************************************************************************************

   //**Assignment operator**********************************************************************************************
   /*!\name Assignment operator */
   //@{
   GenericSendBuffer& operator=( const GenericSendBuffer& sb );
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
   typename boost::enable_if< boost::mpl::or_< boost::is_arithmetic<V>, boost::is_enum<V> >,
                              GenericSendBuffer&  >::type
   operator<<( V value );

   //@}
   //*******************************************************************************************************************

   //**Repositioning ***************************************************************************************************
   /*!\name Repositioning */
   inline void rewind(const size_t & size);
   //@}
   //*******************************************************************************************************************

   //**Utility functions************************************************************************************************
   /*!\name Utility functions */
   //@{
   inline T*             ptr    () const;
   inline std::ptrdiff_t getOffset() const;
   inline T*             getMemoryLocation( const std::ptrdiff_t offset);
   inline void           reserve( size_t newCapacity );
   inline void           clear  ();
   inline void           reset  ();
   inline void           addDebugMarker( const char * marker );
   //@}
   //*******************************************************************************************************************

private:
   //**Utility functions************************************************************************************************
   /*!\name Utility functions */
   //@{
   void extendMemory( size_t newCapacity );

   template< typename V >
   typename boost::enable_if< boost::mpl::or_< boost::is_arithmetic<V>, boost::is_enum<V> >,
                              GenericSendBuffer&  >::type
   put( V value );
   //@}
   //*******************************************************************************************************************

   //**Member variables*************************************************************************************************
   /*!\name Member variables */
   //@{
   T* begin_;  //!< Pointer to the first element of the send buffer.
   T* cur_;    //!< Pointer to the current/last element of the send buffer.
   T* end_;    //!< Pointer to the end of the storage of the send buffer.
   //@}
   //*******************************************************************************************************************

   //**Compile time checks**********************************************************************************************
   static_assert( boost::is_arithmetic<T>::value, "SendBuffer<T>: T has to be native datatype" ) ;
   //*******************************************************************************************************************

   template< typename U >
   friend class GenericRecvBuffer;
};
//**********************************************************************************************************************

typedef GenericSendBuffer<> SendBuffer;



//======================================================================================================================
//
//  CONSTRUCTORS
//
//======================================================================================================================

//**********************************************************************************************************************
/*!\brief Standard constructor for SendBuffer.
//
// \param initCapacity The initial capacity of the send buffer.
//
// The default initial capacity of the send buffer is specified by the selected growth policy.
*/
template< typename T    // Element type
        , typename G >  // Growth policy
inline GenericSendBuffer<T,G>::GenericSendBuffer( size_t initCapacity )
   : begin_( new T[initCapacity] )  // Pointer to the first element
   , cur_  ( begin_ )               // Pointer to the current/last element
   , end_  ( begin_ )               // Pointer to the end of the storage
{}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\brief Copy constructor for SendBuffer.
//
// \param sb The send buffer to be copied.
*/
template< typename T    // Element type
        , typename G >  // Growth policy
inline GenericSendBuffer<T,G>::GenericSendBuffer( const GenericSendBuffer& sb )
   : begin_( new T[sb.size()] )  // Pointer to the first element
   , cur_  ( begin_+sb.size() )  // Pointer to the current/last element
   , end_  ( cur_ )              // Pointer to the end of the storage
{
   for( size_t i=0; i<sb.size(); ++i )
      begin_[i] = sb.begin_[i];
}
//**********************************************************************************************************************




//======================================================================================================================
//
//  DESTRUCTOR
//
//======================================================================================================================

//**********************************************************************************************************************
/*!\brief Destructor for SendBuffer.
*/
template< typename T    // Element type
        , typename G >  // Growth policy
inline GenericSendBuffer<T,G>::~GenericSendBuffer()
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
/*!\brief Copy assignment operator for SendBuffer.
//
// \param sb The send buffer to be copied.
// \return Reference to the assigned send buffer.
*/
template< typename T    // Element type
        , typename G >  // Growth policy
GenericSendBuffer<T,G>& GenericSendBuffer<T,G>::operator=( const GenericSendBuffer& sb )
{
   if( &sb == this ) return *this;

   if( sb.size() > capacity() ) {
      T* newBegin( new T[sb.size()] );
      end_ = std::copy( sb.begin_, sb.cur_, newBegin );
      std::swap( begin_, newBegin );
      delete [] newBegin;

      cur_ = end_;
   }
   else {
      cur_ = std::copy( sb.begin_, sb.end_, begin_ );
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
/*!\brief Returns the maximum possible size of the send buffer.
//
// \return The maximum possible size.
*/
template< typename T    // Element type
        , typename G >  // Growth policy
inline size_t GenericSendBuffer<T,G>::maxSize() const
{
   return size_t(-1) / sizeof(T);
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\brief Returns the current size of the send buffer.
//
// \return The current size.
*/
template< typename T    // Element type
        , typename G >  // Growth policy
inline size_t GenericSendBuffer<T,G>::size() const
{
   return numeric_cast< size_t >( cur_ - begin_ );
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\brief Returns the capacity of the send buffer.
//
// \return The capacity.
*/
template< typename T    // Element type
        , typename G >  // Growth policy
inline size_t GenericSendBuffer<T,G>::capacity() const
{
   return numeric_cast< size_t >( end_ - begin_ );
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\brief Returns \a true if the send buffer is empty.
//
// \return \a true if the send buffer is empty, \a false if it is not.
*/
template< typename T    // Element type
        , typename G >  // Growth policy
inline bool GenericSendBuffer<T,G>::isEmpty() const
{
   return begin_ == cur_;
}
//**********************************************************************************************************************




//======================================================================================================================
//
//  OPERATORS
//
//======================================================================================================================

//**********************************************************************************************************************
/*!\brief Implements operator<< without debugging system
//
// \return Reference to the send buffer.
//
*/
template< typename T    // Element type
        , typename G >  // Growth policy
template< typename V >  // Type of the built-in data value
typename boost::enable_if< boost::mpl::or_< boost::is_arithmetic<V>, boost::is_enum<V> >,
                           GenericSendBuffer<T,G>& >::type
GenericSendBuffer<T,G>::put( V value )
{
   // Compile time check that V is built-in data type
   static_assert( boost::is_arithmetic<V>::value || boost::is_enum<V>::value,
                            "SendBuffer accepts only built-in data types");

   static_assert( sizeof(V) >= sizeof(T), "Type that is stored has to be bigger than T" );
   static_assert( sizeof(V)  % sizeof(T) == 0, "V has to be divisible by T ");

   size_t count =  sizeof(V) / sizeof(T) ;
   const size_t rest = numeric_cast< size_t >( end_ - cur_ );

   // Checking the size of the remaining memory
   if( rest < count ) {
      extendMemory( size() + count );
   }

   // Adding the data value
   V* const tmp( reinterpret_cast<V*>( cur_ ) );
   *tmp  = value;
   cur_ += count;

   // Invariants check
   WALBERLA_ASSERT_LESS_EQUAL( cur_, end_ );

   return *this;
}
//**********************************************************************************************************************



//**********************************************************************************************************************
/*!\brief Adds a built-in data value to the send buffer.
//
// \return Reference to the send buffer.
//
// This function adds one data value of built-in data type to the send buffer.
//
// \b Note: This operator may only be used for built-in data types. The attempt to use
// a user-defined data type results in a compile time error!
*/
template< typename T    // Element type
        , typename G >  // Growth policy
template< typename V >  // Type of the built-in data value
typename boost::enable_if< boost::mpl::or_< boost::is_arithmetic<V>, boost::is_enum<V> >,
                           GenericSendBuffer<T,G>& >::type
GenericSendBuffer<T,G>::operator<<( V value )
{
   addDebugMarker( typeid(V).name() );
   return put( value );
}
//**********************************************************************************************************************







//======================================================================================================================
//
//  REPOSITIONING
//
//======================================================================================================================

//**********************************************************************************************************************
/*!\brief Rewinds the stream to a previous position
//
   \code
   // Filling a send buffer with a double and an integer value
   SendBuffer<byte> buffer;
   double d =  2.1;
   buffer << d;
   size_t memorizedSize = buffer.size();
   int i=3;
   buffer << i;
   buffer.rewind(memorizedSize);
   i=5;
   buffer << i; //overwrites the 3 in the buffer with the new value 5
   \endcode
*/
template< typename T    // Element type
        , typename G >  // Growth policy
inline void GenericSendBuffer<T,G>::rewind(const size_t & s)
{
   WALBERLA_ASSERT_LESS(s, size());
   cur_ = begin_ + s;
   WALBERLA_ASSERT_EQUAL(size(), s);
}


//======================================================================================================================
//
//  UTILITY FUNCTIONS
//
//======================================================================================================================

//**********************************************************************************************************************
/*!\brief Returns a pointer to the first element of the send buffer.
//
// \return Pointer to the first element of the send buffer.
//
// This utility function enables the SendBuffer to be used directly as send buffer in MPI
// send functions:

   \code
   using namespace communication;

   // Filling a send buffer with a double and an integer value
   SendBuffer<byte> buffer;
   double d =  2.1;
   int i= 2;
   buffer << d << i;

   // Sending a MPI message to process 0 with the contents of the send buffer
   MPI_Send( buffer.ptr(), buffer.size(), MPI_BYTE, 0, 0, MPI_COMM_WORLD );
   \endcode
*/
template< typename T    // Element type
        , typename G >  // Growth policy
inline T* GenericSendBuffer<T,G>::ptr() const
{
   return begin_;
}
//**********************************************************************************************************************

/**
 * Returns the offset from the beginning to the current position inside the buffer in bytes.
 *
 * Example:
 * \snippet BufferTest.cpp SendBuffer Overwrite Test
 * The buffer now contains 3, 2, 3
 * \attention This is a low level function. Use with care!
 * \see getMemoryLocation()
 */
template< typename T    // Element type
        , typename G >  // Growth policy
inline std::ptrdiff_t GenericSendBuffer<T,G>::getOffset() const
{
   return cur_ - begin_;
}


/**
 * Returns the memory address corresponding to the offset. Offset is measured in bytes from the beginning of the buffer.
 *
 * Example:
 * \snippet BufferTest.cpp SendBuffer Overwrite Test
 * The buffer now contains 3, 2, 3
 * \attention This is a low level function. Use with care!
 * \see getOffset()
 */
template< typename T    // Element type
        , typename G >  // Growth policy
inline T*   GenericSendBuffer<T,G>::getMemoryLocation( const std::ptrdiff_t offset)
{
   return begin_ + offset;
}


//**********************************************************************************************************************
/*!\brief Setting the minimum capacity of the send buffer.
//
// \param newCapacity The new minimum capacity of the send buffer.
// \return void
//
// This function reserves at least \a newCapacity elements of data type \a T for the send
// buffer.
*/
template< typename T    // Element type
        , typename G >  // Growth policy
inline void GenericSendBuffer<T,G>::reserve( size_t newCapacity )
{
   if( newCapacity > capacity() )
      extendMemory( newCapacity );
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\brief Extending the internal memory of the send buffer.
//
// \param newCapacity The new minimum capacity of the send buffer.
// \return void
*/
template< typename T    // Element type
        , typename G >  // Growth policy
void GenericSendBuffer<T,G>::extendMemory( size_t newCapacity )
{
   // Calculating the new capacity
   WALBERLA_ASSERT_GREATER( newCapacity, capacity() );
   newCapacity = G()( capacity(), newCapacity );
   WALBERLA_ASSERT_GREATER( newCapacity, capacity() );

   // Allocating a new array
   T* tmp = new T[newCapacity];

   // Replacing the old array
   cur_ = std::copy( begin_, cur_, tmp );
   std::swap( tmp, begin_ );
   end_ = begin_ + newCapacity;
   delete [] tmp;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\brief Clearing the send buffer.
//
// \return void
//
// This function performs a complete reset of the send buffer.
*/
template< typename T    // Element type
        , typename G >  // Growth policy
inline void GenericSendBuffer<T,G>::clear()
{
   cur_ = begin_;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\brief Clearing the send buffer.
//
// \return void
//
// This function performs a complete reset of the send buffer - including the deletion of allocated memory!
*/
template< typename T    // Element type
        , typename G >  // Growth policy
inline void GenericSendBuffer<T,G>::reset()
{
   delete [] begin_;
   begin_ = new T[0];
   cur_ = begin_;
   end_ = begin_;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\brief Adds a debug marker if buffer debugging is enabled
//
// See corresponding function GenericRecvBuffer::readDebugMarker()
*/
#ifdef WALBERLA_BUFFER_DEBUG

template< typename T    // Element type
        , typename G >  // Growth policy
inline void GenericSendBuffer<T,G>::addDebugMarker( const char * marker )
{
   uint_t len = std::strlen( marker );
   // Push the first BUFFER_DEBUG_OVERHEAD chars of the marker
   for( uint_t i = 0; i < len && i < BUFFER_DEBUG_OVERHEAD; ++i )
      put( T( marker[i] ) );
   // If marker was shorter then BUFFER_DEBUG_OVERHEAD fill the rest with '-'
   // to always have the same length
   for( uint_t i = len; i < BUFFER_DEBUG_OVERHEAD; ++i )
      put( T( '-' ) );
}
#else
template< typename T    // Element type
        , typename G >  // Growth policy
inline void GenericSendBuffer<T,G>::addDebugMarker( const char *  )
{}
#endif

//**********************************************************************************************************************



} // namespace mpi
} // namespace walberla
