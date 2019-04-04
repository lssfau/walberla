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
//! \file BufferDataTypeExtensions.h
//! \ingroup core
//! \author Martin Bauer <martin.bauer@fau.de>
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//! \brief Streaming Operators, to store/extract compound data types into Buffers
//
//======================================================================================================================

#pragma once

#include "RecvBuffer.h"
#include "SendBuffer.h"
#include "core/Conversion.h"
#include "core/DataTypes.h"
#include "core/Optional.h"
#include "core/RandomUUID.h"

#include <array>
#include <deque>
#include <limits>
#include <list>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>


namespace walberla {
namespace mpi {

// ---------------------------------------------------------------------------------------------------------------------
// -----------------------------   Standard Pair Support   -------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------

template< typename T,    // Element type of SendBuffer
          typename G,    // Growth policy of SendBuffer
          typename T1,   // first type
          typename T2>   // second type
GenericSendBuffer<T,G>& operator<<( GenericSendBuffer<T,G> & buf, const std::pair<T1, T2> & pair )
{
   buf.addDebugMarker( "pa" );
   buf << pair.first << pair.second;
   return buf;
}

template< typename T,    // Element type of RecvBuffer
          typename T1,   // first type
          typename T2>   // second type
GenericRecvBuffer<T>& operator>>( GenericRecvBuffer<T> & buf, std::pair<T1, T2> & pair )
{
   buf.readDebugMarker( "pa" );
   buf >> pair.first >> pair.second;
   return buf;
}

template<typename T1, typename T2>
struct BufferSizeTrait< std::pair<T1,T2> > {
   static const bool constantSize = BufferSizeTrait<T1>::constantSize &&
                                    BufferSizeTrait<T2>::constantSize;
   static const uint_t size = BufferSizeTrait<T1>::size + BufferSizeTrait<T2>::size + BUFFER_DEBUG_OVERHEAD;
};

// ---------------------------------------------------------------------------------------------------------------------
// ----------------------------- Standard Container Support  -----------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------

template< typename T,     // Element type of SendBuffer
          typename G,     // Growth policy of SendBuffer
          typename Cont > // Container
void sendNonResizableContainer( GenericSendBuffer<T,G> & buf, const Cont & container )
{
  buf.addDebugMarker( "ar" );
  for( const auto & it : container )
    buf << it;
}

template< typename T,     // Element type of RecvBuffer
          typename Cont > // Container
void recvNonResizableContainer( GenericRecvBuffer<T> & buf, Cont & container )
{
  buf.readDebugMarker( "ar" );
  for( auto & it : container )
    buf >> it;
}

template< typename T,     // Element type of SendBuffer
          typename G,     // Growth policy of SendBuffer
          typename CT,    // Element type
          std::size_t N > // Array size
GenericSendBuffer<T,G>& operator<<( GenericSendBuffer<T,G> & buf, const std::array<CT, N> & array )
{
  sendNonResizableContainer(buf, array);
  return buf;
}

template< typename T,     // Element type of RecvBuffer
          typename CT,    // Element type
          std::size_t N > // Array size
GenericRecvBuffer<T>& operator>>( GenericRecvBuffer<T> & buf, std::array<CT, N> & array )
{
  recvNonResizableContainer(buf, array);
  return buf;
}

template<typename T, std::size_t N>
struct BufferSizeTrait< std::array< T, N > > {
    static const bool constantSize = true;
    static const uint_t size = N * sizeof(T) + BUFFER_DEBUG_OVERHEAD;
};


template< typename T,    // Element type of SendBuffer
          typename G,    // Growth policy of SendBuffer
          typename Cont> // Container
void sendContainer( GenericSendBuffer<T,G> & buf, const Cont & container )
{
   buf.addDebugMarker( "co" );
   size_t size( container.size() );
   buf << size;
   for( typename Cont::const_iterator it = container.begin(); it != container.end(); ++it )
      buf << *it;
}


#ifdef WALBERLA_CXX_COMPILER_IS_GNU
#if __GNUC__ == 4 && __GNUC_MINOR__ == 9 || __GNUC__ == 6
#pragma GCC push_options
#pragma GCC optimize(2)
#endif
#endif

template< typename T,    // Element type of RecvBuffer
          typename Cont> // Container
void recvContainer( GenericRecvBuffer<T> & buf, Cont & container )
{
   buf.readDebugMarker( "co" );
   size_t size;
   buf >> size;
   container.resize(size);
   for( typename Cont::iterator it = container.begin(); it != container.end(); ++it )
      buf >> *it;
}

#ifdef WALBERLA_CXX_COMPILER_IS_GNU
#if __GNUC__ == 4 && __GNUC_MINOR__ == 9 || __GNUC__ == 6
#pragma GCC pop_options
#endif
#endif



template< typename T,    // Element type of SendBuffer
          typename G,    // Growth policy of SendBuffer
          typename CT,   // Element type
          typename CA>   // Allocator type
GenericSendBuffer<T,G>& operator<<( GenericSendBuffer<T,G> & buf, const std::vector<CT, CA> & c )
{
   sendContainer(buf, c);
   return buf;
}

template< typename T,    // Element type of RecvBuffer
          typename CT,   // Element type
          typename CA>   // Allocator type
GenericRecvBuffer<T>& operator>>( GenericRecvBuffer<T> & buf, std::vector<CT, CA> & c )
{
   recvContainer(buf, c);
   return buf;
}

template<typename T, typename A>
struct BufferSizeTrait< std::vector<T,A> > { static const bool constantSize = false;  };





template< typename T,    // Element type of SendBuffer
          typename G,    // Growth policy of SendBuffer
          typename CT >  // Character type
GenericSendBuffer<T,G>& operator<<( GenericSendBuffer<T,G> & buf, const std::basic_string<CT> & c )
{
   sendContainer(buf, c);
   return buf;
}

template< typename T,    // Element type of RecvBuffer
          typename CT >  // Character type
GenericRecvBuffer<T>& operator>>( GenericRecvBuffer<T> & buf, std::basic_string<CT> & c )
{
   recvContainer(buf, c);
   return buf;
}

template<typename T>
struct BufferSizeTrait< std::basic_string<T> > { static const bool constantSize = false;  };






template< typename T,    // Element type of SendBuffer
          typename G,    // Growth policy of SendBuffer
          typename CT,   // Element type
          typename CA>   // Allocator type
GenericSendBuffer<T,G>& operator<<( GenericSendBuffer<T,G> & buf, const std::deque<CT, CA> & c )
{
   sendContainer(buf, c);
   return buf;
}

template< typename T,    // Element type of RecvBuffer
          typename CT,   // Element type
          typename CA>   // Allocator type
GenericRecvBuffer<T>& operator>>( GenericRecvBuffer<T> & buf, std::deque<CT, CA> & c )
{
   recvContainer(buf, c);
   return buf;
}

template<typename T, typename A>
struct BufferSizeTrait< std::deque<T,A> > { static const bool constantSize = false;  };




template< typename T,    // Element type of SendBuffer
          typename G,    // Growth policy of SendBuffer
          typename CT,   // Element type
          typename CA>   // Allocator type
GenericSendBuffer<T,G>& operator<<( GenericSendBuffer<T,G> & buf, const std::list<CT, CA> & c )
{
   sendContainer(buf, c);
   return buf;
}

template< typename T,    // Element type of RecvBuffer
          typename CT,   // Element type
          typename CA>   // Allocator type
GenericRecvBuffer<T>& operator>>( GenericRecvBuffer<T> & buf, std::list<CT, CA> & c )
{
   recvContainer(buf, c);
   return buf;
}

template<typename T, typename A>
struct BufferSizeTrait< std::list<T,A> > { static const bool constantSize = false;  };



// ---------------------------------------------------------------------------------------------------------------------
// ----------------------------- std::vector<bool> Support -------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------

template< typename T,    // Element type of SendBuffer
          typename G >   // Growth policy of SendBuffer
GenericSendBuffer<T,G>& packBoolVectorWithoutSize(GenericSendBuffer<T,G> & buf, const std::vector<bool> & bools )
{
   // Use an unsigned type at least as large as the SendBuffer base type as container for the bools
   typedef typename leastUnsignedInteger< std::numeric_limits<T>::digits >::type ContainerType;
   static const size_t NUM_BITS = std::numeric_limits<ContainerType>::digits;

   auto it = bools.begin();
   while( it != bools.end() )
   {
      ContainerType container = 0;

      for( size_t bitPos = 0; bitPos < NUM_BITS && it != bools.end(); ++bitPos, ++it )
         container = numeric_cast<ContainerType>( container | numeric_cast<ContainerType>( ( static_cast<bool>( *it ) ? ( 1U << bitPos ) : 0U ) ) );

      buf << container;
   }

   return buf;
}

template< typename T >    // Element type  of RecvBuffer
GenericRecvBuffer<T>& unpackBoolVectorWithoutSize(GenericRecvBuffer<T> & buf, std::vector<bool> & bools, size_t size )
{
   // Use an unsigned type at least as large as the RecvBuffer base type as container for the bools
   typedef typename leastUnsignedInteger<std::numeric_limits<T>::digits>::type ContainerType;
   static const size_t NUM_BITS = std::numeric_limits<ContainerType>::digits;

   bools.resize(size);

   auto it = bools.begin();
   while( it != bools.end() )
   {
      ContainerType container;
      buf >> container;

      for( size_t bitPos = 0; bitPos < NUM_BITS && it != bools.end(); ++bitPos, ++it )
         *it = ( container & ( static_cast<ContainerType>(1) << bitPos ) ) > 0U;
   }
   return buf;
}

template< typename T,    // Element type of SendBuffer
          typename G >   // Growth policy of SendBuffer
GenericSendBuffer<T,G>& operator<<(GenericSendBuffer<T,G> & buf, const std::vector<bool> & bools )
{
   buf.addDebugMarker( "bv" );
   buf << bools.size();
   packBoolVectorWithoutSize( buf, bools );

   return buf;
}

template< typename T >    // Element type  of RecvBuffer
GenericRecvBuffer<T>& operator>>(GenericRecvBuffer<T> & buf, std::vector<bool> & bools )
{
   buf.readDebugMarker( "bv" );
   size_t numBools;
   buf >> numBools;

   unpackBoolVectorWithoutSize( buf,bools, numBools );

   return buf;
}

// ---------------------------------------------------------------------------------------------------------------------
// ----------------------------- Standard Associative Container Support ------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------

template< typename T,    // Element type of SendBuffer
          typename G,    // Growth policy of SendBuffer
          typename Cont> // Container
void sendAssocContainer( GenericSendBuffer<T,G> & buf, const Cont & container )
{
   buf << container.size();
   for( typename Cont::const_iterator it = container.begin(); it != container.end(); ++it )
      buf << *it;
}


template< typename T,    // Element type of RecvBuffer
          typename Cont> // Container
void recvAssocContainer( GenericRecvBuffer<T> & buf, Cont & container )
{
   container.clear();

   size_t size;
   buf >> size;

   for( size_t i = 0; i < size; ++i)
   {
      typename Cont::value_type value;
      buf >> value;
      container.insert(container.end(), value);
   }
}

template< typename T,    // Element type of RecvBuffer
          typename Cont> // Container
void recvMap( GenericRecvBuffer<T> & buf, Cont & container )
{
   container.clear();

   size_t size;
   buf >> size;

   for( size_t i = 0; i < size; ++i)
   {
      std::pair<typename Cont::key_type, typename Cont::mapped_type> value;
      buf >> value;
      container.insert(container.end(), value);
   }
}


template< typename T,    // Element type of SendBuffer
          typename G,    // Growth policy of SendBuffer
          typename CK,   // Key type
          typename CC,   // Comparison type
          typename CA>   // Allocator type
GenericSendBuffer<T,G>& operator<<( GenericSendBuffer<T,G> & buf, const std::set<CK, CC, CA> & c )
{
   buf.addDebugMarker( "se" );
   sendAssocContainer(buf, c);
   return buf;
}

template< typename T,    // Element type of RecvBuffer
          typename CK,   // Key type
          typename CC,   // Comparison type
          typename CA>   // Allocator type
GenericRecvBuffer<T>& operator>>( GenericRecvBuffer<T> & buf, std::set<CK, CC, CA> & c )
{
   buf.readDebugMarker( "se" );
   recvAssocContainer(buf, c);
   return buf;
}

template<typename T, typename C, typename A>
struct BufferSizeTrait< std::set<T,C,A> > { static const bool constantSize = false;  };





template< typename T,    // Element type of SendBuffer
          typename G,    // Growth policy of SendBuffer
          typename CK,   // Key type
          typename CC,   // Comparison type
          typename CA>   // Allocator type
GenericSendBuffer<T,G>& operator<<( GenericSendBuffer<T,G> & buf, const std::multiset<CK, CC, CA> & c )
{
   buf.addDebugMarker( "ms" );
   sendAssocContainer(buf, c);
   return buf;
}

template< typename T,    // Element type of RecvBuffer
          typename CK,   // Key type
          typename CC,   // Comparison type
          typename CA>   // Allocator type
GenericRecvBuffer<T>& operator>>( GenericRecvBuffer<T> & buf, std::multiset<CK, CC, CA> & c )
{
   buf.readDebugMarker( "ms" );
   recvAssocContainer(buf, c);
   return buf;
}


template<typename T, typename C, typename A>
struct BufferSizeTrait< std::multiset<T,C,A> > { static const bool constantSize = false;  };




template< typename T,    // Element type of SendBuffer
          typename G,    // Growth policy of SendBuffer
          typename CK,   // Key type
          typename CT,   // Element type
          typename CC,   // Comparison type
          typename CA>   // Allocator type
GenericSendBuffer<T,G>& operator<<( GenericSendBuffer<T,G> & buf, const std::map<CK, CT, CC, CA> & c )
{
   buf.addDebugMarker( "ma" );
   sendAssocContainer(buf, c);
   return buf;
}

template< typename T,    // Element type of RecvBuffer
          typename CK,   // Key type
          typename CT,   // Element type
          typename CC,   // Comparison type
          typename CA>   // Allocator type
GenericRecvBuffer<T>& operator>>( GenericRecvBuffer<T> & buf, std::map<CK, CT, CC, CA> & c )
{
   buf.readDebugMarker( "ma" );
   recvMap(buf, c);
   return buf;
}

template<typename T, typename K, typename C, typename A>
struct BufferSizeTrait< std::map<K,T,C,A> > { static const bool constantSize = false;  };





template< typename T,    // Element type of SendBuffer
          typename G,    // Growth policy of SendBuffer
          typename CK,   // Key type
          typename CT,   // Element type
          typename CC,   // Comparison type
          typename CA>   // Allocator type
GenericSendBuffer<T,G>& operator<<( GenericSendBuffer<T,G> & buf, const std::multimap<CK, CT, CC, CA> & c )
{
   buf.addDebugMarker( "mm" );
   sendAssocContainer(buf, c);
   return buf;
}

template< typename T,    // Element type of RecvBuffer
          typename CK,   // Key type
          typename CT,   // Element type
          typename CC,   // Comparison type
          typename CA>   // Allocator type
GenericRecvBuffer<T>& operator>>( GenericRecvBuffer<T> & buf, std::multimap<CK, CT, CC, CA> & c )
{
   buf.readDebugMarker( "mm" );
   recvMap(buf, c);
   return buf;
}


template<typename T, typename K, typename C, typename A>
struct BufferSizeTrait< std::multimap<K,T,C,A> > { static const bool constantSize = false;  };


// ---------------------------------------------------------------------------------------------------------------------
// ------------------------------------------- optional Support --------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------

template< typename T,    // Element type of SendBuffer
          typename G,    // Growth policy of SendBuffer
          typename OT>   // Optional type
GenericSendBuffer<T,G>& operator<<( GenericSendBuffer<T,G> & buf, const walberla::optional<OT> & o )
{
   buf.addDebugMarker( "op" );

   bool hasContent = true;
   if (o == walberla::nullopt)
      hasContent = false;

   buf << hasContent;

   if( hasContent )
      buf << *o;

   return buf;
}

template< typename T,    // Element type of RecvBuffer
          typename OT>   // Optional type
GenericRecvBuffer<T>& operator>>( GenericRecvBuffer<T> & buf, walberla::optional<OT> & o )
{
   buf.readDebugMarker( "op" );

   bool hasContent;

   buf >> hasContent;

   if( hasContent )
   {
      if( o )
         buf >> *o;
      else
      {
         o = OT();
         buf >> *o;
      }
   }
   else
   {
      o = walberla::nullopt;
   }

   return buf;
}

// ---------------------------------------------------------------------------------------------------------------------
// --------------------------------------- RandomUUID Support ----------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------

template<>
struct BufferSizeTrait< RandomUUID > {
   static const bool constantSize = true;
   static const uint_t size = 16 + BUFFER_DEBUG_OVERHEAD;
};

inline SendBuffer & operator<<( SendBuffer & buf, const RandomUUID& uuid )
{
   buf.addDebugMarker( "uu" );
   buf << uuid.getFirstUInt();
   buf << uuid.getSecondUInt();

   return buf;
}


inline RecvBuffer & operator>>( RecvBuffer & buf, RandomUUID& uuid )
{
   buf.readDebugMarker( "uu" );
   RandomUUID::UIntType a;
   RandomUUID::UIntType b;
   buf >> a >> b;
   uuid = RandomUUID(a, b);

   return buf;
}


} //namespace mpi
} //namespace walberla
