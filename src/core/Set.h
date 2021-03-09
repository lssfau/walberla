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
//! \file Set.h
//! \ingroup core
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/mpi/BufferDataTypeExtensions.h"

#include <algorithm>
#include <iostream>
#include <iterator>
#include <set>
#include <sstream>


namespace walberla {



template< typename T > class Set;

template< typename T > inline Set<T> setIntersection( const Set<T>& a, const Set<T>& b );
template< typename T > inline Set<T> setUnion       ( const Set<T>& a, const Set<T>& b );
template< typename T > inline Set<T> setDifference  ( const Set<T>& a, const Set<T>& b );

template< typename T > inline bool setIsEqual( const Set<T>& a, const Set<T>& b );

/// \cond internal
namespace set {
template< class InputIterator1, class InputIterator2 >
bool setsIntersect( InputIterator1 first1, InputIterator1 last1, InputIterator2 first2, InputIterator2 last2 );
} // namespace set
/// \endcond


//**********************************************************************************************************************
/*!
*   \brief Wrapper class for managing sets that store elements of type T
*
*   If the equality of two sets must be tested, the operators "==" and "!=" and the member function "isEqual" can be
*   used. If two sets must be compared in terms of size, the operators "<", ">", "<=", and ">=" and the member function
*   "equalSize" can be used:
*
    \code
      Set<int> A = Set<int>(1) + Set<int>(2) + Set<int>(3);
      Set<int> B = Set<int>(1) + Set<int>(2);
      Set<int> C = Set<int>(1) + Set<int>(2) + Set<int>(4);

      bool boolean = ( A == B );     // -> false
           boolean = ( A != B );     // -> true
           boolean = ( A >  B );     // -> true
           boolean = ( A >= B );     // -> true
           boolean = ( A <  B );     // -> false
           boolean = A.equalSize(B); // -> false
           boolean = ( A == C );     // -> false
           boolean = A.equalSize(C); // -> true
      B.insert(3);
           boolean = ( A == B );     // -> true
           boolean = A.equalSize(B); // -> true
           boolean = ( A >  B );     // -> false
           boolean = ( A >= B );     // -> true
           boolean = ( A <= B );     // -> true
    \endcode
*/
//**********************************************************************************************************************

template< typename T >
class Set {

public:
   using value_type = typename std::set<T>::value_type;
   using const_iterator = typename std::set<T>::const_iterator;
   using iterator = typename std::set<T>::iterator;

   friend inline Set<T> operator&( const Set& a, const Set& b ) { return setIntersection(a,b); } ///< intersection
   friend inline Set<T> operator+( const Set& a, const Set& b ) { return setUnion       (a,b); } ///< union
   friend inline Set<T> operator-( const Set& a, const Set& b ) { return setDifference  (a,b); } ///< difference / relative complement

   friend inline bool operator==( const Set& a, const Set& b ) { return  setIsEqual(a,b); } ///< compares the content of two sets
   friend inline bool operator!=( const Set& a, const Set& b ) { return !setIsEqual(a,b); } ///< compares the content of two sets

   inline Set() = default;
   inline Set( const T& element ) { set_.insert( element ); }

   inline virtual ~Set() = default;

   static const Set<T> emptySet() { return {}; }

   inline std::pair<iterator,bool> insert( const T& element )                    { return set_.insert( element ); }
   inline iterator                 insert( iterator position, const T& element ) { return set_.insert( position, element ); }
   template <class InputIterator>
   inline void insert( InputIterator first, InputIterator last ) { set_.insert( first, last); }

   inline void clear() { set_.clear(); } ///< removes all elements from this set

   inline const Set<T>& operator&=( const Set<T>& set ); ///< intersection
   inline const Set<T>& operator+=( const Set<T>& set ); ///< union
   inline const Set<T>& operator-=( const Set<T>& set ); ///< difference / relative complement

   inline bool operator< ( const Set<T>& set ) const { return set_.size() < set.set_.size(); }  ///< compares the size (not the content!) of two sets
   inline bool operator> ( const Set<T>& set ) const { return set_.size() > set.set_.size(); }  ///< compares the size (not the content!) of two sets
   inline bool operator<=( const Set<T>& set ) const { return !(operator>( set )); }            ///< compares the size (not the content!) of two sets
   inline bool operator>=( const Set<T>& set ) const { return !(operator<( set )); }            ///< compares the size (not the content!) of two sets
   inline bool equalSize ( const Set<T>& set ) const { return set_.size() == set.set_.size(); } ///< compares the size (not the content!) of two sets

   bool intersects( const Set<T>& set ) const { return set::setsIntersect( begin(), end(), set.begin(), set.end() ); } ///< true if both sets intersect
   bool contains  ( const Set<T>& set ) const { return std::includes( begin(), end(), set.begin(), set.end() ); }      ///< true if "set" is completely contained within this set
   bool contains  ( const T&  element ) const { return set_.find( element ) != set_.end(); }                           ///< true if "element" is contained within this set
   bool isEqual   ( const Set<T>& set ) const { return set_ == set.set_; }                                             ///< true if both sets contain the same elements

   inline bool   empty() const { return set_.empty(); } ///< true if this set is empty
   inline bool isEmpty() const { return empty(); }      ///< true if this set is empty

   inline size_t size() const { return set_.size(); }

   inline void swap( Set<T>& set ) { set_.swap( set.set_ ); }

          void        toStream( std::ostream& os ) const;
   inline std::string toString() const;

   inline const_iterator begin() const { return set_.begin(); }
   inline iterator       begin()       { return set_.begin(); }

   inline const_iterator end() const { return set_.end(); }
   inline iterator       end()       { return set_.end(); }
   
   inline const std::set<T> & get() const { return set_; }
   inline       std::set<T> & get()       { return set_; }

private:

   std::set<T> set_;

}; // class Set



//**********************************************************************************************************************
/*!
*   \brief Calculates the intersection of "this" and "set", only the resulting set is kept.
*/
//**********************************************************************************************************************
template< typename T >
inline const Set<T>& Set<T>::operator&=( const Set<T>& set ) {  // intersection

   std::set<T> result;
   std::set_intersection( begin(), end(), set.begin(), set.end(), std::inserter( result, result.end() ) );
   set_ = result;

   return *this;
}



//**********************************************************************************************************************
/*!
*   \brief Calculates the union of "this" and "set", only the resulting set is kept.
*/
//**********************************************************************************************************************
template< typename T >
inline const Set<T>&  Set<T>::operator+=( const Set<T>& set ) { // union

   set_.insert( set.begin(), set.end() );

   return *this;
}



//**********************************************************************************************************************
/*!
*   \brief Calculates the difference of "this" and "set", only the resulting set (result = this - set) is kept.
*/
//**********************************************************************************************************************
template< typename T >
inline const Set<T>& Set<T>::operator-=( const Set<T>& set ) { // difference / relative complement

   std::set<T> result;
   std::set_difference( begin(), end(), set.begin(), set.end(), std::inserter( result, result.end() ) );
   set_ = result;

   return *this;
}



template< typename T >
void Set<T>::toStream( std::ostream& os ) const {

   os << "{ ";

   for( const_iterator it = begin(); it != end(); ++it ) {
      const_iterator next = it;
      os << (*it) << ( ( ++next == end() ) ? " " : ", " );
   }

   os << "}";
}



template< typename T >
inline std::string Set<T>::toString() const {

   std::ostringstream oss;
   toStream( oss );

   return oss.str();
}



//////////////////////
// Global Functions //
//////////////////////



template< typename T >
inline Set<T> setIntersection( const Set<T>& a, const Set<T>& b ) {

   Set<T> result(a);
   result &= b;
   return result;
}



template< typename T >
inline Set<T> setUnion( const Set<T>& a, const Set<T>& b ) {

   Set<T> result(a);
   result += b;
   return result;
}



template< typename T >
inline Set<T> setDifference( const Set<T>& a, const Set<T>& b ) {

   Set<T> result(a);
   result -= b;
   return result;
}



template< typename T >
inline bool setIsEqual( const Set<T>& a, const Set<T>& b ) {

   return a.isEqual(b);
}



template< typename T >
inline std::ostream& operator<<( std::ostream& os, const Set<T>& set ) {

   set.toStream( os );
   return os;
}


/// \cond internal
namespace set {

template< class InputIterator1, class InputIterator2 >
bool setsIntersect( InputIterator1 first1, InputIterator1 last1, InputIterator2 first2, InputIterator2 last2 ) {

   while( first1 != last1 && first2 != last2) {
     if( *first1 < *first2 ) ++first1;
     else if( *first2 < *first1 ) ++first2;
     else return true;
   }

   return false;
}

} // namespace set
/// \endcond


} // namespace walberla



//======================================================================================================================
//
//  Send/Recv Buffer Serialization Specialization
//
//======================================================================================================================

namespace walberla {
namespace mpi {

template< typename T, // Element type of SendBuffer
          typename G, // Growth policy of SendBuffer
          typename E >
inline mpi::GenericSendBuffer<T,G> & operator<<( mpi::GenericSendBuffer<T,G> & buffer, const walberla::Set<E> & set )
{
   buffer.addDebugMarker( "se" );
   buffer << set.get();
   return buffer;
}

template< typename T, // Element type  of RecvBuffer
          typename E >
inline mpi::GenericRecvBuffer<T>& operator>>( mpi::GenericRecvBuffer<T> & buffer, walberla::Set<E> & set )
{
   buffer.readDebugMarker( "se" );
   buffer >> set.get();
   return buffer;
}

template< typename T >
struct BufferSizeTrait< walberla::Set<T> > { static const bool constantSize = false;  };

}
}
