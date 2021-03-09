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
//! \file AllSet.h
//! \ingroup core
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//! \author Christian Feichtinger
//
//======================================================================================================================

#pragma once

#include <algorithm>
#include <iostream>
#include <set>
#include <sstream>


namespace walberla{



template< typename T > class AllSet;

template< typename T > inline AllSet<T> setIntersection( const AllSet<T>& a, const AllSet<T>& b );
template< typename T > inline AllSet<T> setUnion       ( const AllSet<T>& a, const AllSet<T>& b );
template< typename T > inline AllSet<T> setDifference  ( const AllSet<T>& a, const AllSet<T>& b );

template< typename T > inline bool setIsEqual( const AllSet<T>& a, const AllSet<T>& b );



//**********************************************************************************************************************
/*!
*   \brief A class for managing sets that supports the ability to define a set that represents all possible elements
*
*   What is special about this set is the possibility to define a set that contains all possible elements:
*
    \code
      AllSet<int> universe   = AllSet<int>::all();                 // "universe" contains all integers
      bool        contains   = universe.contains( [any integer] ); // -> true
                  universe  -= 42;                                 // "universe" contains all integers except 42
                  contains   = universe.contains( 42 );            // -> false
      AllSet<int> set        = AllSet<int>(5) + AllSet<int>(23);   // a set that consists of two integers: 5 and 23
                  contains   = universe.contains( set );           // -> true
                  set       += 42;
                  contains   = universe.contains( set );           // -> false
      bool        intersects = universe.intersects( set );         // -> true
   \endcode
*
*   If the equality of two sets must be tested, the operators "==" and "!=" and the member function "isEqual" can be
*   used. If two sets must be compared in terms of size, the operators "<", ">", "<=", and ">=" and the member function
*   "equalSize" can be used:
*
    \code
      AllSet<int> aAll    = AllSet<int>::all() - 42;
      AllSet<int> bAll    = AllSet<int>::all();
      bool        boolean = ( aAll == bAll );       // -> false
                  boolean = aAll.equalSize( bAll ); // -> false
                  boolean = ( aAll < bAll );        // -> true
                  bAll   -= 23;
                  boolean = ( aAll == bAll );       // -> false
                  boolean = aAll.equalSize( bAll ); // -> true
                  boolean = ( aAll < bAll );        // -> false
                  aAll   -= 23;
                  bAll   -= 42;
                  boolean = ( aAll == bAll );       // -> true
                  boolean = aAll.equalSize( bAll ); // -> true
    \endcode
*/
//**********************************************************************************************************************

template< typename T >
class AllSet {

public:

   friend inline AllSet<T> operator&( const AllSet& a, const AllSet& b ) { return setIntersection(a,b); } ///< intersection
   friend inline AllSet<T> operator+( const AllSet& a, const AllSet& b ) { return setUnion       (a,b); } ///< union
   friend inline AllSet<T> operator-( const AllSet& a, const AllSet& b ) { return setDifference  (a,b); } ///< difference / relative complement

   friend inline bool operator==( const AllSet& a, const AllSet& b ) { return  setIsEqual(a,b); } ///< compares the content of two sets
   friend inline bool operator!=( const AllSet& a, const AllSet& b ) { return !setIsEqual(a,b); } ///< compares the content of two sets

   using ConstIter = typename std::set<T>::const_iterator;



   inline explicit AllSet( const bool _all = false ) : all_( _all ) {}
   inline          AllSet( const T& element )        : all_( false ) { set_.insert( element ); }

   inline static const AllSet<T>& all()      { static AllSet set( true  ); return set; }
   inline static const AllSet<T>& emptySet() { static AllSet set( false ); return set; }

   void   invert() { all_ = !all_; }
   AllSet<T> getComplement() const { AllSet<T> set( *this ); set.all_ = !all_; return set; }

   void insert( const T& element ) { if( isAll() ) set_.erase( element ); else set_.insert( element ); }
   void clear() { all_ = false; set_.clear(); }

   const AllSet<T>& operator&=( const AllSet<T>& set ); ///< intersection
   const AllSet<T>& operator+=( const AllSet<T>& set ); ///< union
   const AllSet<T>& operator-=( const AllSet<T>& set ); ///< difference / relative complement

          bool operator< ( const AllSet<T>& set ) const;                                ///< compares the size (not the content!) of two sets
          bool operator> ( const AllSet<T>& set ) const;                                ///< compares the size (not the content!) of two sets
   inline bool operator<=( const AllSet<T>& set ) const { return !(operator>( set )); } ///< compares the size (not the content!) of two sets
   inline bool operator>=( const AllSet<T>& set ) const { return !(operator<( set )); } ///< compares the size (not the content!) of two sets
   inline bool equalSize ( const AllSet<T>& set ) const;                                ///< compares the size (not the content!) of two sets

          bool intersects( const AllSet<T>& set ) const;
          bool contains  ( const AllSet<T>& set ) const;
   inline bool contains  ( const T&     element ) const;
          bool isEqual   ( const AllSet<T>& set ) const { return all_ == set.all_ && set_ == set.set_; }

   inline bool isUniverse()  const { return  all_ && set_.empty(); }
   inline bool isEmpty()     const { return !all_ && set_.empty(); }
   inline bool isAll()       const { return  all_; }
   inline bool isCountable() const { return !all_; }

          void        toStream( std::ostream& os ) const;
   inline std::string toString() const;

   inline ConstIter begin() const { return set_.begin(); }
   inline ConstIter end()   const { return set_.end(); }

private:

   bool all_; ///< true if the set represents a set of all possible elements (= "the universe" / the complement of the empty set)
   std::set<T> set_;

}; // class AllSet



//**********************************************************************************************************************
/*!
*   \brief Calculates the intersection of "this" and "set", only the resulting set is kept.
*/
//**********************************************************************************************************************
template< typename T >
const AllSet<T>& AllSet<T>::operator&=( const AllSet<T>& set ) {  // intersection

   if( isAll() ) {
      if( set.isAll() ) {
         set_.insert( set.begin(), set.end() );
      }
      else {
         std::set<T> result;
         std::set_difference( set.begin(), set.end(), begin(), end(), std::inserter( result, result.end() ) );
         set_ = result;
         all_ = false;
      }
   }
   else {
      if( set.isAll() ) {
         std::set<T> result;
         std::set_difference( begin(), end(), set.begin(), set.end(), std::inserter( result, result.end() ) );
         set_ = result;
      }
      else {
         std::set<T> result;
         std::set_intersection( begin(), end(), set.begin(), set.end(), std::inserter( result, result.end() ) );
         set_ = result;
      }
   }

   return *this;
}



//**********************************************************************************************************************
/*!
*   \brief Calculates the union of "this" and "set", only the resulting set is kept.
*/
//**********************************************************************************************************************
template< typename T >
const AllSet<T>&  AllSet<T>::operator+=( const AllSet<T>& set ) { // union

   if( isAll() ) {
      if( set.isAll() ) {
         std::set<T> result;
         std::set_intersection( begin(), end(), set.begin(), set.end(), std::inserter( result, result.end() ) );
         set_ = result;
      }
      else {
         std::set<T> result;
         std::set_difference( begin(), end(), set.begin(), set.end(), std::inserter( result, result.end() ) );
         set_ = result;
      }
   }
   else {
      if( set.isAll() ) {
         std::set<T> result;
         std::set_difference( set.begin(), set.end(), begin(), end(), std::inserter( result, result.end() ) );
         set_ = result;
         all_ = true;
      }
      else {
         set_.insert( set.begin(), set.end() );
      }
   }

   return *this;
}



//**********************************************************************************************************************
/*!
*   \brief Calculates the difference of "this" and "set", only the resulting set (result = this - set) is kept.
*/
//**********************************************************************************************************************
template< typename T >
const AllSet<T>& AllSet<T>::operator-=( const AllSet<T>& set ) { // difference / relative complement

   if( isAll() ) {
      if( set.isAll() ) {
         std::set<T> result;
         std::set_difference( set.begin(), set.end(), begin(), end(), std::inserter( result, result.end() ) );
         set_ = result;
         all_ = false;
      }
      else {
         set_.insert( set.begin(), set.end() );
      }
   }
   else {
      if( set.isAll() ) {
         std::set<T> result;
         std::set_intersection( begin(), end(), set.begin(), set.end(), std::inserter( result, result.end() ) );
         set_ = result;
      }
      else {
         std::set<T> result;
         std::set_difference( begin(), end(), set.begin(), set.end(), std::inserter( result, result.end() ) );
         set_ = result;
      }
   }

   return *this;
}



//**********************************************************************************************************************
/*!
*   \brief Returns true only if there are less elements in "this" set than in the set "set".
*/
//**********************************************************************************************************************
template< typename T >
bool AllSet<T>::operator<( const AllSet<T>& set ) const {

   if( isAll() ) {
      if( set.isAll() )
         return set_.size() > set.set_.size();
      else
         return false;
   }
   else if( set.isAll() )
      return true;

   return set_.size() < set.set_.size();
}



//**********************************************************************************************************************
/*!
*   \brief Returns true only if there are more elements in "this" set than in the set "set".
*/
//**********************************************************************************************************************
template< typename T >
bool AllSet<T>::operator>( const AllSet<T>& set ) const {

   if( isAll() ) {
      if( set.isAll() )
         return set_.size() < set.set_.size();
      else
         return true;
   }
   else if( set.isAll() )
      return false;

   return set_.size() > set.set_.size();
}



//**********************************************************************************************************************
/*!
*   \brief Returns true only if there are as many elements in "this" set as there are elements in the set "set".
*/
//**********************************************************************************************************************
template< typename T >
inline bool AllSet<T>::equalSize( const AllSet<T>& set ) const {

   return all_ == set.all_ && set_.size() == set.set_.size();
}


/// \cond internal
namespace allset {
template< class InputIterator1, class InputIterator2 >
bool setsIntersect( InputIterator1 first1, InputIterator1 last1, InputIterator2 first2, InputIterator2 last2 );
} // namespace allset
/// \endcond

//**********************************************************************************************************************
/*!
*   \brief Returns true only if there is an intersection between "this" and "set".
*/
//**********************************************************************************************************************
template< typename T >
bool AllSet<T>::intersects( const AllSet<T>& set ) const {

   if( isAll() ) {
      if( set.isAll() ) {
         return true;
      }
      else {
         return !std::includes( begin(), end(), set.begin(), set.end() );
      }
   }
   else if( set.isAll() ) {
      return !std::includes( set.begin(), set.end(), begin(), end() );
   }

   return allset::setsIntersect( begin(), end(), set.begin(), set.end() );
}



//**********************************************************************************************************************
/*!
*   \brief Returns true only if "set" is a strict subset of "this".
*/
//**********************************************************************************************************************
template< typename T >
bool AllSet<T>::contains( const AllSet<T>& set ) const {

   if( isAll() ) {
      if( set.isAll() ) {
         return std::includes( set.begin(), set.end(), begin(), end() );
      }
      else {
         return !allset::setsIntersect( begin(), end(), set.begin(), set.end() );
      }
   }
   else if( set.isAll() ) {
      return false;
   }

   return std::includes( begin(), end(), set.begin(), set.end() );
}



//**********************************************************************************************************************
/*!
*   \brief Returns true only if the element "element" is included in this set.
*/
//**********************************************************************************************************************
template< typename T >
inline bool AllSet<T>::contains( const T& element ) const {

   if( isAll() )
      return set_.find( element ) == set_.end();

   return set_.find( element ) != set_.end();
}



template< typename T >
void AllSet<T>::toStream( std::ostream& os ) const {

   os << "{ ";

   if( isAll() ) {
      os << "ALL ";
      if( !(set_.empty()) ) os << "- ";
   }

   for( ConstIter it = begin(); it != end(); ++it ) {
      ConstIter next = it;
      os << (*it) << ( ( ++next == end() ) ? " " : ", " );
   }

   os << "}";
}



template< typename T >
inline std::string AllSet<T>::toString() const {

   std::ostringstream oss;
   toStream( oss );

   return oss.str();
}



//////////////////////
// Global Functions //
//////////////////////



template< typename T >
inline AllSet<T> setIntersection( const AllSet<T>& a, const AllSet<T>& b ) {

   AllSet<T> result(a);
   result &= b;
   return result;
}



template< typename T >
inline AllSet<T> setUnion( const AllSet<T>& a, const AllSet<T>& b ) {

   AllSet<T> result(a);
   result += b;
   return result;
}



template< typename T >
inline AllSet<T> setDifference( const AllSet<T>& a, const AllSet<T>& b ) {

   AllSet<T> result(a);
   result -= b;
   return result;
}



template< typename T >
inline bool setIsEqual( const AllSet<T>& a, const AllSet<T>& b ) {

   return a.isEqual(b);
}



template< typename T >
inline std::ostream& operator<<( std::ostream& os, const AllSet<T>& set ) {

   set.toStream( os );
   return os;
}


/// \cond internal
namespace allset {

template< class InputIterator1, class InputIterator2 >
bool setsIntersect( InputIterator1 first1, InputIterator1 last1, InputIterator2 first2, InputIterator2 last2 ) {

   while( first1 != last1 && first2 != last2) {
     if( *first1 < *first2 ) ++first1;
     else if( *first2 < *first1 ) ++first2;
     else return true;
   }

   return false;
}

} // namespace allset
/// \endcond


} // namespace walberla



