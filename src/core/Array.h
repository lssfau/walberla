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
//! \file Array.h
//! \ingroup core
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "DataTypes.h"

#include "core/debug/Debug.h"

#include <iostream>
#include <sstream>
#include <vector>


namespace walberla {



/// Fixed size, dynamically allocated array

template< typename T >
class Array {

public:

   inline Array() : array_( nullptr ), size_( uint_c(0) ) {}
   inline Array( const uint_t n, const T& t = T() );
   inline Array( const std::vector<T>& vector );
   inline Array( const Array& array );

   ~Array() { if( array_ != nullptr ) delete[] array_; }

   uint_t size() const { return size_; }
   bool  empty() const { return size_ == 0; }

   inline bool operator==( const Array& array ) const;
   inline bool operator!=( const Array& array ) const { return !operator==( array );}

   const T& operator[]( const uint_t index ) const { WALBERLA_ASSERT_LESS( index, size_ ); return array_[ index ]; }
         T& operator[]( const uint_t index )       { WALBERLA_ASSERT_LESS( index, size_ ); return array_[ index ]; }

   inline const T* begin() const { return array_; }
   inline T*       begin()       { return array_; }

   inline const T* end() const { return array_ + size_; }
   inline T*       end()       { return array_ + size_; }

   inline void swap( Array& array );

          void        toStream( std::ostream& os ) const;
   inline std::string toString() const;

protected:

   T*     array_;
   uint_t size_;

}; // class Array



template< typename T >
inline Array<T>::Array( const uint_t n, const T& t ) : array_( n == 0 ? NULL : new T[n] ), size_( n )
{
   for( uint_t i = 0; i != n; ++i )
      array_[i] = t;
}



template< typename T >
inline Array<T>::Array( const std::vector<T>& vector ) :

   array_( vector.size() == 0 ? NULL : new T[ vector.size() ] ), size_( vector.size() )
{
   for( uint_t i = 0; i != size_; ++i )
      array_[i] = vector[i];
}



template< typename T >
inline Array<T>::Array( const Array& array ) :

   array_( array.size_ == 0 ? NULL : new T[ array.size_ ] ), size_( array.size_ )
{
   for( uint_t i = 0; i != size_; ++i )
      array_[i] = array.array_[i];
}



template< typename T >
inline bool Array<T>::operator==( const Array& array ) const
{
   if( size_ != array.size_ )
      return false;

   for( uint_t i = 0; i != size_; ++i )
      if( array_[i] != array.array_[i] ) return false;

   return true;
}



template< typename T >
inline void Array<T>::swap( Array& array )
{
   T* ptr = array_;
   array_ = array.array_;
   array.array_ = ptr;

   uint_t s = size_;
   size_ = array.size_;
   array.size_ = s;
}



template< typename T >
void Array<T>::toStream( std::ostream& os ) const {

   os << "{ ";

   for( auto it = begin(); it != end(); ++it ) {
      auto next = it;
      os << (*it) << ( ( ++next == end() ) ? " " : ", " );
   }

   os << "}";
}



template< typename T >
inline std::string Array<T>::toString() const {

   std::ostringstream oss;
   toStream( oss );

   return oss.str();
}



//////////////////////
// Global Functions //
//////////////////////



template< typename T >
inline std::ostream& operator<<( std::ostream& os, const Array<T>& array ) {

   array.toStream( os );
   return os;
}



} // namespace walberla


