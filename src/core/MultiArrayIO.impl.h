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
//! \file MultiArrayIO.h
//! \ingroup config
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/DataTypes.h"
#include "core/StringUtility.h"

#include <cstddef>
#include <sstream>

namespace boost {

//===================================================================================================================
//
//  Helper Functions
//
//===================================================================================================================


template<typename T>
bool parseArray1D( std::vector<T> & arr, std::istream & is,
                   const char openingBracket='[',
                   const char closingBracket=']',
                   const std::string & delimiter = ", \t\n" )
{

   is >> std::skipws;

   char bracket1;
   if( !(is >> bracket1 ) || bracket1 != openingBracket )
      return false;

   std::string line;
   if ( ! std::getline( is, line, closingBracket ) )
      return false;


   std::vector<std::string> stringArr = walberla::string_split( line, delimiter );

   arr.clear();
   arr.reserve( stringArr.size() );
   for( auto sArrIt =stringArr.begin(); sArrIt != stringArr.end(); ++sArrIt )
   {
      if ( *sArrIt == "") continue;

      std::stringstream ss ( *sArrIt );
      T value;
      ss >> value;
      arr.push_back( value );
   }
   return true;
}



template<typename T >
bool parseArray2D( std::vector< std::vector<T> > & arr, std::istream & is,
                   const char openingBracket='[',
                   const char closingBracket=']',
                   const std::string & delimiter = ", \t\n"  )
{
   is >> std::skipws;
   char bracket1;
   if( !(is >> bracket1 ) || bracket1 != openingBracket )
      return false;

   do
   {
      arr.push_back( std::vector<T>() );
      if ( ! parseArray1D( arr.back(), is, openingBracket, closingBracket, delimiter ) )
         return false;

      while( delimiter.find( (char)( is.peek()) ) != std::string::npos )
         is.get();

      if ( is.peek() == closingBracket ) {
         is.get();
         return true;
      }

   } while ( true );

}


//===================================================================================================================
//
//  IO Operators
//
//===================================================================================================================


template<typename T>
std::istream & operator>> ( std::istream & is, boost::multi_array<T,1> & arr )
{
   if ( !is ) return is;

   const std::istream::pos_type pos( is.tellg() );
   const std::istream::fmtflags oldFlags( is.flags() );


   std::vector< T > vec;
   if ( ! parseArray1D( vec, is ) ) {
      is.clear();
      is.seekg( pos );
      is.setstate( std::istream::failbit );
      is.flags( oldFlags );
      return is;
   }

   auto rows = vec.size();
   if ( rows == 0 )
      return is;

   arr.resize( boost::extents[walberla::numeric_cast< boost::multi_array_types::index >(rows)] );

   for( std::size_t r = 0; r < rows; ++r )
      arr[walberla::numeric_cast< boost::multi_array_types::index >(r)] = vec[r];

   return is;
}

template<typename T>
std::ostream & operator<< ( std::ostream & os, const boost::multi_array<T,1> & arr )
{
   os << "[ ";
   for( std::size_t c = 0; c < arr.size(); ++c )
      os << arr[walberla::numeric_cast< boost::multi_array_types::index >(c)] << ",";
   os << "]";

   return os;
}






template<typename T>
std::istream & operator>> ( std::istream & is, boost::multi_array<T,2> & arr )
{
   if ( !is ) return is;

   const std::istream::pos_type pos( is.tellg() );
   const std::istream::fmtflags oldFlags( is.flags() );


   std::vector< std::vector<T> > vec2D;
   if ( ! parseArray2D( vec2D, is ) ) {
      is.clear();
      is.seekg( pos );
      is.setstate( std::istream::failbit );
      is.flags( oldFlags );
      return is;
   }

   std::size_t rows = vec2D.size();
   if ( rows == 0 )
      return is;

   std::size_t cols = vec2D[0].size();
   for( std::size_t r = 0; r < rows; ++r )
   {
      if ( vec2D[r].size() != cols  )
      {
         // non square vector
         is.clear();
         is.seekg( pos );
         is.setstate( std::istream::failbit );
         is.flags( oldFlags );
         return is;
      }
   }

   arr.resize( boost::extents[ walberla::numeric_cast< boost::multi_array_types::index >(rows) ][ walberla::numeric_cast< boost::multi_array_types::index >(cols) ] );

   for( std::size_t r = 0; r < rows; ++r )
      for( std::size_t c = 0; c < cols; ++c )
         arr[walberla::numeric_cast< boost::multi_array_types::index >(r)][walberla::numeric_cast< boost::multi_array_types::index >(c)] = vec2D[r][c];


   return is;
}

template<typename T>
std::ostream & operator<< ( std::ostream & os, const boost::multi_array<T,2> & arr )
{
   os << "[\n";

   for( std::size_t r = 0; r < arr.size(); ++r )
   {
      os << " [ ";
      for( std::size_t c = 0; c < arr[walberla::numeric_cast< boost::multi_array_types::index >(r)].size(); ++c ) {
         os << arr[walberla::numeric_cast< boost::multi_array_types::index >(r)][walberla::numeric_cast< boost::multi_array_types::index >(c)] << "\t";
      }
      os << "] \n";
   }
   os << "]";

   return os;
}








} // namespace walberla


