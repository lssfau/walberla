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
//! \file STLIO.h
//! \ingroup core
//! \author Michael Kuron <mkuron@icp.uni-stuttgart.de>
//
//======================================================================================================================

#pragma once

#include <vector>
#include <iostream>

namespace walberla {

template<typename T>
class ParameterList: public std::vector<T> {
public:
	ParameterList(): std::vector<T>() {};
	template<typename C> ParameterList(C v): std::vector<T>(v.begin(), v.end()) {}
	template<typename C> operator C() { return C(std::vector<T>::begin(), std::vector<T>::end()); }
};

//**********************************************************************************************************************
/*!\fn std::ostream& operator<<( std::ostream& os, const vector<Type>& v )
// \brief Global output operator for vectors.
//
// \param os Reference to the output stream.
// \param v Reference to a constant vector object.
// \return Reference to the output stream.
*/
template< typename Type >
std::ostream& operator<<( std::ostream& os, const ParameterList<Type>& v )
{
   os << "<";
   for (auto it = v.begin(); it != v.end(); ++it)
   {
      os << *it;
      if (it != v.end()-1)
         os << ",";
   }
   os << ">";
   return os;
}
//**********************************************************************************************************************

//**********************************************************************************************************************
/*! RAII Object to read an STL container (e.g. vector or list) from a stream
*/
//**********************************************************************************************************************

template< typename Container >
class ContainerStreamReader
{
public:
   ContainerStreamReader( std::istream& is ) : is_( is ), success_(false)
   {
   }
   
   ~ContainerStreamReader()
   {
      if (!success_)
         setFailed();
      restoreFlags();
   }
   
   void read( Container & v )
   {
      if( !is_ )
         return;

      storeState();

      is_ >> std::skipws;

      char character;
      if( !( is_ >> character ) || character != '<' )
         return;

      typename Container::value_type x;

      do {
         if( !( is_ >> x ) )
            return;
         
         v.push_back( x );
         
         if( !( is_ >> character ) )
            return;
         
      } while( character == ',' );

      if( character != '>' )
         return;
      
      success_ = true;
   }

private:

   void storeState()
   {
      oldPos_   = is_.tellg();
      oldFlags_ = is_.flags();
   }
   
   void restoreFlags()
   {
      is_.flags( oldFlags_ );
   }
   
   void setFailed() 
   {
      is_.clear();
      is_.seekg( oldPos_ );
      is_.setstate( std::istream::failbit );
   }
   
   std::istream& is_;
   std::istream::pos_type oldPos_;
   std::istream::fmtflags oldFlags_;
   bool success_;   
};

//**********************************************************************************************************************
/*!\fn std::istream& operator>>( std::istream& is, vector<Type>& v )
// \brief Global input operator for vectors.
//
// \param is Reference to the input stream.
// \param v Reference to a vector object.
// \return The input stream.
*/
template< typename Type >
std::istream& operator>>( std::istream& is, ParameterList<Type>& v )
{
   ContainerStreamReader< ParameterList<Type> > reader( is );
   reader.read( v );
   return is;
}
//**********************************************************************************************************************

}
