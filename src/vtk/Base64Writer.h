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
//! \file Base64Writer.h
//! \ingroup vtk
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include <ostream>
#include <vector>


namespace walberla {
namespace vtk {



//**********************************************************************************************************************
/*!
*   \brief Class for encoding VTK compatible base64 data streams
*
*   This class is typically used as follows: First, an instance of type Base64Writer is created. It is then filled with
*   data using the left shift operator. Finally, the base64 encoded data is written to a VTK file by calling the member
*   function toStream( std::ostream& os ) where "os" is a reference to an output stream which is attached to the VTK
*   output file:
*   \code
*      std::ofstream vtkFile( [filename] );
*      [...]
*      Base64Writer base64;
*      for( ... )
*         base64 << [data];
*      base64.toStream( vtkFile ); // In accordance with the VTK specification, the beginning of the base64 encoded
*                                  // stream will contain a 32bit unsigned int value specifying the size (in bytes) of
*                                  // the data that follows.
*      [...]
*      vtkFile.close();
*   \endcode
*/
//**********************************************************************************************************************

class Base64Writer {

public:

   //*******************************************************************************************************************
   /*!
   *   This function can be used in order to add data to a Base64Writer instance which then can be converted to a base64
   *   encoded stream. For more information on how a Base64Writer is typically used refer to the documentation/general
   *   description of this class.
   */
   //*******************************************************************************************************************
   template< typename T > Base64Writer& operator<<( const T& data )
   {
      const char * bytePointer = reinterpret_cast<const char*>( &data );
      buffer_.insert( buffer_.end(), bytePointer, bytePointer + sizeof( T ) );
      return *this;
   }

   void toStream( std::ostream& os );

private:

   void encodeblock( unsigned char in[3], unsigned char out[4], int len )
   {
      static const unsigned char cb64[]="ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";

      switch( len )
      {
      case 1:
         out[0] = cb64[ in[0] >> 2 ];
         out[1] = cb64[ ((in[0] & 0x03) << 4) ];
         out[2] = '=';
         out[3] = '=';
         break;
      case 2:
         out[0] = cb64[ in[0] >> 2 ];
         out[1] = cb64[ ((in[0] & 0x03) << 4) | ((in[1] & 0xf0) >> 4) ];
         out[2] = cb64[ ((in[1] & 0x0f) << 2) ];
         out[3] = '=';
         break;
      default:
         out[0] = cb64[ in[0] >> 2 ];
         out[1] = cb64[ ((in[0] & 0x03) << 4) | ((in[1] & 0xf0) >> 4) ];
         out[2] = cb64[ ((in[1] & 0x0f) << 2) | ((in[2] & 0xc0) >> 6) ];
         out[3] = cb64[ in[2] & 0x3f ];
      }
   }

   std::vector<char> buffer_;

}; // class Base64Writer



} // namespace vtk
} // namespace walberla


