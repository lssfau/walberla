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
//! \file Base64Writer.cpp
//! \ingroup vtk
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#include "Base64Writer.h"
#include "core/DataTypes.h"


namespace walberla {
namespace vtk {



//**********************************************************************************************************************
/*!
*   This function can be used in order to convert all the data that was passed to this Base64Writer object into a base64
*   encoded stream. In accordance with the VTK specification, the beginning of the base64 encoded stream will contain
*   a 32bit unsigned int value specifying the size (in bytes) of the data that follows. For more information on how a
*   Base64Writer is typically used refer to the documentation/general description of this class.
*/
//**********************************************************************************************************************
void Base64Writer::toStream( std::ostream& os )
{
   unsigned char input[3];
   unsigned char output[4];

   const uint32_t bytes = uint32_c( buffer_.size() );

   buffer_.insert( buffer_.begin(), 
                   reinterpret_cast<const char*>( &bytes ),
                   reinterpret_cast<const char*>( &bytes ) + sizeof( uint32_t ) );

   for ( uint_t i = 0; i < buffer_.size(); i += 3 ) {
      if ( buffer_.size() - i < 3 ) {
         const uint_t length = buffer_.size() - i;
         for ( uint_t j = 0; j < length; ++j )
            input[j] = static_cast<unsigned char>( buffer_[i+j] );
         for ( uint_t j = length; j < 3; ++j )
            input[j] = static_cast<unsigned char>(0);
         encodeblock( input, output, int_c( length ) );
         os << output[0] << output[1] << output[2] << output[3];
      }
      else {
         input[0] = static_cast<unsigned char>( buffer_[ i ] );
         input[1] = static_cast<unsigned char>( buffer_[i+1] );
         input[2] = static_cast<unsigned char>( buffer_[i+2] );
         encodeblock( input, output, 3 );
         os << output[0] << output[1] << output[2] << output[3];
      }
   }

   os << "\n";
   buffer_.clear();
}



} // namespace vtk
} // namespace walberla
