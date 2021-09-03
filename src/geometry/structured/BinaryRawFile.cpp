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
//! \file BinaryRawFile.cpp
//! \ingroup geometry
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//
//======================================================================================================================

#include "BinaryRawFile.h"


namespace walberla {
namespace geometry {

BinaryRawFile::BinaryRawFile( const std::string & filename, const Vector3< uint_t > & size, const std::string & datatype )
   : size_( size )
{
   init( filename, datatype );
}

BinaryRawFile::BinaryRawFile( const Config::BlockHandle & configBlock )
{
   size_ = configBlock.getParameter< Vector3< uint_t > >( "size" );
   const std::string filename = configBlock.getParameter< std::string >( "filename" );
   const std::string datatype = configBlock.getParameter< std::string >( "datatype" );

   init( filename, datatype );
}

void BinaryRawFile::init( const std::string & filename, const std::string & datatype )
{
   if (datatype == "uint8")
      init<uint8_t>( filename );
   else if (datatype == "int8")
      init<int8_t>( filename );
   else if (datatype == "uint16")
      init<uint16_t>( filename );
   else if (datatype == "int16")
      init<int16_t>( filename );
   else if (datatype == "uint32")
      init<uint32_t>( filename );
   else if (datatype == "int32")
      init<int32_t>( filename );
   else if (datatype == "uint64")
      init<uint64_t>( filename );
   else if (datatype == "int64")
      init<int64_t>( filename );
   else if (datatype == "float")
      init<float>( filename );
   else if (datatype == "double")
      init<double>( filename );
   else
      WALBERLA_ABORT( "Unknown datatype!" );
}


} // namespace geometry
} // namespace walberla