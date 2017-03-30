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
//! \file CheckFunctions.cpp
//! \ingroup core
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//
//======================================================================================================================

#include "CheckFunctions.h"
#include "core/Abort.h"

#include <locale>


/// \cond internal

namespace walberla {
namespace debug {
namespace check_functions_detail {

void ExitHandler::operator()( const std::string & checkErrorMessage )
{
   std::ostringstream oss;
   oss << checkErrorMessage;
   if( !message_.empty() )
      oss << "\nMessage:\n" << message_;

   WALBERLA_ABORT( oss.str() );
}



template <typename T>
static std::ostream & printCharValue( std::ostream & os, const T value )
{
   os << int_c( value ) << " (0x";

   const std::istream::fmtflags oldFlags( os.flags() );
   os << std::hex << int_c( value ) << ", ";
   os.flags( oldFlags );

   if( std::isprint( value, std::locale() ) )
      os << "'" << value << "'";
   else
     os << "[not printable]";

   return os << ")";
}

std::ostream & printValue( std::ostream & os, const unsigned char value )
{
   return printCharValue( os, char(value) );
}

std::ostream & printValue( std::ostream & os, const char value )
{
   return printCharValue( os, value );
}

template< typename T >
static std::ostream & printFPValue( std::ostream & os, T value )
{
   return os << std::scientific << std::setprecision( std::numeric_limits<T>::digits10 + 2 ) << value;
}

std::ostream & printValue( std::ostream & os, float value )
{
   return printFPValue( os, value );
}

std::ostream & printValue( std::ostream & os, double value )
{
   return printFPValue( os, value );
}

std::ostream & printValue( std::ostream & os, long double value )
{
   return printFPValue( os, value );
}

} // namespace check_functions_detail
} // namespace debug
} // namespace walberla

/// \endcond
