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
//! \file Hostname.h
//! \ingroup core
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#include "Hostname.h"

#include "waLBerlaDefinitions.h"

#ifdef WALBERLA_CXX_COMPILER_IS_MSVC
#   include <Winsock2.h>
#   pragma comment(lib, "Ws2_32.lib")
#elif defined(_SX)
#  include <sys/socket.h>
#else
#   include <unistd.h>
#endif

namespace walberla {

std::string getHostName()
{
   char hostname[255];
   gethostname( hostname, 255 );
   return std::string( hostname );
}

} // namespace walberla
