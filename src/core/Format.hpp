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
//! \file Format.h
//! \ingroup core
//! \author Dominik Thoennes <dominik.thoennes@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/Abort.h"

#include <cstdio>
#include <string>

namespace walberla
{
/// format uses the printf syntax to format a given formatString and return a std::string
///\param formatString the format string
///\param args all arguments which will be inserted into the string, these CANNOT be std::string but have to be
/// converted using .c_str()
template< typename... Args >
std::string format(const std::string& formatString, Args&&... args)
{
   /// this size is arbitrary
   const size_t maxBufferSize = 4096;
   char buffer[maxBufferSize];
   int check = snprintf(buffer, maxBufferSize, formatString.c_str(), args...);
   if (check <= 0 || check > int(maxBufferSize)) { WALBERLA_ABORT("snprintf failed"); }
   return std::string(buffer);
}

} // namespace walberla