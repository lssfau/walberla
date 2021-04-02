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
//! \file
//! \author Igor Ostanin <i.ostanin@skoltech.ru>
//! \author Grigorii Drozdov <drozd013@umn.edu>
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#pragma once

#include <string>

namespace walberla {

#ifdef PLAIN_OUTPUT
//----Terminal color manipulators-----
const std::string RED("");      // errors and problems
const std::string GREEN("");    // loading operations
const std::string YELLOW("");   // information (non-Walberla)
const std::string PURPLE("");   // stages of program implementation
const std::string CYAN("");     // save/write operations
const std::string RESET("");
//====================================
#else
//----Terminal color manipulators-----
const std::string RED("\033[1;31m");      // errors and problems
const std::string GREEN("\033[1;32m");    // loading operations
const std::string YELLOW("\033[1;33m");   // information (non-Walberla)
const std::string PURPLE("\033[1;35m");   // stages of program implementation
const std::string CYAN("\033[1;36m");     // save/write operations
const std::string RESET("\033[0m");
//====================================
#endif

} //namespace walberla