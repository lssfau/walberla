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
//! \file Regex.h
//! \ingroup core
//! \author Dominik Thoennes <dominik.thoennes@fau.de>
//
//======================================================================================================================

#pragma once


#include <regex>


namespace walberla {

using std::regex;
using std::regex_match;
using std::regex_error;
using std::regex_search;
using std::regex_replace;

}