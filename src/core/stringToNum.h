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
//! \file stringToNum.h
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include <string>

namespace walberla {

template<typename S>
inline S stringToNum( const std::string & t );
template <> inline float              stringToNum( const std::string & t ) { return std::stof(t); }
template <> inline double             stringToNum( const std::string & t ) { return std::stod(t); }
template <> inline long double        stringToNum( const std::string & t ) { return std::stold(t); }
template <> inline int                stringToNum( const std::string & t ) { return std::stoi(t); }
template <> inline long               stringToNum( const std::string & t ) { return std::stol(t); }
template <> inline long long          stringToNum( const std::string & t ) { return std::stoll(t); }
template <> inline unsigned long      stringToNum( const std::string & t ) { return std::stoul(t); }
template <> inline unsigned long long stringToNum( const std::string & t ) { return std::stoull(t); }

} // namespace walberla


