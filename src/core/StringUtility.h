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
//! \file StringUtility.h
//! \ingroup core
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//
//======================================================================================================================

#pragma once

#include <algorithm>
#include <cctype>
#include <string>
#include <vector>

namespace walberla
{

// Convert (in place) every character in string to uppercase.
inline void string_to_upper(std::string &s);

// Convert (copy) every character in string to uppercase.
inline std::string string_to_upper_copy(const std::string &s);

// Convert (in place) every character in string to lowercase.
inline void string_to_lower(std::string &s);

// Convert (in place) every character in string to lowercase.
inline std::string string_to_lower_copy(const std::string &s);

// Remove (in place) all whitespaces at the beginning of a string.
inline void string_trim_left(std::string &s);

// Remove (in place) all whitespaces at the end of a string.
inline void string_trim_right(std::string &s);

// Remove (in place) all whitespaces at the beginning and at the end of a string.
inline void string_trim(std::string &s);

// Remove (copy) all whitespaces at the beginning of a string.
inline std::string string_trim_left_copy(const std::string &s);

// Remove (copy) all whitespaces at the end of a string.
inline std::string string_trim_right_copy(const std::string &s);

// Remove (copy) all whitespaces at the beginning and at the end of a string.
inline std::string string_trim_copy(const std::string &s);

// Split a string at the given delimiters into a vector of substrings.
// E.g. specify std::string(" |,") in order to split at characters ' ' and ','.
inline std::vector<std::string> string_split(std::string s, const std::string &delimiters);

// Replace (in place) all occurrences of substring "old" with substring "new".
inline void string_replace_all(std::string &s, const std::string &oldSubstr, const std::string &newSubstr);

// Replace (copy) all occurrences of substring "old" with substring "new".
inline std::string string_replace_all_copy(const std::string &s, const std::string &oldSubstr, const std::string &newSubstr);

// Check whether a string ends with a certain substring.
inline bool string_ends_with(const std::string &s, const std::string &substr);

// Case-insensitive std::string::compare.
inline int string_icompare(const std::string &s1, const std::string &s2);

} // namespace walberla

#include "core/StringUtility.impl.h"