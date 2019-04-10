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
//! \file StringUtility.impl.h
//! \ingroup core
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/StringUtility.h"

#include <algorithm>
#include <cctype>
#include <string>
#include <vector>

namespace walberla
{
// Convert (in place) every character in string to uppercase.
inline void string_to_upper(std::string &s) {
   std::transform(s.begin(), s.end(), s.begin(), [](char c){ return static_cast<char>(std::toupper(static_cast<unsigned char>(c))); });
}

// Convert (copy) every character in string to uppercase.
inline std::string string_to_upper_copy(const std::string &s) {
   std::string result = s;
   string_to_upper(result);
   return result;
}

// Convert (in place) every character in string to lowercase.
inline void string_to_lower(std::string &s) {
   std::transform(s.begin(), s.end(), s.begin(), [](char c){ return static_cast<char>(std::tolower(static_cast<unsigned char>(c))); });
}

// Convert (copy) every character in string to lowercase.
inline std::string string_to_lower_copy(const std::string &s) {
   std::string result = s;
   string_to_lower(result);
   return result;
}

// Remove (in place) all whitespaces at the beginning of a string.
inline void string_trim_left(std::string &s) {
   s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](char c) { return ! std::isspace(static_cast<unsigned char>(c)); }));
}

// Remove (in place) all whitespaces at the end of a string.
inline void string_trim_right(std::string &s) {
   s.erase(std::find_if(s.rbegin(), s.rend(), [](char c) { return ! std::isspace(static_cast<unsigned char>(c)); }).base(), s.end());
}

// Remove (in place) all whitespaces at the beginning and at the end of a string.
inline void string_trim(std::string &s) {
   string_trim_left(s);
   string_trim_right(s);
}

// Remove (copy) all whitespaces at the beginning of a string.
inline std::string string_trim_left_copy(const std::string &s) {
   std::string result = s;
   string_trim_left(result);
   return result;
}

// Remove (copy) all whitespaces at the end of a string.
inline std::string string_trim_right_copy(const std::string &s) {
   std::string result = s;
   string_trim_left(result);
   return result;
}

// Remove (copy) all whitespaces at the beginning and at the end of a string.
inline std::string string_trim_copy(const std::string &s) {
   std::string result = s;
   string_trim_left(result);
   return result;
}

// Split a string at the given delimiters into a vector of substrings.
// E.g. specify std::string(" ,") in order to split at characters ' ' and ','.
inline std::vector<std::string> string_split(std::string s, const std::string &delimiters) {
   std::vector<std::string> substrings;

   auto sub_begin = s.begin();   // iterator to the begin and end of a substring
   auto sub_end = sub_begin;

   for (auto it = s.begin(); it != s.end(); ++it) {
      for (auto d : delimiters) {
         if (*it == d) {   // current character in s is a delimiter
            sub_end = it;
            if (sub_begin < sub_end) { // make sure that the substring is not empty
               substrings.push_back(std::string(sub_begin, sub_end));
            }
            sub_begin = ++sub_end;
            continue;
         }
      }
   }

   // add substring from last delimiter to the end of s
   if (sub_begin < s.end()) {
      substrings.push_back(std::string(sub_begin, s.end()));
   }

   return substrings;
}

// Replace (in place) all occurrences of substring "old" with substring "new".
inline void string_replace_all(std::string &s, const std::string &oldSubstr, const std::string &newSubstr) {
   // loop written to avoid infinite-loops when newSubstr contains oldSubstr
   for (size_t pos = s.find(oldSubstr); pos != std::string::npos;) {
      s.replace(pos, oldSubstr.length(), newSubstr);
      pos = s.find(oldSubstr, pos + newSubstr.length());
   }
}

// Replace (copy) all occurrences of substring "old" with substring "new".
inline std::string string_replace_all_copy(const std::string &s, const std::string &oldSubstr, const std::string &newSubstr) {
   std::string result = s;
   string_replace_all(result, oldSubstr, newSubstr);
   return result;
}

// Check whether a string ends with a certain substring.
inline bool string_ends_with(const std::string &s, const std::string &substr) {
   return s.rfind(substr) == (s.length() - substr.length());
}

// Case-insensitive wrapper for std::string::compare (return values as for std::string::compare).
inline int string_icompare(const std::string &s1, const std::string &s2) {
   // function std::string::compare returns 0 in case of equality => invert result to obtain expected bool behavior
   return string_to_lower_copy(s1).compare(string_to_lower_copy(s2));
}

} // namespace walberla
