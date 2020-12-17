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
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#pragma once

#include <mesa_pd/data/DataTypes.h>

#include <map>
#include <set>
#include <unordered_set>
#include <unordered_map>
#include <vector>

namespace walberla {
namespace mesa_pd {
namespace data {

template< typename Type >
std::ostream& operator<<( std::ostream& os, const std::vector<Type>& v )
{
   os << "<";
   for (auto it = v.begin(); it != v.end(); ++it)
   {
      os << *it;
      if (it != v.end()-1)
         os << ",";
   }
   os << ">";
   return os;
}

template< typename Key, typename T, typename Compare, typename Allocator >
std::ostream& operator<<( std::ostream& os, const std::map<Key, T, Compare, Allocator>& m )
{
   os << "{";
   for (auto& v : m)
   {
      os << v.first << ":" << v.second << ", ";
   }
   os << "}";
   return os;
}

template< typename Key, typename Compare, typename Allocator >
std::ostream& operator<<( std::ostream& os, const std::set<Key, Compare, Allocator>& m )
{
   os << "{";
   for (auto& v : m)
   {
      os << v << ", ";
   }
   os << "}";
   return os;
}

template< typename Key, typename T, typename Compare, typename Allocator >
std::ostream& operator<<( std::ostream& os, const std::unordered_map<Key, T, Compare, Allocator>& m )
{
   os << "{";
   for (auto& v : m)
   {
      os << v.first << ":" << v.second << ", ";
   }
   os << "}";
   return os;
}

template< typename Key, typename Compare, typename Allocator >
std::ostream& operator<<( std::ostream& os, const std::unordered_set<Key, Compare, Allocator>& m )
{
   os << "{";
   for (auto& v : m)
   {
      os << v << ", ";
   }
   os << "}";
   return os;
}

} //namespace data
} //namespace mesa_pd
} //namespace walberla
