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
//! \file RandomUUID.cpp
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include "RandomUUID.h"

#include <iomanip>
#include <random>
#include <sstream>

namespace walberla {

RandomUUID::RandomUUID()
{
   std::random_device rd;
   std::mt19937 gen( rd() );
   std::uniform_int_distribution<uint64_t> dis(std::numeric_limits<uint64_t>::min(), std::numeric_limits<uint64_t>::max());
   a_ = dis(gen);
   std::random_device rd2;
   std::mt19937 gen2( rd2() );
   b_ = dis(gen2);
   //unsure how independent these two random numbers really are

   //setting version numbers
   //https://en.m.wikipedia.org/wiki/Universally_unique_identifier
   a_ = ((a_ & 0xFFFFFFFFFFFF0FFF) | 0x0000000000004000);
   b_ = ((b_ & 0x3FFFFFFFFFFFFFFF) | 0x8000000000000000);
}

RandomUUID::RandomUUID( const UIntType a, const UIntType b)
   : a_(a)
   , b_(b)
{}

std::string RandomUUID::toString() const
{
   std::stringstream ss;
   ss << std::hex
      << std::setfill('0')
      << std::setw(8) << (a_ >> 32) << "-"
      << std::setw(4) << ((a_&0xFFFFFFFF)>>16) << "-"
      << std::setw(4) << (a_&0xFFFF) << "-"
      << std::setw(4) << (b_ >> 48) << "-"
      << std::setw(12) << (b_&0xFFFFFFFFFFFF);
   return ss.str();
}

bool operator==(const RandomUUID& lhs, const RandomUUID& rhs)
{
   return (lhs.a_ == rhs.a_) && (lhs.b_ == rhs.b_);
}

bool operator!=(const RandomUUID& lhs, const RandomUUID& rhs)
{
   return (lhs.a_ != rhs.a_) || (lhs.b_ != rhs.b_);
}

std::ostream& operator<<(std::ostream& os, const RandomUUID& uuid)
{
   os << uuid.toString();
   return os;
}

}
