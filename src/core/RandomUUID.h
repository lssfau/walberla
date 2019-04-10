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
//! \file RandomUUID.h
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#pragma once

#include <ostream>
#include <string>

namespace walberla {

/**
 * Replacement for boost::uuids::uuid and boost::uuids::random_generator
 *
 * Uses two 64 bit random numbers to create a 128 bit uuid.
 */
class RandomUUID
{
   friend bool operator==(const RandomUUID& lhs, const RandomUUID& rhs);
   friend bool operator!=(const RandomUUID& lhs, const RandomUUID& rhs);
public:
   using UIntType = uint64_t;

   RandomUUID();
   RandomUUID(const UIntType a, const UIntType b);

   /**
    * returns a string representation of the uuid
    *
    * format: hhhhhhhh-hhhh-hhhh-hhhh-hhhhhhhhhhhh
    */
   std::string toString() const;

   UIntType getFirstUInt() const {return a_;}
   UIntType getSecondUInt() const {return b_;}
private:
   UIntType a_; ///< first part of the uuid
   UIntType b_; ///< second part of the uuid
};

bool operator==(const RandomUUID& lhs, const RandomUUID& rhs);
bool operator!=(const RandomUUID& lhs, const RandomUUID& rhs);
std::ostream& operator<<(std::ostream& os, const RandomUUID& uuid);

}
