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
//! \file IFCD.h
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#pragma once

#include "pe/Types.h"
#include "pe/contact/Contact.h"
#include "core/NonCopyable.h"

#include <vector>

namespace walberla {
namespace pe {

class BodyStorage;

namespace fcd {

class IFCD : private NonCopyable{
public:
   virtual ~IFCD() {};

   virtual Contacts& generateContacts(PossibleContacts& possibleContacts) = 0;
   Contacts& getContacts() {return contacts_;}
protected:
   Contacts contacts_;
};

//*************************************************************************************************
/*!\brief Compare if two fine collision detectors are equal.
 *
 * Since collision detectors are uncopyable two collision detectors are considered equal if their adresses are equal.
 */
inline bool operator==(const IFCD& lhs, const IFCD& rhs) {return &lhs == &rhs;}

}  // namespace fcd
}  // namespace pe
}  // namespace walberla
