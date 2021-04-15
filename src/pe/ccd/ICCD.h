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
//! \file ICCD.h
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#pragma once

#include "pe/Types.h"
#include "core/NonCopyable.h"
#include "core/timing/TimingTree.h"

namespace walberla{
namespace pe{

class BodyStorage;

namespace ccd {

class ICCD : private NonCopyable {
public:
   virtual ~ICCD() = default;

   /// Generates a list of possible contact pairs.
   /// This list is also stored in the member variable contacts_ for reuse lateron.
   virtual PossibleContacts& generatePossibleContacts( WcTimingTree* tt = nullptr ) = 0;
   PossibleContacts& getPossibleContacts() {return contacts_;}

   virtual void reloadBodies() {}
   virtual int  getObservedBodyCount() const = 0;
protected:
   PossibleContacts contacts_;
};

//*************************************************************************************************
/*!\brief Compare if two coarse collision detectors are equal.
 *
 * Since collision detectors are uncopyable two collision detectors are considered equal if their addresses are equal.
 */
inline bool operator==(const ICCD& lhs, const ICCD& rhs) {return &lhs == &rhs;}

}  // namespace ccd
}  // namespace pe
}  // namespace walberla
