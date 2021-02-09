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
//! \file GenericFCD.h
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#pragma once

#include "IFCD.h"

#include "pe/utility/BodyCast.h"

#include "blockforest/BlockDataHandling.h"

namespace walberla{
namespace pe{
namespace fcd {

///
/// \brief Uses CollideFunctor to call collide function without additional namespace inclusion.
///
template <typename BodyTypeTuple, template <typename Container> class CollisionFunctor >
class GenericFCD : public IFCD{
public:
   Contacts& generateContacts(PossibleContacts& possibleContacts) override
   {
      contacts_.clear();
      CollisionFunctor<decltype(contacts_)> func(contacts_);
      for (auto it = possibleContacts.begin(); it != possibleContacts.end(); ++it)
      {
         DoubleCast<BodyTypeTuple, BodyTypeTuple, CollisionFunctor<decltype(contacts_)>, bool>::execute(it->first, it->second, func);
      }
      return contacts_;
   }
};

template <typename BodyTypeTuple, template <typename Container> class CollisionFunctor>
shared_ptr< blockforest::AlwaysCreateBlockDataHandling<GenericFCD<BodyTypeTuple, CollisionFunctor> > > createGenericFCDDataHandling()
{
   return make_shared< blockforest::AlwaysCreateBlockDataHandling<GenericFCD<BodyTypeTuple, CollisionFunctor> > >( );
}

}
}
}
