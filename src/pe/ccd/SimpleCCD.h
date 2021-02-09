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

#include "ICCD.h"

namespace walberla{
namespace pe{
namespace ccd {

class SimpleCCD : public ICCD{
public:
   explicit SimpleCCD(BodyStorage& globalStorage, Storage& storage);
   ~SimpleCCD() override;

   PossibleContacts& generatePossibleContacts( WcTimingTree* tt = nullptr ) override;

   int getObservedBodyCount() const override;
private:
   //**Add/remove functions*********************************************************************
   /*!\name Add/remove functions */
   //@{
   void add   ( BodyID body );
   void remove( BodyID body );
   //@}
   //*******************************************************************************************

   std::string identifier_;

   BodyStorage& globalStorage_;
   Storage& storage_;

   std::vector<BodyID> bodies_;
};

}  // namespace ccd
}  // namespace pe
}  // namespace walberla
