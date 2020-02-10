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

#include <cstdlib>

#include "SimpleCCD.h"

#include "pe/rigidbody/BodyStorage.h"

#include "core/logging/Logging.h"

#include <vector>
#include <utility>

namespace walberla {
namespace pe {
namespace ccd {

SimpleCCD::SimpleCCD(BodyStorage& globalStorage, Storage& storage) : globalStorage_(globalStorage), storage_(storage)
{
   storage_[0].registerAddCallback( "SimpleCCD", std::bind(&SimpleCCD::add, this, std::placeholders::_1) );
   storage_[0].registerRemoveCallback( "SimpleCCD", std::bind(&SimpleCCD::remove, this, std::placeholders::_1) );

   storage_[1].registerAddCallback( "SimpleCCD", std::bind(&SimpleCCD::add, this, std::placeholders::_1) );
   storage_[1].registerRemoveCallback( "SimpleCCD", std::bind(&SimpleCCD::remove, this, std::placeholders::_1) );
}

SimpleCCD::~SimpleCCD()
{
   storage_[0].deregisterAddCallback( "SimpleCCD" );
   storage_[0].deregisterRemoveCallback( "SimpleCCD" );

   storage_[1].deregisterAddCallback( "SimpleCCD" );
   storage_[1].deregisterRemoveCallback( "SimpleCCD" );
}


PossibleContacts& SimpleCCD::generatePossibleContacts( WcTimingTree* tt ){
   contacts_.clear();

   if (tt != nullptr) tt->start("SimpleCCD");
   for (auto it1 = bodies_.begin(); it1 != bodies_.end(); ++it1){
      for (auto it2 = it1 + 1; it2 !=bodies_.end(); ++it2)
      {
         if (!((*it1)->hasInfiniteMass() && (*it2)->hasInfiniteMass()))
         {
            if ( (*it1)->getSystemID() > (*it2)->getSystemID() )
               contacts_.push_back(std::make_pair(*it2, *it1));
            else
               contacts_.push_back(std::make_pair(*it1, *it2));
         }
      }

      for (auto it2 = globalStorage_.begin(); it2 != globalStorage_.end(); ++it2)
      {
         if (!((*it1)->hasInfiniteMass() && it2->hasInfiniteMass()))
         {
            if ( (*it1)->getSystemID() > it2->getSystemID() )
               contacts_.push_back(std::make_pair(it2.getBodyID(), *it1));
            else
               contacts_.push_back(std::make_pair(*it1, it2.getBodyID()));
         }
      }
   }
   if (tt != nullptr) tt->stop("SimpleCCD");

   return contacts_;
}

int SimpleCCD::getObservedBodyCount() const
{
   return static_cast<int> (globalStorage_.size() + bodies_.size());
}

void SimpleCCD::add   ( BodyID body )
{
   bodies_.push_back( body );
}

void SimpleCCD::remove( BodyID body )
{
   WALBERLA_LOG_DETAIL( "Removing body " << body->getSystemID() << " from CCD." );
   auto bodyIt = std::find(bodies_.begin(), bodies_.end(), body);
   if (bodyIt != bodies_.end())
   {
      std::swap(*bodyIt, bodies_.back());
      bodies_.pop_back();
   }
}

}  // namespace ccd
}  // namespace pe
}  // namespace walberla
