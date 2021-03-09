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
//! \file BodiesForceTorqueContainer.h
//! \ingroup pe_coupling
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//
//======================================================================================================================

#pragma once

#include "blockforest/StructuredBlockForest.h"
#include "BodySelectorFunctions.h"

#include <map>
#include <array>

namespace walberla {
namespace pe_coupling {

class BodiesForceTorqueContainer
{  
public:

   using ForceTorqueStorage_T = std::map<walberla::id_t, std::array<real_t, 6>>;

   BodiesForceTorqueContainer( const shared_ptr<StructuredBlockForest> & blockForest, const BlockDataID & bodyStorageID, const std::function<bool(
            pe::BodyID)> &bodySelectorFct = selectRegularBodies)
   : blockForest_( blockForest ), bodyStorageID_( bodyStorageID ), bodySelectorFct_( bodySelectorFct )
   {
      // has to be added to the forest (not the storage) to register correctly
      bodyForceTorqueStorageID_ = blockForest->addBlockData(make_shared<blockforest::AlwaysCreateBlockDataHandling<ForceTorqueStorage_T> >(), "BodiesForceTorqueContainer");
   }

   void operator()()
   {
      store();
   }

   void store();

   void setOnBodies();

   void clear();

   void swap( BodiesForceTorqueContainer & other );

private:

   shared_ptr<StructuredBlockStorage> blockForest_;
   const BlockDataID bodyStorageID_;
   BlockDataID bodyForceTorqueStorageID_;
   std::function<bool(pe::BodyID)> bodySelectorFct_;
};


class BodyContainerSwapper
{
public:
   BodyContainerSwapper( const shared_ptr<BodiesForceTorqueContainer> & cont1, const shared_ptr<BodiesForceTorqueContainer> & cont2 )
   : cont1_( cont1 ), cont2_( cont2 )
   { }

   void operator()()
   {
      cont1_->swap( *cont2_ );
   }

private:
   shared_ptr<BodiesForceTorqueContainer> cont1_;
   shared_ptr<BodiesForceTorqueContainer> cont2_;
};

} // namespace pe_coupling
} // namespace walberla
