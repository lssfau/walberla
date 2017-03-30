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
//! \file SetBodyTypeIDs.h
//! \ingroup RigidBody
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/UniqueID.h"
#include "core/logging/Logging.h"

#include <boost/tuple/tuple.hpp>

namespace walberla {
namespace pe {

template < typename BodyTypeTuple >
struct SetBodyTypeIDs{
   static void execute(){
      auto typeID = UniqueID<SetBodyTypeIDs<int> >::createGlobal();
      BodyTypeTuple::head_type::setStaticTypeID( typeID );
      WALBERLA_LOG_DETAIL_ON_ROOT("SetBodyTypeID " << typeID << " set.");
      SetBodyTypeIDs<typename BodyTypeTuple::tail_type>::execute();
   }
};

template < >
struct SetBodyTypeIDs< boost::tuples::null_type>{
   static void execute(){
   }
};

}
}
