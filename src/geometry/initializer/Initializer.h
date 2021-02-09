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
//! \file Initializer.h
//! \ingroup geometry
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/DataTypes.h"
#include "core/config/Config.h"
#include "core/uid/UID.h"
#include "core/uid/UIDGenerators.h"

#include <vector>


namespace walberla {

namespace domain_decomposition {
class BlockStorage;
}

namespace geometry {
namespace initializer {

   class InitializerUIDGenerator : public uid::IndexGenerator< InitializerUIDGenerator, uint_t > { };



   //*******************************************************************************************************************
   /*! Abstract base class for all Initializers
   *
   * An initializer takes one configuration block, parses it, and sets up the domain accordingly
   */
   //*******************************************************************************************************************
   class Initializer
   {
   public:
      virtual ~Initializer() = default;

      virtual void init( domain_decomposition::BlockStorage & blockStorage, const Config::BlockHandle & blockHandle ) = 0;
   };




} // namespace initializer
} // namespace geometry
} // namespace walberla
