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
//! \file InitializationManager.h
//! \ingroup geometry
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/DataTypes.h"
#include "core/config/Config.h"

#include <map>
#include <string>
#include <vector>


namespace walberla {

namespace domain_decomposition {
class BlockStorage;
}
namespace uid {
template< typename T > class UID;
}

namespace geometry {
namespace initializer {


   class Initializer;
   class InitializerUIDGenerator;
   using InitializerUID = uid::UID<InitializerUIDGenerator>;


   //*******************************************************************************************************************
   /*! Manages domain initialization from configuration file
   *
   * The class parses one configuration block. The sub-block names have to match the
   * InitializerUID's. The sub-block is then processed by the Initializer which was registered for this UID.
   *
   */
   //*******************************************************************************************************************
   class InitializationManager
   {
   public:
      InitializationManager( domain_decomposition::BlockStorage & blockStorage );

      void registerInitializer( const InitializerUID & uid, const shared_ptr<Initializer> & geometry );

      void init( const Config::BlockHandle & blockHandle );

   protected:
      domain_decomposition::BlockStorage & blockStorage_;
      std::map< InitializerUID, shared_ptr<Initializer> > geometryRegistry_;
   };



} // namespace initializer
} // namespace geometry
} // namespace walberla
