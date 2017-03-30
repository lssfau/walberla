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
//! \file InitializationManager.cpp
//! \ingroup geometry
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//
//======================================================================================================================

#include "InitializationManager.h"
#include "Initializer.h"
#include "core/Abort.h"
#include "domain_decomposition/BlockStorage.h"
#include "domain_decomposition/IBlock.h"


namespace walberla {
namespace geometry {
namespace initializer {


InitializationManager::InitializationManager( BlockStorage & blockStorage ) : blockStorage_(blockStorage)
{
}

void InitializationManager::registerInitializer( const InitializerUID & uid, const shared_ptr<Initializer> & geometry )
{
   auto result = geometryRegistry_.insert( std::make_pair(uid, geometry) );
   if( !result.second )
      WALBERLA_ABORT("You are trying to register GeometryFactory '" << uid.getIdentifier() << "', but a GeometryFactory with that"
                     " name has already been registered!");
}

void InitializationManager::init( const Config::BlockHandle & blockHandle )
{
   Config::Blocks blocks;
   blockHandle.getBlocks(blocks);

   for( auto blockIt = blocks.begin(); blockIt != blocks.end(); ++blockIt)
   {
      auto foundGeometryIt = geometryRegistry_.find( blockIt->getKey() );

      if( foundGeometryIt == geometryRegistry_.end() )
      {
         std::ostringstream oss;
         oss << "No Geometry for block " << blockIt->getKey() << " registered at GeometryMaster!\n"
             << geometryRegistry_.size() << " Geometry instances are registered" << ( geometryRegistry_.empty() ? "." : ":" );

         for( auto geometryIt = geometryRegistry_.begin(); geometryIt != geometryRegistry_.end(); ++geometryIt )
            oss <<  "\n   " << geometryIt->first.getIdentifier();

         WALBERLA_ABORT( oss.str() );
      }

      foundGeometryIt->second->init( blockStorage_, *blockIt );
   }
}


} // namespace initializer
} // namespace geometry
} // namespace walberla
