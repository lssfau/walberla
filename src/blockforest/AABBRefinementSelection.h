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
//! \file Block.h
//! \ingroup blockforest
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "SetupBlockForest.h"

#include "core/config/Config.h"
#include "core/math/AABB.h"

#include <vector>



namespace walberla {
namespace blockforest {



class AABBRefinementSelection
{
public:

   AABBRefinementSelection( const Config::BlockHandle & configBlock )
   {
      if( configBlock )
      {
         auto refinementBlock = configBlock.getBlock( "AABBRefinementSelection" );

         if( refinementBlock )
         {
            Config::Blocks aabbBlocks;
            refinementBlock.getBlocks( "AABB", aabbBlocks );
            for( auto block = aabbBlocks.begin(); block != aabbBlocks.end(); ++block )
            {
               const math::AABB aabb  = block->getParameter< math::AABB >( "AABB" );
               const uint_t     level = block->getParameter< uint_t >( "level" );
               addAABB( aabb, level );
            }

            // For regions, the simulation space is assumed to be [ <0,0,0>, <1,1,1> ]
            Config::Blocks regionBlocks;
            refinementBlock.getBlocks( "Region", regionBlocks );
            for( auto block = regionBlocks.begin(); block != regionBlocks.end(); ++block )
            {
               const math::AABB region = block->getParameter< math::AABB >( "region" );
               const uint_t     level  = block->getParameter< uint_t >( "level" );
               addRegion( region, level );
            }
         }
      }
   }

   void addAABB( const math::AABB & aabb, const uint_t level )
   {
      aabbs_.push_back( std::make_pair( aabb, level ) );
   }

   void addRegion( const math::AABB & region, const uint_t level )
   {
      regions_.push_back( std::make_pair( region, level ) );
   }

   std::vector< std::pair< math::AABB, uint_t > > transformRegionsToAABBs( SetupBlockForest & forest ) const
   {
      std::vector< std::pair< math::AABB, uint_t > > aabbs;
      math::AABB aabb = forest.getDomain();
      for( auto region = regions_.begin(); region != regions_.end(); ++region )
      {
         aabbs.push_back( std::make_pair( math::AABB( aabb.xMin() + region->first.xMin() * aabb.xSize(),
                                                      aabb.yMin() + region->first.yMin() * aabb.ySize(),
                                                      aabb.zMin() + region->first.zMin() * aabb.zSize(),
                                                      aabb.xMin() + region->first.xMax() * aabb.xSize(),
                                                      aabb.yMin() + region->first.yMax() * aabb.ySize(),
                                                      aabb.zMin() + region->first.zMax() * aabb.zSize() ), region->second ) );
      }
      return aabbs;
   }

   void operator()( SetupBlockForest & forest )
   {
      std::vector< std::pair< math::AABB, uint_t > > aabbs = transformRegionsToAABBs( forest );
      aabbs.insert( aabbs.end(), aabbs_.begin(), aabbs_.end() );

      if( aabbs.empty() )
         return;

      for( auto block = forest.begin(); block != forest.end(); ++block )
      {
         for( auto aabb = aabbs.begin(); aabb != aabbs.end(); ++aabb )
         {
            if( block->getAABB().intersects( aabb->first ) && block->getLevel() < aabb->second )
               block->setMarker( true );
         }
      }
   }

private:

   std::vector< std::pair< math::AABB, uint_t > > aabbs_;
   std::vector< std::pair< math::AABB, uint_t > > regions_;

}; // class AABBRefinementSelection



} // namespace blockforest
} // namespace walberla
