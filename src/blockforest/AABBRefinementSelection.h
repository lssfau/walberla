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
#include "BlockForest.h"

#include "core/config/Config.h"
#include "core/math/AABB.h"

#include <vector>



namespace walberla {
namespace blockforest {



class AABBRefinementSelection
{
public:

   AABBRefinementSelection()= default;

   AABBRefinementSelection( const Config::BlockHandle & configBlock )
   {
      if( configBlock )
      {
         auto refinementBlock = configBlock.getKey() == "AABBRefinementSelection" ? configBlock : configBlock.getBlock( "AABBRefinementSelection" );

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
      aabbs_.emplace_back( aabb, level );
   }

   void addRegion( const math::AABB & region, const uint_t level )
   {
      regions_.emplace_back( region, level );
   }

   // for static refinement
   void operator()( SetupBlockForest & forest )
   {
      std::vector< std::pair< math::AABB, uint_t > > aabbs = transformRegionsToAABBs( forest.getDomain() );
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

   // for dynamic refinement
   void operator()(std::vector< std::pair< const Block *, uint_t > > & minTargetLevels,
                   std::vector< const Block * > &, const BlockForest & forest )
   {
      std::vector< std::pair< math::AABB, uint_t > > aabbs = transformRegionsToAABBs( forest.getDomain() );
      aabbs.insert( aabbs.end(), aabbs_.begin(), aabbs_.end() );

      for( auto it = minTargetLevels.begin(); it != minTargetLevels.end(); ++it )
      {
         uint_t currentLevelOfBlock = it->first->getLevel();
         uint_t targetLevelOfBlock = currentLevelOfBlock;

         for( auto aabb = aabbs.begin(); aabb != aabbs.end(); ++aabb )
         {
            if( it->first->getAABB().intersects( aabb->first ) )
            {
               uint_t targetLevelOfAABB = aabb->second;
               if( currentLevelOfBlock > targetLevelOfAABB )
               {
                  targetLevelOfBlock = currentLevelOfBlock - uint_t(1);
               }
               else if ( currentLevelOfBlock < targetLevelOfBlock )
               {
                  targetLevelOfBlock = currentLevelOfBlock + uint_t(1);
               }
               // only the first found intersecting AABB is taken into account
               break;
            }
         }

         WALBERLA_CHECK_LESS_EQUAL(std::abs(int_c(targetLevelOfBlock) - int_c(currentLevelOfBlock)), uint_t(1), "Only level difference of maximum 1 allowed!");
         it->second = targetLevelOfBlock;
      }
   }


private:

   std::vector< std::pair< math::AABB, uint_t > > transformRegionsToAABBs( const math::AABB & simulationDomain ) const
   {
      std::vector< std::pair< math::AABB, uint_t > > aabbs;
      for( auto region = regions_.begin(); region != regions_.end(); ++region )
      {
         aabbs.emplace_back( math::AABB( simulationDomain.xMin() + region->first.xMin() * simulationDomain.xSize(),
                                         simulationDomain.yMin() + region->first.yMin() * simulationDomain.ySize(),
                                         simulationDomain.zMin() + region->first.zMin() * simulationDomain.zSize(),
                                         simulationDomain.xMin() + region->first.xMax() * simulationDomain.xSize(),
                                         simulationDomain.yMin() + region->first.yMax() * simulationDomain.ySize(),
                                         simulationDomain.zMin() + region->first.zMax() * simulationDomain.zSize() ), region->second );
      }
      return aabbs;
   }


   std::vector< std::pair< math::AABB, uint_t > > aabbs_;
   std::vector< std::pair< math::AABB, uint_t > > regions_;

}; // class AABBRefinementSelection



} // namespace blockforest
} // namespace walberla
