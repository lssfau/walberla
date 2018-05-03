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
//! \file StaticCurve.cpp
//! \ingroup blockforest
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#include "StaticCurve.h"

#include <cmath>



namespace walberla {
namespace blockforest {



uint_t StaticLevelwiseCurveBalance::operator()( SetupBlockForest & forest, const uint_t numberOfProcesses, const memory_t /*perProcessMemoryLimit*/ )
{
   // TODO: take per process memory limit into account?

   std::vector< SetupBlock * > blocks;
   if( hilbert_ )
      forest.getHilbertOrder( blocks );
   else
      forest.getMortonOrder( blocks );

   uint_t border = uint_t(0);

   for( uint_t level = forest.getNumberOfLevels(); level-- > uint_t(0); )
   {
      std::vector< SetupBlock * > blocksOnLevel;

      for( auto block = blocks.begin(); block != blocks.end(); ++block )
         if( (*block)->getLevel() == level )
            blocksOnLevel.push_back( *block );

      const uint_t nBlocks = blocksOnLevel.size();

      if( nBlocks <= ( numberOfProcesses - border ) )
      {
         for( auto block = blocksOnLevel.begin(); block != blocksOnLevel.end(); ++block )
            (*block)->assignTargetProcess( border++ );

         WALBERLA_ASSERT_LESS_EQUAL( border, numberOfProcesses );

         if( border == numberOfProcesses )
            border = uint_t(0);
      }
      else
      {
         const uint_t reducedNBlocks = nBlocks - ( numberOfProcesses - border);
         const uint_t div = reducedNBlocks / numberOfProcesses;
         const uint_t mod = reducedNBlocks % numberOfProcesses;

         uint_t bIndex = uint_t(0);
         for( uint_t p = 0; p != numberOfProcesses; ++p )
         {
            uint_t count = div;
            if( p < mod ) ++count;
            if( p >= border ) ++count;

            WALBERLA_ASSERT_LESS_EQUAL( bIndex + count, blocksOnLevel.size() );

            for( uint_t i = bIndex; i < ( bIndex + count ); ++i )
               blocksOnLevel[i]->assignTargetProcess( p );
            bIndex += count;
         }

         border = mod;
      }
   }

   return std::min( numberOfProcesses, blocks.size() );
}



uint_t StaticLevelwiseCurveBalanceWeighted::operator()( SetupBlockForest & forest, const uint_t numberOfProcesses, const memory_t /*perProcessMemoryLimit*/ )
{
   // TODO: take per process memory limit into account?

   std::vector< SetupBlock * > blocks;
   if( hilbert_ )
      forest.getHilbertOrder( blocks );
   else
      forest.getMortonOrder( blocks );

   uint_t usedProcesses( uint_t(0) );

   for( uint_t level = uint_t(0); level < forest.getNumberOfLevels(); ++level )
   {
      std::vector< SetupBlock * > blocksOnLevel;

      for( auto block = blocks.begin(); block != blocks.end(); ++block )
         if( (*block)->getLevel() == level )
            blocksOnLevel.push_back( *block );

      workload_t totalWeight( 0 );
      for( auto block = blocksOnLevel.begin(); block != blocksOnLevel.end(); ++block )
      {
         WALBERLA_ASSERT( !( (*block)->getWorkload() < workload_t(0) ) );
         totalWeight += (*block)->getWorkload();
      }

      uint_t c( uint_t(0) );
      for( uint_t p = uint_t(0); p != numberOfProcesses; ++p )
      {
         const workload_t pWeight = totalWeight / workload_c( numberOfProcesses - p );
         workload_t weight( 0 );
         while( c < blocksOnLevel.size() && ( isIdentical(weight, workload_t(0)) ||
                std::abs( pWeight - weight - blocksOnLevel[c]->getWorkload() ) <=
                std::abs( pWeight - weight ) ) )
         {
            blocksOnLevel[c]->assignTargetProcess(p);

            WALBERLA_ASSERT_LESS_EQUAL( p, usedProcesses );
            usedProcesses = p + uint_t(1);

            const workload_t addedWeight = blocksOnLevel[c]->getWorkload();
            weight += addedWeight;
            totalWeight -= addedWeight;
            ++c;
         }
      }
      while( c < blocksOnLevel.size() )
      {
         blocksOnLevel[c]->assignTargetProcess( numberOfProcesses - uint_t(1) );

         WALBERLA_ASSERT_LESS_EQUAL( numberOfProcesses - uint_t(1), usedProcesses );
         usedProcesses = numberOfProcesses;

         ++c;
      }


   }

   return usedProcesses;
}



} // namespace blockforest
} // namespace walberla
