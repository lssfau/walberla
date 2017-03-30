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
//! \file BlockForestEvaluation.cpp
//! \ingroup blockforest
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#include "BlockForestEvaluation.h"
#include "core/math/DistributedSample.h"
#include "core/math/FPClassify.h"
#include "core/mpi/MPIManager.h"

#include <cmath>


namespace walberla {
namespace blockforest {

namespace internal {
std::string boolToText( const bool b )
{
   std::string result( b ? "yes" : "no" );
   return result;
}
}



BlockForestEvaluation::BlockForestEvaluation( const BlockForest & forest ) :
   forest_( forest ), blockStatistics_( forest.getNumberOfLevels() + uint_t(1) )
{
   for( uint_t i = uint_t(0); i <= forest.getNumberOfLevels(); ++i )
   {
      math::DistributedSample sample;

      const uint_t blocks = (i == forest.getNumberOfLevels()) ? forest.getNumberOfBlocks() : forest.getNumberOfBlocks(i);

      sample.insert( real_c( blocks ) );
      sample.mpiGatherRoot();

      WALBERLA_ROOT_SECTION()
      {
         blockStatistics_[i].sum = sample.sum();
         blockStatistics_[i].min = sample.min();
         blockStatistics_[i].max = sample.max();
         blockStatistics_[i].avg = sample.mean();
         blockStatistics_[i].stdDev = sample.stdDeviation();
         const auto relStdDev = sample.relativeStdDeviation();
         blockStatistics_[i].relStdDev = math::isnan( relStdDev ) ? real_t(0) : relStdDev;
      }
   }
}



void BlockForestEvaluation::toStream( std::ostream & os ) const
{
   WALBERLA_NON_ROOT_SECTION()
   {
      os << "WARNING/ERROR: information about the global block forest structure are only available on the root process!";
   }
   WALBERLA_ROOT_SECTION()
   {
      os << "- AABB:             " << forest_.getDomain() << "\n"
         << "- forest size:      " << forest_.getXSize() << " x " << forest_.getYSize() << " x " << forest_.getZSize() << " blocks (number of octree root blocks)\n"
         << "- periodicity:      " << ( forest_.isXPeriodic() ? "true" : "false" ) << " x " << ( forest_.isYPeriodic() ? "true" : "false" ) << " x " << ( forest_.isZPeriodic() ? "true" : "false" ) << "\n"
         << "- number of levels: " << forest_.getNumberOfLevels() << " (max. number of levels: ";

      if( forest_.limitedLevels() )
         os << forest_.getMaxLevels() << ")\n";
      else
         os << "unlimited)\n";

      os << "- tree ID digits:   " << forest_.getTreeIdDigits() << " (-> block ID bytes = " << forest_.getBlockIdBytes() << ")\n"
         << "- total number of blocks: " << uint_c( blockStatistics_.back().sum + real_c(0.5) );

      if( mpi::MPIManager::instance()->numProcesses() > 1 )
      {
         os << "\n- blocks per process:"
            << "\n   + min       = " << blockStatistics_.back().min
            << "\n   + max       = " << blockStatistics_.back().max
            << "\n   + avg       = " << blockStatistics_.back().avg
            << "\n   + stdDev    = " << blockStatistics_.back().stdDev
            << "\n   + relStdDev = " << blockStatistics_.back().relStdDev;
      }

      if( forest_.getNumberOfLevels() > uint_t(1) )
      {
         os << "\n- distribution to different levels:";
         for( uint_t l = uint_t(0); l != forest_.getNumberOfLevels(); ++l )
         {
            os << "\n   + level " << l
               << "\n      - " << uint_c( blockStatistics_[l].sum + real_c(0.5) ) << " blocks ..."
               << "\n      - ... cover " << space(l) << " % of the total simulation space"
               << "\n      - ... account for " << ( real_c(100) * blockStatistics_[l].sum / blockStatistics_.back().sum ) << " % of all blocks";
            if( mpi::MPIManager::instance()->numProcesses() > 1 )
            {
               os << "\n      - blocks on this level per process:"
                  << "\n         + min       = " << blockStatistics_[l].min
                  << "\n         + max       = " << blockStatistics_[l].max
                  << "\n         + avg       = " << blockStatistics_[l].avg
                  << "\n         + stdDev    = " << blockStatistics_[l].stdDev
                  << "\n         + relStdDev = " << blockStatistics_[l].relStdDev;
            }
         }
      }

      auto blockData = forest_.getBlockDataIdentifiers();
      if( !blockData.empty() )
      {
         os << "\n- block data:";
         for( auto data = blockData.begin(); data != blockData.end(); ++data )
            os << "\n   + " << *data;
      }

      os << "\n- number of already performed restructure cycles: " << forest_.getModificationStamp()
         << "\n- allow multiple restructure cycles per refresh cycle: " << internal::boolToText( forest_.allowMultipleRefreshCycles() )
         << "\n- recalculate block levels in refresh: " << internal::boolToText( forest_.recalculateBlockLevelsInRefresh() );
      if( forest_.recalculateBlockLevelsInRefresh() )
      {
         os << "\n   + allow for changing maximum level/depth: " << internal::boolToText( forest_.allowRefreshChangingDepth() ) << " (requires allreduce of one unsigned integer)"
            << "\n   + check for early out during recalculation of block levels: " << internal::boolToText( forest_.checkForEarlyOutInRefresh() ) << " (requires allreduce of one boolean)"
            << "\n   + check for late out during recalculation of block levels: " <<  internal::boolToText( forest_.checkForLateOutInRefresh() ) << " (requires allreduce of one boolean)";
      }
      os << "\n- always perform redistribution/load balancing in refresh: " << internal::boolToText( forest_.alwaysRebalanceInRefresh() )
         << "\n- redistribution/load balancing function registered: " << internal::boolToText( forest_.loadBalancingFunctionRegistered() )
         << "\n- check for early out after load balancing: " << internal::boolToText( forest_.checkForEarlyOutAfterLoadBalancing() );
   }
}



} // namespace blockforest
} // namespace walberla
