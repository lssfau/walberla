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
//! \file BlockForestEvaluation.h
//! \ingroup blockforest
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "BlockForest.h"
#include "core/logging/Logging.h"


namespace walberla {
namespace blockforest {



class BlockForestEvaluation
{
public:

   struct BlockStatistics
   {
      real_t sum;
      real_t min;
      real_t max;
      real_t avg;
      real_t stdDev;
      real_t relStdDev;
   };

   BlockForestEvaluation( const BlockForest & forest );

   void toStream( std::ostream & os ) const;

   inline std::string toString() const
   {
      std::ostringstream oss;
      toStream( oss );
      return oss.str();
   }

private:

   real_t space( const uint_t level ) const
   {
      return ( real_c(100) * blockStatistics_[level].sum ) /
             ( real_c( forest_.getXSize() * forest_.getYSize() * forest_.getZSize() ) * real_c( math::uintPow8(level) ) );
   }

   const BlockForest & forest_;

   std::vector< BlockStatistics> blockStatistics_;

}; // class BlockForestEvaluation



inline void logDuringRefresh( BlockForest & forest, const PhantomBlockForest & )
{
   BlockForestEvaluation evaluation( forest );
   WALBERLA_LOG_INFO_ON_ROOT( "BlockForest refresh: Current state of the block data structure:\n" << evaluation.toString() );
}



//////////////////////
// Global Functions //
//////////////////////



inline std::ostream & operator<<( std::ostream & os, const BlockForest & forest )
{
   BlockForestEvaluation evaluation( forest );
   evaluation.toStream( os );
   return os;
}



} // namespace blockforest
} // namespace walberla
