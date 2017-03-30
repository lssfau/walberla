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
//! \file PerformanceLogger.h
//! \ingroup lbm
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//
//======================================================================================================================

#pragma once

#include "PerformanceEvaluation.h"

#include "core/DataTypes.h"
#include "core/Set.h"
#include "core/uid/SUID.h"
#include "core/timing/Timer.h"

#include "domain_decomposition/BlockDataID.h"
#include "domain_decomposition/StructuredBlockStorage.h"

#include "field/FlagUID.h"

#include <limits>

namespace walberla {
namespace lbm {

//**********************************************************************************************************************
/*!
*   \brief Class for using the PerformanceEvaluation in a timeloop
*
*   Providing a measurement interval, this class will regularly (every <interval> time steps) measure and report the
*   LBM performance. At the end of the simulation logOverallResults() may be called to output minimum, maximum and
*   average performance during the simulation run.   
*/
//**********************************************************************************************************************
template< typename FlagField_T >
class PerformanceLogger
{
public:
   PerformanceLogger( const shared_ptr< StructuredBlockStorage > & blocks,
                      const ConstBlockDataID & flagFieldId, const Set< FlagUID > & fluid,
                      const uint_t interval,
                      const Set<SUID> & requiredSelectors = Set<SUID>::emptySet(),
                      const Set<SUID> & incompatibleSelectors = Set<SUID>::emptySet() );

   void operator()();
   void logOverallResultsOnRoot() const;

private:
   enum Mode { MIN, MAX, AVG, LAST };
   real_t getTiming( Mode mode ) const;

   PerformanceEvaluation<FlagField_T> performanceEvaluation_;
   uint_t interval_;
   uint_t timestep_;
   WcTimer timer_;

}; // class PerformanceLogger


template< typename FlagField_T >
PerformanceLogger<FlagField_T>::PerformanceLogger( const shared_ptr< StructuredBlockStorage > & blocks,
                   const ConstBlockDataID & flagFieldId, const Set< FlagUID > & fluid,
                   const uint_t interval,
                   const Set<SUID> & requiredSelectors /*= Set<SUID>::emptySet()*/,
                   const Set<SUID> & incompatibleSelectors /*= Set<SUID>::emptySet()*/ )
                   : performanceEvaluation_( blocks, flagFieldId, fluid, requiredSelectors, incompatibleSelectors ),
                     interval_(interval), timestep_(1)
{
}


template< typename FlagField_T >
real_t PerformanceLogger<FlagField_T>::getTiming( Mode mode ) const
{
   switch( mode )
   {
   case MIN:
      return timer_.max();
   case MAX:
      return timer_.min();
   case AVG:
      return timer_.average();
   case LAST:
      return timer_.last();
   default:
      WALBERLA_ASSERT( false );
      return std::numeric_limits< real_t >::signaling_NaN();
   }
}

template< typename FlagField_T >
void PerformanceLogger<FlagField_T>::operator()()
{
   if( timestep_ == interval_ )
   {
      WALBERLA_MPI_BARRIER();
      timer_.end();
      performanceEvaluation_.logResultOnRoot( interval_, getTiming( LAST ) );
      timer_.start();
      timestep_ = 0;
   }

   ++timestep_;
}

template< typename FlagField_T >
void PerformanceLogger<FlagField_T>::logOverallResultsOnRoot() const
{
   WALBERLA_LOG_RESULT_ON_ROOT( "Min Performance:\n" << performanceEvaluation_.loggingString( interval_, getTiming( MIN ) ); );
   WALBERLA_LOG_RESULT_ON_ROOT( "Max Performance:\n" << performanceEvaluation_.loggingString( interval_, getTiming( MAX ) ); );
   WALBERLA_LOG_RESULT_ON_ROOT( "Avg Performance:\n" << performanceEvaluation_.loggingString( interval_, getTiming( AVG ) ); );
}


} // namespace lbm
} // namespace walberla
