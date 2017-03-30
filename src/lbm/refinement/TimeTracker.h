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
//! \file BlockSweepWrapper.h
//! \ingroup lbm
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/DataTypes.h"
#include "core/debug/Debug.h"


namespace walberla {
namespace lbm {



/// Keeps track of the simulation time (One step on the coarsest level counts as '1', one step on the
/// next finer level counts as '0.5', one step on the next level as '0.25' etc.)
class TimeTracker
{
public:

   TimeTracker( const uint_t initialTime = uint_t(0) )
   {
      fractionCounter_.push_back( initialTime );
      fractionBase_.push_back( uint_t(1) );
      time_.push_back( real_c( initialTime ) );
   }

   // for use without refinement time step
   void operator()()
   {
      WALBERLA_ASSERT_GREATER( fractionCounter_.size(), 0 );
      WALBERLA_ASSERT_GREATER( fractionBase_.size(), 0 );
      WALBERLA_ASSERT_GREATER( time_.size(), 0 );
      fractionCounter_[0] += uint_t(1);
      time_[0] = real_c( fractionCounter_[0] );
   }

   // for use with refinement time step (addPost*VoidFunction)
#ifdef NDEBUG
   void operator()( const uint_t level, const uint_t )
#else
   void operator()( const uint_t level, const uint_t executionCount )
#endif
   {
      checkLevel( level );

      WALBERLA_ASSERT_EQUAL( fractionCounter_[ level ] % fractionBase_[ level ], executionCount );
      fractionCounter_[ level ] += uint_t(1);

      time_[ level ] = real_c( fractionCounter_[ level ] ) / real_c( fractionBase_[ level ] );
   }

   real_t getTime( const uint_t level = uint_t(0) )
   {
      checkLevel( level );

      return time_[ level ];
   }

private:

   void checkLevel( const uint_t level )
   {
      if( level >= time_.size() )
      {
         const uint_t previousSize = time_.size();

         fractionCounter_.resize( level + uint_t(1) );
         fractionBase_.resize( level + uint_t(1) );
         time_.resize( level + uint_t(1) );

         for( uint_t i = previousSize; i <= level; ++i )
         {
            fractionCounter_[ i ] = fractionCounter_[0] * math::uintPow2( i );
            fractionBase_[ i ] = math::uintPow2( i );
            time_[ i ] = time_[0];
         }
      }
   }

   std::vector< uint_t > fractionCounter_;
   std::vector< uint_t > fractionBase_;

   std::vector< real_t > time_;

}; // class TimeTracker



} // namespace lbm
} // namespace walberla
