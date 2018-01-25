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
//! \file RBGSIteration.cpp
//! \ingroup pde
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#include "RBGSIteration.h"
#include "core/debug/CheckFunctions.h"
#include "core/logging/Logging.h"



namespace walberla {
namespace pde {



void RBGSIteration::operator()()
{
   WALBERLA_LOG_PROGRESS_ON_ROOT( "Starting red-black Gauss-Seidel iteration with a maximum number of " << iterations_ << " iterations" );

   iterationsPerformed_ = uint_t(0);
   thresholdReached_ = false;

   uint_t i( uint_t(0) );
   while( i < iterations_ )
   {
      if( boundary_ )
         boundary_();
      communication_();

      for( auto block = blocks_.begin( requiredSelectors_, incompatibleSelectors_ ); block != blocks_.end(); ++block )
         redUpdate_( block.get() );

      if( boundary_ )
         boundary_();
      communication_();

      for( auto block = blocks_.begin( requiredSelectors_, incompatibleSelectors_ ); block != blocks_.end(); ++block )
         blackUpdate_( block.get() );

      if( residualNormThreshold_ > real_t(0) && residualCheckFrequency_ > uint_t(0) )
      {
         if( (i % residualCheckFrequency_) == uint_t(0) )
         {
            if( boundary_ )
               boundary_();
            const real_t residualNorm = residualNorm_();
            WALBERLA_CHECK( math::finite(residualNorm), "Non-finite residual norm detected during the red-black Gauss-Seidel iteration, "
                                                        "the simulation has probably diverged." );
            WALBERLA_LOG_DETAIL_ON_ROOT( "Residual norm after " << (i+1) << " red-black Gauss-Seidel iterations: " << residualNorm );
            if( residualNorm < residualNormThreshold_ )
            {
               WALBERLA_LOG_PROGRESS_ON_ROOT( "Aborting red-black Gauss-Seidel iteration (residual norm threshold reached):"
                                              "\n  residual norm threshold: " << residualNormThreshold_ <<
                                              "\n  residual norm:           " << residualNorm );
               thresholdReached_ = true;
               break;
            }
         }
      }

      ++i;
   }

   iterationsPerformed_ = i;
   WALBERLA_LOG_PROGRESS_ON_ROOT( "Red-black Gauss-Seidel iteration finished after " << i << " iterations" );
}



} // namespace pde
} // namespace walberla
