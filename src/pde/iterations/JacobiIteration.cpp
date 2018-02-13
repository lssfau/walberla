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
//! \file JacobiIteration.cpp
//! \ingroup pde
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#include "JacobiIteration.h"
#include "core/debug/CheckFunctions.h"
#include "core/logging/Logging.h"



namespace walberla {
namespace pde {



void JacobiIteration::operator()()
{
   WALBERLA_LOG_PROGRESS_ON_ROOT( "Starting Jacobi iteration with a maximum number of " << iterations_ << " iterations" );

   uint_t i( uint_t(0) );
   while( i < iterations_ )
   {
      if( boundary_ )
         boundary_();
      communication_();

      for( auto block = blocks_.begin( requiredSelectors_, incompatibleSelectors_ ); block != blocks_.end(); ++block )
         jacobi_( block.get() );

      if( residualNormThreshold_ > real_t(0) && residualCheckFrequency_ > uint_t(0) )
      {
         if( (i % residualCheckFrequency_) == uint_t(0) )
         {
            if( boundary_ )
               boundary_();
            const real_t residualNorm = residualNorm_();
            WALBERLA_CHECK( math::finite(residualNorm), "Non-finite residual norm detected during the Jacobi iteration, "
                                                        "the simulation has probably diverged." );
            WALBERLA_LOG_DETAIL_ON_ROOT( "Residual norm after " << (i+1) << " Jacobi iterations: " << residualNorm );
            if( residualNorm < residualNormThreshold_ )
            {
               WALBERLA_LOG_PROGRESS_ON_ROOT( "Aborting Jacobi iteration (residual norm threshold reached):"
                                              "\n  residual norm threshold: " << residualNormThreshold_ <<
                                              "\n  residual norm:           " << residualNorm );
               break;
            }
         }
      }

      ++i;
   }

   WALBERLA_LOG_PROGRESS_ON_ROOT( "Jacobi iteration finished after " << i << " iterations" );
}



} // namespace pde
} // namespace walberla
