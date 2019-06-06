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
//! \file   TemperatureIntegration.cpp
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include <mesa_pd/data/DataTypes.h>
#include <mesa_pd/data/ParticleAccessor.h>

#include <mesa_pd/kernel/TemperatureIntegration.h>

#include <core/Environment.h>
#include <core/logging/Logging.h>

#include <iostream>

namespace walberla {

using namespace walberla::mesa_pd;

int main( int argc, char ** argv )
{
   Environment env(argc, argv);
   WALBERLA_UNUSED(env);
   mpi::MPIManager::instance()->useWorldComm();

   //init data structures
   data::SingleParticleAccessor accessor;

   //initialize particle
   accessor.setType(        0, 0 );
   accessor.setHeatFlux(    0, real_t(8) );
   accessor.setTemperature( 0, real_t(5) );

   //init kernels
   const real_t dt = real_t(1);
   kernel::TemperatureIntegration integrator( dt, 1 );
   integrator.setInvHeatCapacity( 0, real_t(2) );

   integrator(0, accessor);

   //check force
   WALBERLA_CHECK_FLOAT_EQUAL(accessor.getHeatFlux(0) + real_t(1), real_t(1));

   //check velocity
   WALBERLA_CHECK_FLOAT_EQUAL(accessor.getTemperature(0), real_t(21));

   return EXIT_SUCCESS;
}

} //namespace walberla

int main( int argc, char ** argv )
{
   return walberla::main(argc, argv);
}
