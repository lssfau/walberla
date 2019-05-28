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
//! \file   HeatConduction.cpp
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include <mesa_pd/data/ParticleAccessor.h>
#include <mesa_pd/data/ParticleStorage.h>
#include <mesa_pd/kernel/HeatConduction.h>

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

   if (std::is_same<real_t, float>::value)
   {
      WALBERLA_LOG_WARNING("waLBerla build in sp mode: skipping test due to low precision");
      return EXIT_SUCCESS;
   }

   //init data structures
   auto ps = std::make_shared<data::ParticleStorage>(100);

   data::ParticleAccessor ac(ps);

   data::Particle&& p1 = *ps->create();
   p1.setTemperature(10);
   p1.setType( 0 );

   data::Particle&& p2 = *ps->create();
   p2.setTemperature(20);
   p2.setType( 0 );

   // Init kernels
   kernel::HeatConduction heatConduction(1);
   heatConduction.setConductance(0, 0, real_t(0.2));

   // single contact test
   heatConduction(0,
                  1,
                  ac);

   WALBERLA_CHECK_FLOAT_EQUAL( ps->getHeatFlux(0), -ps->getHeatFlux(1) );
   WALBERLA_CHECK_FLOAT_EQUAL( ps->getHeatFlux(0), real_t(2) );

   return EXIT_SUCCESS;
}

} //namespace walberla

int main( int argc, char ** argv )
{
   return walberla::main(argc, argv);
}
