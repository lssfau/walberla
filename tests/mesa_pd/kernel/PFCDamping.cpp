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
//! \file   PFCDamping.cpp
//! \author Igor Ostanin <i.ostanin@skoltech.ru>
//! \author Grigorii Drozdov <drozd013@umn.edu>
//
//======================================================================================================================

#include <mesa_pd/data/DataTypes.h>
#include <mesa_pd/data/ParticleAccessor.h>

#include <mesa_pd/kernel/PFCDamping.h>

#include <core/Environment.h>

namespace walberla {

using namespace walberla::mesa_pd;

int main( int argc, char ** argv )
{
   Environment env(argc, argv);
   WALBERLA_UNUSED(env);
   mpi::MPIManager::instance()->useWorldComm();

   //init data structures
   data::SingleParticleAccessor ac;

   ac.setLinearVelocity( 0, Vec3(+2_r,-2_r,+2_r) );
   ac.setForce         ( 0, Vec3(+2_r,+3_r,-4_r) );

   ac.setAngularVelocity( 0, Vec3(+2_r,-2_r,+2_r) );
   ac.setTorque         ( 0, Vec3(+3_r,+5_r,-2_r) );

   //init kernels
   kernel::PFCDamping damping( 0.1_r );

   damping(0, ac);

   WALBERLA_CHECK_FLOAT_EQUAL(ac.getForce(0),  Vec3(1.8_r, 3.3_r, -4.4_r));
   WALBERLA_CHECK_FLOAT_EQUAL(ac.getTorque(0), Vec3(2.7_r, 5.5_r, -2.2_r));

   return EXIT_SUCCESS;
}

} //namespace walberla

int main( int argc, char ** argv )
{
   return walberla::main(argc, argv);
}
