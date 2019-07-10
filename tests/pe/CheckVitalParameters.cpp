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
//! \file CheckVitalParameters.cpp
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include "CheckVitalParameters.h"
#include "pe/Materials.h"

#include "core/Environment.h"
#include "core/debug/TestSubsystem.h"
#include "core/UniqueID.h"

using namespace walberla::pe;

int main( int argc, char ** argv )
{
   walberla::debug::enterTestMode();
   walberla::Environment env(argc, argv);
   WALBERLA_UNUSED(env);

   MaterialID iron = Material::find("iron");
   Sphere sphere(walberla::UniqueID<RigidBody>::create(), 0, Vec3(15, 15, 15), Quat(), 3, iron, false, true, false );
   sphere.MPITrait.setOwner( Owner(0, walberla::blockforest::BlockID().getID() ) );
   checkVitalParameters( &sphere, &sphere );

   return EXIT_SUCCESS;
}
