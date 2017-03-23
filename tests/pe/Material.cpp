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
//! \file Material.cpp
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include "core/debug/TestSubsystem.h"

#include "pe/Materials.h"

#include <algorithm>

using namespace walberla;
using namespace walberla::pe;

int main( int /*argc*/, char** /*argv*/ )
{
   walberla::debug::enterTestMode();

   createMaterial("test", real_c(0.2), real_c(0.2), real_c(0.2), real_c(0.2), real_c(0.2), real_c(0.2), real_c(0.2), real_c(0.2), real_c(0.2));

   MaterialID mat2 = Material::find("test");
   WALBERLA_CHECK_UNEQUAL(mat2, invalid_material);
   WALBERLA_CHECK_FLOAT_EQUAL(Material::getDensity(mat2), real_c(0.2));
   WALBERLA_CHECK_FLOAT_EQUAL(Material::getPoissonRatio(mat2), real_c(0.2));
   WALBERLA_CHECK_FLOAT_EQUAL(Material::getStiffness(mat2), real_c(0.2));
   WALBERLA_CHECK_FLOAT_EQUAL(Material::getRestitution(mat2), real_c(0.2));

   MaterialID mat3 = Material::find("test2");
   WALBERLA_CHECK_EQUAL(mat3, invalid_material);

   return EXIT_SUCCESS;
}
