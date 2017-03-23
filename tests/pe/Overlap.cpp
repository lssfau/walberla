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
//! \file Overlap.cpp
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include "core/debug/TestSubsystem.h"

#include "pe/utility/Overlap.h"

#include "core/math/Constants.h"

using namespace walberla;
using namespace walberla::pe;

int main( int /*argc*/, char** /*argv*/ )
{
    walberla::debug::enterTestMode();

    WALBERLA_CHECK_FLOAT_EQUAL(getSphereSphereOverlap(real_c(3), real_c(1), real_c(1)), real_c(0));
    WALBERLA_CHECK_FLOAT_EQUAL(getSphereSphereOverlap(real_c(1.1), real_c(1), real_c(1)), real_c(1.08149327099828607));

    return EXIT_SUCCESS;
}
