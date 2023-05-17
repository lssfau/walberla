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
//! \file FieldTransferTest.h
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================


#include "core/debug/TestSubsystem.h"
#include "core/Environment.h"

#include "field/Field.h"

#include "gpu/GPUField.h"
#include "gpu/FieldCopy.h"
#include "core/math/Random.h"


using namespace walberla;

void simpleTransfer()
{
   Field<double, 4> h_f1( 16, 20, 30, 42.0, field::fzyx );
   Field<double, 4> h_f2( 16, 20, 30, 0.0, field::fzyx );

   WALBERLA_FOR_ALL_CELLS_XYZ(&h_f1,
      h_f1(x, y, z, 0) = math::realRandom<double>();
   )

   gpu::GPUField<double> d_f( 16, 20, 30, 4, 0, field::fzyx );

   WALBERLA_CHECK_EQUAL( h_f1.xSize(), d_f.xSize())
   WALBERLA_CHECK_EQUAL( h_f1.ySize(), d_f.ySize())
   WALBERLA_CHECK_EQUAL( h_f1.zSize(), d_f.zSize())
   WALBERLA_CHECK_EQUAL( h_f1.fSize(), d_f.fSize())
   WALBERLA_CHECK_EQUAL( h_f1.layout(), d_f.layout())

   gpu::fieldCpy( d_f, h_f1 );
   gpu::fieldCpy( h_f2, d_f );

   WALBERLA_CHECK_EQUAL( h_f1, h_f2 )
}


int main( int argc, char **argv )
{
   debug::enterTestMode();
   walberla::Environment walberlaEnv( argc, argv );

   simpleTransfer();

   return EXIT_SUCCESS;
}
