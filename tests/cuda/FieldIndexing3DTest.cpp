
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
//! \file
//! \author Paulo Carvalho <prcjunior@inf.ufpr.br>
//
//======================================================================================================================


#include "core/debug/TestSubsystem.h"
#include "core/Environment.h"
#include "core/mpi/Datatype.h"

#include "field/GhostLayerField.h"

#include "cuda/GPUField.h"
#include "cuda/FieldCopy.h"
#include "cuda/Kernel.h"
#include "cuda/FieldIndexing3D.h"

#include "FieldIndexing3DTest.h"

using namespace walberla;


typedef cuda::FieldIndexing3D<int> FieldIdx3D_T;
typedef GhostLayerField<int , F_SIZE> HostField_T;
typedef cuda::GPUField<int> GPUField_T;



void xyzTest()
{
   const HostField_T emptyField( X_SIZE, Y_SIZE, Z_SIZE, GL_SIZE, -1, LAYOUT );
   GPUField_T deviceField( X_SIZE, Y_SIZE, Z_SIZE, F_SIZE, 1, LAYOUT );
   cuda::fieldCpy( deviceField, emptyField );

   auto setValue = cuda::make_kernel( &setValueKernel );
   setValue.addFieldIndexingParam( FieldIdx3D_T::xyz( deviceField ) );
   setValue();

   HostField_T resultField( X_SIZE, Y_SIZE, Z_SIZE, GL_SIZE, -1, LAYOUT );
   cuda::fieldCpy( resultField, deviceField );

   HostField_T expectedField( X_SIZE, Y_SIZE, Z_SIZE, GL_SIZE, -1, LAYOUT );
   WALBERLA_FOR_ALL_CELLS_XYZ( &expectedField,
      for ( uint_t f = 0; f < expectedField.fSize(); ++f )
      {
         expectedField.get( x, y, z, f ) = IDX4D( x, y, z, f );
      }
   )

   WALBERLA_ASSERT( resultField == expectedField );
}


void sliceBeforeGhostLayerXYZTest()
{
   const HostField_T emptyField( X_SIZE, Y_SIZE, Z_SIZE, GL_SIZE, -1, LAYOUT );
   GPUField_T deviceField( X_SIZE, Y_SIZE, Z_SIZE, F_SIZE, 1, LAYOUT );
   cuda::fieldCpy( deviceField, emptyField );

   auto setValue = cuda::make_kernel( &setValueKernel );
   setValue.addFieldIndexingParam( FieldIdx3D_T::sliceBeforeGhostLayerXYZ( deviceField, 1, stencil::B, true ) );
   setValue();

   HostField_T resultField( X_SIZE, Y_SIZE, Z_SIZE, GL_SIZE, -1, LAYOUT );
   cuda::fieldCpy( resultField, deviceField );

   HostField_T expectedField( X_SIZE, Y_SIZE, Z_SIZE, GL_SIZE, -1, LAYOUT );
   CellInterval ci;
   expectedField.getSliceBeforeGhostLayer( stencil::B, ci, 1, true );
   WALBERLA_FOR_ALL_CELLS_IN_INTERVAL_XYZ( ci,
      for ( uint_t f = 0; f < expectedField.fSize(); ++f )
      {
         expectedField.get( x, y, z, f ) = IDX4D( x - ci.xMin(), y - ci.yMin(), z - ci.zMin(), f );
      }
   )
   WALBERLA_ASSERT( resultField == expectedField );
}


int main( int argc, char ** argv )
{
   debug::enterTestMode();
   walberla::Environment walberlaEnv( argc, argv );

   xyzTest();
   sliceBeforeGhostLayerXYZTest();

   return 0;
}
