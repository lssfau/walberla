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
//! \file InterpolationTest.cpp
//! \ingroup field
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#include "field/GhostLayerField.h"
#include "field/interpolators/NearestNeighborInterpolator.h"
#include "field/interpolators/TrilinearInterpolator.h"
#include "core/debug/TestSubsystem.h"
#include "core/mpi/Environment.h"


using namespace walberla;



void testNearestNeighbor1D()
{
   GhostLayerField<real_t,1> field ( 2,2,2,2, 0);

   for( auto it = field.beginWithGhostLayer(); it != field.end(); ++it )
      *it = real_c(it.x());

   field::NearestNeighborInterpolator<GhostLayerField<real_t,1> > ip ( field );


   WALBERLA_CHECK_FLOAT_EQUAL( ip(   real_c(0), real_c(0), real_c(0) ),  real_c(0)   );
   WALBERLA_CHECK_FLOAT_EQUAL( ip( real_c(0.5), real_c(0), real_c(0) ),  real_c(0)   );
   WALBERLA_CHECK_FLOAT_EQUAL( ip(-real_c(0.5), real_c(0), real_c(0) ), -real_c(1) );
   WALBERLA_CHECK_FLOAT_EQUAL( ip( real_c(1.9), real_c(0), real_c(0) ),  real_c(1) );
}

void testTrilinear1D()
{
   GhostLayerField<real_t,1> field ( 2,2,2,2, 0);

   for( auto it = field.beginWithGhostLayer(); it != field.end(); ++it )
      *it = real_c(it.x());

   field::TrilinearInterpolator<GhostLayerField<real_t,1> > ip ( field );

   WALBERLA_CHECK_FLOAT_EQUAL( ip( real_c(0.5), real_c(0.5), real_c(0.5) ), real_c(0) );
   // Left Hit
   WALBERLA_CHECK_FLOAT_EQUAL( ip( real_c(0), real_c(0), real_c(0) ), -real_c(0.5) );
   // Right Hit
   WALBERLA_CHECK_FLOAT_EQUAL( ip( real_c(1), real_c(0), real_c(0) ),  real_c(0.5) );

   // Left boundary included
   WALBERLA_CHECK_FLOAT_EQUAL( ip( -real_c(0.5), real_c(0), real_c(0) ), -real_c(1) );
}


int main( int argc, char**argv )
{
   mpi::Environment mpiEnv( argc, argv );
   debug::enterTestMode();

   testTrilinear1D();
   testNearestNeighbor1D();
   return 0;
}
