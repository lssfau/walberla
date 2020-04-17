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
//! \file MatrixVectorOperationsTest.cpp
//! \ingroup mesh
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//
//======================================================================================================================

#include "core/debug/TestSubsystem.h"
#include "core/mpi/Environment.h"

#include "mesh_common/MatrixVectorOperations.h"

namespace walberla {
namespace mesh {


int main( int /*argc*/, char * /*argv*/[] )
{
   debug::enterTestMode();

   const math::Matrix3<real_t> m( real_t(2), real_t(3), real_t(4),
                                  real_t(5), real_t(6), real_t(7),
                                  real_t(8), real_t(9), real_t(10) );

   const math::Vector3<real_t> vwb( real_t(2), real_t(3), real_t(4) );
   const OpenMesh::VectorT<real_t, 3> vom( real_t(2), real_t(3), real_t(4) );

   WALBERLA_CHECK_EQUAL( vwb, toWalberla( vom ) );
   WALBERLA_CHECK_EQUAL( toOpenMesh(vwb), vom );

   WALBERLA_CHECK_EQUAL( toWalberla( vom ), vwb );
   WALBERLA_CHECK_EQUAL( vom, toOpenMesh(vwb) );

   WALBERLA_CHECK_EQUAL( m * vwb, toWalberla(m * vom) );


   const math::Vector3<double> vwbd( double(2), double(3), double(4) );
   const OpenMesh::VectorT<double, 3> vomd( double(2), double(3), double(4) );
   const math::Vector3<float> vwbf( float(2), float(3), float(4) );
   const OpenMesh::VectorT<float, 3> vomf( float(2), float(3), float(4) );

   WALBERLA_CHECK_LESS( ( vwbd - toWalberlaNumericCast<double>(vomf) ).length(), real_comparison::Epsilon<float>::value );
   WALBERLA_CHECK_LESS( ( vomd - toOpenMeshNumericCast<double>(vwbf) ).length(), real_comparison::Epsilon<float>::value );
   WALBERLA_CHECK_LESS( ( vwbf - toWalberlaNumericCast<float>(vomd)  ).length(), real_comparison::Epsilon<float>::value );
   WALBERLA_CHECK_LESS( ( vomf - toOpenMeshNumericCast<float>(vwbd)  ).length(), real_comparison::Epsilon<float>::value );

   return EXIT_SUCCESS;
}


} // namespace mesh
} // namespace walberla

int main( int argc, char * argv[] )
{
   return walberla::mesh::main( argc, argv );
}