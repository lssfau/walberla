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

#include "FieldIndexing3DTest.h"


namespace walberla {

__global__ void setValueKernel( FieldAccessor3D_T fa )
{
   unsigned int const x = blockIdx.x * blockDim.x + threadIdx.x;
   unsigned int const y = blockIdx.y * blockDim.y + threadIdx.y;
   unsigned int const z = blockIdx.z * blockDim.z + threadIdx.z;
   fa.set( blockIdx, threadIdx );

   if ( fa.isValidPosition() )
   {
      for ( int f = 0; f < F_SIZE; ++f )
      {
         fa.get(f) = IDX4D( x, y, z, f );
      }
   }
}


} // namespace walberla
