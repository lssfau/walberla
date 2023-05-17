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

#pragma once

#include "gpu/FieldAccessor3D.h"

#define X_SIZE    (64-2)
#define Y_SIZE    (64-2)
#define Z_SIZE    (64-2)
#define F_SIZE    19
#define LAYOUT    field::fzyx
#define GL_SIZE   1

#define YOFFSET               ( int( X_SIZE ) )
#define ZOFFSET               ( int( Y_SIZE ) * int( YOFFSET ) )
#define FOFFSET               ( int( Z_SIZE ) * int( ZOFFSET ) )
#define IDX4D( x, y, z, f )   ( (int)( int(f) * int(FOFFSET) + int(z) * int(Z_SIZE) + int(y) * int(YOFFSET) + int(x) ) )


namespace walberla {

using FieldAccessor3D_T = gpu::FieldAccessor3D<int>;

__global__ void setValueKernel( FieldAccessor3D_T fa );


} // namespace walberla
