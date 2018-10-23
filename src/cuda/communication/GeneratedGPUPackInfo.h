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
//! \file GeneratedGPUPackInfo.h
//! \ingroup core
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================


#pragma once
#include "stencil/Directions.h"
#include "domain_decomposition/IBlock.h"
#include <cuda_runtime.h>


namespace walberla {
namespace cuda {


class GeneratedGPUPackInfo
{
public:
   virtual void pack  ( stencil::Direction dir, unsigned char *buffer, IBlock *block, cudaStream_t stream ) = 0;
   virtual void unpack( stencil::Direction dir, unsigned char *buffer, IBlock *block, cudaStream_t stream ) = 0;
   virtual uint_t size( stencil::Direction dir, IBlock *block ) = 0;
};



} //namespace cuda
} //namespace walberla