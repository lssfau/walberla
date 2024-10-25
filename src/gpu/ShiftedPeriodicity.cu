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
//! \file ShiftedPeriodicity.cu
//! \ingroup gpu
//! \author Helen Schottenhamml <helen.schottenhamml@fau.de>
//
//======================================================================================================================

#include "gpu/ShiftedPeriodicity.h"

namespace walberla {
namespace gpu
{

namespace internal {
#ifdef WALBERLA_BUILD_WITH_GPU_SUPPORT
   __global__ void packBufferGPU( gpu::FieldAccessor<real_t> fa, real_t * const buffer ) {
      fa.set(blockIdx, threadIdx);
      if(fa.isValidPosition()) {
         buffer[fa.getLinearIndex(blockIdx, threadIdx, gridDim, blockDim)] = fa.get();
      }
   }
   __global__ void unpackBufferGPU( gpu::FieldAccessor<real_t> fa, const real_t * const buffer ) {
      fa.set(blockIdx, threadIdx);
      if(fa.isValidPosition()) {
         fa.get() = buffer[fa.getLinearIndex(blockIdx, threadIdx, gridDim, blockDim)];
      }
   }
#else
   __global__ void packBufferGPU( gpu::FieldAccessor<real_t>, real_t * const ) {
      WALBERLA_ABORT("gpu/ShiftedPeriodicity only supported when built with GPU support. Please use boundary/ShiftedPeriodicity on CPUs")
   }
   __global__ void unpackBufferGPU( gpu::FieldAccessor<real_t>, const real_t * const ) {
      WALBERLA_ABORT("gpu/ShiftedPeriodicity only supported when built with GPU support. Please use boundary/ShiftedPeriodicity on CPUs")
   }
#endif
}

} // namespace gpu
} // namespace walberla
