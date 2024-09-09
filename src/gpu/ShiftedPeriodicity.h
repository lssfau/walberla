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
//! \file ShiftedPeriodicity.h
//! \ingroup gpu
//! \author Helen Schottenhamml <helen.schottenhamml@fau.de>
//
//======================================================================================================================

#include "blockforest/StructuredBlockForest.h"

#include "boundary/ShiftedPeriodicity.h"

#include "core/DataTypes.h"
#include "core/cell/CellInterval.h"
#include "core/debug/Debug.h"
#include "core/math/Vector3.h"

#include "domain_decomposition/BlockDataID.h"
#include "domain_decomposition/IBlock.h"

#include "gpu/DeviceWrapper.h"
#include "gpu/FieldAccessor.h"
#include "gpu/FieldIndexing.h"
#include "gpu/GPUField.h"
#include "gpu/GPUWrapper.h"
#include "gpu/Kernel.h"

#include <cstdlib>
#include <device_launch_parameters.h>
#include <memory>
#include <vector>

#include "ErrorChecking.h"

namespace walberla {
namespace gpu
{

namespace internal {
   // GPU kernels - can be extended for other data types
   __global__ void packBufferGPU( gpu::FieldAccessor<real_t> fa, real_t * const buffer );
   __global__ void unpackBufferGPU( gpu::FieldAccessor<real_t> fa, const real_t * const buffer );
}

//*******************************************************************************************************************
/*!
 * A periodicity boundary condition that adds a user-defined spatial shift to the field when applied.
 * This shift can prevent the locking of large-scale turbulent features in the flow direction, see e.g.,
 * Munters et al. (https://doi.org/10.1063/1.4941912).
 *
 * Periodicity defined in the blockforest must be turned off in the normal-direction.
 *
 * This class handles the GPU-specific packing and unpacking of the communication buffers.
 *
 * @tparam GhostLayerField_T Type of the ghost-layer field that is shifted periodically
 */
//*******************************************************************************************************************
template< typename GPUField_T >
class ShiftedPeriodicityGPU : public boundary::ShiftedPeriodicityBase<ShiftedPeriodicityGPU<GPUField_T>, GPUField_T> {

   using Base = boundary::ShiftedPeriodicityBase<ShiftedPeriodicityGPU<GPUField_T>, GPUField_T>;
   friend Base;

 public:

   using ValueType = typename GPUField_T::value_type;
   using ShiftType = typename Base::ShiftType;
   using FieldIdx_T = gpu::FieldIndexing<ValueType>;

   ShiftedPeriodicityGPU(const std::weak_ptr< StructuredBlockForest >& blockForest,
                         const BlockDataID& fieldID, const uint_t fieldGhostLayers,
                         const uint_t normalDir, const uint_t shiftDir, const ShiftType shiftValue)
      : Base(blockForest, fieldID, fieldGhostLayers, normalDir, shiftDir, shiftValue)
   {}


 private:

   void packBuffer(IBlock* const block, const CellInterval& ci, std::vector< ValueType >& h_buffer) {

      // get field
      auto d_field = block->getData< GPUField_T >(this->fieldID_);
      WALBERLA_ASSERT_NOT_NULLPTR(d_field)

      const uint_t nValues = ci.numCells() * uint_c(this->fSize_);

      // create GPU buffer
      ValueType * d_buffer{};
      WALBERLA_GPU_CHECK(gpuMalloc(&d_buffer, nValues * sizeof(ValueType)))

      // fill buffer on GPU
      auto packKernel = gpu::make_kernel( &internal::packBufferGPU );
      packKernel.addFieldIndexingParam( FieldIdx_T::interval( *d_field, ci, 0, this->fSize_ ) );
      packKernel.addParam<real_t*>(d_buffer);
      packKernel();

      // copy from device to host buffer
      WALBERLA_GPU_CHECK(gpuMemcpy(h_buffer.data(), d_buffer, nValues * sizeof(ValueType), gpuMemcpyDeviceToHost))

      WALBERLA_GPU_CHECK(gpuFree(d_buffer))

   }

   void unpackBuffer(IBlock* const block, const CellInterval& ci, const std::vector< ValueType >& h_buffer) {

      // get field
      auto d_field = block->getData< GPUField_T >(this->fieldID_);
      WALBERLA_ASSERT_NOT_NULLPTR(d_field)

      const uint_t nValues = ci.numCells() * uint_c(this->fSize_);

      // create GPU buffer
      ValueType * d_buffer{};
      WALBERLA_GPU_CHECK(gpuMalloc(&d_buffer, nValues * sizeof(ValueType)))

      // copy from host to device buffer
      WALBERLA_GPU_CHECK(gpuMemcpy(d_buffer, h_buffer.data(), nValues * sizeof(ValueType), gpuMemcpyHostToDevice))

      // unpack buffer on GPU
      auto unpackKernel = gpu::make_kernel( &internal::unpackBufferGPU );
      unpackKernel.addFieldIndexingParam( FieldIdx_T::interval( *d_field, ci, 0, this->fSize_ ) );
      unpackKernel.addParam<const real_t*>(d_buffer);
      unpackKernel();

      WALBERLA_GPU_CHECK(gpuFree(d_buffer))
   }

}; // class ShiftedPeriodicity

} // namespace gpu
} // namespace walberla
