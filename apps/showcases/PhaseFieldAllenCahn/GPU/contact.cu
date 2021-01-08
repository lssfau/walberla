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
//! \file contact.cu
//! \author Markus Holzer <markus.holzer@fau.de>
//
//======================================================================================================================

#include "core/DataTypes.h"
#include "core/Macros.h"

#include <cmath>

#include "contact.h"

#define FUNC_PREFIX __global__

namespace walberla
{
namespace lbm
{
#ifdef __GNUC__
#   pragma GCC diagnostic push
#endif

#ifdef __CUDACC__
#   pragma push
#endif

namespace internal_boundary_contact
{
static FUNC_PREFIX void contact_angle_treatment(uint8_t* WALBERLA_RESTRICT const _data_indexVector, double* WALBERLA_RESTRICT _data_phase,
                                                int64_t const _stride_phase_0, int64_t const _stride_phase_1,
                                                int64_t const _stride_phase_2, int64_t indexVectorSize, double alpha)
{
   if (blockDim.x * blockIdx.x + threadIdx.x < indexVectorSize)
   {
      uint8_t* WALBERLA_RESTRICT _data_indexVector_10 = _data_indexVector;
      const int32_t x = *((int32_t*) (&_data_indexVector_10[24 * blockDim.x * blockIdx.x + 24 * threadIdx.x]));
      uint8_t* WALBERLA_RESTRICT _data_indexVector_14 = _data_indexVector + 4;
      const int32_t y = *((int32_t*) (&_data_indexVector_14[24 * blockDim.x * blockIdx.x + 24 * threadIdx.x]));
      uint8_t* WALBERLA_RESTRICT _data_indexVector_18 = _data_indexVector + 8;
      const int32_t z = *((int32_t*) (&_data_indexVector_18[24 * blockDim.x * blockIdx.x + 24 * threadIdx.x]));
      uint8_t* WALBERLA_RESTRICT _data_indexVector_112 = _data_indexVector + 12;
      const int32_t nx = *((int32_t*) (&_data_indexVector_112[24 * blockDim.x * blockIdx.x + 24 * threadIdx.x]));
      const int32_t x1 = x + nx;
      uint8_t* WALBERLA_RESTRICT _data_indexVector_116 = _data_indexVector + 16;
      const int32_t ny = *((int32_t*) (&_data_indexVector_116[24 * blockDim.x * blockIdx.x + 24 * threadIdx.x]));
      const int32_t y1 = y + ny;
      uint8_t* WALBERLA_RESTRICT _data_indexVector_200 = _data_indexVector + 20;
      const int32_t nz = *((int32_t*) (&_data_indexVector_200[24 * blockDim.x * blockIdx.x + 24 * threadIdx.x]));
      const int32_t z1 = z + nz;

      const double h = 0.5 * sqrt((float) (nx * nx + ny * ny + nz * nz));
      const double a = cos(alpha);
      const double W = 5;

      double* WALBERLA_RESTRICT _phase_wall     = _data_phase + _stride_phase_1 * y + _stride_phase_2 * z;
      double* WALBERLA_RESTRICT _phase_interior = _data_phase + _stride_phase_1 * y1 + _stride_phase_2 * z1;
      if (h < 0.001) { _phase_wall[_stride_phase_0 * x] = 1.0; }
      else if (a > 1e-8 || a < -1e-8)
      {
         const double var = -h * (4.0 / W) * a;
         _phase_wall[_stride_phase_0 * x] =
            (1 + var - sqrt((1 + var) * (1 + var) - 4 * var * _phase_interior[_stride_phase_0 * x1])) / (var + 1e-12) -
            _phase_interior[_stride_phase_0 * x1];
      }
      else
      {
         _phase_wall[_stride_phase_0 * x] = _phase_interior[_stride_phase_0 * x1];
      }
   }
}
} // namespace internal_boundary_contact

#ifdef __GNUC__
#   pragma GCC diagnostic pop
#endif

#ifdef __CUDACC__
#   pragma pop
#endif

void contact::run(IBlock* block, IndexVectors::Type type, cudaStream_t stream)
{
   auto* indexVectors      = block->getData< IndexVectors >(indexVectorID);
   int64_t indexVectorSize = int64_c(indexVectors->indexVector(type).size());
   if (indexVectorSize == 0) return;

   auto pointer = indexVectors->pointerGpu(type);

   auto* _data_indexVector = reinterpret_cast< uint8_t* >(pointer);

   auto phaseField = block->getData< cuda::GPUField< double > >(phaseFieldID);

   auto& alpha = this->alpha_;
   WALBERLA_ASSERT_GREATER_EQUAL(0, -int_c(phaseField->nrOfGhostLayers()))
   double* WALBERLA_RESTRICT _data_phase = phaseField->dataAt(0, 0, 0, 0);
   const auto _stride_pdfs_0    = int64_t(phaseField->xStride());
   const auto _stride_pdfs_1    = int64_t(phaseField->yStride());
   const auto _stride_pdfs_2    = int64_t(phaseField->zStride());
   dim3 _block(int(((256 < indexVectorSize) ? 256 : indexVectorSize)), int(1), int(1));
   dim3 _grid(int(((indexVectorSize) % (((256 < indexVectorSize) ? 256 : indexVectorSize)) == 0 ?
                      (int64_t)(indexVectorSize) / (int64_t)(((256 < indexVectorSize) ? 256 : indexVectorSize)) :
                      ((int64_t)(indexVectorSize) / (int64_t)(((256 < indexVectorSize) ? 256 : indexVectorSize))) + 1)),
              int(1), int(1));
   internal_boundary_contact::contact_angle_treatment<<< _grid, _block, 0, stream >>>(
      _data_indexVector, _data_phase, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, indexVectorSize, alpha);
}

void contact::operator()(IBlock* block, cudaStream_t stream) { run(block, IndexVectors::ALL, stream); }

void contact::inner(IBlock* block, cudaStream_t stream) { run(block, IndexVectors::INNER, stream); }

void contact::outer(IBlock* block, cudaStream_t stream) { run(block, IndexVectors::OUTER, stream); }

} // namespace lbm
} // namespace walberla
