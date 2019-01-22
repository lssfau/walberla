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
//! \\file UniformGridGPU_NoSlip.cpp
//! \\ingroup lbm
//! \\author lbmpy
//======================================================================================================================

#include <cmath>

#include "core/DataTypes.h"
#include "core/Macros.h"
#include "UniformGridGPU_NoSlip.h"
#include "cuda/ErrorChecking.h"


#define FUNC_PREFIX __global__

using namespace std;

namespace walberla {
namespace lbm {

#ifdef __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wstrict-aliasing"
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wconversion"
#endif

#ifdef __CUDACC__
#pragma push
#pragma diag_suppress = declared_but_not_referenced
#endif


namespace internal_boundary_UniformGridGPU_NoSlip {
static FUNC_PREFIX void boundary_UniformGridGPU_NoSlip(uint8_t * const _data_indexVector, double * _data_pdfs, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2, int64_t const _stride_pdfs_3, int64_t indexVectorSize)
{
   if (blockDim.x*blockIdx.x + threadIdx.x < indexVectorSize)
   {
      uint8_t * const _data_indexVector_10 = _data_indexVector;
      const int32_t x = *((int32_t *)(& _data_indexVector_10[16*blockDim.x*blockIdx.x + 16*threadIdx.x]));
      uint8_t * const _data_indexVector_14 = _data_indexVector + 4;
      const int32_t y = *((int32_t *)(& _data_indexVector_14[16*blockDim.x*blockIdx.x + 16*threadIdx.x]));
      uint8_t * const _data_indexVector_18 = _data_indexVector + 8;
      const int32_t z = *((int32_t *)(& _data_indexVector_18[16*blockDim.x*blockIdx.x + 16*threadIdx.x]));
      
      
      const int64_t cx [] = { 0, 0, 0, -1, 1, 0, 0, -1, 1, -1, 1, 0, 0, -1, 1, 0, 0, -1, 1 };
      const int64_t cy [] = { 0, 1, -1, 0, 0, 0, 0, 1, 1, -1, -1, 1, -1, 0, 0, 1, -1, 0, 0 };
      const int64_t cz [] = { 0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 1, 1, 1, 1, -1, -1, -1, -1 };
      const int invdir [] = { 0, 2, 1, 4, 3, 6, 5, 10, 9, 8, 7, 16, 15, 18, 17, 12, 11, 14, 13 };
      
      
      const double weights [] = { 0.333333333333333,0.0555555555555556,0.0555555555555556,0.0555555555555556,0.0555555555555556,0.0555555555555556,0.0555555555555556,0.0277777777777778,0.0277777777777778,0.0277777777777778,0.0277777777777778,0.0277777777777778,0.0277777777777778,0.0277777777777778,0.0277777777777778,0.0277777777777778,0.0277777777777778,0.0277777777777778,0.0277777777777778 };
      
      uint8_t * const _data_indexVector_112 = _data_indexVector + 12;
      const int32_t dir = *((int32_t *)(& _data_indexVector_112[16*blockDim.x*blockIdx.x + 16*threadIdx.x]));
      double * _data_pdfsf9cc34cc4e2b6261 = _data_pdfs + _stride_pdfs_1*y + _stride_pdfs_1*cy[dir] + _stride_pdfs_2*z + _stride_pdfs_2*cz[dir] + _stride_pdfs_3*invdir[dir];
      double * _data_pdfs_10_2011ac6bf6446d4afa = _data_pdfs + _stride_pdfs_1*y + _stride_pdfs_2*z + _stride_pdfs_3*dir;
      _data_pdfsf9cc34cc4e2b6261[_stride_pdfs_0*x + _stride_pdfs_0*cx[dir]] = _data_pdfs_10_2011ac6bf6446d4afa[_stride_pdfs_0*x];
   } 
}
}

#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif

#ifdef __CUDACC__
#pragma pop
#endif


void UniformGridGPU_NoSlip::run( IBlock * block, IndexVectors::Type type , cudaStream_t stream )
{
    auto * indexVectors = block->getData<IndexVectors>(indexVectorID);

    auto pointer = indexVectors->pointerGpu(type);
    

    int64_t indexVectorSize = int64_c( indexVectors->indexVector(type).size() );
    if( indexVectorSize == 0)
        return;

    uint8_t * _data_indexVector = reinterpret_cast<uint8_t*>(pointer);

    auto pdfs = block->getData< cuda::GPUField<double> >(pdfsID);

    WALBERLA_ASSERT_GREATER_EQUAL(0, -int_c(pdfs->nrOfGhostLayers()));
    double * _data_pdfs = pdfs->dataAt(0, 0, 0, 0);
    const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
    const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
    const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
    const int64_t _stride_pdfs_3 = int64_t(pdfs->fStride());
    dim3 _block(int(((256 < indexVectorSize) ? 256 : indexVectorSize)), int(1), int(1));
    dim3 _grid(int(( (indexVectorSize) % (((256 < indexVectorSize) ? 256 : indexVectorSize)) == 0 ? (int64_t)(indexVectorSize) / (int64_t)(((256 < indexVectorSize) ? 256 : indexVectorSize)) : ( (int64_t)(indexVectorSize) / (int64_t)(((256 < indexVectorSize) ? 256 : indexVectorSize)) ) +1 )), int(1), int(1));
    internal_boundary_UniformGridGPU_NoSlip::boundary_UniformGridGPU_NoSlip<<<_grid, _block, 0, stream>>>(_data_indexVector, _data_pdfs, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3, indexVectorSize);
}

void UniformGridGPU_NoSlip::operator() ( IBlock * block, cudaStream_t stream  )
{
    run( block, IndexVectors::ALL, stream );
}

void UniformGridGPU_NoSlip::inner( IBlock * block, cudaStream_t stream  )
{
    run( block, IndexVectors::INNER, stream  );
}

void UniformGridGPU_NoSlip::outer( IBlock * block, cudaStream_t stream  )
{
    run( block, IndexVectors::OUTER, stream  );
}


} // namespace lbm
} // namespace walberla

