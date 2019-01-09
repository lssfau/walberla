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
//! \\file .cpp
//! \\ingroup lbm
//! \\author lbmpy
//======================================================================================================================

#include <cmath>

#include "core/DataTypes.h"
#include "core/Macros.h"
#include "UniformGridGPU_LbKernel.h"


#define FUNC_PREFIX __global__

#if ( defined WALBERLA_CXX_COMPILER_IS_GNU ) || ( defined WALBERLA_CXX_COMPILER_IS_CLANG )
#   pragma GCC diagnostic push
#   pragma GCC diagnostic ignored "-Wfloat-equal"
#   pragma GCC diagnostic ignored "-Wshadow"
#   pragma GCC diagnostic ignored "-Wconversion"
#endif

using namespace std;

namespace walberla {
namespace pystencils {

namespace internal_UniformGridGPU_LbKernel {
static FUNC_PREFIX void UniformGridGPU_LbKernel(double * const _data_pdfs, double * _data_pdfs_tmp, int64_t const _size_pdfs_0, int64_t const _size_pdfs_1, int64_t const _size_pdfs_2, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2, int64_t const _stride_pdfs_3, double omega)
{
   if (blockDim.x*blockIdx.x + threadIdx.x + 1 < _size_pdfs_0 - 1 && blockDim.y*blockIdx.y + threadIdx.y + 1 < _size_pdfs_1 - 1 && blockDim.z*blockIdx.z + threadIdx.z + 1 < _size_pdfs_2 - 1)
   {
      const int64_t ctr_0 = blockDim.x*blockIdx.x + threadIdx.x + 1;
      const int64_t ctr_1 = blockDim.y*blockIdx.y + threadIdx.y + 1;
      const int64_t ctr_2 = blockDim.z*blockIdx.z + threadIdx.z + 1;
      double * const _data_pdfs_10_21_317 = _data_pdfs + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 17*_stride_pdfs_3;
      const double xi_18 = -_data_pdfs_10_21_317[_stride_pdfs_0*ctr_0 + _stride_pdfs_0];
      double * const _data_pdfs_11_20_39 = _data_pdfs + _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _stride_pdfs_2*ctr_2 + 9*_stride_pdfs_3;
      const double xi_19 = -_data_pdfs_11_20_39[_stride_pdfs_0*ctr_0 + _stride_pdfs_0];
      double * const _data_pdfs_11_21_316 = _data_pdfs + _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 16*_stride_pdfs_3;
      const double xi_20 = -_data_pdfs_11_21_316[_stride_pdfs_0*ctr_0];
      double * const _data_pdfs_10_2m1_314 = _data_pdfs + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 14*_stride_pdfs_3;
      double * const _data_pdfs_11_20_310 = _data_pdfs + _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _stride_pdfs_2*ctr_2 + 10*_stride_pdfs_3;
      double * const _data_pdfs_1m1_20_38 = _data_pdfs + _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _stride_pdfs_2*ctr_2 + 8*_stride_pdfs_3;
      double * const _data_pdfs_10_21_318 = _data_pdfs + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 18*_stride_pdfs_3;
      double * const _data_pdfs_10_20_34 = _data_pdfs + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 4*_stride_pdfs_3;
      const double vel0Term = _data_pdfs_10_20_34[_stride_pdfs_0*ctr_0 - _stride_pdfs_0] + _data_pdfs_10_21_318[_stride_pdfs_0*ctr_0 - _stride_pdfs_0] + _data_pdfs_10_2m1_314[_stride_pdfs_0*ctr_0 - _stride_pdfs_0] + _data_pdfs_11_20_310[_stride_pdfs_0*ctr_0 - _stride_pdfs_0] + _data_pdfs_1m1_20_38[_stride_pdfs_0*ctr_0 - _stride_pdfs_0];
      double * const _data_pdfs_1m1_2m1_311 = _data_pdfs + _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 11*_stride_pdfs_3;
      double * const _data_pdfs_1m1_20_37 = _data_pdfs + _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _stride_pdfs_2*ctr_2 + 7*_stride_pdfs_3;
      double * const _data_pdfs_1m1_20_31 = _data_pdfs + _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _stride_pdfs_2*ctr_2 + _stride_pdfs_3;
      double * const _data_pdfs_1m1_21_315 = _data_pdfs + _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 15*_stride_pdfs_3;
      const double vel1Term = _data_pdfs_1m1_20_31[_stride_pdfs_0*ctr_0] + _data_pdfs_1m1_20_37[_stride_pdfs_0*ctr_0 + _stride_pdfs_0] + _data_pdfs_1m1_21_315[_stride_pdfs_0*ctr_0] + _data_pdfs_1m1_2m1_311[_stride_pdfs_0*ctr_0];
      double * const _data_pdfs_10_2m1_313 = _data_pdfs + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 13*_stride_pdfs_3;
      double * const _data_pdfs_11_2m1_312 = _data_pdfs + _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 12*_stride_pdfs_3;
      double * const _data_pdfs_10_2m1_35 = _data_pdfs + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 5*_stride_pdfs_3;
      const double vel2Term = _data_pdfs_10_2m1_313[_stride_pdfs_0*ctr_0 + _stride_pdfs_0] + _data_pdfs_10_2m1_35[_stride_pdfs_0*ctr_0] + _data_pdfs_11_2m1_312[_stride_pdfs_0*ctr_0];
      double * const _data_pdfs_10_20_30 = _data_pdfs + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2;
      double * const _data_pdfs_10_20_33 = _data_pdfs + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 3*_stride_pdfs_3;
      double * const _data_pdfs_11_20_32 = _data_pdfs + _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _stride_pdfs_2*ctr_2 + 2*_stride_pdfs_3;
      double * const _data_pdfs_10_21_36 = _data_pdfs + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 6*_stride_pdfs_3;
      const double rho = vel0Term + vel1Term + vel2Term + _data_pdfs_10_20_30[_stride_pdfs_0*ctr_0] + _data_pdfs_10_20_33[_stride_pdfs_0*ctr_0 + _stride_pdfs_0] + _data_pdfs_10_21_317[_stride_pdfs_0*ctr_0 + _stride_pdfs_0] + _data_pdfs_10_21_36[_stride_pdfs_0*ctr_0] + _data_pdfs_11_20_32[_stride_pdfs_0*ctr_0] + _data_pdfs_11_20_39[_stride_pdfs_0*ctr_0 + _stride_pdfs_0] + _data_pdfs_11_21_316[_stride_pdfs_0*ctr_0];
      const double u_0 = vel0Term + xi_18 + xi_19 - _data_pdfs_10_20_33[_stride_pdfs_0*ctr_0 + _stride_pdfs_0] - _data_pdfs_10_2m1_313[_stride_pdfs_0*ctr_0 + _stride_pdfs_0] - _data_pdfs_1m1_20_37[_stride_pdfs_0*ctr_0 + _stride_pdfs_0];
      const double xi_23 = (u_0*u_0);
      const double u_1 = vel1Term + xi_19 + xi_20 - _data_pdfs_11_20_310[_stride_pdfs_0*ctr_0 - _stride_pdfs_0] - _data_pdfs_11_20_32[_stride_pdfs_0*ctr_0] - _data_pdfs_11_2m1_312[_stride_pdfs_0*ctr_0] + _data_pdfs_1m1_20_38[_stride_pdfs_0*ctr_0 - _stride_pdfs_0];
      const double xi_21 = -u_1;
      const double xi_24 = (u_1*u_1);
      const double u_2 = vel2Term + xi_18 + xi_20 - _data_pdfs_10_21_318[_stride_pdfs_0*ctr_0 - _stride_pdfs_0] - _data_pdfs_10_21_36[_stride_pdfs_0*ctr_0] + _data_pdfs_10_2m1_314[_stride_pdfs_0*ctr_0 - _stride_pdfs_0] - _data_pdfs_1m1_21_315[_stride_pdfs_0*ctr_0] + _data_pdfs_1m1_2m1_311[_stride_pdfs_0*ctr_0];
      const double xi_22 = -u_2;
      const double xi_25 = (u_2*u_2);
      const double u0Mu1 = u_0 + xi_21;
      const double u0Pu1 = u_0 + u_1;
      const double u1Pu2 = u_1 + u_2;
      const double u1Mu2 = u_1 + xi_22;
      const double u0Mu2 = u_0 + xi_22;
      const double u0Pu2 = u_0 + u_2;
      const double f_eq_common = rho - xi_23 - xi_24 - xi_25;
      const double xi_26 = f_eq_common + rho*-0.666666666666667;
      const double xi_27 = f_eq_common + rho*-0.333333333333333;
      const double xi_28 = xi_25 + xi_27;
      const double xi_29 = xi_23 + xi_27;
      const double xi_30 = xi_24 + xi_27;
      const double xi_2 = xi_24*2 + xi_26;
      const double xi_3 = xi_23*2 + xi_26;
      const double xi_4 = xi_25*2 + xi_26;
      const double xi_6 = u0Mu1*2;
      const double xi_7 = (u0Mu1*u0Mu1)*3 + xi_28;
      const double xi_8 = u0Pu1*2;
      const double xi_9 = (u0Pu1*u0Pu1)*3 + xi_28;
      const double xi_10 = u1Pu2*2;
      const double xi_11 = (u1Pu2*u1Pu2)*3 + xi_29;
      const double xi_12 = u1Mu2*2;
      const double xi_13 = (u1Mu2*u1Mu2)*3 + xi_29;
      const double xi_14 = u0Mu2*2;
      const double xi_15 = (u0Mu2*u0Mu2)*3 + xi_30;
      const double xi_16 = u0Pu2*2;
      const double xi_17 = (u0Pu2*u0Pu2)*3 + xi_30;
      const double xi_1 = omega*0.166666666666667;
      const double xi_5 = omega*0.0416666666666667;
      double * _data_pdfs_tmp_10_20_30 = _data_pdfs_tmp + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2;
      _data_pdfs_tmp_10_20_30[_stride_pdfs_0*ctr_0] = omega*(f_eq_common*0.333333333333333 - _data_pdfs_10_20_30[_stride_pdfs_0*ctr_0]) + _data_pdfs_10_20_30[_stride_pdfs_0*ctr_0];
      double * _data_pdfs_tmp_10_20_31 = _data_pdfs_tmp + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + _stride_pdfs_3;
      _data_pdfs_tmp_10_20_31[_stride_pdfs_0*ctr_0] = xi_1*(u_1 + xi_2 - 6*_data_pdfs_1m1_20_31[_stride_pdfs_0*ctr_0]) + _data_pdfs_1m1_20_31[_stride_pdfs_0*ctr_0];
      double * _data_pdfs_tmp_10_20_32 = _data_pdfs_tmp + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 2*_stride_pdfs_3;
      _data_pdfs_tmp_10_20_32[_stride_pdfs_0*ctr_0] = xi_1*(xi_2 + xi_21 - 6*_data_pdfs_11_20_32[_stride_pdfs_0*ctr_0]) + _data_pdfs_11_20_32[_stride_pdfs_0*ctr_0];
      double * _data_pdfs_tmp_10_20_33 = _data_pdfs_tmp + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 3*_stride_pdfs_3;
      _data_pdfs_tmp_10_20_33[_stride_pdfs_0*ctr_0] = xi_1*(-u_0 + xi_3 - 6*_data_pdfs_10_20_33[_stride_pdfs_0*ctr_0 + _stride_pdfs_0]) + _data_pdfs_10_20_33[_stride_pdfs_0*ctr_0 + _stride_pdfs_0];
      double * _data_pdfs_tmp_10_20_34 = _data_pdfs_tmp + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 4*_stride_pdfs_3;
      _data_pdfs_tmp_10_20_34[_stride_pdfs_0*ctr_0] = xi_1*(u_0 + xi_3 - 6*_data_pdfs_10_20_34[_stride_pdfs_0*ctr_0 - _stride_pdfs_0]) + _data_pdfs_10_20_34[_stride_pdfs_0*ctr_0 - _stride_pdfs_0];
      double * _data_pdfs_tmp_10_20_35 = _data_pdfs_tmp + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 5*_stride_pdfs_3;
      _data_pdfs_tmp_10_20_35[_stride_pdfs_0*ctr_0] = xi_1*(u_2 + xi_4 - 6*_data_pdfs_10_2m1_35[_stride_pdfs_0*ctr_0]) + _data_pdfs_10_2m1_35[_stride_pdfs_0*ctr_0];
      double * _data_pdfs_tmp_10_20_36 = _data_pdfs_tmp + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 6*_stride_pdfs_3;
      _data_pdfs_tmp_10_20_36[_stride_pdfs_0*ctr_0] = xi_1*(xi_22 + xi_4 - 6*_data_pdfs_10_21_36[_stride_pdfs_0*ctr_0]) + _data_pdfs_10_21_36[_stride_pdfs_0*ctr_0];
      double * _data_pdfs_tmp_10_20_37 = _data_pdfs_tmp + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 7*_stride_pdfs_3;
      _data_pdfs_tmp_10_20_37[_stride_pdfs_0*ctr_0] = xi_5*(-xi_6 + xi_7 - 24*_data_pdfs_1m1_20_37[_stride_pdfs_0*ctr_0 + _stride_pdfs_0]) + _data_pdfs_1m1_20_37[_stride_pdfs_0*ctr_0 + _stride_pdfs_0];
      double * _data_pdfs_tmp_10_20_38 = _data_pdfs_tmp + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 8*_stride_pdfs_3;
      _data_pdfs_tmp_10_20_38[_stride_pdfs_0*ctr_0] = xi_5*(xi_8 + xi_9 - 24*_data_pdfs_1m1_20_38[_stride_pdfs_0*ctr_0 - _stride_pdfs_0]) + _data_pdfs_1m1_20_38[_stride_pdfs_0*ctr_0 - _stride_pdfs_0];
      double * _data_pdfs_tmp_10_20_39 = _data_pdfs_tmp + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 9*_stride_pdfs_3;
      _data_pdfs_tmp_10_20_39[_stride_pdfs_0*ctr_0] = xi_5*(-xi_8 + xi_9 - 24*_data_pdfs_11_20_39[_stride_pdfs_0*ctr_0 + _stride_pdfs_0]) + _data_pdfs_11_20_39[_stride_pdfs_0*ctr_0 + _stride_pdfs_0];
      double * _data_pdfs_tmp_10_20_310 = _data_pdfs_tmp + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 10*_stride_pdfs_3;
      _data_pdfs_tmp_10_20_310[_stride_pdfs_0*ctr_0] = xi_5*(xi_6 + xi_7 - 24*_data_pdfs_11_20_310[_stride_pdfs_0*ctr_0 - _stride_pdfs_0]) + _data_pdfs_11_20_310[_stride_pdfs_0*ctr_0 - _stride_pdfs_0];
      double * _data_pdfs_tmp_10_20_311 = _data_pdfs_tmp + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 11*_stride_pdfs_3;
      _data_pdfs_tmp_10_20_311[_stride_pdfs_0*ctr_0] = xi_5*(xi_10 + xi_11 - 24*_data_pdfs_1m1_2m1_311[_stride_pdfs_0*ctr_0]) + _data_pdfs_1m1_2m1_311[_stride_pdfs_0*ctr_0];
      double * _data_pdfs_tmp_10_20_312 = _data_pdfs_tmp + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 12*_stride_pdfs_3;
      _data_pdfs_tmp_10_20_312[_stride_pdfs_0*ctr_0] = xi_5*(-xi_12 + xi_13 - 24*_data_pdfs_11_2m1_312[_stride_pdfs_0*ctr_0]) + _data_pdfs_11_2m1_312[_stride_pdfs_0*ctr_0];
      double * _data_pdfs_tmp_10_20_313 = _data_pdfs_tmp + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 13*_stride_pdfs_3;
      _data_pdfs_tmp_10_20_313[_stride_pdfs_0*ctr_0] = xi_5*(-xi_14 + xi_15 - 24*_data_pdfs_10_2m1_313[_stride_pdfs_0*ctr_0 + _stride_pdfs_0]) + _data_pdfs_10_2m1_313[_stride_pdfs_0*ctr_0 + _stride_pdfs_0];
      double * _data_pdfs_tmp_10_20_314 = _data_pdfs_tmp + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 14*_stride_pdfs_3;
      _data_pdfs_tmp_10_20_314[_stride_pdfs_0*ctr_0] = xi_5*(xi_16 + xi_17 - 24*_data_pdfs_10_2m1_314[_stride_pdfs_0*ctr_0 - _stride_pdfs_0]) + _data_pdfs_10_2m1_314[_stride_pdfs_0*ctr_0 - _stride_pdfs_0];
      double * _data_pdfs_tmp_10_20_315 = _data_pdfs_tmp + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 15*_stride_pdfs_3;
      _data_pdfs_tmp_10_20_315[_stride_pdfs_0*ctr_0] = xi_5*(xi_12 + xi_13 - 24*_data_pdfs_1m1_21_315[_stride_pdfs_0*ctr_0]) + _data_pdfs_1m1_21_315[_stride_pdfs_0*ctr_0];
      double * _data_pdfs_tmp_10_20_316 = _data_pdfs_tmp + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 16*_stride_pdfs_3;
      _data_pdfs_tmp_10_20_316[_stride_pdfs_0*ctr_0] = xi_5*(-xi_10 + xi_11 - 24*_data_pdfs_11_21_316[_stride_pdfs_0*ctr_0]) + _data_pdfs_11_21_316[_stride_pdfs_0*ctr_0];
      double * _data_pdfs_tmp_10_20_317 = _data_pdfs_tmp + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 17*_stride_pdfs_3;
      _data_pdfs_tmp_10_20_317[_stride_pdfs_0*ctr_0] = xi_5*(-xi_16 + xi_17 - 24*_data_pdfs_10_21_317[_stride_pdfs_0*ctr_0 + _stride_pdfs_0]) + _data_pdfs_10_21_317[_stride_pdfs_0*ctr_0 + _stride_pdfs_0];
      double * _data_pdfs_tmp_10_20_318 = _data_pdfs_tmp + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 18*_stride_pdfs_3;
      _data_pdfs_tmp_10_20_318[_stride_pdfs_0*ctr_0] = xi_5*(xi_14 + xi_15 - 24*_data_pdfs_10_21_318[_stride_pdfs_0*ctr_0 - _stride_pdfs_0]) + _data_pdfs_10_21_318[_stride_pdfs_0*ctr_0 - _stride_pdfs_0];
   } 
}
}

void UniformGridGPU_LbKernel::operator() ( IBlock * block , cudaStream_t stream )
{
    auto pdfs = block->getData< cuda::GPUField<double> >(pdfsID);
    cuda::GPUField<double> * pdfs_tmp;
    // Getting temporary field pdfs_tmp
    auto it = cache_pdfs_.find( pdfs );
    if( it != cache_pdfs_.end() )
    {
        pdfs_tmp = *it;
    }
    else 
    {
        pdfs_tmp = pdfs->cloneUninitialized();
        cache_pdfs_.insert(pdfs_tmp);
    }

    WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(pdfs->nrOfGhostLayers()));
    double * const _data_pdfs = pdfs->dataAt(-1, -1, -1, 0);
    WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(pdfs_tmp->nrOfGhostLayers()));
    double * _data_pdfs_tmp = pdfs_tmp->dataAt(-1, -1, -1, 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(), int64_t(pdfs->xSize() + 2));
    const int64_t _size_pdfs_0 = int64_t(pdfs->xSize() + 2);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(), int64_t(pdfs->ySize() + 2));
    const int64_t _size_pdfs_1 = int64_t(pdfs->ySize() + 2);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(), int64_t(pdfs->zSize() + 2));
    const int64_t _size_pdfs_2 = int64_t(pdfs->zSize() + 2);
    const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
    const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
    const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
    const int64_t _stride_pdfs_3 = int64_t(pdfs->fStride());
    dim3 _block(int(((128 < _size_pdfs_0 - 2) ? 128 : _size_pdfs_0 - 2)), int(((1 < _size_pdfs_1 - 2) ? 1 : _size_pdfs_1 - 2)), int(((1 < _size_pdfs_2 - 2) ? 1 : _size_pdfs_2 - 2)));
    dim3 _grid(int(( (_size_pdfs_0 - 2) % (((128 < _size_pdfs_0 - 2) ? 128 : _size_pdfs_0 - 2)) == 0 ? (int64_t)(_size_pdfs_0 - 2) / (int64_t)(((128 < _size_pdfs_0 - 2) ? 128 : _size_pdfs_0 - 2)) : ( (int64_t)(_size_pdfs_0 - 2) / (int64_t)(((128 < _size_pdfs_0 - 2) ? 128 : _size_pdfs_0 - 2)) ) +1 )), int(( (_size_pdfs_1 - 2) % (((1 < _size_pdfs_1 - 2) ? 1 : _size_pdfs_1 - 2)) == 0 ? (int64_t)(_size_pdfs_1 - 2) / (int64_t)(((1 < _size_pdfs_1 - 2) ? 1 : _size_pdfs_1 - 2)) : ( (int64_t)(_size_pdfs_1 - 2) / (int64_t)(((1 < _size_pdfs_1 - 2) ? 1 : _size_pdfs_1 - 2)) ) +1 )), int(( (_size_pdfs_2 - 2) % (((1 < _size_pdfs_2 - 2) ? 1 : _size_pdfs_2 - 2)) == 0 ? (int64_t)(_size_pdfs_2 - 2) / (int64_t)(((1 < _size_pdfs_2 - 2) ? 1 : _size_pdfs_2 - 2)) : ( (int64_t)(_size_pdfs_2 - 2) / (int64_t)(((1 < _size_pdfs_2 - 2) ? 1 : _size_pdfs_2 - 2)) ) +1 )));
    internal_UniformGridGPU_LbKernel::UniformGridGPU_LbKernel<<<_grid, _block, 0, stream>>>(_data_pdfs, _data_pdfs_tmp, _size_pdfs_0, _size_pdfs_1, _size_pdfs_2, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3, omega);
    pdfs->swapDataPointers(pdfs_tmp);

}



void UniformGridGPU_LbKernel::inner( IBlock * block , cudaStream_t stream )
{
    auto pdfs = block->getData< cuda::GPUField<double> >(pdfsID);
    cuda::GPUField<double> * pdfs_tmp;
    // Getting temporary field pdfs_tmp
    auto it = cache_pdfs_.find( pdfs );
    if( it != cache_pdfs_.end() )
    {
        pdfs_tmp = *it;
    }
    else 
    {
        pdfs_tmp = pdfs->cloneUninitialized();
        cache_pdfs_.insert(pdfs_tmp);
    }


    CellInterval inner = pdfs->xyzSize();
    inner.expand(-1);

    WALBERLA_ASSERT_GREATER_EQUAL(inner.xMin() - 1, -int_c(pdfs->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(inner.yMin() - 1, -int_c(pdfs->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(inner.zMin() - 1, -int_c(pdfs->nrOfGhostLayers()));
    double * const _data_pdfs = pdfs->dataAt(inner.xMin() - 1, inner.yMin() - 1, inner.zMin() - 1, 0);
    WALBERLA_ASSERT_GREATER_EQUAL(inner.xMin() - 1, -int_c(pdfs_tmp->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(inner.yMin() - 1, -int_c(pdfs_tmp->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(inner.zMin() - 1, -int_c(pdfs_tmp->nrOfGhostLayers()));
    double * _data_pdfs_tmp = pdfs_tmp->dataAt(inner.xMin() - 1, inner.yMin() - 1, inner.zMin() - 1, 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(), int64_t(inner.xSize() + 2));
    const int64_t _size_pdfs_0 = int64_t(inner.xSize() + 2);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(), int64_t(inner.ySize() + 2));
    const int64_t _size_pdfs_1 = int64_t(inner.ySize() + 2);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(), int64_t(inner.zSize() + 2));
    const int64_t _size_pdfs_2 = int64_t(inner.zSize() + 2);
    const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
    const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
    const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
    const int64_t _stride_pdfs_3 = int64_t(pdfs->fStride());
    dim3 _block(int(((128 < _size_pdfs_0 - 2) ? 128 : _size_pdfs_0 - 2)), int(((1 < _size_pdfs_1 - 2) ? 1 : _size_pdfs_1 - 2)), int(((1 < _size_pdfs_2 - 2) ? 1 : _size_pdfs_2 - 2)));
    dim3 _grid(int(( (_size_pdfs_0 - 2) % (((128 < _size_pdfs_0 - 2) ? 128 : _size_pdfs_0 - 2)) == 0 ? (int64_t)(_size_pdfs_0 - 2) / (int64_t)(((128 < _size_pdfs_0 - 2) ? 128 : _size_pdfs_0 - 2)) : ( (int64_t)(_size_pdfs_0 - 2) / (int64_t)(((128 < _size_pdfs_0 - 2) ? 128 : _size_pdfs_0 - 2)) ) +1 )), int(( (_size_pdfs_1 - 2) % (((1 < _size_pdfs_1 - 2) ? 1 : _size_pdfs_1 - 2)) == 0 ? (int64_t)(_size_pdfs_1 - 2) / (int64_t)(((1 < _size_pdfs_1 - 2) ? 1 : _size_pdfs_1 - 2)) : ( (int64_t)(_size_pdfs_1 - 2) / (int64_t)(((1 < _size_pdfs_1 - 2) ? 1 : _size_pdfs_1 - 2)) ) +1 )), int(( (_size_pdfs_2 - 2) % (((1 < _size_pdfs_2 - 2) ? 1 : _size_pdfs_2 - 2)) == 0 ? (int64_t)(_size_pdfs_2 - 2) / (int64_t)(((1 < _size_pdfs_2 - 2) ? 1 : _size_pdfs_2 - 2)) : ( (int64_t)(_size_pdfs_2 - 2) / (int64_t)(((1 < _size_pdfs_2 - 2) ? 1 : _size_pdfs_2 - 2)) ) +1 )));
    internal_UniformGridGPU_LbKernel::UniformGridGPU_LbKernel<<<_grid, _block, 0, stream>>>(_data_pdfs, _data_pdfs_tmp, _size_pdfs_0, _size_pdfs_1, _size_pdfs_2, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3, omega);
}


void UniformGridGPU_LbKernel::outer( IBlock * block , cudaStream_t stream  )
{
    static std::vector<CellInterval> layers;

    auto pdfs = block->getData< cuda::GPUField<double> >(pdfsID);
    cuda::GPUField<double> * pdfs_tmp;
    // Getting temporary field pdfs_tmp
    auto it = cache_pdfs_.find( pdfs );
    if( it != cache_pdfs_.end() )
    {
        pdfs_tmp = *it;
    }
    else 
    {
        pdfs_tmp = pdfs->cloneUninitialized();
        cache_pdfs_.insert(pdfs_tmp);
    }


    if( layers.size() == 0 )
    {
        CellInterval ci;

        pdfs->getSliceBeforeGhostLayer(stencil::T, ci, 1, false);
        layers.push_back(ci);
        pdfs->getSliceBeforeGhostLayer(stencil::B, ci, 1, false);
        layers.push_back(ci);

        pdfs->getSliceBeforeGhostLayer(stencil::N, ci, 1, false);
        ci.expand(Cell(0, 0, -1));
        layers.push_back(ci);
        pdfs->getSliceBeforeGhostLayer(stencil::S, ci, 1, false);
        ci.expand(Cell(0, 0, -1));
        layers.push_back(ci);

        pdfs->getSliceBeforeGhostLayer(stencil::E, ci, 1, false);
        ci.expand(Cell(0, -1, -1));
        layers.push_back(ci);
        pdfs->getSliceBeforeGhostLayer(stencil::W, ci, 1, false);
        ci.expand(Cell(0, -1, -1));
        layers.push_back(ci);
    }

    
    {
        auto parallelSection_ = parallelStreams_.parallelSection( stream );
        for( auto & ci: layers )
        {
            parallelSection_.run([&]( auto s ) {
                WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(pdfs->nrOfGhostLayers()));
                WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(pdfs->nrOfGhostLayers()));
                WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(pdfs->nrOfGhostLayers()));
                double * const _data_pdfs = pdfs->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
                WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(pdfs_tmp->nrOfGhostLayers()));
                WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(pdfs_tmp->nrOfGhostLayers()));
                WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(pdfs_tmp->nrOfGhostLayers()));
                double * _data_pdfs_tmp = pdfs_tmp->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
                WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(), int64_t(ci.xSize() + 2));
                const int64_t _size_pdfs_0 = int64_t(ci.xSize() + 2);
                WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(), int64_t(ci.ySize() + 2));
                const int64_t _size_pdfs_1 = int64_t(ci.ySize() + 2);
                WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(), int64_t(ci.zSize() + 2));
                const int64_t _size_pdfs_2 = int64_t(ci.zSize() + 2);
                const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
                const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
                const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
                const int64_t _stride_pdfs_3 = int64_t(pdfs->fStride());
                dim3 _block(int(((128 < _size_pdfs_0 - 2) ? 128 : _size_pdfs_0 - 2)), int(((1 < _size_pdfs_1 - 2) ? 1 : _size_pdfs_1 - 2)), int(((1 < _size_pdfs_2 - 2) ? 1 : _size_pdfs_2 - 2)));
                dim3 _grid(int(( (_size_pdfs_0 - 2) % (((128 < _size_pdfs_0 - 2) ? 128 : _size_pdfs_0 - 2)) == 0 ? (int64_t)(_size_pdfs_0 - 2) / (int64_t)(((128 < _size_pdfs_0 - 2) ? 128 : _size_pdfs_0 - 2)) : ( (int64_t)(_size_pdfs_0 - 2) / (int64_t)(((128 < _size_pdfs_0 - 2) ? 128 : _size_pdfs_0 - 2)) ) +1 )), int(( (_size_pdfs_1 - 2) % (((1 < _size_pdfs_1 - 2) ? 1 : _size_pdfs_1 - 2)) == 0 ? (int64_t)(_size_pdfs_1 - 2) / (int64_t)(((1 < _size_pdfs_1 - 2) ? 1 : _size_pdfs_1 - 2)) : ( (int64_t)(_size_pdfs_1 - 2) / (int64_t)(((1 < _size_pdfs_1 - 2) ? 1 : _size_pdfs_1 - 2)) ) +1 )), int(( (_size_pdfs_2 - 2) % (((1 < _size_pdfs_2 - 2) ? 1 : _size_pdfs_2 - 2)) == 0 ? (int64_t)(_size_pdfs_2 - 2) / (int64_t)(((1 < _size_pdfs_2 - 2) ? 1 : _size_pdfs_2 - 2)) : ( (int64_t)(_size_pdfs_2 - 2) / (int64_t)(((1 < _size_pdfs_2 - 2) ? 1 : _size_pdfs_2 - 2)) ) +1 )));
                internal_UniformGridGPU_LbKernel::UniformGridGPU_LbKernel<<<_grid, _block, 0, s>>>(_data_pdfs, _data_pdfs_tmp, _size_pdfs_0, _size_pdfs_1, _size_pdfs_2, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3, omega);
            });
        }
    }
    

    pdfs->swapDataPointers(pdfs_tmp);

}


} // namespace pystencils
} // namespace walberla


#if ( defined WALBERLA_CXX_COMPILER_IS_GNU ) || ( defined WALBERLA_CXX_COMPILER_IS_CLANG )
#   pragma GCC diagnostic pop
#endif