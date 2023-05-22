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
//! \\file FixedDensityD3Q19.cpp
//! \\author pystencils
//======================================================================================================================

#include "core/DataTypes.h"
#include "core/Macros.h"
#include "FixedDensityD3Q19.h"



#define FUNC_PREFIX

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
#ifdef __NVCC_DIAG_PRAGMA_SUPPORT__
#pragma nv_diag_suppress 177
#else
#pragma diag_suppress 177
#endif
#endif

namespace internal_fixeddensityd3q19_even {
static FUNC_PREFIX void fixeddensityd3q19_even(const uint8_t * RESTRICT const _data_indexVector, double * RESTRICT  _data_pdfs, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2, int64_t const _stride_pdfs_3, double density, int32_t indexVectorSize)
{
   
   const int32_t f_in_inv_dir_idx [] = { 0,2,1,4,3,6,5,10,9,8,7,16,15,18,17,12,11,14,13 }; 
   const int32_t f_in_inv_offsets_x [] = { 0,0,0,-1,1,0,0,-1,1,-1,1,0,0,-1,1,0,0,-1,1 }; 
   const int32_t f_in_inv_offsets_y [] = { 0,1,-1,0,0,0,0,1,1,-1,-1,1,-1,0,0,1,-1,0,0 }; 
   const int32_t f_in_inv_offsets_z [] = { 0,0,0,0,0,1,-1,0,0,0,0,1,1,1,1,-1,-1,-1,-1 }; 
   
   const double rho = density;
   const double delta_rho = rho - 1.0;
   for (int64_t ctr_0 = 0; ctr_0 < indexVectorSize; ctr_0 += 1)
   {
      const int32_t x = *((int32_t * )(& _data_indexVector[16*ctr_0]));
      const int32_t y = *((int32_t * )(& _data_indexVector[16*ctr_0 + 4]));
      const int32_t z = *((int32_t * )(& _data_indexVector[16*ctr_0 + 8]));
      const int32_t dir = *((int32_t * )(& _data_indexVector[16*ctr_0 + 12]));
      const double vel0Term = _data_pdfs[_stride_pdfs_0*x + _stride_pdfs_1*y + _stride_pdfs_2*z + 10*_stride_pdfs_3] + _data_pdfs[_stride_pdfs_0*x + _stride_pdfs_1*y + _stride_pdfs_2*z + 14*_stride_pdfs_3] + _data_pdfs[_stride_pdfs_0*x + _stride_pdfs_1*y + _stride_pdfs_2*z + 18*_stride_pdfs_3] + _data_pdfs[_stride_pdfs_0*x + _stride_pdfs_1*y + _stride_pdfs_2*z + 4*_stride_pdfs_3] + _data_pdfs[_stride_pdfs_0*x + _stride_pdfs_1*y + _stride_pdfs_2*z + 8*_stride_pdfs_3];
      const double vel1Term = _data_pdfs[_stride_pdfs_0*x + _stride_pdfs_1*y + _stride_pdfs_2*z + 11*_stride_pdfs_3] + _data_pdfs[_stride_pdfs_0*x + _stride_pdfs_1*y + _stride_pdfs_2*z + 15*_stride_pdfs_3] + _data_pdfs[_stride_pdfs_0*x + _stride_pdfs_1*y + _stride_pdfs_2*z + 7*_stride_pdfs_3] + _data_pdfs[_stride_pdfs_0*x + _stride_pdfs_1*y + _stride_pdfs_2*z + _stride_pdfs_3];
      const double vel2Term = _data_pdfs[_stride_pdfs_0*x + _stride_pdfs_1*y + _stride_pdfs_2*z + 12*_stride_pdfs_3] + _data_pdfs[_stride_pdfs_0*x + _stride_pdfs_1*y + _stride_pdfs_2*z + 13*_stride_pdfs_3] + _data_pdfs[_stride_pdfs_0*x + _stride_pdfs_1*y + _stride_pdfs_2*z + 5*_stride_pdfs_3];
      const double u_0 = vel0Term - 1.0*_data_pdfs[_stride_pdfs_0*x + _stride_pdfs_1*y + _stride_pdfs_2*z + 13*_stride_pdfs_3] - 1.0*_data_pdfs[_stride_pdfs_0*x + _stride_pdfs_1*y + _stride_pdfs_2*z + 17*_stride_pdfs_3] - 1.0*_data_pdfs[_stride_pdfs_0*x + _stride_pdfs_1*y + _stride_pdfs_2*z + 3*_stride_pdfs_3] - 1.0*_data_pdfs[_stride_pdfs_0*x + _stride_pdfs_1*y + _stride_pdfs_2*z + 7*_stride_pdfs_3] - 1.0*_data_pdfs[_stride_pdfs_0*x + _stride_pdfs_1*y + _stride_pdfs_2*z + 9*_stride_pdfs_3];
      const double u_1 = vel1Term - 1.0*_data_pdfs[_stride_pdfs_0*x + _stride_pdfs_1*y + _stride_pdfs_2*z + 10*_stride_pdfs_3] - 1.0*_data_pdfs[_stride_pdfs_0*x + _stride_pdfs_1*y + _stride_pdfs_2*z + 12*_stride_pdfs_3] - 1.0*_data_pdfs[_stride_pdfs_0*x + _stride_pdfs_1*y + _stride_pdfs_2*z + 16*_stride_pdfs_3] - 1.0*_data_pdfs[_stride_pdfs_0*x + _stride_pdfs_1*y + _stride_pdfs_2*z + 2*_stride_pdfs_3] - 1.0*_data_pdfs[_stride_pdfs_0*x + _stride_pdfs_1*y + _stride_pdfs_2*z + 9*_stride_pdfs_3] + _data_pdfs[_stride_pdfs_0*x + _stride_pdfs_1*y + _stride_pdfs_2*z + 8*_stride_pdfs_3];
      const double u_2 = vel2Term - 1.0*_data_pdfs[_stride_pdfs_0*x + _stride_pdfs_1*y + _stride_pdfs_2*z + 15*_stride_pdfs_3] - 1.0*_data_pdfs[_stride_pdfs_0*x + _stride_pdfs_1*y + _stride_pdfs_2*z + 16*_stride_pdfs_3] - 1.0*_data_pdfs[_stride_pdfs_0*x + _stride_pdfs_1*y + _stride_pdfs_2*z + 17*_stride_pdfs_3] - 1.0*_data_pdfs[_stride_pdfs_0*x + _stride_pdfs_1*y + _stride_pdfs_2*z + 18*_stride_pdfs_3] - 1.0*_data_pdfs[_stride_pdfs_0*x + _stride_pdfs_1*y + _stride_pdfs_2*z + 6*_stride_pdfs_3] + _data_pdfs[_stride_pdfs_0*x + _stride_pdfs_1*y + _stride_pdfs_2*z + 11*_stride_pdfs_3] + _data_pdfs[_stride_pdfs_0*x + _stride_pdfs_1*y + _stride_pdfs_2*z + 14*_stride_pdfs_3];
      const double u0Mu1 = u_0 + u_1*-1.0;
      const double u0Pu1 = u_0 + u_1;
      const double u1Pu2 = u_1 + u_2;
      const double u1Mu2 = u_1 + u_2*-1.0;
      const double u0Mu2 = u_0 + u_2*-1.0;
      const double u0Pu2 = u_0 + u_2;
      const double f_eq_common = delta_rho - 1.0*(u_0*u_0) - 1.0*(u_1*u_1) - 1.0*(u_2*u_2);
      _data_pdfs[_stride_pdfs_0*x + _stride_pdfs_0*f_in_inv_offsets_x[dir] + _stride_pdfs_1*y + _stride_pdfs_1*f_in_inv_offsets_y[dir] + _stride_pdfs_2*z + _stride_pdfs_2*f_in_inv_offsets_z[dir] + _stride_pdfs_3*f_in_inv_dir_idx[dir]] = -1.0*_data_pdfs[_stride_pdfs_0*x + _stride_pdfs_1*y + _stride_pdfs_2*z + _stride_pdfs_3*dir] + 2.0*((((dir) == (0))) ? (f_eq_common*0.33333333333333331): ((((dir) == (1)) || ((dir) == (2))) ? (delta_rho*-0.1111111111111111 + f_eq_common*0.16666666666666666 + 0.33333333333333331*(u_1*u_1)): ((((dir) == (3)) || ((dir) == (4))) ? (delta_rho*-0.1111111111111111 + f_eq_common*0.16666666666666666 + 0.33333333333333331*(u_0*u_0)): ((((dir) == (5)) || ((dir) == (6))) ? (delta_rho*-0.1111111111111111 + f_eq_common*0.16666666666666666 + 0.33333333333333331*(u_2*u_2)): ((((dir) == (7))) ? (delta_rho*-0.013888888888888888 + f_eq_common*0.041666666666666664 + 0.041666666666666664*(u_2*u_2) + 0.125*(u0Mu1*u0Mu1)): ((((dir) == (8)) || ((dir) == (9))) ? (delta_rho*-0.013888888888888888 + f_eq_common*0.041666666666666664 + 0.041666666666666664*(u_2*u_2) + 0.125*(u0Pu1*u0Pu1)): ((((dir) == (10))) ? (delta_rho*-0.013888888888888888 + f_eq_common*0.041666666666666664 + 0.041666666666666664*(u_2*u_2) + 0.125*(u0Mu1*u0Mu1)): ((((dir) == (11))) ? (delta_rho*-0.013888888888888888 + f_eq_common*0.041666666666666664 + 0.041666666666666664*(u_0*u_0) + 0.125*(u1Pu2*u1Pu2)): ((((dir) == (12))) ? (delta_rho*-0.013888888888888888 + f_eq_common*0.041666666666666664 + 0.041666666666666664*(u_0*u_0) + 0.125*(u1Mu2*u1Mu2)): ((((dir) == (13))) ? (delta_rho*-0.013888888888888888 + f_eq_common*0.041666666666666664 + 0.041666666666666664*(u_1*u_1) + 0.125*(u0Mu2*u0Mu2)): ((((dir) == (14))) ? (delta_rho*-0.013888888888888888 + f_eq_common*0.041666666666666664 + 0.041666666666666664*(u_1*u_1) + 0.125*(u0Pu2*u0Pu2)): ((((dir) == (15))) ? (delta_rho*-0.013888888888888888 + f_eq_common*0.041666666666666664 + 0.041666666666666664*(u_0*u_0) + 0.125*(u1Mu2*u1Mu2)): ((((dir) == (16))) ? (delta_rho*-0.013888888888888888 + f_eq_common*0.041666666666666664 + 0.041666666666666664*(u_0*u_0) + 0.125*(u1Pu2*u1Pu2)): ((((dir) == (17))) ? (delta_rho*-0.013888888888888888 + f_eq_common*0.041666666666666664 + 0.041666666666666664*(u_1*u_1) + 0.125*(u0Pu2*u0Pu2)): ((((dir) == (18))) ? (delta_rho*-0.013888888888888888 + f_eq_common*0.041666666666666664 + 0.041666666666666664*(u_1*u_1) + 0.125*(u0Mu2*u0Mu2)): (0.0))))))))))))))));
   }
}
}


#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif

#ifdef __CUDACC__
#pragma pop
#endif


void FixedDensityD3Q19::run_impl(IBlock * block, IndexVectors::Type type)
{
   auto * indexVectors = block->getData<IndexVectors>(indexVectorID);
   int32_t indexVectorSize = int32_c( indexVectors->indexVector(type).size() );
   if( indexVectorSize == 0)
      return;

   
   auto pointer = indexVectors->pointerCpu(type);
   

   uint8_t * _data_indexVector = reinterpret_cast<uint8_t*>(pointer);

   auto pdfs = block->getData< field::GhostLayerField<double, 19> >(pdfsID);

   uint8_t timestep = pdfs->getTimestep();
   auto & density = density_;
   WALBERLA_ASSERT_GREATER_EQUAL(0, -int_c(pdfs->nrOfGhostLayers()))
    double * RESTRICT  _data_pdfs = pdfs->dataAt(0, 0, 0, 0);
    const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
    const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
    const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
    const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
    if(((timestep & 1) ^ 1)) {
        internal_fixeddensityd3q19_even::fixeddensityd3q19_even(_data_indexVector, _data_pdfs, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3, density, indexVectorSize);
    } else {
        internal_fixeddensityd3q19_even::fixeddensityd3q19_even(_data_indexVector, _data_pdfs, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3, density, indexVectorSize);
    }
}

void FixedDensityD3Q19::run(IBlock * block)
{
   run_impl(block, IndexVectors::ALL);
}

void FixedDensityD3Q19::inner(IBlock * block)
{
   run_impl(block, IndexVectors::INNER);
}

void FixedDensityD3Q19::outer(IBlock * block)
{
   run_impl(block, IndexVectors::OUTER);
}

} // namespace lbm
} // namespace walberla

