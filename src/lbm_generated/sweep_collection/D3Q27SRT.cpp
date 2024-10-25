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
//! \\file D3Q27SRT.cpp
//! \\author pystencils
//======================================================================================================================
#include "D3Q27SRT.h"

#define FUNC_PREFIX

#if ( defined WALBERLA_CXX_COMPILER_IS_GNU ) || ( defined WALBERLA_CXX_COMPILER_IS_CLANG )
#   pragma GCC diagnostic push
#   pragma GCC diagnostic ignored "-Wunused-variable"
#endif

#if ( defined WALBERLA_CXX_COMPILER_IS_INTEL )
#pragma warning push
#pragma warning( disable :  1599 )
#endif

using namespace std;

namespace walberla {
namespace lbm {


namespace internal_d3q27srt_kernel_streamCollide {
static FUNC_PREFIX void d3q27srt_kernel_streamCollide(double * RESTRICT const _data_pdfs, double * RESTRICT  _data_pdfs_tmp, int64_t const _size_pdfs_0, int64_t const _size_pdfs_1, int64_t const _size_pdfs_2, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2, int64_t const _stride_pdfs_3, int64_t const _stride_pdfs_tmp_0, int64_t const _stride_pdfs_tmp_1, int64_t const _stride_pdfs_tmp_2, int64_t const _stride_pdfs_tmp_3, double omega)
{
   for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_2; ctr_2 += 1)
   {
      for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_1; ctr_1 += 1)
      {
         for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_0; ctr_0 += 1)
         {
            const double vel0Term = _data_pdfs[_stride_pdfs_0*ctr_0 - _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _stride_pdfs_2*ctr_2 + 10*_stride_pdfs_3] + _data_pdfs[_stride_pdfs_0*ctr_0 - _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 25*_stride_pdfs_3] + _data_pdfs[_stride_pdfs_0*ctr_0 - _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 21*_stride_pdfs_3] + _data_pdfs[_stride_pdfs_0*ctr_0 - _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 4*_stride_pdfs_3] + _data_pdfs[_stride_pdfs_0*ctr_0 - _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 18*_stride_pdfs_3] + _data_pdfs[_stride_pdfs_0*ctr_0 - _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 14*_stride_pdfs_3] + _data_pdfs[_stride_pdfs_0*ctr_0 - _stride_pdfs_0 + _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _stride_pdfs_2*ctr_2 + 8*_stride_pdfs_3] + _data_pdfs[_stride_pdfs_0*ctr_0 - _stride_pdfs_0 + _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 23*_stride_pdfs_3] + _data_pdfs[_stride_pdfs_0*ctr_0 - _stride_pdfs_0 + _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 19*_stride_pdfs_3];
            const double vel1Term = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_0 + _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _stride_pdfs_2*ctr_2 + 7*_stride_pdfs_3] + _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_0 + _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 24*_stride_pdfs_3] + _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_0 + _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 20*_stride_pdfs_3] + _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 15*_stride_pdfs_3] + _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _stride_pdfs_2*ctr_2 + _stride_pdfs_3] + _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 11*_stride_pdfs_3];
            const double vel2Term = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 22*_stride_pdfs_3] + _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 13*_stride_pdfs_3] + _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 12*_stride_pdfs_3] + _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 5*_stride_pdfs_3];
            const double delta_rho = vel0Term + vel1Term + vel2Term + _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _stride_pdfs_2*ctr_2 + 9*_stride_pdfs_3] + _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 26*_stride_pdfs_3] + _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 3*_stride_pdfs_3] + _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 17*_stride_pdfs_3] + _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _stride_pdfs_2*ctr_2 + 2*_stride_pdfs_3] + _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 16*_stride_pdfs_3] + _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 6*_stride_pdfs_3] + _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2];
            const double u_0 = vel0Term - _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _stride_pdfs_2*ctr_2 + 9*_stride_pdfs_3] - _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 26*_stride_pdfs_3] - _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 22*_stride_pdfs_3] - _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 3*_stride_pdfs_3] - _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 17*_stride_pdfs_3] - _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 13*_stride_pdfs_3] - _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_0 + _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _stride_pdfs_2*ctr_2 + 7*_stride_pdfs_3] - _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_0 + _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 24*_stride_pdfs_3] - _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_0 + _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 20*_stride_pdfs_3];
            const double u_1 = vel1Term - _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _stride_pdfs_2*ctr_2 + 9*_stride_pdfs_3] - _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 26*_stride_pdfs_3] - _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 22*_stride_pdfs_3] - _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _stride_pdfs_2*ctr_2 + 2*_stride_pdfs_3] - _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 16*_stride_pdfs_3] - _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 12*_stride_pdfs_3] - _data_pdfs[_stride_pdfs_0*ctr_0 - _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _stride_pdfs_2*ctr_2 + 10*_stride_pdfs_3] - _data_pdfs[_stride_pdfs_0*ctr_0 - _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 25*_stride_pdfs_3] - _data_pdfs[_stride_pdfs_0*ctr_0 - _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 21*_stride_pdfs_3] + _data_pdfs[_stride_pdfs_0*ctr_0 - _stride_pdfs_0 + _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _stride_pdfs_2*ctr_2 + 8*_stride_pdfs_3] + _data_pdfs[_stride_pdfs_0*ctr_0 - _stride_pdfs_0 + _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 23*_stride_pdfs_3] + _data_pdfs[_stride_pdfs_0*ctr_0 - _stride_pdfs_0 + _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 19*_stride_pdfs_3];
            const double u_2 = vel2Term - _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 26*_stride_pdfs_3] - _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 17*_stride_pdfs_3] - _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_0 + _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 24*_stride_pdfs_3] + _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_0 + _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 20*_stride_pdfs_3] - _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 16*_stride_pdfs_3] - _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 6*_stride_pdfs_3] - _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 15*_stride_pdfs_3] + _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 11*_stride_pdfs_3] - _data_pdfs[_stride_pdfs_0*ctr_0 - _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 25*_stride_pdfs_3] + _data_pdfs[_stride_pdfs_0*ctr_0 - _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 21*_stride_pdfs_3] - _data_pdfs[_stride_pdfs_0*ctr_0 - _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 18*_stride_pdfs_3] + _data_pdfs[_stride_pdfs_0*ctr_0 - _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 14*_stride_pdfs_3] - _data_pdfs[_stride_pdfs_0*ctr_0 - _stride_pdfs_0 + _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 23*_stride_pdfs_3] + _data_pdfs[_stride_pdfs_0*ctr_0 - _stride_pdfs_0 + _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 19*_stride_pdfs_3];
            const double u0Mu1 = u_0 - u_1;
            const double u0Pu1 = u_0 + u_1;
            const double u1Pu2 = u_1 + u_2;
            const double u1Mu2 = u_1 - u_2;
            const double u0Mu2 = u_0 - u_2;
            const double u0Pu2 = u_0 + u_2;
            const double f_eq_common = delta_rho - 1.5*(u_0*u_0) - 1.5*(u_1*u_1) - 1.5*(u_2*u_2);
            _data_pdfs_tmp[_stride_pdfs_tmp_0*ctr_0 + _stride_pdfs_tmp_1*ctr_1 + _stride_pdfs_tmp_2*ctr_2] = omega*(f_eq_common*0.29629629629629628 - _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2]) + _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2];
            _data_pdfs_tmp[_stride_pdfs_tmp_0*ctr_0 + _stride_pdfs_tmp_1*ctr_1 + _stride_pdfs_tmp_2*ctr_2 + _stride_pdfs_tmp_3] = omega*(f_eq_common*0.07407407407407407 + u_1*0.22222222222222221 + 0.33333333333333331*(u_1*u_1) - _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _stride_pdfs_2*ctr_2 + _stride_pdfs_3]) + _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _stride_pdfs_2*ctr_2 + _stride_pdfs_3];
            _data_pdfs_tmp[_stride_pdfs_tmp_0*ctr_0 + _stride_pdfs_tmp_1*ctr_1 + _stride_pdfs_tmp_2*ctr_2 + 2*_stride_pdfs_tmp_3] = omega*(f_eq_common*0.07407407407407407 + u_1*-0.22222222222222221 + 0.33333333333333331*(u_1*u_1) - _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _stride_pdfs_2*ctr_2 + 2*_stride_pdfs_3]) + _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _stride_pdfs_2*ctr_2 + 2*_stride_pdfs_3];
            _data_pdfs_tmp[_stride_pdfs_tmp_0*ctr_0 + _stride_pdfs_tmp_1*ctr_1 + _stride_pdfs_tmp_2*ctr_2 + 3*_stride_pdfs_tmp_3] = omega*(f_eq_common*0.07407407407407407 + u_0*-0.22222222222222221 + 0.33333333333333331*(u_0*u_0) - _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 3*_stride_pdfs_3]) + _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 3*_stride_pdfs_3];
            _data_pdfs_tmp[_stride_pdfs_tmp_0*ctr_0 + _stride_pdfs_tmp_1*ctr_1 + _stride_pdfs_tmp_2*ctr_2 + 4*_stride_pdfs_tmp_3] = omega*(f_eq_common*0.07407407407407407 + u_0*0.22222222222222221 + 0.33333333333333331*(u_0*u_0) - _data_pdfs[_stride_pdfs_0*ctr_0 - _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 4*_stride_pdfs_3]) + _data_pdfs[_stride_pdfs_0*ctr_0 - _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 4*_stride_pdfs_3];
            _data_pdfs_tmp[_stride_pdfs_tmp_0*ctr_0 + _stride_pdfs_tmp_1*ctr_1 + _stride_pdfs_tmp_2*ctr_2 + 5*_stride_pdfs_tmp_3] = omega*(f_eq_common*0.07407407407407407 + u_2*0.22222222222222221 + 0.33333333333333331*(u_2*u_2) - _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 5*_stride_pdfs_3]) + _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 5*_stride_pdfs_3];
            _data_pdfs_tmp[_stride_pdfs_tmp_0*ctr_0 + _stride_pdfs_tmp_1*ctr_1 + _stride_pdfs_tmp_2*ctr_2 + 6*_stride_pdfs_tmp_3] = omega*(f_eq_common*0.07407407407407407 + u_2*-0.22222222222222221 + 0.33333333333333331*(u_2*u_2) - _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 6*_stride_pdfs_3]) + _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 6*_stride_pdfs_3];
            _data_pdfs_tmp[_stride_pdfs_tmp_0*ctr_0 + _stride_pdfs_tmp_1*ctr_1 + _stride_pdfs_tmp_2*ctr_2 + 7*_stride_pdfs_tmp_3] = omega*(f_eq_common*0.018518518518518517 + u0Mu1*-0.055555555555555552 + 0.083333333333333329*(u0Mu1*u0Mu1) - _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_0 + _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _stride_pdfs_2*ctr_2 + 7*_stride_pdfs_3]) + _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_0 + _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _stride_pdfs_2*ctr_2 + 7*_stride_pdfs_3];
            _data_pdfs_tmp[_stride_pdfs_tmp_0*ctr_0 + _stride_pdfs_tmp_1*ctr_1 + _stride_pdfs_tmp_2*ctr_2 + 8*_stride_pdfs_tmp_3] = omega*(f_eq_common*0.018518518518518517 + u0Pu1*0.055555555555555552 + 0.083333333333333329*(u0Pu1*u0Pu1) - _data_pdfs[_stride_pdfs_0*ctr_0 - _stride_pdfs_0 + _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _stride_pdfs_2*ctr_2 + 8*_stride_pdfs_3]) + _data_pdfs[_stride_pdfs_0*ctr_0 - _stride_pdfs_0 + _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _stride_pdfs_2*ctr_2 + 8*_stride_pdfs_3];
            _data_pdfs_tmp[_stride_pdfs_tmp_0*ctr_0 + _stride_pdfs_tmp_1*ctr_1 + _stride_pdfs_tmp_2*ctr_2 + 9*_stride_pdfs_tmp_3] = omega*(f_eq_common*0.018518518518518517 + u0Pu1*-0.055555555555555552 + 0.083333333333333329*(u0Pu1*u0Pu1) - _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _stride_pdfs_2*ctr_2 + 9*_stride_pdfs_3]) + _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _stride_pdfs_2*ctr_2 + 9*_stride_pdfs_3];
            _data_pdfs_tmp[_stride_pdfs_tmp_0*ctr_0 + _stride_pdfs_tmp_1*ctr_1 + _stride_pdfs_tmp_2*ctr_2 + 10*_stride_pdfs_tmp_3] = omega*(f_eq_common*0.018518518518518517 + u0Mu1*0.055555555555555552 + 0.083333333333333329*(u0Mu1*u0Mu1) - _data_pdfs[_stride_pdfs_0*ctr_0 - _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _stride_pdfs_2*ctr_2 + 10*_stride_pdfs_3]) + _data_pdfs[_stride_pdfs_0*ctr_0 - _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _stride_pdfs_2*ctr_2 + 10*_stride_pdfs_3];
            _data_pdfs_tmp[_stride_pdfs_tmp_0*ctr_0 + _stride_pdfs_tmp_1*ctr_1 + _stride_pdfs_tmp_2*ctr_2 + 11*_stride_pdfs_tmp_3] = omega*(f_eq_common*0.018518518518518517 + u1Pu2*0.055555555555555552 + 0.083333333333333329*(u1Pu2*u1Pu2) - _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 11*_stride_pdfs_3]) + _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 11*_stride_pdfs_3];
            _data_pdfs_tmp[_stride_pdfs_tmp_0*ctr_0 + _stride_pdfs_tmp_1*ctr_1 + _stride_pdfs_tmp_2*ctr_2 + 12*_stride_pdfs_tmp_3] = omega*(f_eq_common*0.018518518518518517 + u1Mu2*-0.055555555555555552 + 0.083333333333333329*(u1Mu2*u1Mu2) - _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 12*_stride_pdfs_3]) + _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 12*_stride_pdfs_3];
            _data_pdfs_tmp[_stride_pdfs_tmp_0*ctr_0 + _stride_pdfs_tmp_1*ctr_1 + _stride_pdfs_tmp_2*ctr_2 + 13*_stride_pdfs_tmp_3] = omega*(f_eq_common*0.018518518518518517 + u0Mu2*-0.055555555555555552 + 0.083333333333333329*(u0Mu2*u0Mu2) - _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 13*_stride_pdfs_3]) + _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 13*_stride_pdfs_3];
            _data_pdfs_tmp[_stride_pdfs_tmp_0*ctr_0 + _stride_pdfs_tmp_1*ctr_1 + _stride_pdfs_tmp_2*ctr_2 + 14*_stride_pdfs_tmp_3] = omega*(f_eq_common*0.018518518518518517 + u0Pu2*0.055555555555555552 + 0.083333333333333329*(u0Pu2*u0Pu2) - _data_pdfs[_stride_pdfs_0*ctr_0 - _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 14*_stride_pdfs_3]) + _data_pdfs[_stride_pdfs_0*ctr_0 - _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 14*_stride_pdfs_3];
            _data_pdfs_tmp[_stride_pdfs_tmp_0*ctr_0 + _stride_pdfs_tmp_1*ctr_1 + _stride_pdfs_tmp_2*ctr_2 + 15*_stride_pdfs_tmp_3] = omega*(f_eq_common*0.018518518518518517 + u1Mu2*0.055555555555555552 + 0.083333333333333329*(u1Mu2*u1Mu2) - _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 15*_stride_pdfs_3]) + _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 15*_stride_pdfs_3];
            _data_pdfs_tmp[_stride_pdfs_tmp_0*ctr_0 + _stride_pdfs_tmp_1*ctr_1 + _stride_pdfs_tmp_2*ctr_2 + 16*_stride_pdfs_tmp_3] = omega*(f_eq_common*0.018518518518518517 + u1Pu2*-0.055555555555555552 + 0.083333333333333329*(u1Pu2*u1Pu2) - _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 16*_stride_pdfs_3]) + _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 16*_stride_pdfs_3];
            _data_pdfs_tmp[_stride_pdfs_tmp_0*ctr_0 + _stride_pdfs_tmp_1*ctr_1 + _stride_pdfs_tmp_2*ctr_2 + 17*_stride_pdfs_tmp_3] = omega*(f_eq_common*0.018518518518518517 + u0Pu2*-0.055555555555555552 + 0.083333333333333329*(u0Pu2*u0Pu2) - _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 17*_stride_pdfs_3]) + _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 17*_stride_pdfs_3];
            _data_pdfs_tmp[_stride_pdfs_tmp_0*ctr_0 + _stride_pdfs_tmp_1*ctr_1 + _stride_pdfs_tmp_2*ctr_2 + 18*_stride_pdfs_tmp_3] = omega*(f_eq_common*0.018518518518518517 + u0Mu2*0.055555555555555552 + 0.083333333333333329*(u0Mu2*u0Mu2) - _data_pdfs[_stride_pdfs_0*ctr_0 - _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 18*_stride_pdfs_3]) + _data_pdfs[_stride_pdfs_0*ctr_0 - _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 18*_stride_pdfs_3];
            _data_pdfs_tmp[_stride_pdfs_tmp_0*ctr_0 + _stride_pdfs_tmp_1*ctr_1 + _stride_pdfs_tmp_2*ctr_2 + 19*_stride_pdfs_tmp_3] = omega*(delta_rho*-0.013888888888888888 + f_eq_common*0.018518518518518517 + u0Pu1*0.013888888888888888 + u_2*0.013888888888888888 + 0.020833333333333332*(u0Pu1*u0Pu1) + 0.020833333333333332*(u0Pu2*u0Pu2) + 0.020833333333333332*(u1Pu2*u1Pu2) - _data_pdfs[_stride_pdfs_0*ctr_0 - _stride_pdfs_0 + _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 19*_stride_pdfs_3]) + _data_pdfs[_stride_pdfs_0*ctr_0 - _stride_pdfs_0 + _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 19*_stride_pdfs_3];
            _data_pdfs_tmp[_stride_pdfs_tmp_0*ctr_0 + _stride_pdfs_tmp_1*ctr_1 + _stride_pdfs_tmp_2*ctr_2 + 20*_stride_pdfs_tmp_3] = omega*(delta_rho*-0.013888888888888888 + f_eq_common*0.018518518518518517 + u0Mu1*-0.013888888888888888 + u_2*0.013888888888888888 + 0.020833333333333332*(u0Mu1*u0Mu1) + 0.020833333333333332*(u0Mu2*u0Mu2) + 0.020833333333333332*(u1Pu2*u1Pu2) - _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_0 + _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 20*_stride_pdfs_3]) + _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_0 + _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 20*_stride_pdfs_3];
            _data_pdfs_tmp[_stride_pdfs_tmp_0*ctr_0 + _stride_pdfs_tmp_1*ctr_1 + _stride_pdfs_tmp_2*ctr_2 + 21*_stride_pdfs_tmp_3] = omega*(delta_rho*-0.013888888888888888 + f_eq_common*0.018518518518518517 + u0Mu1*0.013888888888888888 + u_2*0.013888888888888888 + 0.020833333333333332*(u0Mu1*u0Mu1) + 0.020833333333333332*(u0Pu2*u0Pu2) + 0.020833333333333332*(u1Mu2*u1Mu2) - _data_pdfs[_stride_pdfs_0*ctr_0 - _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 21*_stride_pdfs_3]) + _data_pdfs[_stride_pdfs_0*ctr_0 - _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 21*_stride_pdfs_3];
            _data_pdfs_tmp[_stride_pdfs_tmp_0*ctr_0 + _stride_pdfs_tmp_1*ctr_1 + _stride_pdfs_tmp_2*ctr_2 + 22*_stride_pdfs_tmp_3] = omega*(delta_rho*-0.013888888888888888 + f_eq_common*0.018518518518518517 + u0Pu1*-0.013888888888888888 + u_2*0.013888888888888888 + 0.020833333333333332*(u0Mu2*u0Mu2) + 0.020833333333333332*(u0Pu1*u0Pu1) + 0.020833333333333332*(u1Mu2*u1Mu2) - _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 22*_stride_pdfs_3]) + _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 22*_stride_pdfs_3];
            _data_pdfs_tmp[_stride_pdfs_tmp_0*ctr_0 + _stride_pdfs_tmp_1*ctr_1 + _stride_pdfs_tmp_2*ctr_2 + 23*_stride_pdfs_tmp_3] = omega*(delta_rho*-0.013888888888888888 + f_eq_common*0.018518518518518517 + u0Pu1*0.013888888888888888 + u_2*-0.013888888888888888 + 0.020833333333333332*(u0Mu2*u0Mu2) + 0.020833333333333332*(u0Pu1*u0Pu1) + 0.020833333333333332*(u1Mu2*u1Mu2) - _data_pdfs[_stride_pdfs_0*ctr_0 - _stride_pdfs_0 + _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 23*_stride_pdfs_3]) + _data_pdfs[_stride_pdfs_0*ctr_0 - _stride_pdfs_0 + _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 23*_stride_pdfs_3];
            _data_pdfs_tmp[_stride_pdfs_tmp_0*ctr_0 + _stride_pdfs_tmp_1*ctr_1 + _stride_pdfs_tmp_2*ctr_2 + 24*_stride_pdfs_tmp_3] = omega*(delta_rho*-0.013888888888888888 + f_eq_common*0.018518518518518517 + u0Mu1*-0.013888888888888888 + u_2*-0.013888888888888888 + 0.020833333333333332*(u0Mu1*u0Mu1) + 0.020833333333333332*(u0Pu2*u0Pu2) + 0.020833333333333332*(u1Mu2*u1Mu2) - _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_0 + _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 24*_stride_pdfs_3]) + _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_0 + _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 24*_stride_pdfs_3];
            _data_pdfs_tmp[_stride_pdfs_tmp_0*ctr_0 + _stride_pdfs_tmp_1*ctr_1 + _stride_pdfs_tmp_2*ctr_2 + 25*_stride_pdfs_tmp_3] = omega*(delta_rho*-0.013888888888888888 + f_eq_common*0.018518518518518517 + u0Mu1*0.013888888888888888 + u_2*-0.013888888888888888 + 0.020833333333333332*(u0Mu1*u0Mu1) + 0.020833333333333332*(u0Mu2*u0Mu2) + 0.020833333333333332*(u1Pu2*u1Pu2) - _data_pdfs[_stride_pdfs_0*ctr_0 - _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 25*_stride_pdfs_3]) + _data_pdfs[_stride_pdfs_0*ctr_0 - _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 25*_stride_pdfs_3];
            _data_pdfs_tmp[_stride_pdfs_tmp_0*ctr_0 + _stride_pdfs_tmp_1*ctr_1 + _stride_pdfs_tmp_2*ctr_2 + 26*_stride_pdfs_tmp_3] = omega*(delta_rho*-0.013888888888888888 + f_eq_common*0.018518518518518517 + u0Pu1*-0.013888888888888888 + u_2*-0.013888888888888888 + 0.020833333333333332*(u0Pu1*u0Pu1) + 0.020833333333333332*(u0Pu2*u0Pu2) + 0.020833333333333332*(u1Pu2*u1Pu2) - _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 26*_stride_pdfs_3]) + _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 26*_stride_pdfs_3];
         }
      }
   }
}
}


namespace internal_d3q27srt_kernel_collide {
static FUNC_PREFIX void d3q27srt_kernel_collide(double * RESTRICT  _data_pdfs, int64_t const _size_pdfs_0, int64_t const _size_pdfs_1, int64_t const _size_pdfs_2, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2, int64_t const _stride_pdfs_3, double omega)
{
   for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_2; ctr_2 += 1)
   {
      for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_1; ctr_1 += 1)
      {
         for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_0; ctr_0 += 1)
         {
            const double xi_1 = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 23*_stride_pdfs_3];
            const double xi_2 = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 16*_stride_pdfs_3];
            const double xi_3 = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2];
            const double xi_4 = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + _stride_pdfs_3];
            const double xi_5 = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 13*_stride_pdfs_3];
            const double xi_6 = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 5*_stride_pdfs_3];
            const double xi_7 = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 12*_stride_pdfs_3];
            const double xi_8 = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 26*_stride_pdfs_3];
            const double xi_9 = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 11*_stride_pdfs_3];
            const double xi_10 = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 19*_stride_pdfs_3];
            const double xi_11 = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 17*_stride_pdfs_3];
            const double xi_12 = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 20*_stride_pdfs_3];
            const double xi_13 = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 18*_stride_pdfs_3];
            const double xi_14 = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 4*_stride_pdfs_3];
            const double xi_15 = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 22*_stride_pdfs_3];
            const double xi_16 = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 25*_stride_pdfs_3];
            const double xi_17 = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 14*_stride_pdfs_3];
            const double xi_18 = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 6*_stride_pdfs_3];
            const double xi_19 = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 8*_stride_pdfs_3];
            const double xi_20 = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 9*_stride_pdfs_3];
            const double xi_21 = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 2*_stride_pdfs_3];
            const double xi_22 = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 15*_stride_pdfs_3];
            const double xi_23 = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 24*_stride_pdfs_3];
            const double xi_24 = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 21*_stride_pdfs_3];
            const double xi_25 = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 10*_stride_pdfs_3];
            const double xi_26 = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 3*_stride_pdfs_3];
            const double xi_27 = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 7*_stride_pdfs_3];
            const double vel0Term = xi_1 + xi_10 + xi_13 + xi_14 + xi_16 + xi_17 + xi_19 + xi_24 + xi_25;
            const double vel1Term = xi_12 + xi_22 + xi_23 + xi_27 + xi_4 + xi_9;
            const double vel2Term = xi_15 + xi_5 + xi_6 + xi_7;
            const double delta_rho = vel0Term + vel1Term + vel2Term + xi_11 + xi_18 + xi_2 + xi_20 + xi_21 + xi_26 + xi_3 + xi_8;
            const double u_0 = vel0Term - xi_11 - xi_12 - xi_15 - xi_20 - xi_23 - xi_26 - xi_27 - xi_5 - xi_8;
            const double u_1 = vel1Term + xi_1 + xi_10 - xi_15 - xi_16 + xi_19 - xi_2 - xi_20 - xi_21 - xi_24 - xi_25 - xi_7 - xi_8;
            const double u_2 = vel2Term - xi_1 + xi_10 - xi_11 + xi_12 - xi_13 - xi_16 + xi_17 - xi_18 - xi_2 - xi_22 - xi_23 + xi_24 - xi_8 + xi_9;
            const double u0Mu1 = u_0 - u_1;
            const double u0Pu1 = u_0 + u_1;
            const double u1Pu2 = u_1 + u_2;
            const double u1Mu2 = u_1 - u_2;
            const double u0Mu2 = u_0 - u_2;
            const double u0Pu2 = u_0 + u_2;
            const double f_eq_common = delta_rho - 1.5*(u_0*u_0) - 1.5*(u_1*u_1) - 1.5*(u_2*u_2);
            _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2] = omega*(f_eq_common*0.29629629629629628 - xi_3) + xi_3;
            _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + _stride_pdfs_3] = omega*(f_eq_common*0.07407407407407407 + u_1*0.22222222222222221 - xi_4 + 0.33333333333333331*(u_1*u_1)) + xi_4;
            _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 2*_stride_pdfs_3] = omega*(f_eq_common*0.07407407407407407 + u_1*-0.22222222222222221 - xi_21 + 0.33333333333333331*(u_1*u_1)) + xi_21;
            _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 3*_stride_pdfs_3] = omega*(f_eq_common*0.07407407407407407 + u_0*-0.22222222222222221 - xi_26 + 0.33333333333333331*(u_0*u_0)) + xi_26;
            _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 4*_stride_pdfs_3] = omega*(f_eq_common*0.07407407407407407 + u_0*0.22222222222222221 - xi_14 + 0.33333333333333331*(u_0*u_0)) + xi_14;
            _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 5*_stride_pdfs_3] = omega*(f_eq_common*0.07407407407407407 + u_2*0.22222222222222221 - xi_6 + 0.33333333333333331*(u_2*u_2)) + xi_6;
            _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 6*_stride_pdfs_3] = omega*(f_eq_common*0.07407407407407407 + u_2*-0.22222222222222221 - xi_18 + 0.33333333333333331*(u_2*u_2)) + xi_18;
            _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 7*_stride_pdfs_3] = omega*(f_eq_common*0.018518518518518517 + u0Mu1*-0.055555555555555552 - xi_27 + 0.083333333333333329*(u0Mu1*u0Mu1)) + xi_27;
            _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 8*_stride_pdfs_3] = omega*(f_eq_common*0.018518518518518517 + u0Pu1*0.055555555555555552 - xi_19 + 0.083333333333333329*(u0Pu1*u0Pu1)) + xi_19;
            _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 9*_stride_pdfs_3] = omega*(f_eq_common*0.018518518518518517 + u0Pu1*-0.055555555555555552 - xi_20 + 0.083333333333333329*(u0Pu1*u0Pu1)) + xi_20;
            _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 10*_stride_pdfs_3] = omega*(f_eq_common*0.018518518518518517 + u0Mu1*0.055555555555555552 - xi_25 + 0.083333333333333329*(u0Mu1*u0Mu1)) + xi_25;
            _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 11*_stride_pdfs_3] = omega*(f_eq_common*0.018518518518518517 + u1Pu2*0.055555555555555552 - xi_9 + 0.083333333333333329*(u1Pu2*u1Pu2)) + xi_9;
            _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 12*_stride_pdfs_3] = omega*(f_eq_common*0.018518518518518517 + u1Mu2*-0.055555555555555552 - xi_7 + 0.083333333333333329*(u1Mu2*u1Mu2)) + xi_7;
            _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 13*_stride_pdfs_3] = omega*(f_eq_common*0.018518518518518517 + u0Mu2*-0.055555555555555552 - xi_5 + 0.083333333333333329*(u0Mu2*u0Mu2)) + xi_5;
            _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 14*_stride_pdfs_3] = omega*(f_eq_common*0.018518518518518517 + u0Pu2*0.055555555555555552 - xi_17 + 0.083333333333333329*(u0Pu2*u0Pu2)) + xi_17;
            _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 15*_stride_pdfs_3] = omega*(f_eq_common*0.018518518518518517 + u1Mu2*0.055555555555555552 - xi_22 + 0.083333333333333329*(u1Mu2*u1Mu2)) + xi_22;
            _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 16*_stride_pdfs_3] = omega*(f_eq_common*0.018518518518518517 + u1Pu2*-0.055555555555555552 - xi_2 + 0.083333333333333329*(u1Pu2*u1Pu2)) + xi_2;
            _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 17*_stride_pdfs_3] = omega*(f_eq_common*0.018518518518518517 + u0Pu2*-0.055555555555555552 - xi_11 + 0.083333333333333329*(u0Pu2*u0Pu2)) + xi_11;
            _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 18*_stride_pdfs_3] = omega*(f_eq_common*0.018518518518518517 + u0Mu2*0.055555555555555552 - xi_13 + 0.083333333333333329*(u0Mu2*u0Mu2)) + xi_13;
            _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 19*_stride_pdfs_3] = omega*(delta_rho*-0.013888888888888888 + f_eq_common*0.018518518518518517 + u0Pu1*0.013888888888888888 + u_2*0.013888888888888888 - xi_10 + 0.020833333333333332*(u0Pu1*u0Pu1) + 0.020833333333333332*(u0Pu2*u0Pu2) + 0.020833333333333332*(u1Pu2*u1Pu2)) + xi_10;
            _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 20*_stride_pdfs_3] = omega*(delta_rho*-0.013888888888888888 + f_eq_common*0.018518518518518517 + u0Mu1*-0.013888888888888888 + u_2*0.013888888888888888 - xi_12 + 0.020833333333333332*(u0Mu1*u0Mu1) + 0.020833333333333332*(u0Mu2*u0Mu2) + 0.020833333333333332*(u1Pu2*u1Pu2)) + xi_12;
            _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 21*_stride_pdfs_3] = omega*(delta_rho*-0.013888888888888888 + f_eq_common*0.018518518518518517 + u0Mu1*0.013888888888888888 + u_2*0.013888888888888888 - xi_24 + 0.020833333333333332*(u0Mu1*u0Mu1) + 0.020833333333333332*(u0Pu2*u0Pu2) + 0.020833333333333332*(u1Mu2*u1Mu2)) + xi_24;
            _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 22*_stride_pdfs_3] = omega*(delta_rho*-0.013888888888888888 + f_eq_common*0.018518518518518517 + u0Pu1*-0.013888888888888888 + u_2*0.013888888888888888 - xi_15 + 0.020833333333333332*(u0Mu2*u0Mu2) + 0.020833333333333332*(u0Pu1*u0Pu1) + 0.020833333333333332*(u1Mu2*u1Mu2)) + xi_15;
            _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 23*_stride_pdfs_3] = omega*(delta_rho*-0.013888888888888888 + f_eq_common*0.018518518518518517 + u0Pu1*0.013888888888888888 + u_2*-0.013888888888888888 - xi_1 + 0.020833333333333332*(u0Mu2*u0Mu2) + 0.020833333333333332*(u0Pu1*u0Pu1) + 0.020833333333333332*(u1Mu2*u1Mu2)) + xi_1;
            _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 24*_stride_pdfs_3] = omega*(delta_rho*-0.013888888888888888 + f_eq_common*0.018518518518518517 + u0Mu1*-0.013888888888888888 + u_2*-0.013888888888888888 - xi_23 + 0.020833333333333332*(u0Mu1*u0Mu1) + 0.020833333333333332*(u0Pu2*u0Pu2) + 0.020833333333333332*(u1Mu2*u1Mu2)) + xi_23;
            _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 25*_stride_pdfs_3] = omega*(delta_rho*-0.013888888888888888 + f_eq_common*0.018518518518518517 + u0Mu1*0.013888888888888888 + u_2*-0.013888888888888888 - xi_16 + 0.020833333333333332*(u0Mu1*u0Mu1) + 0.020833333333333332*(u0Mu2*u0Mu2) + 0.020833333333333332*(u1Pu2*u1Pu2)) + xi_16;
            _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 26*_stride_pdfs_3] = omega*(delta_rho*-0.013888888888888888 + f_eq_common*0.018518518518518517 + u0Pu1*-0.013888888888888888 + u_2*-0.013888888888888888 - xi_8 + 0.020833333333333332*(u0Pu1*u0Pu1) + 0.020833333333333332*(u0Pu2*u0Pu2) + 0.020833333333333332*(u1Pu2*u1Pu2)) + xi_8;
         }
      }
   }
}
}


namespace internal_d3q27srt_kernel_stream {
static FUNC_PREFIX void d3q27srt_kernel_stream(double * RESTRICT const _data_pdfs, double * RESTRICT  _data_pdfs_tmp, int64_t const _size_pdfs_0, int64_t const _size_pdfs_1, int64_t const _size_pdfs_2, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2, int64_t const _stride_pdfs_3, int64_t const _stride_pdfs_tmp_0, int64_t const _stride_pdfs_tmp_1, int64_t const _stride_pdfs_tmp_2, int64_t const _stride_pdfs_tmp_3)
{
   for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_2; ctr_2 += 1)
   {
      for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_1; ctr_1 += 1)
      {
         for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_0; ctr_0 += 1)
         {
            const double streamed_0 = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2];
            const double streamed_1 = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _stride_pdfs_2*ctr_2 + _stride_pdfs_3];
            const double streamed_2 = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _stride_pdfs_2*ctr_2 + 2*_stride_pdfs_3];
            const double streamed_3 = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 3*_stride_pdfs_3];
            const double streamed_4 = _data_pdfs[_stride_pdfs_0*ctr_0 - _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 4*_stride_pdfs_3];
            const double streamed_5 = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 5*_stride_pdfs_3];
            const double streamed_6 = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 6*_stride_pdfs_3];
            const double streamed_7 = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_0 + _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _stride_pdfs_2*ctr_2 + 7*_stride_pdfs_3];
            const double streamed_8 = _data_pdfs[_stride_pdfs_0*ctr_0 - _stride_pdfs_0 + _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _stride_pdfs_2*ctr_2 + 8*_stride_pdfs_3];
            const double streamed_9 = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _stride_pdfs_2*ctr_2 + 9*_stride_pdfs_3];
            const double streamed_10 = _data_pdfs[_stride_pdfs_0*ctr_0 - _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _stride_pdfs_2*ctr_2 + 10*_stride_pdfs_3];
            const double streamed_11 = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 11*_stride_pdfs_3];
            const double streamed_12 = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 12*_stride_pdfs_3];
            const double streamed_13 = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 13*_stride_pdfs_3];
            const double streamed_14 = _data_pdfs[_stride_pdfs_0*ctr_0 - _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 14*_stride_pdfs_3];
            const double streamed_15 = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 15*_stride_pdfs_3];
            const double streamed_16 = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 16*_stride_pdfs_3];
            const double streamed_17 = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 17*_stride_pdfs_3];
            const double streamed_18 = _data_pdfs[_stride_pdfs_0*ctr_0 - _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 18*_stride_pdfs_3];
            const double streamed_19 = _data_pdfs[_stride_pdfs_0*ctr_0 - _stride_pdfs_0 + _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 19*_stride_pdfs_3];
            const double streamed_20 = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_0 + _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 20*_stride_pdfs_3];
            const double streamed_21 = _data_pdfs[_stride_pdfs_0*ctr_0 - _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 21*_stride_pdfs_3];
            const double streamed_22 = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 22*_stride_pdfs_3];
            const double streamed_23 = _data_pdfs[_stride_pdfs_0*ctr_0 - _stride_pdfs_0 + _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 23*_stride_pdfs_3];
            const double streamed_24 = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_0 + _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 24*_stride_pdfs_3];
            const double streamed_25 = _data_pdfs[_stride_pdfs_0*ctr_0 - _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 25*_stride_pdfs_3];
            const double streamed_26 = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 26*_stride_pdfs_3];
            _data_pdfs_tmp[_stride_pdfs_tmp_0*ctr_0 + _stride_pdfs_tmp_1*ctr_1 + _stride_pdfs_tmp_2*ctr_2] = streamed_0;
            _data_pdfs_tmp[_stride_pdfs_tmp_0*ctr_0 + _stride_pdfs_tmp_1*ctr_1 + _stride_pdfs_tmp_2*ctr_2 + _stride_pdfs_tmp_3] = streamed_1;
            _data_pdfs_tmp[_stride_pdfs_tmp_0*ctr_0 + _stride_pdfs_tmp_1*ctr_1 + _stride_pdfs_tmp_2*ctr_2 + 2*_stride_pdfs_tmp_3] = streamed_2;
            _data_pdfs_tmp[_stride_pdfs_tmp_0*ctr_0 + _stride_pdfs_tmp_1*ctr_1 + _stride_pdfs_tmp_2*ctr_2 + 3*_stride_pdfs_tmp_3] = streamed_3;
            _data_pdfs_tmp[_stride_pdfs_tmp_0*ctr_0 + _stride_pdfs_tmp_1*ctr_1 + _stride_pdfs_tmp_2*ctr_2 + 4*_stride_pdfs_tmp_3] = streamed_4;
            _data_pdfs_tmp[_stride_pdfs_tmp_0*ctr_0 + _stride_pdfs_tmp_1*ctr_1 + _stride_pdfs_tmp_2*ctr_2 + 5*_stride_pdfs_tmp_3] = streamed_5;
            _data_pdfs_tmp[_stride_pdfs_tmp_0*ctr_0 + _stride_pdfs_tmp_1*ctr_1 + _stride_pdfs_tmp_2*ctr_2 + 6*_stride_pdfs_tmp_3] = streamed_6;
            _data_pdfs_tmp[_stride_pdfs_tmp_0*ctr_0 + _stride_pdfs_tmp_1*ctr_1 + _stride_pdfs_tmp_2*ctr_2 + 7*_stride_pdfs_tmp_3] = streamed_7;
            _data_pdfs_tmp[_stride_pdfs_tmp_0*ctr_0 + _stride_pdfs_tmp_1*ctr_1 + _stride_pdfs_tmp_2*ctr_2 + 8*_stride_pdfs_tmp_3] = streamed_8;
            _data_pdfs_tmp[_stride_pdfs_tmp_0*ctr_0 + _stride_pdfs_tmp_1*ctr_1 + _stride_pdfs_tmp_2*ctr_2 + 9*_stride_pdfs_tmp_3] = streamed_9;
            _data_pdfs_tmp[_stride_pdfs_tmp_0*ctr_0 + _stride_pdfs_tmp_1*ctr_1 + _stride_pdfs_tmp_2*ctr_2 + 10*_stride_pdfs_tmp_3] = streamed_10;
            _data_pdfs_tmp[_stride_pdfs_tmp_0*ctr_0 + _stride_pdfs_tmp_1*ctr_1 + _stride_pdfs_tmp_2*ctr_2 + 11*_stride_pdfs_tmp_3] = streamed_11;
            _data_pdfs_tmp[_stride_pdfs_tmp_0*ctr_0 + _stride_pdfs_tmp_1*ctr_1 + _stride_pdfs_tmp_2*ctr_2 + 12*_stride_pdfs_tmp_3] = streamed_12;
            _data_pdfs_tmp[_stride_pdfs_tmp_0*ctr_0 + _stride_pdfs_tmp_1*ctr_1 + _stride_pdfs_tmp_2*ctr_2 + 13*_stride_pdfs_tmp_3] = streamed_13;
            _data_pdfs_tmp[_stride_pdfs_tmp_0*ctr_0 + _stride_pdfs_tmp_1*ctr_1 + _stride_pdfs_tmp_2*ctr_2 + 14*_stride_pdfs_tmp_3] = streamed_14;
            _data_pdfs_tmp[_stride_pdfs_tmp_0*ctr_0 + _stride_pdfs_tmp_1*ctr_1 + _stride_pdfs_tmp_2*ctr_2 + 15*_stride_pdfs_tmp_3] = streamed_15;
            _data_pdfs_tmp[_stride_pdfs_tmp_0*ctr_0 + _stride_pdfs_tmp_1*ctr_1 + _stride_pdfs_tmp_2*ctr_2 + 16*_stride_pdfs_tmp_3] = streamed_16;
            _data_pdfs_tmp[_stride_pdfs_tmp_0*ctr_0 + _stride_pdfs_tmp_1*ctr_1 + _stride_pdfs_tmp_2*ctr_2 + 17*_stride_pdfs_tmp_3] = streamed_17;
            _data_pdfs_tmp[_stride_pdfs_tmp_0*ctr_0 + _stride_pdfs_tmp_1*ctr_1 + _stride_pdfs_tmp_2*ctr_2 + 18*_stride_pdfs_tmp_3] = streamed_18;
            _data_pdfs_tmp[_stride_pdfs_tmp_0*ctr_0 + _stride_pdfs_tmp_1*ctr_1 + _stride_pdfs_tmp_2*ctr_2 + 19*_stride_pdfs_tmp_3] = streamed_19;
            _data_pdfs_tmp[_stride_pdfs_tmp_0*ctr_0 + _stride_pdfs_tmp_1*ctr_1 + _stride_pdfs_tmp_2*ctr_2 + 20*_stride_pdfs_tmp_3] = streamed_20;
            _data_pdfs_tmp[_stride_pdfs_tmp_0*ctr_0 + _stride_pdfs_tmp_1*ctr_1 + _stride_pdfs_tmp_2*ctr_2 + 21*_stride_pdfs_tmp_3] = streamed_21;
            _data_pdfs_tmp[_stride_pdfs_tmp_0*ctr_0 + _stride_pdfs_tmp_1*ctr_1 + _stride_pdfs_tmp_2*ctr_2 + 22*_stride_pdfs_tmp_3] = streamed_22;
            _data_pdfs_tmp[_stride_pdfs_tmp_0*ctr_0 + _stride_pdfs_tmp_1*ctr_1 + _stride_pdfs_tmp_2*ctr_2 + 23*_stride_pdfs_tmp_3] = streamed_23;
            _data_pdfs_tmp[_stride_pdfs_tmp_0*ctr_0 + _stride_pdfs_tmp_1*ctr_1 + _stride_pdfs_tmp_2*ctr_2 + 24*_stride_pdfs_tmp_3] = streamed_24;
            _data_pdfs_tmp[_stride_pdfs_tmp_0*ctr_0 + _stride_pdfs_tmp_1*ctr_1 + _stride_pdfs_tmp_2*ctr_2 + 25*_stride_pdfs_tmp_3] = streamed_25;
            _data_pdfs_tmp[_stride_pdfs_tmp_0*ctr_0 + _stride_pdfs_tmp_1*ctr_1 + _stride_pdfs_tmp_2*ctr_2 + 26*_stride_pdfs_tmp_3] = streamed_26;
         }
      }
   }
}
}


namespace internal_d3q27srt_kernel_streamOnlyNoAdvancement {
static FUNC_PREFIX void d3q27srt_kernel_streamOnlyNoAdvancement(double * RESTRICT const _data_pdfs, double * RESTRICT  _data_pdfs_tmp, int64_t const _size_pdfs_0, int64_t const _size_pdfs_1, int64_t const _size_pdfs_2, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2, int64_t const _stride_pdfs_3, int64_t const _stride_pdfs_tmp_0, int64_t const _stride_pdfs_tmp_1, int64_t const _stride_pdfs_tmp_2, int64_t const _stride_pdfs_tmp_3)
{
   for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_2; ctr_2 += 1)
   {
      for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_1; ctr_1 += 1)
      {
         for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_0; ctr_0 += 1)
         {
            const double streamed_0 = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2];
            const double streamed_1 = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _stride_pdfs_2*ctr_2 + _stride_pdfs_3];
            const double streamed_2 = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _stride_pdfs_2*ctr_2 + 2*_stride_pdfs_3];
            const double streamed_3 = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 3*_stride_pdfs_3];
            const double streamed_4 = _data_pdfs[_stride_pdfs_0*ctr_0 - _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 4*_stride_pdfs_3];
            const double streamed_5 = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 5*_stride_pdfs_3];
            const double streamed_6 = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 6*_stride_pdfs_3];
            const double streamed_7 = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_0 + _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _stride_pdfs_2*ctr_2 + 7*_stride_pdfs_3];
            const double streamed_8 = _data_pdfs[_stride_pdfs_0*ctr_0 - _stride_pdfs_0 + _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _stride_pdfs_2*ctr_2 + 8*_stride_pdfs_3];
            const double streamed_9 = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _stride_pdfs_2*ctr_2 + 9*_stride_pdfs_3];
            const double streamed_10 = _data_pdfs[_stride_pdfs_0*ctr_0 - _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _stride_pdfs_2*ctr_2 + 10*_stride_pdfs_3];
            const double streamed_11 = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 11*_stride_pdfs_3];
            const double streamed_12 = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 12*_stride_pdfs_3];
            const double streamed_13 = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 13*_stride_pdfs_3];
            const double streamed_14 = _data_pdfs[_stride_pdfs_0*ctr_0 - _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 14*_stride_pdfs_3];
            const double streamed_15 = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 15*_stride_pdfs_3];
            const double streamed_16 = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 16*_stride_pdfs_3];
            const double streamed_17 = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 17*_stride_pdfs_3];
            const double streamed_18 = _data_pdfs[_stride_pdfs_0*ctr_0 - _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 18*_stride_pdfs_3];
            const double streamed_19 = _data_pdfs[_stride_pdfs_0*ctr_0 - _stride_pdfs_0 + _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 19*_stride_pdfs_3];
            const double streamed_20 = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_0 + _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 20*_stride_pdfs_3];
            const double streamed_21 = _data_pdfs[_stride_pdfs_0*ctr_0 - _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 21*_stride_pdfs_3];
            const double streamed_22 = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 22*_stride_pdfs_3];
            const double streamed_23 = _data_pdfs[_stride_pdfs_0*ctr_0 - _stride_pdfs_0 + _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 23*_stride_pdfs_3];
            const double streamed_24 = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_0 + _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 24*_stride_pdfs_3];
            const double streamed_25 = _data_pdfs[_stride_pdfs_0*ctr_0 - _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 25*_stride_pdfs_3];
            const double streamed_26 = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 26*_stride_pdfs_3];
            _data_pdfs_tmp[_stride_pdfs_tmp_0*ctr_0 + _stride_pdfs_tmp_1*ctr_1 + _stride_pdfs_tmp_2*ctr_2] = streamed_0;
            _data_pdfs_tmp[_stride_pdfs_tmp_0*ctr_0 + _stride_pdfs_tmp_1*ctr_1 + _stride_pdfs_tmp_2*ctr_2 + _stride_pdfs_tmp_3] = streamed_1;
            _data_pdfs_tmp[_stride_pdfs_tmp_0*ctr_0 + _stride_pdfs_tmp_1*ctr_1 + _stride_pdfs_tmp_2*ctr_2 + 2*_stride_pdfs_tmp_3] = streamed_2;
            _data_pdfs_tmp[_stride_pdfs_tmp_0*ctr_0 + _stride_pdfs_tmp_1*ctr_1 + _stride_pdfs_tmp_2*ctr_2 + 3*_stride_pdfs_tmp_3] = streamed_3;
            _data_pdfs_tmp[_stride_pdfs_tmp_0*ctr_0 + _stride_pdfs_tmp_1*ctr_1 + _stride_pdfs_tmp_2*ctr_2 + 4*_stride_pdfs_tmp_3] = streamed_4;
            _data_pdfs_tmp[_stride_pdfs_tmp_0*ctr_0 + _stride_pdfs_tmp_1*ctr_1 + _stride_pdfs_tmp_2*ctr_2 + 5*_stride_pdfs_tmp_3] = streamed_5;
            _data_pdfs_tmp[_stride_pdfs_tmp_0*ctr_0 + _stride_pdfs_tmp_1*ctr_1 + _stride_pdfs_tmp_2*ctr_2 + 6*_stride_pdfs_tmp_3] = streamed_6;
            _data_pdfs_tmp[_stride_pdfs_tmp_0*ctr_0 + _stride_pdfs_tmp_1*ctr_1 + _stride_pdfs_tmp_2*ctr_2 + 7*_stride_pdfs_tmp_3] = streamed_7;
            _data_pdfs_tmp[_stride_pdfs_tmp_0*ctr_0 + _stride_pdfs_tmp_1*ctr_1 + _stride_pdfs_tmp_2*ctr_2 + 8*_stride_pdfs_tmp_3] = streamed_8;
            _data_pdfs_tmp[_stride_pdfs_tmp_0*ctr_0 + _stride_pdfs_tmp_1*ctr_1 + _stride_pdfs_tmp_2*ctr_2 + 9*_stride_pdfs_tmp_3] = streamed_9;
            _data_pdfs_tmp[_stride_pdfs_tmp_0*ctr_0 + _stride_pdfs_tmp_1*ctr_1 + _stride_pdfs_tmp_2*ctr_2 + 10*_stride_pdfs_tmp_3] = streamed_10;
            _data_pdfs_tmp[_stride_pdfs_tmp_0*ctr_0 + _stride_pdfs_tmp_1*ctr_1 + _stride_pdfs_tmp_2*ctr_2 + 11*_stride_pdfs_tmp_3] = streamed_11;
            _data_pdfs_tmp[_stride_pdfs_tmp_0*ctr_0 + _stride_pdfs_tmp_1*ctr_1 + _stride_pdfs_tmp_2*ctr_2 + 12*_stride_pdfs_tmp_3] = streamed_12;
            _data_pdfs_tmp[_stride_pdfs_tmp_0*ctr_0 + _stride_pdfs_tmp_1*ctr_1 + _stride_pdfs_tmp_2*ctr_2 + 13*_stride_pdfs_tmp_3] = streamed_13;
            _data_pdfs_tmp[_stride_pdfs_tmp_0*ctr_0 + _stride_pdfs_tmp_1*ctr_1 + _stride_pdfs_tmp_2*ctr_2 + 14*_stride_pdfs_tmp_3] = streamed_14;
            _data_pdfs_tmp[_stride_pdfs_tmp_0*ctr_0 + _stride_pdfs_tmp_1*ctr_1 + _stride_pdfs_tmp_2*ctr_2 + 15*_stride_pdfs_tmp_3] = streamed_15;
            _data_pdfs_tmp[_stride_pdfs_tmp_0*ctr_0 + _stride_pdfs_tmp_1*ctr_1 + _stride_pdfs_tmp_2*ctr_2 + 16*_stride_pdfs_tmp_3] = streamed_16;
            _data_pdfs_tmp[_stride_pdfs_tmp_0*ctr_0 + _stride_pdfs_tmp_1*ctr_1 + _stride_pdfs_tmp_2*ctr_2 + 17*_stride_pdfs_tmp_3] = streamed_17;
            _data_pdfs_tmp[_stride_pdfs_tmp_0*ctr_0 + _stride_pdfs_tmp_1*ctr_1 + _stride_pdfs_tmp_2*ctr_2 + 18*_stride_pdfs_tmp_3] = streamed_18;
            _data_pdfs_tmp[_stride_pdfs_tmp_0*ctr_0 + _stride_pdfs_tmp_1*ctr_1 + _stride_pdfs_tmp_2*ctr_2 + 19*_stride_pdfs_tmp_3] = streamed_19;
            _data_pdfs_tmp[_stride_pdfs_tmp_0*ctr_0 + _stride_pdfs_tmp_1*ctr_1 + _stride_pdfs_tmp_2*ctr_2 + 20*_stride_pdfs_tmp_3] = streamed_20;
            _data_pdfs_tmp[_stride_pdfs_tmp_0*ctr_0 + _stride_pdfs_tmp_1*ctr_1 + _stride_pdfs_tmp_2*ctr_2 + 21*_stride_pdfs_tmp_3] = streamed_21;
            _data_pdfs_tmp[_stride_pdfs_tmp_0*ctr_0 + _stride_pdfs_tmp_1*ctr_1 + _stride_pdfs_tmp_2*ctr_2 + 22*_stride_pdfs_tmp_3] = streamed_22;
            _data_pdfs_tmp[_stride_pdfs_tmp_0*ctr_0 + _stride_pdfs_tmp_1*ctr_1 + _stride_pdfs_tmp_2*ctr_2 + 23*_stride_pdfs_tmp_3] = streamed_23;
            _data_pdfs_tmp[_stride_pdfs_tmp_0*ctr_0 + _stride_pdfs_tmp_1*ctr_1 + _stride_pdfs_tmp_2*ctr_2 + 24*_stride_pdfs_tmp_3] = streamed_24;
            _data_pdfs_tmp[_stride_pdfs_tmp_0*ctr_0 + _stride_pdfs_tmp_1*ctr_1 + _stride_pdfs_tmp_2*ctr_2 + 25*_stride_pdfs_tmp_3] = streamed_25;
            _data_pdfs_tmp[_stride_pdfs_tmp_0*ctr_0 + _stride_pdfs_tmp_1*ctr_1 + _stride_pdfs_tmp_2*ctr_2 + 26*_stride_pdfs_tmp_3] = streamed_26;
         }
      }
   }
}
}


namespace internal_d3q27srt_kernel_initialise {
static FUNC_PREFIX void d3q27srt_kernel_initialise(double * RESTRICT const _data_density, double * RESTRICT  _data_pdfs, double * RESTRICT const _data_velocity, int64_t const _size_density_0, int64_t const _size_density_1, int64_t const _size_density_2, int64_t const _stride_density_0, int64_t const _stride_density_1, int64_t const _stride_density_2, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2, int64_t const _stride_pdfs_3, int64_t const _stride_velocity_0, int64_t const _stride_velocity_1, int64_t const _stride_velocity_2, int64_t const _stride_velocity_3)
{
   for (int64_t ctr_2 = 0; ctr_2 < _size_density_2; ctr_2 += 1)
   {
      for (int64_t ctr_1 = 0; ctr_1 < _size_density_1; ctr_1 += 1)
      {
         for (int64_t ctr_0 = 0; ctr_0 < _size_density_0; ctr_0 += 1)
         {
            const double rho = _data_density[_stride_density_0*ctr_0 + _stride_density_1*ctr_1 + _stride_density_2*ctr_2];
            const double delta_rho = rho - 1.0;
            const double u_0 = _data_velocity[_stride_velocity_0*ctr_0 + _stride_velocity_1*ctr_1 + _stride_velocity_2*ctr_2];
            const double u_1 = _data_velocity[_stride_velocity_0*ctr_0 + _stride_velocity_1*ctr_1 + _stride_velocity_2*ctr_2 + _stride_velocity_3];
            const double u_2 = _data_velocity[_stride_velocity_0*ctr_0 + _stride_velocity_1*ctr_1 + _stride_velocity_2*ctr_2 + 2*_stride_velocity_3];
            _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2] = delta_rho*0.29629629629629628 - 0.44444444444444442*(u_0*u_0) - 0.44444444444444442*(u_1*u_1) - 0.44444444444444442*(u_2*u_2);
            _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _stride_pdfs_2*ctr_2 + _stride_pdfs_3] = delta_rho*0.07407407407407407 + u_1*0.22222222222222221 - 0.1111111111111111*(u_0*u_0) - 0.1111111111111111*(u_2*u_2) + 0.22222222222222221*(u_1*u_1);
            _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _stride_pdfs_2*ctr_2 + 2*_stride_pdfs_3] = delta_rho*0.07407407407407407 + u_1*-0.22222222222222221 - 0.1111111111111111*(u_0*u_0) - 0.1111111111111111*(u_2*u_2) + 0.22222222222222221*(u_1*u_1);
            _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 3*_stride_pdfs_3] = delta_rho*0.07407407407407407 + u_0*-0.22222222222222221 - 0.1111111111111111*(u_1*u_1) - 0.1111111111111111*(u_2*u_2) + 0.22222222222222221*(u_0*u_0);
            _data_pdfs[_stride_pdfs_0*ctr_0 - _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 4*_stride_pdfs_3] = delta_rho*0.07407407407407407 + u_0*0.22222222222222221 - 0.1111111111111111*(u_1*u_1) - 0.1111111111111111*(u_2*u_2) + 0.22222222222222221*(u_0*u_0);
            _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 5*_stride_pdfs_3] = delta_rho*0.07407407407407407 + u_2*0.22222222222222221 - 0.1111111111111111*(u_0*u_0) - 0.1111111111111111*(u_1*u_1) + 0.22222222222222221*(u_2*u_2);
            _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 6*_stride_pdfs_3] = delta_rho*0.07407407407407407 + u_2*-0.22222222222222221 - 0.1111111111111111*(u_0*u_0) - 0.1111111111111111*(u_1*u_1) + 0.22222222222222221*(u_2*u_2);
            _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_0 + _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _stride_pdfs_2*ctr_2 + 7*_stride_pdfs_3] = delta_rho*0.018518518518518517 + u_0*u_1*-0.16666666666666666 + u_0*-0.055555555555555552 + u_1*0.055555555555555552 - 0.027777777777777776*(u_2*u_2) + 0.055555555555555552*(u_0*u_0) + 0.055555555555555552*(u_1*u_1);
            _data_pdfs[_stride_pdfs_0*ctr_0 - _stride_pdfs_0 + _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _stride_pdfs_2*ctr_2 + 8*_stride_pdfs_3] = delta_rho*0.018518518518518517 + u_0*u_1*0.16666666666666666 + u_0*0.055555555555555552 + u_1*0.055555555555555552 - 0.027777777777777776*(u_2*u_2) + 0.055555555555555552*(u_0*u_0) + 0.055555555555555552*(u_1*u_1);
            _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _stride_pdfs_2*ctr_2 + 9*_stride_pdfs_3] = delta_rho*0.018518518518518517 + u_0*u_1*0.16666666666666666 + u_0*-0.055555555555555552 + u_1*-0.055555555555555552 - 0.027777777777777776*(u_2*u_2) + 0.055555555555555552*(u_0*u_0) + 0.055555555555555552*(u_1*u_1);
            _data_pdfs[_stride_pdfs_0*ctr_0 - _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _stride_pdfs_2*ctr_2 + 10*_stride_pdfs_3] = delta_rho*0.018518518518518517 + u_0*u_1*-0.16666666666666666 + u_0*0.055555555555555552 + u_1*-0.055555555555555552 - 0.027777777777777776*(u_2*u_2) + 0.055555555555555552*(u_0*u_0) + 0.055555555555555552*(u_1*u_1);
            _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 11*_stride_pdfs_3] = delta_rho*0.018518518518518517 + u_1*u_2*0.16666666666666666 + u_1*0.055555555555555552 + u_2*0.055555555555555552 - 0.027777777777777776*(u_0*u_0) + 0.055555555555555552*(u_1*u_1) + 0.055555555555555552*(u_2*u_2);
            _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 12*_stride_pdfs_3] = delta_rho*0.018518518518518517 + u_1*u_2*-0.16666666666666666 + u_1*-0.055555555555555552 + u_2*0.055555555555555552 - 0.027777777777777776*(u_0*u_0) + 0.055555555555555552*(u_1*u_1) + 0.055555555555555552*(u_2*u_2);
            _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 13*_stride_pdfs_3] = delta_rho*0.018518518518518517 + u_0*u_2*-0.16666666666666666 + u_0*-0.055555555555555552 + u_2*0.055555555555555552 - 0.027777777777777776*(u_1*u_1) + 0.055555555555555552*(u_0*u_0) + 0.055555555555555552*(u_2*u_2);
            _data_pdfs[_stride_pdfs_0*ctr_0 - _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 14*_stride_pdfs_3] = delta_rho*0.018518518518518517 + u_0*u_2*0.16666666666666666 + u_0*0.055555555555555552 + u_2*0.055555555555555552 - 0.027777777777777776*(u_1*u_1) + 0.055555555555555552*(u_0*u_0) + 0.055555555555555552*(u_2*u_2);
            _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 15*_stride_pdfs_3] = delta_rho*0.018518518518518517 + u_1*u_2*-0.16666666666666666 + u_1*0.055555555555555552 + u_2*-0.055555555555555552 - 0.027777777777777776*(u_0*u_0) + 0.055555555555555552*(u_1*u_1) + 0.055555555555555552*(u_2*u_2);
            _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 16*_stride_pdfs_3] = delta_rho*0.018518518518518517 + u_1*u_2*0.16666666666666666 + u_1*-0.055555555555555552 + u_2*-0.055555555555555552 - 0.027777777777777776*(u_0*u_0) + 0.055555555555555552*(u_1*u_1) + 0.055555555555555552*(u_2*u_2);
            _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 17*_stride_pdfs_3] = delta_rho*0.018518518518518517 + u_0*u_2*0.16666666666666666 + u_0*-0.055555555555555552 + u_2*-0.055555555555555552 - 0.027777777777777776*(u_1*u_1) + 0.055555555555555552*(u_0*u_0) + 0.055555555555555552*(u_2*u_2);
            _data_pdfs[_stride_pdfs_0*ctr_0 - _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 18*_stride_pdfs_3] = delta_rho*0.018518518518518517 + u_0*u_2*-0.16666666666666666 + u_0*0.055555555555555552 + u_2*-0.055555555555555552 - 0.027777777777777776*(u_1*u_1) + 0.055555555555555552*(u_0*u_0) + 0.055555555555555552*(u_2*u_2);
            _data_pdfs[_stride_pdfs_0*ctr_0 - _stride_pdfs_0 + _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 19*_stride_pdfs_3] = delta_rho*0.0046296296296296294 + u_0*u_1*0.041666666666666664 + u_0*u_2*0.041666666666666664 + u_0*0.013888888888888888 + u_1*u_2*0.041666666666666664 + u_1*0.013888888888888888 + u_2*0.013888888888888888 + 0.013888888888888888*(u_0*u_0) + 0.013888888888888888*(u_1*u_1) + 0.013888888888888888*(u_2*u_2);
            _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_0 + _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 20*_stride_pdfs_3] = delta_rho*0.0046296296296296294 + u_0*u_1*-0.041666666666666664 + u_0*u_2*-0.041666666666666664 + u_0*-0.013888888888888888 + u_1*u_2*0.041666666666666664 + u_1*0.013888888888888888 + u_2*0.013888888888888888 + 0.013888888888888888*(u_0*u_0) + 0.013888888888888888*(u_1*u_1) + 0.013888888888888888*(u_2*u_2);
            _data_pdfs[_stride_pdfs_0*ctr_0 - _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 21*_stride_pdfs_3] = delta_rho*0.0046296296296296294 + u_0*u_1*-0.041666666666666664 + u_0*u_2*0.041666666666666664 + u_0*0.013888888888888888 + u_1*u_2*-0.041666666666666664 + u_1*-0.013888888888888888 + u_2*0.013888888888888888 + 0.013888888888888888*(u_0*u_0) + 0.013888888888888888*(u_1*u_1) + 0.013888888888888888*(u_2*u_2);
            _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 22*_stride_pdfs_3] = delta_rho*0.0046296296296296294 + u_0*u_1*0.041666666666666664 + u_0*u_2*-0.041666666666666664 + u_0*-0.013888888888888888 + u_1*u_2*-0.041666666666666664 + u_1*-0.013888888888888888 + u_2*0.013888888888888888 + 0.013888888888888888*(u_0*u_0) + 0.013888888888888888*(u_1*u_1) + 0.013888888888888888*(u_2*u_2);
            _data_pdfs[_stride_pdfs_0*ctr_0 - _stride_pdfs_0 + _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 23*_stride_pdfs_3] = delta_rho*0.0046296296296296294 + u_0*u_1*0.041666666666666664 + u_0*u_2*-0.041666666666666664 + u_0*0.013888888888888888 + u_1*u_2*-0.041666666666666664 + u_1*0.013888888888888888 + u_2*-0.013888888888888888 + 0.013888888888888888*(u_0*u_0) + 0.013888888888888888*(u_1*u_1) + 0.013888888888888888*(u_2*u_2);
            _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_0 + _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 24*_stride_pdfs_3] = delta_rho*0.0046296296296296294 + u_0*u_1*-0.041666666666666664 + u_0*u_2*0.041666666666666664 + u_0*-0.013888888888888888 + u_1*u_2*-0.041666666666666664 + u_1*0.013888888888888888 + u_2*-0.013888888888888888 + 0.013888888888888888*(u_0*u_0) + 0.013888888888888888*(u_1*u_1) + 0.013888888888888888*(u_2*u_2);
            _data_pdfs[_stride_pdfs_0*ctr_0 - _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 25*_stride_pdfs_3] = delta_rho*0.0046296296296296294 + u_0*u_1*-0.041666666666666664 + u_0*u_2*-0.041666666666666664 + u_0*0.013888888888888888 + u_1*u_2*0.041666666666666664 + u_1*-0.013888888888888888 + u_2*-0.013888888888888888 + 0.013888888888888888*(u_0*u_0) + 0.013888888888888888*(u_1*u_1) + 0.013888888888888888*(u_2*u_2);
            _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 26*_stride_pdfs_3] = delta_rho*0.0046296296296296294 + u_0*u_1*0.041666666666666664 + u_0*u_2*0.041666666666666664 + u_0*-0.013888888888888888 + u_1*u_2*0.041666666666666664 + u_1*-0.013888888888888888 + u_2*-0.013888888888888888 + 0.013888888888888888*(u_0*u_0) + 0.013888888888888888*(u_1*u_1) + 0.013888888888888888*(u_2*u_2);
         }
      }
   }
}
}


namespace internal_d3q27srt_kernel_getter {
static FUNC_PREFIX void d3q27srt_kernel_getter(double * RESTRICT  _data_density, double * RESTRICT const _data_pdfs, double * RESTRICT  _data_velocity, int64_t const _size_density_0, int64_t const _size_density_1, int64_t const _size_density_2, int64_t const _stride_density_0, int64_t const _stride_density_1, int64_t const _stride_density_2, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2, int64_t const _stride_pdfs_3, int64_t const _stride_velocity_0, int64_t const _stride_velocity_1, int64_t const _stride_velocity_2, int64_t const _stride_velocity_3)
{
   for (int64_t ctr_2 = 0; ctr_2 < _size_density_2; ctr_2 += 1)
   {
      for (int64_t ctr_1 = 0; ctr_1 < _size_density_1; ctr_1 += 1)
      {
         for (int64_t ctr_0 = 0; ctr_0 < _size_density_0; ctr_0 += 1)
         {
            const double vel0Term = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 10*_stride_pdfs_3] + _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 14*_stride_pdfs_3] + _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 18*_stride_pdfs_3] + _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 19*_stride_pdfs_3] + _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 21*_stride_pdfs_3] + _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 23*_stride_pdfs_3] + _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 25*_stride_pdfs_3] + _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 4*_stride_pdfs_3] + _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 8*_stride_pdfs_3];
            const double momdensity_0 = vel0Term - _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 13*_stride_pdfs_3] - _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 17*_stride_pdfs_3] - _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 20*_stride_pdfs_3] - _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 22*_stride_pdfs_3] - _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 24*_stride_pdfs_3] - _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 26*_stride_pdfs_3] - _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 3*_stride_pdfs_3] - _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 7*_stride_pdfs_3] - _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 9*_stride_pdfs_3];
            const double vel1Term = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 11*_stride_pdfs_3] + _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 15*_stride_pdfs_3] + _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 20*_stride_pdfs_3] + _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 24*_stride_pdfs_3] + _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 7*_stride_pdfs_3] + _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + _stride_pdfs_3];
            const double momdensity_1 = vel1Term - _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 10*_stride_pdfs_3] - _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 12*_stride_pdfs_3] - _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 16*_stride_pdfs_3] + _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 19*_stride_pdfs_3] - _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 2*_stride_pdfs_3] - _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 21*_stride_pdfs_3] - _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 22*_stride_pdfs_3] + _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 23*_stride_pdfs_3] - _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 25*_stride_pdfs_3] - _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 26*_stride_pdfs_3] + _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 8*_stride_pdfs_3] - _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 9*_stride_pdfs_3];
            const double vel2Term = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 12*_stride_pdfs_3] + _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 13*_stride_pdfs_3] + _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 22*_stride_pdfs_3] + _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 5*_stride_pdfs_3];
            const double delta_rho = vel0Term + vel1Term + vel2Term + _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 16*_stride_pdfs_3] + _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 17*_stride_pdfs_3] + _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 2*_stride_pdfs_3] + _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 26*_stride_pdfs_3] + _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 3*_stride_pdfs_3] + _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 6*_stride_pdfs_3] + _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 9*_stride_pdfs_3] + _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2];
            const double momdensity_2 = vel2Term + _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 11*_stride_pdfs_3] + _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 14*_stride_pdfs_3] - _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 15*_stride_pdfs_3] - _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 16*_stride_pdfs_3] - _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 17*_stride_pdfs_3] - _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 18*_stride_pdfs_3] + _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 19*_stride_pdfs_3] + _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 20*_stride_pdfs_3] + _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 21*_stride_pdfs_3] - _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 23*_stride_pdfs_3] - _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 24*_stride_pdfs_3] - _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 25*_stride_pdfs_3] - _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 26*_stride_pdfs_3] - _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 6*_stride_pdfs_3];
            const double rho = delta_rho + 1.0;
            const double u_0 = momdensity_0;
            const double u_1 = momdensity_1;
            const double u_2 = momdensity_2;
            _data_density[_stride_density_0*ctr_0 + _stride_density_1*ctr_1 + _stride_density_2*ctr_2] = rho;
            _data_velocity[_stride_velocity_0*ctr_0 + _stride_velocity_1*ctr_1 + _stride_velocity_2*ctr_2] = u_0;
            _data_velocity[_stride_velocity_0*ctr_0 + _stride_velocity_1*ctr_1 + _stride_velocity_2*ctr_2 + _stride_velocity_3] = u_1;
            _data_velocity[_stride_velocity_0*ctr_0 + _stride_velocity_1*ctr_1 + _stride_velocity_2*ctr_2 + 2*_stride_velocity_3] = u_2;
         }
      }
   }
}
}





void D3Q27SRT::streamCollide( field::GhostLayerField<double, 27> * pdfs, field::GhostLayerField<double, 27> * pdfs_tmp, double omega, const cell_idx_t ghost_layers )
{
   WALBERLA_ASSERT_GREATER_EQUAL(-ghost_layers, -int_c(pdfs->nrOfGhostLayers()))
   double * RESTRICT const _data_pdfs = pdfs->dataAt(-ghost_layers, -ghost_layers, -ghost_layers, 0);
   WALBERLA_ASSERT_GREATER_EQUAL(-ghost_layers, -int_c(pdfs_tmp->nrOfGhostLayers()))
   double * RESTRICT  _data_pdfs_tmp = pdfs_tmp->dataAt(-ghost_layers, -ghost_layers, -ghost_layers, 0);
   WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(), int64_t(int64_c(pdfs->xSize()) + 2*ghost_layers))
   const int64_t _size_pdfs_0 = int64_t(int64_c(pdfs->xSize()) + 2*ghost_layers);
   WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(), int64_t(int64_c(pdfs->ySize()) + 2*ghost_layers))
   const int64_t _size_pdfs_1 = int64_t(int64_c(pdfs->ySize()) + 2*ghost_layers);
   WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(), int64_t(int64_c(pdfs->zSize()) + 2*ghost_layers))
   const int64_t _size_pdfs_2 = int64_t(int64_c(pdfs->zSize()) + 2*ghost_layers);
   const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
   const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
   const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
   const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
   const int64_t _stride_pdfs_tmp_0 = int64_t(pdfs_tmp->xStride());
   const int64_t _stride_pdfs_tmp_1 = int64_t(pdfs_tmp->yStride());
   const int64_t _stride_pdfs_tmp_2 = int64_t(pdfs_tmp->zStride());
   const int64_t _stride_pdfs_tmp_3 = int64_t(1 * int64_t(pdfs_tmp->fStride()));
   internal_d3q27srt_kernel_streamCollide::d3q27srt_kernel_streamCollide(_data_pdfs, _data_pdfs_tmp, _size_pdfs_0, _size_pdfs_1, _size_pdfs_2, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3, _stride_pdfs_tmp_0, _stride_pdfs_tmp_1, _stride_pdfs_tmp_2, _stride_pdfs_tmp_3, omega);
}
void D3Q27SRT::streamCollideCellInterval( field::GhostLayerField<double, 27> * pdfs, field::GhostLayerField<double, 27> * pdfs_tmp, double omega, const CellInterval & ci)
{
   WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()))
   WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()))
   WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()))
   double * RESTRICT const _data_pdfs = pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
   WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs_tmp->nrOfGhostLayers()))
   WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs_tmp->nrOfGhostLayers()))
   WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs_tmp->nrOfGhostLayers()))
   double * RESTRICT  _data_pdfs_tmp = pdfs_tmp->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
   WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(), int64_t(int64_c(ci.xSize()) + 0))
   const int64_t _size_pdfs_0 = int64_t(int64_c(ci.xSize()) + 0);
   WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(), int64_t(int64_c(ci.ySize()) + 0))
   const int64_t _size_pdfs_1 = int64_t(int64_c(ci.ySize()) + 0);
   WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(), int64_t(int64_c(ci.zSize()) + 0))
   const int64_t _size_pdfs_2 = int64_t(int64_c(ci.zSize()) + 0);
   const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
   const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
   const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
   const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
   const int64_t _stride_pdfs_tmp_0 = int64_t(pdfs_tmp->xStride());
   const int64_t _stride_pdfs_tmp_1 = int64_t(pdfs_tmp->yStride());
   const int64_t _stride_pdfs_tmp_2 = int64_t(pdfs_tmp->zStride());
   const int64_t _stride_pdfs_tmp_3 = int64_t(1 * int64_t(pdfs_tmp->fStride()));
   internal_d3q27srt_kernel_streamCollide::d3q27srt_kernel_streamCollide(_data_pdfs, _data_pdfs_tmp, _size_pdfs_0, _size_pdfs_1, _size_pdfs_2, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3, _stride_pdfs_tmp_0, _stride_pdfs_tmp_1, _stride_pdfs_tmp_2, _stride_pdfs_tmp_3, omega);
}

void D3Q27SRT::collide( field::GhostLayerField<double, 27> * pdfs, double omega, const cell_idx_t ghost_layers )
{
   WALBERLA_ASSERT_GREATER_EQUAL(-ghost_layers, -int_c(pdfs->nrOfGhostLayers()))
   double * RESTRICT  _data_pdfs = pdfs->dataAt(-ghost_layers, -ghost_layers, -ghost_layers, 0);
   WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(), int64_t(int64_c(pdfs->xSize()) + 2*ghost_layers))
   const int64_t _size_pdfs_0 = int64_t(int64_c(pdfs->xSize()) + 2*ghost_layers);
   WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(), int64_t(int64_c(pdfs->ySize()) + 2*ghost_layers))
   const int64_t _size_pdfs_1 = int64_t(int64_c(pdfs->ySize()) + 2*ghost_layers);
   WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(), int64_t(int64_c(pdfs->zSize()) + 2*ghost_layers))
   const int64_t _size_pdfs_2 = int64_t(int64_c(pdfs->zSize()) + 2*ghost_layers);
   const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
   const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
   const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
   const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
   internal_d3q27srt_kernel_collide::d3q27srt_kernel_collide(_data_pdfs, _size_pdfs_0, _size_pdfs_1, _size_pdfs_2, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3, omega);
}
void D3Q27SRT::collideCellInterval( field::GhostLayerField<double, 27> * pdfs, double omega, const CellInterval & ci)
{
   WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()))
   WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()))
   WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()))
   double * RESTRICT  _data_pdfs = pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
   WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(), int64_t(int64_c(ci.xSize()) + 0))
   const int64_t _size_pdfs_0 = int64_t(int64_c(ci.xSize()) + 0);
   WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(), int64_t(int64_c(ci.ySize()) + 0))
   const int64_t _size_pdfs_1 = int64_t(int64_c(ci.ySize()) + 0);
   WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(), int64_t(int64_c(ci.zSize()) + 0))
   const int64_t _size_pdfs_2 = int64_t(int64_c(ci.zSize()) + 0);
   const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
   const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
   const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
   const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
   internal_d3q27srt_kernel_collide::d3q27srt_kernel_collide(_data_pdfs, _size_pdfs_0, _size_pdfs_1, _size_pdfs_2, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3, omega);
}

void D3Q27SRT::stream( field::GhostLayerField<double, 27> * pdfs, field::GhostLayerField<double, 27> * pdfs_tmp, const cell_idx_t ghost_layers )
{
   WALBERLA_ASSERT_GREATER_EQUAL(-ghost_layers, -int_c(pdfs->nrOfGhostLayers()))
   double * RESTRICT const _data_pdfs = pdfs->dataAt(-ghost_layers, -ghost_layers, -ghost_layers, 0);
   WALBERLA_ASSERT_GREATER_EQUAL(-ghost_layers, -int_c(pdfs_tmp->nrOfGhostLayers()))
   double * RESTRICT  _data_pdfs_tmp = pdfs_tmp->dataAt(-ghost_layers, -ghost_layers, -ghost_layers, 0);
   WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(), int64_t(int64_c(pdfs->xSize()) + 2*ghost_layers))
   const int64_t _size_pdfs_0 = int64_t(int64_c(pdfs->xSize()) + 2*ghost_layers);
   WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(), int64_t(int64_c(pdfs->ySize()) + 2*ghost_layers))
   const int64_t _size_pdfs_1 = int64_t(int64_c(pdfs->ySize()) + 2*ghost_layers);
   WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(), int64_t(int64_c(pdfs->zSize()) + 2*ghost_layers))
   const int64_t _size_pdfs_2 = int64_t(int64_c(pdfs->zSize()) + 2*ghost_layers);
   const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
   const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
   const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
   const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
   const int64_t _stride_pdfs_tmp_0 = int64_t(pdfs_tmp->xStride());
   const int64_t _stride_pdfs_tmp_1 = int64_t(pdfs_tmp->yStride());
   const int64_t _stride_pdfs_tmp_2 = int64_t(pdfs_tmp->zStride());
   const int64_t _stride_pdfs_tmp_3 = int64_t(1 * int64_t(pdfs_tmp->fStride()));
   internal_d3q27srt_kernel_stream::d3q27srt_kernel_stream(_data_pdfs, _data_pdfs_tmp, _size_pdfs_0, _size_pdfs_1, _size_pdfs_2, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3, _stride_pdfs_tmp_0, _stride_pdfs_tmp_1, _stride_pdfs_tmp_2, _stride_pdfs_tmp_3);
}
void D3Q27SRT::streamCellInterval( field::GhostLayerField<double, 27> * pdfs, field::GhostLayerField<double, 27> * pdfs_tmp, const CellInterval & ci)
{
   WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()))
   WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()))
   WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()))
   double * RESTRICT const _data_pdfs = pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
   WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs_tmp->nrOfGhostLayers()))
   WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs_tmp->nrOfGhostLayers()))
   WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs_tmp->nrOfGhostLayers()))
   double * RESTRICT  _data_pdfs_tmp = pdfs_tmp->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
   WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(), int64_t(int64_c(ci.xSize()) + 0))
   const int64_t _size_pdfs_0 = int64_t(int64_c(ci.xSize()) + 0);
   WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(), int64_t(int64_c(ci.ySize()) + 0))
   const int64_t _size_pdfs_1 = int64_t(int64_c(ci.ySize()) + 0);
   WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(), int64_t(int64_c(ci.zSize()) + 0))
   const int64_t _size_pdfs_2 = int64_t(int64_c(ci.zSize()) + 0);
   const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
   const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
   const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
   const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
   const int64_t _stride_pdfs_tmp_0 = int64_t(pdfs_tmp->xStride());
   const int64_t _stride_pdfs_tmp_1 = int64_t(pdfs_tmp->yStride());
   const int64_t _stride_pdfs_tmp_2 = int64_t(pdfs_tmp->zStride());
   const int64_t _stride_pdfs_tmp_3 = int64_t(1 * int64_t(pdfs_tmp->fStride()));
   internal_d3q27srt_kernel_stream::d3q27srt_kernel_stream(_data_pdfs, _data_pdfs_tmp, _size_pdfs_0, _size_pdfs_1, _size_pdfs_2, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3, _stride_pdfs_tmp_0, _stride_pdfs_tmp_1, _stride_pdfs_tmp_2, _stride_pdfs_tmp_3);
}

void D3Q27SRT::streamOnlyNoAdvancement( field::GhostLayerField<double, 27> * pdfs, field::GhostLayerField<double, 27> * pdfs_tmp, const cell_idx_t ghost_layers )
{
   WALBERLA_ASSERT_GREATER_EQUAL(-ghost_layers, -int_c(pdfs->nrOfGhostLayers()))
   double * RESTRICT const _data_pdfs = pdfs->dataAt(-ghost_layers, -ghost_layers, -ghost_layers, 0);
   WALBERLA_ASSERT_GREATER_EQUAL(-ghost_layers, -int_c(pdfs_tmp->nrOfGhostLayers()))
   double * RESTRICT  _data_pdfs_tmp = pdfs_tmp->dataAt(-ghost_layers, -ghost_layers, -ghost_layers, 0);
   WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(), int64_t(int64_c(pdfs->xSize()) + 2*ghost_layers))
   const int64_t _size_pdfs_0 = int64_t(int64_c(pdfs->xSize()) + 2*ghost_layers);
   WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(), int64_t(int64_c(pdfs->ySize()) + 2*ghost_layers))
   const int64_t _size_pdfs_1 = int64_t(int64_c(pdfs->ySize()) + 2*ghost_layers);
   WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(), int64_t(int64_c(pdfs->zSize()) + 2*ghost_layers))
   const int64_t _size_pdfs_2 = int64_t(int64_c(pdfs->zSize()) + 2*ghost_layers);
   const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
   const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
   const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
   const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
   const int64_t _stride_pdfs_tmp_0 = int64_t(pdfs_tmp->xStride());
   const int64_t _stride_pdfs_tmp_1 = int64_t(pdfs_tmp->yStride());
   const int64_t _stride_pdfs_tmp_2 = int64_t(pdfs_tmp->zStride());
   const int64_t _stride_pdfs_tmp_3 = int64_t(1 * int64_t(pdfs_tmp->fStride()));
   internal_d3q27srt_kernel_streamOnlyNoAdvancement::d3q27srt_kernel_streamOnlyNoAdvancement(_data_pdfs, _data_pdfs_tmp, _size_pdfs_0, _size_pdfs_1, _size_pdfs_2, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3, _stride_pdfs_tmp_0, _stride_pdfs_tmp_1, _stride_pdfs_tmp_2, _stride_pdfs_tmp_3);
}
void D3Q27SRT::streamOnlyNoAdvancementCellInterval( field::GhostLayerField<double, 27> * pdfs, field::GhostLayerField<double, 27> * pdfs_tmp, const CellInterval & ci)
{
   WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()))
   WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()))
   WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()))
   double * RESTRICT const _data_pdfs = pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
   WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs_tmp->nrOfGhostLayers()))
   WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs_tmp->nrOfGhostLayers()))
   WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs_tmp->nrOfGhostLayers()))
   double * RESTRICT  _data_pdfs_tmp = pdfs_tmp->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
   WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(), int64_t(int64_c(ci.xSize()) + 0))
   const int64_t _size_pdfs_0 = int64_t(int64_c(ci.xSize()) + 0);
   WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(), int64_t(int64_c(ci.ySize()) + 0))
   const int64_t _size_pdfs_1 = int64_t(int64_c(ci.ySize()) + 0);
   WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(), int64_t(int64_c(ci.zSize()) + 0))
   const int64_t _size_pdfs_2 = int64_t(int64_c(ci.zSize()) + 0);
   const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
   const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
   const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
   const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
   const int64_t _stride_pdfs_tmp_0 = int64_t(pdfs_tmp->xStride());
   const int64_t _stride_pdfs_tmp_1 = int64_t(pdfs_tmp->yStride());
   const int64_t _stride_pdfs_tmp_2 = int64_t(pdfs_tmp->zStride());
   const int64_t _stride_pdfs_tmp_3 = int64_t(1 * int64_t(pdfs_tmp->fStride()));
   internal_d3q27srt_kernel_streamOnlyNoAdvancement::d3q27srt_kernel_streamOnlyNoAdvancement(_data_pdfs, _data_pdfs_tmp, _size_pdfs_0, _size_pdfs_1, _size_pdfs_2, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3, _stride_pdfs_tmp_0, _stride_pdfs_tmp_1, _stride_pdfs_tmp_2, _stride_pdfs_tmp_3);
}

void D3Q27SRT::initialise( field::GhostLayerField<double, 1> * density, field::GhostLayerField<double, 27> * pdfs, field::GhostLayerField<double, 3> * velocity, const cell_idx_t ghost_layers )
{
   WALBERLA_ASSERT_GREATER_EQUAL(-ghost_layers, -int_c(density->nrOfGhostLayers()))
   double * RESTRICT const _data_density = density->dataAt(-ghost_layers, -ghost_layers, -ghost_layers, 0);
   WALBERLA_ASSERT_GREATER_EQUAL(-ghost_layers, -int_c(pdfs->nrOfGhostLayers()))
   double * RESTRICT  _data_pdfs = pdfs->dataAt(-ghost_layers, -ghost_layers, -ghost_layers, 0);
   WALBERLA_ASSERT_GREATER_EQUAL(-ghost_layers, -int_c(velocity->nrOfGhostLayers()))
   double * RESTRICT const _data_velocity = velocity->dataAt(-ghost_layers, -ghost_layers, -ghost_layers, 0);
   WALBERLA_ASSERT_GREATER_EQUAL(density->xSizeWithGhostLayer(), int64_t(int64_c(density->xSize()) + 2*ghost_layers))
   const int64_t _size_density_0 = int64_t(int64_c(density->xSize()) + 2*ghost_layers);
   WALBERLA_ASSERT_GREATER_EQUAL(density->ySizeWithGhostLayer(), int64_t(int64_c(density->ySize()) + 2*ghost_layers))
   const int64_t _size_density_1 = int64_t(int64_c(density->ySize()) + 2*ghost_layers);
   WALBERLA_ASSERT_GREATER_EQUAL(density->zSizeWithGhostLayer(), int64_t(int64_c(density->zSize()) + 2*ghost_layers))
   const int64_t _size_density_2 = int64_t(int64_c(density->zSize()) + 2*ghost_layers);
   const int64_t _stride_density_0 = int64_t(density->xStride());
   const int64_t _stride_density_1 = int64_t(density->yStride());
   const int64_t _stride_density_2 = int64_t(density->zStride());
   const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
   const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
   const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
   const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
   const int64_t _stride_velocity_0 = int64_t(velocity->xStride());
   const int64_t _stride_velocity_1 = int64_t(velocity->yStride());
   const int64_t _stride_velocity_2 = int64_t(velocity->zStride());
   const int64_t _stride_velocity_3 = int64_t(1 * int64_t(velocity->fStride()));
   internal_d3q27srt_kernel_initialise::d3q27srt_kernel_initialise(_data_density, _data_pdfs, _data_velocity, _size_density_0, _size_density_1, _size_density_2, _stride_density_0, _stride_density_1, _stride_density_2, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3, _stride_velocity_0, _stride_velocity_1, _stride_velocity_2, _stride_velocity_3);
}
void D3Q27SRT::initialiseCellInterval( field::GhostLayerField<double, 1> * density, field::GhostLayerField<double, 27> * pdfs, field::GhostLayerField<double, 3> * velocity, const CellInterval & ci)
{
   WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(density->nrOfGhostLayers()))
   WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(density->nrOfGhostLayers()))
   WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(density->nrOfGhostLayers()))
   double * RESTRICT const _data_density = density->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
   WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()))
   WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()))
   WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()))
   double * RESTRICT  _data_pdfs = pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
   WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(velocity->nrOfGhostLayers()))
   WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(velocity->nrOfGhostLayers()))
   WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(velocity->nrOfGhostLayers()))
   double * RESTRICT const _data_velocity = velocity->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
   WALBERLA_ASSERT_GREATER_EQUAL(density->xSizeWithGhostLayer(), int64_t(int64_c(ci.xSize()) + 0))
   const int64_t _size_density_0 = int64_t(int64_c(ci.xSize()) + 0);
   WALBERLA_ASSERT_GREATER_EQUAL(density->ySizeWithGhostLayer(), int64_t(int64_c(ci.ySize()) + 0))
   const int64_t _size_density_1 = int64_t(int64_c(ci.ySize()) + 0);
   WALBERLA_ASSERT_GREATER_EQUAL(density->zSizeWithGhostLayer(), int64_t(int64_c(ci.zSize()) + 0))
   const int64_t _size_density_2 = int64_t(int64_c(ci.zSize()) + 0);
   const int64_t _stride_density_0 = int64_t(density->xStride());
   const int64_t _stride_density_1 = int64_t(density->yStride());
   const int64_t _stride_density_2 = int64_t(density->zStride());
   const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
   const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
   const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
   const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
   const int64_t _stride_velocity_0 = int64_t(velocity->xStride());
   const int64_t _stride_velocity_1 = int64_t(velocity->yStride());
   const int64_t _stride_velocity_2 = int64_t(velocity->zStride());
   const int64_t _stride_velocity_3 = int64_t(1 * int64_t(velocity->fStride()));
   internal_d3q27srt_kernel_initialise::d3q27srt_kernel_initialise(_data_density, _data_pdfs, _data_velocity, _size_density_0, _size_density_1, _size_density_2, _stride_density_0, _stride_density_1, _stride_density_2, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3, _stride_velocity_0, _stride_velocity_1, _stride_velocity_2, _stride_velocity_3);
}

void D3Q27SRT::calculateMacroscopicParameters( field::GhostLayerField<double, 1> * density, field::GhostLayerField<double, 27> * pdfs, field::GhostLayerField<double, 3> * velocity, const cell_idx_t ghost_layers )
{
   WALBERLA_ASSERT_GREATER_EQUAL(-ghost_layers, -int_c(density->nrOfGhostLayers()))
   double * RESTRICT  _data_density = density->dataAt(-ghost_layers, -ghost_layers, -ghost_layers, 0);
   WALBERLA_ASSERT_GREATER_EQUAL(-ghost_layers, -int_c(pdfs->nrOfGhostLayers()))
   double * RESTRICT const _data_pdfs = pdfs->dataAt(-ghost_layers, -ghost_layers, -ghost_layers, 0);
   WALBERLA_ASSERT_GREATER_EQUAL(-ghost_layers, -int_c(velocity->nrOfGhostLayers()))
   double * RESTRICT  _data_velocity = velocity->dataAt(-ghost_layers, -ghost_layers, -ghost_layers, 0);
   WALBERLA_ASSERT_GREATER_EQUAL(density->xSizeWithGhostLayer(), int64_t(int64_c(density->xSize()) + 2*ghost_layers))
   const int64_t _size_density_0 = int64_t(int64_c(density->xSize()) + 2*ghost_layers);
   WALBERLA_ASSERT_GREATER_EQUAL(density->ySizeWithGhostLayer(), int64_t(int64_c(density->ySize()) + 2*ghost_layers))
   const int64_t _size_density_1 = int64_t(int64_c(density->ySize()) + 2*ghost_layers);
   WALBERLA_ASSERT_GREATER_EQUAL(density->zSizeWithGhostLayer(), int64_t(int64_c(density->zSize()) + 2*ghost_layers))
   const int64_t _size_density_2 = int64_t(int64_c(density->zSize()) + 2*ghost_layers);
   const int64_t _stride_density_0 = int64_t(density->xStride());
   const int64_t _stride_density_1 = int64_t(density->yStride());
   const int64_t _stride_density_2 = int64_t(density->zStride());
   const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
   const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
   const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
   const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
   const int64_t _stride_velocity_0 = int64_t(velocity->xStride());
   const int64_t _stride_velocity_1 = int64_t(velocity->yStride());
   const int64_t _stride_velocity_2 = int64_t(velocity->zStride());
   const int64_t _stride_velocity_3 = int64_t(1 * int64_t(velocity->fStride()));
   internal_d3q27srt_kernel_getter::d3q27srt_kernel_getter(_data_density, _data_pdfs, _data_velocity, _size_density_0, _size_density_1, _size_density_2, _stride_density_0, _stride_density_1, _stride_density_2, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3, _stride_velocity_0, _stride_velocity_1, _stride_velocity_2, _stride_velocity_3);
}
void D3Q27SRT::calculateMacroscopicParametersCellInterval( field::GhostLayerField<double, 1> * density, field::GhostLayerField<double, 27> * pdfs, field::GhostLayerField<double, 3> * velocity, const CellInterval & ci)
{
   WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(density->nrOfGhostLayers()))
   WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(density->nrOfGhostLayers()))
   WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(density->nrOfGhostLayers()))
   double * RESTRICT  _data_density = density->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
   WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()))
   WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()))
   WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()))
   double * RESTRICT const _data_pdfs = pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
   WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(velocity->nrOfGhostLayers()))
   WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(velocity->nrOfGhostLayers()))
   WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(velocity->nrOfGhostLayers()))
   double * RESTRICT  _data_velocity = velocity->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
   WALBERLA_ASSERT_GREATER_EQUAL(density->xSizeWithGhostLayer(), int64_t(int64_c(ci.xSize()) + 0))
   const int64_t _size_density_0 = int64_t(int64_c(ci.xSize()) + 0);
   WALBERLA_ASSERT_GREATER_EQUAL(density->ySizeWithGhostLayer(), int64_t(int64_c(ci.ySize()) + 0))
   const int64_t _size_density_1 = int64_t(int64_c(ci.ySize()) + 0);
   WALBERLA_ASSERT_GREATER_EQUAL(density->zSizeWithGhostLayer(), int64_t(int64_c(ci.zSize()) + 0))
   const int64_t _size_density_2 = int64_t(int64_c(ci.zSize()) + 0);
   const int64_t _stride_density_0 = int64_t(density->xStride());
   const int64_t _stride_density_1 = int64_t(density->yStride());
   const int64_t _stride_density_2 = int64_t(density->zStride());
   const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
   const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
   const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
   const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
   const int64_t _stride_velocity_0 = int64_t(velocity->xStride());
   const int64_t _stride_velocity_1 = int64_t(velocity->yStride());
   const int64_t _stride_velocity_2 = int64_t(velocity->zStride());
   const int64_t _stride_velocity_3 = int64_t(1 * int64_t(velocity->fStride()));
   internal_d3q27srt_kernel_getter::d3q27srt_kernel_getter(_data_density, _data_pdfs, _data_velocity, _size_density_0, _size_density_1, _size_density_2, _stride_density_0, _stride_density_1, _stride_density_2, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3, _stride_velocity_0, _stride_velocity_1, _stride_velocity_2, _stride_velocity_3);
}



} // namespace lbm
} // namespace walberla


#if ( defined WALBERLA_CXX_COMPILER_IS_GNU ) || ( defined WALBERLA_CXX_COMPILER_IS_CLANG )
#   pragma GCC diagnostic pop
#endif

#if ( defined WALBERLA_CXX_COMPILER_IS_INTEL )
#pragma warning pop
#endif