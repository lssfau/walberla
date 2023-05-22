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
      double * RESTRICT _data_pdfs_2m1_321 = _data_pdfs + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 21*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_2m1_319 = _data_pdfs + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 19*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_2m1_314 = _data_pdfs + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 14*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_20_310 = _data_pdfs + _stride_pdfs_2*ctr_2 + 10*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_20_38 = _data_pdfs + _stride_pdfs_2*ctr_2 + 8*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_21_325 = _data_pdfs + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 25*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_21_323 = _data_pdfs + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 23*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_21_318 = _data_pdfs + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 18*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_20_34 = _data_pdfs + _stride_pdfs_2*ctr_2 + 4*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_2m1_320 = _data_pdfs + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 20*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_2m1_311 = _data_pdfs + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 11*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_20_37 = _data_pdfs + _stride_pdfs_2*ctr_2 + 7*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_20_31 = _data_pdfs + _stride_pdfs_2*ctr_2 + _stride_pdfs_3;
      double * RESTRICT _data_pdfs_21_324 = _data_pdfs + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 24*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_21_315 = _data_pdfs + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 15*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_2m1_313 = _data_pdfs + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 13*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_2m1_322 = _data_pdfs + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 22*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_2m1_312 = _data_pdfs + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 12*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_2m1_35 = _data_pdfs + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 5*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_20_30 = _data_pdfs + _stride_pdfs_2*ctr_2;
      double * RESTRICT _data_pdfs_20_33 = _data_pdfs + _stride_pdfs_2*ctr_2 + 3*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_20_39 = _data_pdfs + _stride_pdfs_2*ctr_2 + 9*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_20_32 = _data_pdfs + _stride_pdfs_2*ctr_2 + 2*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_21_317 = _data_pdfs + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 17*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_21_326 = _data_pdfs + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 26*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_21_316 = _data_pdfs + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 16*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_21_36 = _data_pdfs + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 6*_stride_pdfs_3;
      double * RESTRICT  _data_pdfs_tmp_20_30 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2;
      double * RESTRICT  _data_pdfs_tmp_20_31 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + _stride_pdfs_tmp_3;
      double * RESTRICT  _data_pdfs_tmp_20_32 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 2*_stride_pdfs_tmp_3;
      double * RESTRICT  _data_pdfs_tmp_20_33 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 3*_stride_pdfs_tmp_3;
      double * RESTRICT  _data_pdfs_tmp_20_34 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 4*_stride_pdfs_tmp_3;
      double * RESTRICT  _data_pdfs_tmp_20_35 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 5*_stride_pdfs_tmp_3;
      double * RESTRICT  _data_pdfs_tmp_20_36 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 6*_stride_pdfs_tmp_3;
      double * RESTRICT  _data_pdfs_tmp_20_37 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 7*_stride_pdfs_tmp_3;
      double * RESTRICT  _data_pdfs_tmp_20_38 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 8*_stride_pdfs_tmp_3;
      double * RESTRICT  _data_pdfs_tmp_20_39 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 9*_stride_pdfs_tmp_3;
      double * RESTRICT  _data_pdfs_tmp_20_310 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 10*_stride_pdfs_tmp_3;
      double * RESTRICT  _data_pdfs_tmp_20_311 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 11*_stride_pdfs_tmp_3;
      double * RESTRICT  _data_pdfs_tmp_20_312 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 12*_stride_pdfs_tmp_3;
      double * RESTRICT  _data_pdfs_tmp_20_313 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 13*_stride_pdfs_tmp_3;
      double * RESTRICT  _data_pdfs_tmp_20_314 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 14*_stride_pdfs_tmp_3;
      double * RESTRICT  _data_pdfs_tmp_20_315 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 15*_stride_pdfs_tmp_3;
      double * RESTRICT  _data_pdfs_tmp_20_316 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 16*_stride_pdfs_tmp_3;
      double * RESTRICT  _data_pdfs_tmp_20_317 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 17*_stride_pdfs_tmp_3;
      double * RESTRICT  _data_pdfs_tmp_20_318 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 18*_stride_pdfs_tmp_3;
      double * RESTRICT  _data_pdfs_tmp_20_319 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 19*_stride_pdfs_tmp_3;
      double * RESTRICT  _data_pdfs_tmp_20_320 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 20*_stride_pdfs_tmp_3;
      double * RESTRICT  _data_pdfs_tmp_20_321 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 21*_stride_pdfs_tmp_3;
      double * RESTRICT  _data_pdfs_tmp_20_322 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 22*_stride_pdfs_tmp_3;
      double * RESTRICT  _data_pdfs_tmp_20_323 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 23*_stride_pdfs_tmp_3;
      double * RESTRICT  _data_pdfs_tmp_20_324 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 24*_stride_pdfs_tmp_3;
      double * RESTRICT  _data_pdfs_tmp_20_325 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 25*_stride_pdfs_tmp_3;
      double * RESTRICT  _data_pdfs_tmp_20_326 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 26*_stride_pdfs_tmp_3;
      for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_1; ctr_1 += 1)
      {
         double * RESTRICT _data_pdfs_2m1_321_11 = _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _data_pdfs_2m1_321;
         double * RESTRICT _data_pdfs_2m1_319_1m1 = _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _data_pdfs_2m1_319;
         double * RESTRICT _data_pdfs_2m1_314_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_2m1_314;
         double * RESTRICT _data_pdfs_20_310_11 = _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _data_pdfs_20_310;
         double * RESTRICT _data_pdfs_20_38_1m1 = _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _data_pdfs_20_38;
         double * RESTRICT _data_pdfs_21_325_11 = _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _data_pdfs_21_325;
         double * RESTRICT _data_pdfs_21_323_1m1 = _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _data_pdfs_21_323;
         double * RESTRICT _data_pdfs_21_318_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_21_318;
         double * RESTRICT _data_pdfs_20_34_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_34;
         double * RESTRICT _data_pdfs_2m1_320_1m1 = _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _data_pdfs_2m1_320;
         double * RESTRICT _data_pdfs_2m1_311_1m1 = _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _data_pdfs_2m1_311;
         double * RESTRICT _data_pdfs_20_37_1m1 = _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _data_pdfs_20_37;
         double * RESTRICT _data_pdfs_20_31_1m1 = _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _data_pdfs_20_31;
         double * RESTRICT _data_pdfs_21_324_1m1 = _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _data_pdfs_21_324;
         double * RESTRICT _data_pdfs_21_315_1m1 = _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _data_pdfs_21_315;
         double * RESTRICT _data_pdfs_2m1_313_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_2m1_313;
         double * RESTRICT _data_pdfs_2m1_322_11 = _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _data_pdfs_2m1_322;
         double * RESTRICT _data_pdfs_2m1_312_11 = _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _data_pdfs_2m1_312;
         double * RESTRICT _data_pdfs_2m1_35_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_2m1_35;
         double * RESTRICT _data_pdfs_20_30_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_30;
         double * RESTRICT _data_pdfs_20_33_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_33;
         double * RESTRICT _data_pdfs_20_39_11 = _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _data_pdfs_20_39;
         double * RESTRICT _data_pdfs_20_32_11 = _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _data_pdfs_20_32;
         double * RESTRICT _data_pdfs_21_317_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_21_317;
         double * RESTRICT _data_pdfs_21_326_11 = _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _data_pdfs_21_326;
         double * RESTRICT _data_pdfs_21_316_11 = _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _data_pdfs_21_316;
         double * RESTRICT _data_pdfs_21_36_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_21_36;
         double * RESTRICT  _data_pdfs_tmp_20_30_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_30;
         double * RESTRICT  _data_pdfs_tmp_20_31_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_31;
         double * RESTRICT  _data_pdfs_tmp_20_32_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_32;
         double * RESTRICT  _data_pdfs_tmp_20_33_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_33;
         double * RESTRICT  _data_pdfs_tmp_20_34_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_34;
         double * RESTRICT  _data_pdfs_tmp_20_35_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_35;
         double * RESTRICT  _data_pdfs_tmp_20_36_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_36;
         double * RESTRICT  _data_pdfs_tmp_20_37_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_37;
         double * RESTRICT  _data_pdfs_tmp_20_38_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_38;
         double * RESTRICT  _data_pdfs_tmp_20_39_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_39;
         double * RESTRICT  _data_pdfs_tmp_20_310_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_310;
         double * RESTRICT  _data_pdfs_tmp_20_311_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_311;
         double * RESTRICT  _data_pdfs_tmp_20_312_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_312;
         double * RESTRICT  _data_pdfs_tmp_20_313_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_313;
         double * RESTRICT  _data_pdfs_tmp_20_314_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_314;
         double * RESTRICT  _data_pdfs_tmp_20_315_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_315;
         double * RESTRICT  _data_pdfs_tmp_20_316_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_316;
         double * RESTRICT  _data_pdfs_tmp_20_317_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_317;
         double * RESTRICT  _data_pdfs_tmp_20_318_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_318;
         double * RESTRICT  _data_pdfs_tmp_20_319_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_319;
         double * RESTRICT  _data_pdfs_tmp_20_320_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_320;
         double * RESTRICT  _data_pdfs_tmp_20_321_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_321;
         double * RESTRICT  _data_pdfs_tmp_20_322_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_322;
         double * RESTRICT  _data_pdfs_tmp_20_323_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_323;
         double * RESTRICT  _data_pdfs_tmp_20_324_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_324;
         double * RESTRICT  _data_pdfs_tmp_20_325_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_325;
         double * RESTRICT  _data_pdfs_tmp_20_326_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_326;
         for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_0; ctr_0 += 1)
         {
            const double vel0Term = _data_pdfs_20_310_11[_stride_pdfs_0*ctr_0 - _stride_pdfs_0] + _data_pdfs_20_34_10[_stride_pdfs_0*ctr_0 - _stride_pdfs_0] + _data_pdfs_20_38_1m1[_stride_pdfs_0*ctr_0 - _stride_pdfs_0] + _data_pdfs_21_318_10[_stride_pdfs_0*ctr_0 - _stride_pdfs_0] + _data_pdfs_21_323_1m1[_stride_pdfs_0*ctr_0 - _stride_pdfs_0] + _data_pdfs_21_325_11[_stride_pdfs_0*ctr_0 - _stride_pdfs_0] + _data_pdfs_2m1_314_10[_stride_pdfs_0*ctr_0 - _stride_pdfs_0] + _data_pdfs_2m1_319_1m1[_stride_pdfs_0*ctr_0 - _stride_pdfs_0] + _data_pdfs_2m1_321_11[_stride_pdfs_0*ctr_0 - _stride_pdfs_0];
            const double vel1Term = _data_pdfs_20_31_1m1[_stride_pdfs_0*ctr_0] + _data_pdfs_20_37_1m1[_stride_pdfs_0*ctr_0 + _stride_pdfs_0] + _data_pdfs_21_315_1m1[_stride_pdfs_0*ctr_0] + _data_pdfs_21_324_1m1[_stride_pdfs_0*ctr_0 + _stride_pdfs_0] + _data_pdfs_2m1_311_1m1[_stride_pdfs_0*ctr_0] + _data_pdfs_2m1_320_1m1[_stride_pdfs_0*ctr_0 + _stride_pdfs_0];
            const double vel2Term = _data_pdfs_2m1_312_11[_stride_pdfs_0*ctr_0] + _data_pdfs_2m1_313_10[_stride_pdfs_0*ctr_0 + _stride_pdfs_0] + _data_pdfs_2m1_322_11[_stride_pdfs_0*ctr_0 + _stride_pdfs_0] + _data_pdfs_2m1_35_10[_stride_pdfs_0*ctr_0];
            const double delta_rho = vel0Term + vel1Term + vel2Term + _data_pdfs_20_30_10[_stride_pdfs_0*ctr_0] + _data_pdfs_20_32_11[_stride_pdfs_0*ctr_0] + _data_pdfs_20_33_10[_stride_pdfs_0*ctr_0 + _stride_pdfs_0] + _data_pdfs_20_39_11[_stride_pdfs_0*ctr_0 + _stride_pdfs_0] + _data_pdfs_21_316_11[_stride_pdfs_0*ctr_0] + _data_pdfs_21_317_10[_stride_pdfs_0*ctr_0 + _stride_pdfs_0] + _data_pdfs_21_326_11[_stride_pdfs_0*ctr_0 + _stride_pdfs_0] + _data_pdfs_21_36_10[_stride_pdfs_0*ctr_0];
            const double u_0 = vel0Term - 1.0*_data_pdfs_20_33_10[_stride_pdfs_0*ctr_0 + _stride_pdfs_0] - 1.0*_data_pdfs_20_37_1m1[_stride_pdfs_0*ctr_0 + _stride_pdfs_0] - 1.0*_data_pdfs_20_39_11[_stride_pdfs_0*ctr_0 + _stride_pdfs_0] - 1.0*_data_pdfs_21_317_10[_stride_pdfs_0*ctr_0 + _stride_pdfs_0] - 1.0*_data_pdfs_21_324_1m1[_stride_pdfs_0*ctr_0 + _stride_pdfs_0] - 1.0*_data_pdfs_21_326_11[_stride_pdfs_0*ctr_0 + _stride_pdfs_0] - 1.0*_data_pdfs_2m1_313_10[_stride_pdfs_0*ctr_0 + _stride_pdfs_0] - 1.0*_data_pdfs_2m1_320_1m1[_stride_pdfs_0*ctr_0 + _stride_pdfs_0] - 1.0*_data_pdfs_2m1_322_11[_stride_pdfs_0*ctr_0 + _stride_pdfs_0];
            const double u_1 = vel1Term - 1.0*_data_pdfs_20_310_11[_stride_pdfs_0*ctr_0 - _stride_pdfs_0] - 1.0*_data_pdfs_20_32_11[_stride_pdfs_0*ctr_0] - 1.0*_data_pdfs_20_39_11[_stride_pdfs_0*ctr_0 + _stride_pdfs_0] - 1.0*_data_pdfs_21_316_11[_stride_pdfs_0*ctr_0] - 1.0*_data_pdfs_21_325_11[_stride_pdfs_0*ctr_0 - _stride_pdfs_0] - 1.0*_data_pdfs_21_326_11[_stride_pdfs_0*ctr_0 + _stride_pdfs_0] - 1.0*_data_pdfs_2m1_312_11[_stride_pdfs_0*ctr_0] - 1.0*_data_pdfs_2m1_321_11[_stride_pdfs_0*ctr_0 - _stride_pdfs_0] - 1.0*_data_pdfs_2m1_322_11[_stride_pdfs_0*ctr_0 + _stride_pdfs_0] + _data_pdfs_20_38_1m1[_stride_pdfs_0*ctr_0 - _stride_pdfs_0] + _data_pdfs_21_323_1m1[_stride_pdfs_0*ctr_0 - _stride_pdfs_0] + _data_pdfs_2m1_319_1m1[_stride_pdfs_0*ctr_0 - _stride_pdfs_0];
            const double u_2 = vel2Term - 1.0*_data_pdfs_21_315_1m1[_stride_pdfs_0*ctr_0] - 1.0*_data_pdfs_21_316_11[_stride_pdfs_0*ctr_0] - 1.0*_data_pdfs_21_317_10[_stride_pdfs_0*ctr_0 + _stride_pdfs_0] - 1.0*_data_pdfs_21_318_10[_stride_pdfs_0*ctr_0 - _stride_pdfs_0] - 1.0*_data_pdfs_21_323_1m1[_stride_pdfs_0*ctr_0 - _stride_pdfs_0] - 1.0*_data_pdfs_21_324_1m1[_stride_pdfs_0*ctr_0 + _stride_pdfs_0] - 1.0*_data_pdfs_21_325_11[_stride_pdfs_0*ctr_0 - _stride_pdfs_0] - 1.0*_data_pdfs_21_326_11[_stride_pdfs_0*ctr_0 + _stride_pdfs_0] - 1.0*_data_pdfs_21_36_10[_stride_pdfs_0*ctr_0] + _data_pdfs_2m1_311_1m1[_stride_pdfs_0*ctr_0] + _data_pdfs_2m1_314_10[_stride_pdfs_0*ctr_0 - _stride_pdfs_0] + _data_pdfs_2m1_319_1m1[_stride_pdfs_0*ctr_0 - _stride_pdfs_0] + _data_pdfs_2m1_320_1m1[_stride_pdfs_0*ctr_0 + _stride_pdfs_0] + _data_pdfs_2m1_321_11[_stride_pdfs_0*ctr_0 - _stride_pdfs_0];
            const double u0Mu1 = u_0 + u_1*-1.0;
            const double u0Pu1 = u_0 + u_1;
            const double u1Pu2 = u_1 + u_2;
            const double u1Mu2 = u_1 + u_2*-1.0;
            const double u0Mu2 = u_0 + u_2*-1.0;
            const double u0Pu2 = u_0 + u_2;
            const double f_eq_common = delta_rho - 1.5*(u_0*u_0) - 1.5*(u_1*u_1) - 1.5*(u_2*u_2);
            _data_pdfs_tmp_20_30_10[_stride_pdfs_tmp_0*ctr_0] = omega*(f_eq_common*0.29629629629629628 - 1.0*_data_pdfs_20_30_10[_stride_pdfs_0*ctr_0]) + _data_pdfs_20_30_10[_stride_pdfs_0*ctr_0];
            _data_pdfs_tmp_20_31_10[_stride_pdfs_tmp_0*ctr_0] = omega*(f_eq_common*0.07407407407407407 + u_1*0.22222222222222221 - 1.0*_data_pdfs_20_31_1m1[_stride_pdfs_0*ctr_0] + 0.33333333333333331*(u_1*u_1)) + _data_pdfs_20_31_1m1[_stride_pdfs_0*ctr_0];
            _data_pdfs_tmp_20_32_10[_stride_pdfs_tmp_0*ctr_0] = omega*(f_eq_common*0.07407407407407407 + u_1*-0.22222222222222221 - 1.0*_data_pdfs_20_32_11[_stride_pdfs_0*ctr_0] + 0.33333333333333331*(u_1*u_1)) + _data_pdfs_20_32_11[_stride_pdfs_0*ctr_0];
            _data_pdfs_tmp_20_33_10[_stride_pdfs_tmp_0*ctr_0] = omega*(f_eq_common*0.07407407407407407 + u_0*-0.22222222222222221 - 1.0*_data_pdfs_20_33_10[_stride_pdfs_0*ctr_0 + _stride_pdfs_0] + 0.33333333333333331*(u_0*u_0)) + _data_pdfs_20_33_10[_stride_pdfs_0*ctr_0 + _stride_pdfs_0];
            _data_pdfs_tmp_20_34_10[_stride_pdfs_tmp_0*ctr_0] = omega*(f_eq_common*0.07407407407407407 + u_0*0.22222222222222221 - 1.0*_data_pdfs_20_34_10[_stride_pdfs_0*ctr_0 - _stride_pdfs_0] + 0.33333333333333331*(u_0*u_0)) + _data_pdfs_20_34_10[_stride_pdfs_0*ctr_0 - _stride_pdfs_0];
            _data_pdfs_tmp_20_35_10[_stride_pdfs_tmp_0*ctr_0] = omega*(f_eq_common*0.07407407407407407 + u_2*0.22222222222222221 - 1.0*_data_pdfs_2m1_35_10[_stride_pdfs_0*ctr_0] + 0.33333333333333331*(u_2*u_2)) + _data_pdfs_2m1_35_10[_stride_pdfs_0*ctr_0];
            _data_pdfs_tmp_20_36_10[_stride_pdfs_tmp_0*ctr_0] = omega*(f_eq_common*0.07407407407407407 + u_2*-0.22222222222222221 - 1.0*_data_pdfs_21_36_10[_stride_pdfs_0*ctr_0] + 0.33333333333333331*(u_2*u_2)) + _data_pdfs_21_36_10[_stride_pdfs_0*ctr_0];
            _data_pdfs_tmp_20_37_10[_stride_pdfs_tmp_0*ctr_0] = omega*(f_eq_common*0.018518518518518517 + u0Mu1*-0.055555555555555552 - 1.0*_data_pdfs_20_37_1m1[_stride_pdfs_0*ctr_0 + _stride_pdfs_0] + 0.083333333333333329*(u0Mu1*u0Mu1)) + _data_pdfs_20_37_1m1[_stride_pdfs_0*ctr_0 + _stride_pdfs_0];
            _data_pdfs_tmp_20_38_10[_stride_pdfs_tmp_0*ctr_0] = omega*(f_eq_common*0.018518518518518517 + u0Pu1*0.055555555555555552 - 1.0*_data_pdfs_20_38_1m1[_stride_pdfs_0*ctr_0 - _stride_pdfs_0] + 0.083333333333333329*(u0Pu1*u0Pu1)) + _data_pdfs_20_38_1m1[_stride_pdfs_0*ctr_0 - _stride_pdfs_0];
            _data_pdfs_tmp_20_39_10[_stride_pdfs_tmp_0*ctr_0] = omega*(f_eq_common*0.018518518518518517 + u0Pu1*-0.055555555555555552 - 1.0*_data_pdfs_20_39_11[_stride_pdfs_0*ctr_0 + _stride_pdfs_0] + 0.083333333333333329*(u0Pu1*u0Pu1)) + _data_pdfs_20_39_11[_stride_pdfs_0*ctr_0 + _stride_pdfs_0];
            _data_pdfs_tmp_20_310_10[_stride_pdfs_tmp_0*ctr_0] = omega*(f_eq_common*0.018518518518518517 + u0Mu1*0.055555555555555552 - 1.0*_data_pdfs_20_310_11[_stride_pdfs_0*ctr_0 - _stride_pdfs_0] + 0.083333333333333329*(u0Mu1*u0Mu1)) + _data_pdfs_20_310_11[_stride_pdfs_0*ctr_0 - _stride_pdfs_0];
            _data_pdfs_tmp_20_311_10[_stride_pdfs_tmp_0*ctr_0] = omega*(f_eq_common*0.018518518518518517 + u1Pu2*0.055555555555555552 - 1.0*_data_pdfs_2m1_311_1m1[_stride_pdfs_0*ctr_0] + 0.083333333333333329*(u1Pu2*u1Pu2)) + _data_pdfs_2m1_311_1m1[_stride_pdfs_0*ctr_0];
            _data_pdfs_tmp_20_312_10[_stride_pdfs_tmp_0*ctr_0] = omega*(f_eq_common*0.018518518518518517 + u1Mu2*-0.055555555555555552 - 1.0*_data_pdfs_2m1_312_11[_stride_pdfs_0*ctr_0] + 0.083333333333333329*(u1Mu2*u1Mu2)) + _data_pdfs_2m1_312_11[_stride_pdfs_0*ctr_0];
            _data_pdfs_tmp_20_313_10[_stride_pdfs_tmp_0*ctr_0] = omega*(f_eq_common*0.018518518518518517 + u0Mu2*-0.055555555555555552 - 1.0*_data_pdfs_2m1_313_10[_stride_pdfs_0*ctr_0 + _stride_pdfs_0] + 0.083333333333333329*(u0Mu2*u0Mu2)) + _data_pdfs_2m1_313_10[_stride_pdfs_0*ctr_0 + _stride_pdfs_0];
            _data_pdfs_tmp_20_314_10[_stride_pdfs_tmp_0*ctr_0] = omega*(f_eq_common*0.018518518518518517 + u0Pu2*0.055555555555555552 - 1.0*_data_pdfs_2m1_314_10[_stride_pdfs_0*ctr_0 - _stride_pdfs_0] + 0.083333333333333329*(u0Pu2*u0Pu2)) + _data_pdfs_2m1_314_10[_stride_pdfs_0*ctr_0 - _stride_pdfs_0];
            _data_pdfs_tmp_20_315_10[_stride_pdfs_tmp_0*ctr_0] = omega*(f_eq_common*0.018518518518518517 + u1Mu2*0.055555555555555552 - 1.0*_data_pdfs_21_315_1m1[_stride_pdfs_0*ctr_0] + 0.083333333333333329*(u1Mu2*u1Mu2)) + _data_pdfs_21_315_1m1[_stride_pdfs_0*ctr_0];
            _data_pdfs_tmp_20_316_10[_stride_pdfs_tmp_0*ctr_0] = omega*(f_eq_common*0.018518518518518517 + u1Pu2*-0.055555555555555552 - 1.0*_data_pdfs_21_316_11[_stride_pdfs_0*ctr_0] + 0.083333333333333329*(u1Pu2*u1Pu2)) + _data_pdfs_21_316_11[_stride_pdfs_0*ctr_0];
            _data_pdfs_tmp_20_317_10[_stride_pdfs_tmp_0*ctr_0] = omega*(f_eq_common*0.018518518518518517 + u0Pu2*-0.055555555555555552 - 1.0*_data_pdfs_21_317_10[_stride_pdfs_0*ctr_0 + _stride_pdfs_0] + 0.083333333333333329*(u0Pu2*u0Pu2)) + _data_pdfs_21_317_10[_stride_pdfs_0*ctr_0 + _stride_pdfs_0];
            _data_pdfs_tmp_20_318_10[_stride_pdfs_tmp_0*ctr_0] = omega*(f_eq_common*0.018518518518518517 + u0Mu2*0.055555555555555552 - 1.0*_data_pdfs_21_318_10[_stride_pdfs_0*ctr_0 - _stride_pdfs_0] + 0.083333333333333329*(u0Mu2*u0Mu2)) + _data_pdfs_21_318_10[_stride_pdfs_0*ctr_0 - _stride_pdfs_0];
            _data_pdfs_tmp_20_319_10[_stride_pdfs_tmp_0*ctr_0] = omega*(delta_rho*-0.013888888888888888 + f_eq_common*0.018518518518518517 + u0Pu1*0.013888888888888888 + u_2*0.013888888888888888 - 1.0*_data_pdfs_2m1_319_1m1[_stride_pdfs_0*ctr_0 - _stride_pdfs_0] + 0.020833333333333332*(u0Pu1*u0Pu1) + 0.020833333333333332*(u0Pu2*u0Pu2) + 0.020833333333333332*(u1Pu2*u1Pu2)) + _data_pdfs_2m1_319_1m1[_stride_pdfs_0*ctr_0 - _stride_pdfs_0];
            _data_pdfs_tmp_20_320_10[_stride_pdfs_tmp_0*ctr_0] = omega*(delta_rho*-0.013888888888888888 + f_eq_common*0.018518518518518517 + u0Mu1*-0.013888888888888888 + u_2*0.013888888888888888 - 1.0*_data_pdfs_2m1_320_1m1[_stride_pdfs_0*ctr_0 + _stride_pdfs_0] + 0.020833333333333332*(u0Mu1*u0Mu1) + 0.020833333333333332*(u0Mu2*u0Mu2) + 0.020833333333333332*(u1Pu2*u1Pu2)) + _data_pdfs_2m1_320_1m1[_stride_pdfs_0*ctr_0 + _stride_pdfs_0];
            _data_pdfs_tmp_20_321_10[_stride_pdfs_tmp_0*ctr_0] = omega*(delta_rho*-0.013888888888888888 + f_eq_common*0.018518518518518517 + u0Mu1*0.013888888888888888 + u_2*0.013888888888888888 - 1.0*_data_pdfs_2m1_321_11[_stride_pdfs_0*ctr_0 - _stride_pdfs_0] + 0.020833333333333332*(u0Mu1*u0Mu1) + 0.020833333333333332*(u0Pu2*u0Pu2) + 0.020833333333333332*(u1Mu2*u1Mu2)) + _data_pdfs_2m1_321_11[_stride_pdfs_0*ctr_0 - _stride_pdfs_0];
            _data_pdfs_tmp_20_322_10[_stride_pdfs_tmp_0*ctr_0] = omega*(delta_rho*-0.013888888888888888 + f_eq_common*0.018518518518518517 + u0Pu1*-0.013888888888888888 + u_2*0.013888888888888888 - 1.0*_data_pdfs_2m1_322_11[_stride_pdfs_0*ctr_0 + _stride_pdfs_0] + 0.020833333333333332*(u0Mu2*u0Mu2) + 0.020833333333333332*(u0Pu1*u0Pu1) + 0.020833333333333332*(u1Mu2*u1Mu2)) + _data_pdfs_2m1_322_11[_stride_pdfs_0*ctr_0 + _stride_pdfs_0];
            _data_pdfs_tmp_20_323_10[_stride_pdfs_tmp_0*ctr_0] = omega*(delta_rho*-0.013888888888888888 + f_eq_common*0.018518518518518517 + u0Pu1*0.013888888888888888 + u_2*-0.013888888888888888 - 1.0*_data_pdfs_21_323_1m1[_stride_pdfs_0*ctr_0 - _stride_pdfs_0] + 0.020833333333333332*(u0Mu2*u0Mu2) + 0.020833333333333332*(u0Pu1*u0Pu1) + 0.020833333333333332*(u1Mu2*u1Mu2)) + _data_pdfs_21_323_1m1[_stride_pdfs_0*ctr_0 - _stride_pdfs_0];
            _data_pdfs_tmp_20_324_10[_stride_pdfs_tmp_0*ctr_0] = omega*(delta_rho*-0.013888888888888888 + f_eq_common*0.018518518518518517 + u0Mu1*-0.013888888888888888 + u_2*-0.013888888888888888 - 1.0*_data_pdfs_21_324_1m1[_stride_pdfs_0*ctr_0 + _stride_pdfs_0] + 0.020833333333333332*(u0Mu1*u0Mu1) + 0.020833333333333332*(u0Pu2*u0Pu2) + 0.020833333333333332*(u1Mu2*u1Mu2)) + _data_pdfs_21_324_1m1[_stride_pdfs_0*ctr_0 + _stride_pdfs_0];
            _data_pdfs_tmp_20_325_10[_stride_pdfs_tmp_0*ctr_0] = omega*(delta_rho*-0.013888888888888888 + f_eq_common*0.018518518518518517 + u0Mu1*0.013888888888888888 + u_2*-0.013888888888888888 - 1.0*_data_pdfs_21_325_11[_stride_pdfs_0*ctr_0 - _stride_pdfs_0] + 0.020833333333333332*(u0Mu1*u0Mu1) + 0.020833333333333332*(u0Mu2*u0Mu2) + 0.020833333333333332*(u1Pu2*u1Pu2)) + _data_pdfs_21_325_11[_stride_pdfs_0*ctr_0 - _stride_pdfs_0];
            _data_pdfs_tmp_20_326_10[_stride_pdfs_tmp_0*ctr_0] = omega*(delta_rho*-0.013888888888888888 + f_eq_common*0.018518518518518517 + u0Pu1*-0.013888888888888888 + u_2*-0.013888888888888888 - 1.0*_data_pdfs_21_326_11[_stride_pdfs_0*ctr_0 + _stride_pdfs_0] + 0.020833333333333332*(u0Pu1*u0Pu1) + 0.020833333333333332*(u0Pu2*u0Pu2) + 0.020833333333333332*(u1Pu2*u1Pu2)) + _data_pdfs_21_326_11[_stride_pdfs_0*ctr_0 + _stride_pdfs_0];
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
      double * RESTRICT  _data_pdfs_20_313 = _data_pdfs + _stride_pdfs_2*ctr_2 + 13*_stride_pdfs_3;
      double * RESTRICT  _data_pdfs_20_36 = _data_pdfs + _stride_pdfs_2*ctr_2 + 6*_stride_pdfs_3;
      double * RESTRICT  _data_pdfs_20_326 = _data_pdfs + _stride_pdfs_2*ctr_2 + 26*_stride_pdfs_3;
      double * RESTRICT  _data_pdfs_20_323 = _data_pdfs + _stride_pdfs_2*ctr_2 + 23*_stride_pdfs_3;
      double * RESTRICT  _data_pdfs_20_310 = _data_pdfs + _stride_pdfs_2*ctr_2 + 10*_stride_pdfs_3;
      double * RESTRICT  _data_pdfs_20_316 = _data_pdfs + _stride_pdfs_2*ctr_2 + 16*_stride_pdfs_3;
      double * RESTRICT  _data_pdfs_20_312 = _data_pdfs + _stride_pdfs_2*ctr_2 + 12*_stride_pdfs_3;
      double * RESTRICT  _data_pdfs_20_321 = _data_pdfs + _stride_pdfs_2*ctr_2 + 21*_stride_pdfs_3;
      double * RESTRICT  _data_pdfs_20_311 = _data_pdfs + _stride_pdfs_2*ctr_2 + 11*_stride_pdfs_3;
      double * RESTRICT  _data_pdfs_20_324 = _data_pdfs + _stride_pdfs_2*ctr_2 + 24*_stride_pdfs_3;
      double * RESTRICT  _data_pdfs_20_37 = _data_pdfs + _stride_pdfs_2*ctr_2 + 7*_stride_pdfs_3;
      double * RESTRICT  _data_pdfs_20_319 = _data_pdfs + _stride_pdfs_2*ctr_2 + 19*_stride_pdfs_3;
      double * RESTRICT  _data_pdfs_20_31 = _data_pdfs + _stride_pdfs_2*ctr_2 + _stride_pdfs_3;
      double * RESTRICT  _data_pdfs_20_38 = _data_pdfs + _stride_pdfs_2*ctr_2 + 8*_stride_pdfs_3;
      double * RESTRICT  _data_pdfs_20_325 = _data_pdfs + _stride_pdfs_2*ctr_2 + 25*_stride_pdfs_3;
      double * RESTRICT  _data_pdfs_20_32 = _data_pdfs + _stride_pdfs_2*ctr_2 + 2*_stride_pdfs_3;
      double * RESTRICT  _data_pdfs_20_320 = _data_pdfs + _stride_pdfs_2*ctr_2 + 20*_stride_pdfs_3;
      double * RESTRICT  _data_pdfs_20_39 = _data_pdfs + _stride_pdfs_2*ctr_2 + 9*_stride_pdfs_3;
      double * RESTRICT  _data_pdfs_20_35 = _data_pdfs + _stride_pdfs_2*ctr_2 + 5*_stride_pdfs_3;
      double * RESTRICT  _data_pdfs_20_30 = _data_pdfs + _stride_pdfs_2*ctr_2;
      double * RESTRICT  _data_pdfs_20_315 = _data_pdfs + _stride_pdfs_2*ctr_2 + 15*_stride_pdfs_3;
      double * RESTRICT  _data_pdfs_20_322 = _data_pdfs + _stride_pdfs_2*ctr_2 + 22*_stride_pdfs_3;
      double * RESTRICT  _data_pdfs_20_317 = _data_pdfs + _stride_pdfs_2*ctr_2 + 17*_stride_pdfs_3;
      double * RESTRICT  _data_pdfs_20_314 = _data_pdfs + _stride_pdfs_2*ctr_2 + 14*_stride_pdfs_3;
      double * RESTRICT  _data_pdfs_20_318 = _data_pdfs + _stride_pdfs_2*ctr_2 + 18*_stride_pdfs_3;
      double * RESTRICT  _data_pdfs_20_34 = _data_pdfs + _stride_pdfs_2*ctr_2 + 4*_stride_pdfs_3;
      double * RESTRICT  _data_pdfs_20_33 = _data_pdfs + _stride_pdfs_2*ctr_2 + 3*_stride_pdfs_3;
      for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_1; ctr_1 += 1)
      {
         double * RESTRICT  _data_pdfs_20_313_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_313;
         double * RESTRICT  _data_pdfs_20_36_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_36;
         double * RESTRICT  _data_pdfs_20_326_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_326;
         double * RESTRICT  _data_pdfs_20_323_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_323;
         double * RESTRICT  _data_pdfs_20_310_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_310;
         double * RESTRICT  _data_pdfs_20_316_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_316;
         double * RESTRICT  _data_pdfs_20_312_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_312;
         double * RESTRICT  _data_pdfs_20_321_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_321;
         double * RESTRICT  _data_pdfs_20_311_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_311;
         double * RESTRICT  _data_pdfs_20_324_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_324;
         double * RESTRICT  _data_pdfs_20_37_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_37;
         double * RESTRICT  _data_pdfs_20_319_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_319;
         double * RESTRICT  _data_pdfs_20_31_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_31;
         double * RESTRICT  _data_pdfs_20_38_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_38;
         double * RESTRICT  _data_pdfs_20_325_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_325;
         double * RESTRICT  _data_pdfs_20_32_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_32;
         double * RESTRICT  _data_pdfs_20_320_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_320;
         double * RESTRICT  _data_pdfs_20_39_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_39;
         double * RESTRICT  _data_pdfs_20_35_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_35;
         double * RESTRICT  _data_pdfs_20_30_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_30;
         double * RESTRICT  _data_pdfs_20_315_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_315;
         double * RESTRICT  _data_pdfs_20_322_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_322;
         double * RESTRICT  _data_pdfs_20_317_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_317;
         double * RESTRICT  _data_pdfs_20_314_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_314;
         double * RESTRICT  _data_pdfs_20_318_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_318;
         double * RESTRICT  _data_pdfs_20_34_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_34;
         double * RESTRICT  _data_pdfs_20_33_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_33;
         for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_0; ctr_0 += 1)
         {
            const double xi_1 = _data_pdfs_20_313_10[_stride_pdfs_0*ctr_0];
            const double xi_2 = _data_pdfs_20_36_10[_stride_pdfs_0*ctr_0];
            const double xi_3 = _data_pdfs_20_326_10[_stride_pdfs_0*ctr_0];
            const double xi_4 = _data_pdfs_20_323_10[_stride_pdfs_0*ctr_0];
            const double xi_5 = _data_pdfs_20_310_10[_stride_pdfs_0*ctr_0];
            const double xi_6 = _data_pdfs_20_316_10[_stride_pdfs_0*ctr_0];
            const double xi_7 = _data_pdfs_20_312_10[_stride_pdfs_0*ctr_0];
            const double xi_8 = _data_pdfs_20_321_10[_stride_pdfs_0*ctr_0];
            const double xi_9 = _data_pdfs_20_311_10[_stride_pdfs_0*ctr_0];
            const double xi_10 = _data_pdfs_20_324_10[_stride_pdfs_0*ctr_0];
            const double xi_11 = _data_pdfs_20_37_10[_stride_pdfs_0*ctr_0];
            const double xi_12 = _data_pdfs_20_319_10[_stride_pdfs_0*ctr_0];
            const double xi_13 = _data_pdfs_20_31_10[_stride_pdfs_0*ctr_0];
            const double xi_14 = _data_pdfs_20_38_10[_stride_pdfs_0*ctr_0];
            const double xi_15 = _data_pdfs_20_325_10[_stride_pdfs_0*ctr_0];
            const double xi_16 = _data_pdfs_20_32_10[_stride_pdfs_0*ctr_0];
            const double xi_17 = _data_pdfs_20_320_10[_stride_pdfs_0*ctr_0];
            const double xi_18 = _data_pdfs_20_39_10[_stride_pdfs_0*ctr_0];
            const double xi_19 = _data_pdfs_20_35_10[_stride_pdfs_0*ctr_0];
            const double xi_20 = _data_pdfs_20_30_10[_stride_pdfs_0*ctr_0];
            const double xi_21 = _data_pdfs_20_315_10[_stride_pdfs_0*ctr_0];
            const double xi_22 = _data_pdfs_20_322_10[_stride_pdfs_0*ctr_0];
            const double xi_23 = _data_pdfs_20_317_10[_stride_pdfs_0*ctr_0];
            const double xi_24 = _data_pdfs_20_314_10[_stride_pdfs_0*ctr_0];
            const double xi_25 = _data_pdfs_20_318_10[_stride_pdfs_0*ctr_0];
            const double xi_26 = _data_pdfs_20_34_10[_stride_pdfs_0*ctr_0];
            const double xi_27 = _data_pdfs_20_33_10[_stride_pdfs_0*ctr_0];
            const double vel0Term = xi_12 + xi_14 + xi_15 + xi_24 + xi_25 + xi_26 + xi_4 + xi_5 + xi_8;
            const double vel1Term = xi_10 + xi_11 + xi_13 + xi_17 + xi_21 + xi_9;
            const double vel2Term = xi_1 + xi_19 + xi_22 + xi_7;
            const double delta_rho = vel0Term + vel1Term + vel2Term + xi_16 + xi_18 + xi_2 + xi_20 + xi_23 + xi_27 + xi_3 + xi_6;
            const double u_0 = vel0Term + xi_1*-1.0 + xi_10*-1.0 + xi_11*-1.0 + xi_17*-1.0 + xi_18*-1.0 + xi_22*-1.0 + xi_23*-1.0 + xi_27*-1.0 + xi_3*-1.0;
            const double u_1 = vel1Term + xi_12 + xi_14 + xi_15*-1.0 + xi_16*-1.0 + xi_18*-1.0 + xi_22*-1.0 + xi_3*-1.0 + xi_4 + xi_5*-1.0 + xi_6*-1.0 + xi_7*-1.0 + xi_8*-1.0;
            const double u_2 = vel2Term + xi_10*-1.0 + xi_12 + xi_15*-1.0 + xi_17 + xi_2*-1.0 + xi_21*-1.0 + xi_23*-1.0 + xi_24 + xi_25*-1.0 + xi_3*-1.0 + xi_4*-1.0 + xi_6*-1.0 + xi_8 + xi_9;
            const double u0Mu1 = u_0 + u_1*-1.0;
            const double u0Pu1 = u_0 + u_1;
            const double u1Pu2 = u_1 + u_2;
            const double u1Mu2 = u_1 + u_2*-1.0;
            const double u0Mu2 = u_0 + u_2*-1.0;
            const double u0Pu2 = u_0 + u_2;
            const double f_eq_common = delta_rho - 1.5*(u_0*u_0) - 1.5*(u_1*u_1) - 1.5*(u_2*u_2);
            _data_pdfs_20_30_10[_stride_pdfs_0*ctr_0] = omega*(f_eq_common*0.29629629629629628 + xi_20*-1.0) + xi_20;
            _data_pdfs_20_31_10[_stride_pdfs_0*ctr_0] = omega*(f_eq_common*0.07407407407407407 + u_1*0.22222222222222221 + xi_13*-1.0 + 0.33333333333333331*(u_1*u_1)) + xi_13;
            _data_pdfs_20_32_10[_stride_pdfs_0*ctr_0] = omega*(f_eq_common*0.07407407407407407 + u_1*-0.22222222222222221 + xi_16*-1.0 + 0.33333333333333331*(u_1*u_1)) + xi_16;
            _data_pdfs_20_33_10[_stride_pdfs_0*ctr_0] = omega*(f_eq_common*0.07407407407407407 + u_0*-0.22222222222222221 + xi_27*-1.0 + 0.33333333333333331*(u_0*u_0)) + xi_27;
            _data_pdfs_20_34_10[_stride_pdfs_0*ctr_0] = omega*(f_eq_common*0.07407407407407407 + u_0*0.22222222222222221 + xi_26*-1.0 + 0.33333333333333331*(u_0*u_0)) + xi_26;
            _data_pdfs_20_35_10[_stride_pdfs_0*ctr_0] = omega*(f_eq_common*0.07407407407407407 + u_2*0.22222222222222221 + xi_19*-1.0 + 0.33333333333333331*(u_2*u_2)) + xi_19;
            _data_pdfs_20_36_10[_stride_pdfs_0*ctr_0] = omega*(f_eq_common*0.07407407407407407 + u_2*-0.22222222222222221 + xi_2*-1.0 + 0.33333333333333331*(u_2*u_2)) + xi_2;
            _data_pdfs_20_37_10[_stride_pdfs_0*ctr_0] = omega*(f_eq_common*0.018518518518518517 + u0Mu1*-0.055555555555555552 + xi_11*-1.0 + 0.083333333333333329*(u0Mu1*u0Mu1)) + xi_11;
            _data_pdfs_20_38_10[_stride_pdfs_0*ctr_0] = omega*(f_eq_common*0.018518518518518517 + u0Pu1*0.055555555555555552 + xi_14*-1.0 + 0.083333333333333329*(u0Pu1*u0Pu1)) + xi_14;
            _data_pdfs_20_39_10[_stride_pdfs_0*ctr_0] = omega*(f_eq_common*0.018518518518518517 + u0Pu1*-0.055555555555555552 + xi_18*-1.0 + 0.083333333333333329*(u0Pu1*u0Pu1)) + xi_18;
            _data_pdfs_20_310_10[_stride_pdfs_0*ctr_0] = omega*(f_eq_common*0.018518518518518517 + u0Mu1*0.055555555555555552 + xi_5*-1.0 + 0.083333333333333329*(u0Mu1*u0Mu1)) + xi_5;
            _data_pdfs_20_311_10[_stride_pdfs_0*ctr_0] = omega*(f_eq_common*0.018518518518518517 + u1Pu2*0.055555555555555552 + xi_9*-1.0 + 0.083333333333333329*(u1Pu2*u1Pu2)) + xi_9;
            _data_pdfs_20_312_10[_stride_pdfs_0*ctr_0] = omega*(f_eq_common*0.018518518518518517 + u1Mu2*-0.055555555555555552 + xi_7*-1.0 + 0.083333333333333329*(u1Mu2*u1Mu2)) + xi_7;
            _data_pdfs_20_313_10[_stride_pdfs_0*ctr_0] = omega*(f_eq_common*0.018518518518518517 + u0Mu2*-0.055555555555555552 + xi_1*-1.0 + 0.083333333333333329*(u0Mu2*u0Mu2)) + xi_1;
            _data_pdfs_20_314_10[_stride_pdfs_0*ctr_0] = omega*(f_eq_common*0.018518518518518517 + u0Pu2*0.055555555555555552 + xi_24*-1.0 + 0.083333333333333329*(u0Pu2*u0Pu2)) + xi_24;
            _data_pdfs_20_315_10[_stride_pdfs_0*ctr_0] = omega*(f_eq_common*0.018518518518518517 + u1Mu2*0.055555555555555552 + xi_21*-1.0 + 0.083333333333333329*(u1Mu2*u1Mu2)) + xi_21;
            _data_pdfs_20_316_10[_stride_pdfs_0*ctr_0] = omega*(f_eq_common*0.018518518518518517 + u1Pu2*-0.055555555555555552 + xi_6*-1.0 + 0.083333333333333329*(u1Pu2*u1Pu2)) + xi_6;
            _data_pdfs_20_317_10[_stride_pdfs_0*ctr_0] = omega*(f_eq_common*0.018518518518518517 + u0Pu2*-0.055555555555555552 + xi_23*-1.0 + 0.083333333333333329*(u0Pu2*u0Pu2)) + xi_23;
            _data_pdfs_20_318_10[_stride_pdfs_0*ctr_0] = omega*(f_eq_common*0.018518518518518517 + u0Mu2*0.055555555555555552 + xi_25*-1.0 + 0.083333333333333329*(u0Mu2*u0Mu2)) + xi_25;
            _data_pdfs_20_319_10[_stride_pdfs_0*ctr_0] = omega*(delta_rho*-0.013888888888888888 + f_eq_common*0.018518518518518517 + u0Pu1*0.013888888888888888 + u_2*0.013888888888888888 + xi_12*-1.0 + 0.020833333333333332*(u0Pu1*u0Pu1) + 0.020833333333333332*(u0Pu2*u0Pu2) + 0.020833333333333332*(u1Pu2*u1Pu2)) + xi_12;
            _data_pdfs_20_320_10[_stride_pdfs_0*ctr_0] = omega*(delta_rho*-0.013888888888888888 + f_eq_common*0.018518518518518517 + u0Mu1*-0.013888888888888888 + u_2*0.013888888888888888 + xi_17*-1.0 + 0.020833333333333332*(u0Mu1*u0Mu1) + 0.020833333333333332*(u0Mu2*u0Mu2) + 0.020833333333333332*(u1Pu2*u1Pu2)) + xi_17;
            _data_pdfs_20_321_10[_stride_pdfs_0*ctr_0] = omega*(delta_rho*-0.013888888888888888 + f_eq_common*0.018518518518518517 + u0Mu1*0.013888888888888888 + u_2*0.013888888888888888 + xi_8*-1.0 + 0.020833333333333332*(u0Mu1*u0Mu1) + 0.020833333333333332*(u0Pu2*u0Pu2) + 0.020833333333333332*(u1Mu2*u1Mu2)) + xi_8;
            _data_pdfs_20_322_10[_stride_pdfs_0*ctr_0] = omega*(delta_rho*-0.013888888888888888 + f_eq_common*0.018518518518518517 + u0Pu1*-0.013888888888888888 + u_2*0.013888888888888888 + xi_22*-1.0 + 0.020833333333333332*(u0Mu2*u0Mu2) + 0.020833333333333332*(u0Pu1*u0Pu1) + 0.020833333333333332*(u1Mu2*u1Mu2)) + xi_22;
            _data_pdfs_20_323_10[_stride_pdfs_0*ctr_0] = omega*(delta_rho*-0.013888888888888888 + f_eq_common*0.018518518518518517 + u0Pu1*0.013888888888888888 + u_2*-0.013888888888888888 + xi_4*-1.0 + 0.020833333333333332*(u0Mu2*u0Mu2) + 0.020833333333333332*(u0Pu1*u0Pu1) + 0.020833333333333332*(u1Mu2*u1Mu2)) + xi_4;
            _data_pdfs_20_324_10[_stride_pdfs_0*ctr_0] = omega*(delta_rho*-0.013888888888888888 + f_eq_common*0.018518518518518517 + u0Mu1*-0.013888888888888888 + u_2*-0.013888888888888888 + xi_10*-1.0 + 0.020833333333333332*(u0Mu1*u0Mu1) + 0.020833333333333332*(u0Pu2*u0Pu2) + 0.020833333333333332*(u1Mu2*u1Mu2)) + xi_10;
            _data_pdfs_20_325_10[_stride_pdfs_0*ctr_0] = omega*(delta_rho*-0.013888888888888888 + f_eq_common*0.018518518518518517 + u0Mu1*0.013888888888888888 + u_2*-0.013888888888888888 + xi_15*-1.0 + 0.020833333333333332*(u0Mu1*u0Mu1) + 0.020833333333333332*(u0Mu2*u0Mu2) + 0.020833333333333332*(u1Pu2*u1Pu2)) + xi_15;
            _data_pdfs_20_326_10[_stride_pdfs_0*ctr_0] = omega*(delta_rho*-0.013888888888888888 + f_eq_common*0.018518518518518517 + u0Pu1*-0.013888888888888888 + u_2*-0.013888888888888888 + xi_3*-1.0 + 0.020833333333333332*(u0Pu1*u0Pu1) + 0.020833333333333332*(u0Pu2*u0Pu2) + 0.020833333333333332*(u1Pu2*u1Pu2)) + xi_3;
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
      double * RESTRICT _data_pdfs_20_30 = _data_pdfs + _stride_pdfs_2*ctr_2;
      double * RESTRICT _data_pdfs_20_31 = _data_pdfs + _stride_pdfs_2*ctr_2 + _stride_pdfs_3;
      double * RESTRICT _data_pdfs_20_32 = _data_pdfs + _stride_pdfs_2*ctr_2 + 2*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_20_33 = _data_pdfs + _stride_pdfs_2*ctr_2 + 3*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_20_34 = _data_pdfs + _stride_pdfs_2*ctr_2 + 4*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_2m1_35 = _data_pdfs + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 5*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_21_36 = _data_pdfs + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 6*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_20_37 = _data_pdfs + _stride_pdfs_2*ctr_2 + 7*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_20_38 = _data_pdfs + _stride_pdfs_2*ctr_2 + 8*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_20_39 = _data_pdfs + _stride_pdfs_2*ctr_2 + 9*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_20_310 = _data_pdfs + _stride_pdfs_2*ctr_2 + 10*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_2m1_311 = _data_pdfs + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 11*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_2m1_312 = _data_pdfs + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 12*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_2m1_313 = _data_pdfs + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 13*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_2m1_314 = _data_pdfs + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 14*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_21_315 = _data_pdfs + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 15*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_21_316 = _data_pdfs + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 16*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_21_317 = _data_pdfs + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 17*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_21_318 = _data_pdfs + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 18*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_2m1_319 = _data_pdfs + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 19*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_2m1_320 = _data_pdfs + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 20*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_2m1_321 = _data_pdfs + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 21*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_2m1_322 = _data_pdfs + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 22*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_21_323 = _data_pdfs + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 23*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_21_324 = _data_pdfs + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 24*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_21_325 = _data_pdfs + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 25*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_21_326 = _data_pdfs + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 26*_stride_pdfs_3;
      double * RESTRICT  _data_pdfs_tmp_20_30 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2;
      double * RESTRICT  _data_pdfs_tmp_20_31 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + _stride_pdfs_tmp_3;
      double * RESTRICT  _data_pdfs_tmp_20_32 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 2*_stride_pdfs_tmp_3;
      double * RESTRICT  _data_pdfs_tmp_20_33 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 3*_stride_pdfs_tmp_3;
      double * RESTRICT  _data_pdfs_tmp_20_34 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 4*_stride_pdfs_tmp_3;
      double * RESTRICT  _data_pdfs_tmp_20_35 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 5*_stride_pdfs_tmp_3;
      double * RESTRICT  _data_pdfs_tmp_20_36 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 6*_stride_pdfs_tmp_3;
      double * RESTRICT  _data_pdfs_tmp_20_37 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 7*_stride_pdfs_tmp_3;
      double * RESTRICT  _data_pdfs_tmp_20_38 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 8*_stride_pdfs_tmp_3;
      double * RESTRICT  _data_pdfs_tmp_20_39 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 9*_stride_pdfs_tmp_3;
      double * RESTRICT  _data_pdfs_tmp_20_310 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 10*_stride_pdfs_tmp_3;
      double * RESTRICT  _data_pdfs_tmp_20_311 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 11*_stride_pdfs_tmp_3;
      double * RESTRICT  _data_pdfs_tmp_20_312 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 12*_stride_pdfs_tmp_3;
      double * RESTRICT  _data_pdfs_tmp_20_313 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 13*_stride_pdfs_tmp_3;
      double * RESTRICT  _data_pdfs_tmp_20_314 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 14*_stride_pdfs_tmp_3;
      double * RESTRICT  _data_pdfs_tmp_20_315 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 15*_stride_pdfs_tmp_3;
      double * RESTRICT  _data_pdfs_tmp_20_316 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 16*_stride_pdfs_tmp_3;
      double * RESTRICT  _data_pdfs_tmp_20_317 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 17*_stride_pdfs_tmp_3;
      double * RESTRICT  _data_pdfs_tmp_20_318 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 18*_stride_pdfs_tmp_3;
      double * RESTRICT  _data_pdfs_tmp_20_319 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 19*_stride_pdfs_tmp_3;
      double * RESTRICT  _data_pdfs_tmp_20_320 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 20*_stride_pdfs_tmp_3;
      double * RESTRICT  _data_pdfs_tmp_20_321 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 21*_stride_pdfs_tmp_3;
      double * RESTRICT  _data_pdfs_tmp_20_322 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 22*_stride_pdfs_tmp_3;
      double * RESTRICT  _data_pdfs_tmp_20_323 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 23*_stride_pdfs_tmp_3;
      double * RESTRICT  _data_pdfs_tmp_20_324 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 24*_stride_pdfs_tmp_3;
      double * RESTRICT  _data_pdfs_tmp_20_325 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 25*_stride_pdfs_tmp_3;
      double * RESTRICT  _data_pdfs_tmp_20_326 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 26*_stride_pdfs_tmp_3;
      for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_1; ctr_1 += 1)
      {
         double * RESTRICT _data_pdfs_20_30_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_30;
         double * RESTRICT _data_pdfs_20_31_1m1 = _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _data_pdfs_20_31;
         double * RESTRICT _data_pdfs_20_32_11 = _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _data_pdfs_20_32;
         double * RESTRICT _data_pdfs_20_33_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_33;
         double * RESTRICT _data_pdfs_20_34_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_34;
         double * RESTRICT _data_pdfs_2m1_35_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_2m1_35;
         double * RESTRICT _data_pdfs_21_36_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_21_36;
         double * RESTRICT _data_pdfs_20_37_1m1 = _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _data_pdfs_20_37;
         double * RESTRICT _data_pdfs_20_38_1m1 = _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _data_pdfs_20_38;
         double * RESTRICT _data_pdfs_20_39_11 = _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _data_pdfs_20_39;
         double * RESTRICT _data_pdfs_20_310_11 = _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _data_pdfs_20_310;
         double * RESTRICT _data_pdfs_2m1_311_1m1 = _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _data_pdfs_2m1_311;
         double * RESTRICT _data_pdfs_2m1_312_11 = _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _data_pdfs_2m1_312;
         double * RESTRICT _data_pdfs_2m1_313_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_2m1_313;
         double * RESTRICT _data_pdfs_2m1_314_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_2m1_314;
         double * RESTRICT _data_pdfs_21_315_1m1 = _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _data_pdfs_21_315;
         double * RESTRICT _data_pdfs_21_316_11 = _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _data_pdfs_21_316;
         double * RESTRICT _data_pdfs_21_317_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_21_317;
         double * RESTRICT _data_pdfs_21_318_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_21_318;
         double * RESTRICT _data_pdfs_2m1_319_1m1 = _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _data_pdfs_2m1_319;
         double * RESTRICT _data_pdfs_2m1_320_1m1 = _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _data_pdfs_2m1_320;
         double * RESTRICT _data_pdfs_2m1_321_11 = _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _data_pdfs_2m1_321;
         double * RESTRICT _data_pdfs_2m1_322_11 = _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _data_pdfs_2m1_322;
         double * RESTRICT _data_pdfs_21_323_1m1 = _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _data_pdfs_21_323;
         double * RESTRICT _data_pdfs_21_324_1m1 = _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _data_pdfs_21_324;
         double * RESTRICT _data_pdfs_21_325_11 = _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _data_pdfs_21_325;
         double * RESTRICT _data_pdfs_21_326_11 = _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _data_pdfs_21_326;
         double * RESTRICT  _data_pdfs_tmp_20_30_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_30;
         double * RESTRICT  _data_pdfs_tmp_20_31_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_31;
         double * RESTRICT  _data_pdfs_tmp_20_32_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_32;
         double * RESTRICT  _data_pdfs_tmp_20_33_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_33;
         double * RESTRICT  _data_pdfs_tmp_20_34_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_34;
         double * RESTRICT  _data_pdfs_tmp_20_35_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_35;
         double * RESTRICT  _data_pdfs_tmp_20_36_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_36;
         double * RESTRICT  _data_pdfs_tmp_20_37_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_37;
         double * RESTRICT  _data_pdfs_tmp_20_38_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_38;
         double * RESTRICT  _data_pdfs_tmp_20_39_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_39;
         double * RESTRICT  _data_pdfs_tmp_20_310_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_310;
         double * RESTRICT  _data_pdfs_tmp_20_311_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_311;
         double * RESTRICT  _data_pdfs_tmp_20_312_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_312;
         double * RESTRICT  _data_pdfs_tmp_20_313_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_313;
         double * RESTRICT  _data_pdfs_tmp_20_314_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_314;
         double * RESTRICT  _data_pdfs_tmp_20_315_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_315;
         double * RESTRICT  _data_pdfs_tmp_20_316_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_316;
         double * RESTRICT  _data_pdfs_tmp_20_317_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_317;
         double * RESTRICT  _data_pdfs_tmp_20_318_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_318;
         double * RESTRICT  _data_pdfs_tmp_20_319_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_319;
         double * RESTRICT  _data_pdfs_tmp_20_320_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_320;
         double * RESTRICT  _data_pdfs_tmp_20_321_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_321;
         double * RESTRICT  _data_pdfs_tmp_20_322_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_322;
         double * RESTRICT  _data_pdfs_tmp_20_323_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_323;
         double * RESTRICT  _data_pdfs_tmp_20_324_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_324;
         double * RESTRICT  _data_pdfs_tmp_20_325_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_325;
         double * RESTRICT  _data_pdfs_tmp_20_326_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_326;
         for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_0; ctr_0 += 1)
         {
            const double streamed_0 = _data_pdfs_20_30_10[_stride_pdfs_0*ctr_0];
            const double streamed_1 = _data_pdfs_20_31_1m1[_stride_pdfs_0*ctr_0];
            const double streamed_2 = _data_pdfs_20_32_11[_stride_pdfs_0*ctr_0];
            const double streamed_3 = _data_pdfs_20_33_10[_stride_pdfs_0*ctr_0 + _stride_pdfs_0];
            const double streamed_4 = _data_pdfs_20_34_10[_stride_pdfs_0*ctr_0 - _stride_pdfs_0];
            const double streamed_5 = _data_pdfs_2m1_35_10[_stride_pdfs_0*ctr_0];
            const double streamed_6 = _data_pdfs_21_36_10[_stride_pdfs_0*ctr_0];
            const double streamed_7 = _data_pdfs_20_37_1m1[_stride_pdfs_0*ctr_0 + _stride_pdfs_0];
            const double streamed_8 = _data_pdfs_20_38_1m1[_stride_pdfs_0*ctr_0 - _stride_pdfs_0];
            const double streamed_9 = _data_pdfs_20_39_11[_stride_pdfs_0*ctr_0 + _stride_pdfs_0];
            const double streamed_10 = _data_pdfs_20_310_11[_stride_pdfs_0*ctr_0 - _stride_pdfs_0];
            const double streamed_11 = _data_pdfs_2m1_311_1m1[_stride_pdfs_0*ctr_0];
            const double streamed_12 = _data_pdfs_2m1_312_11[_stride_pdfs_0*ctr_0];
            const double streamed_13 = _data_pdfs_2m1_313_10[_stride_pdfs_0*ctr_0 + _stride_pdfs_0];
            const double streamed_14 = _data_pdfs_2m1_314_10[_stride_pdfs_0*ctr_0 - _stride_pdfs_0];
            const double streamed_15 = _data_pdfs_21_315_1m1[_stride_pdfs_0*ctr_0];
            const double streamed_16 = _data_pdfs_21_316_11[_stride_pdfs_0*ctr_0];
            const double streamed_17 = _data_pdfs_21_317_10[_stride_pdfs_0*ctr_0 + _stride_pdfs_0];
            const double streamed_18 = _data_pdfs_21_318_10[_stride_pdfs_0*ctr_0 - _stride_pdfs_0];
            const double streamed_19 = _data_pdfs_2m1_319_1m1[_stride_pdfs_0*ctr_0 - _stride_pdfs_0];
            const double streamed_20 = _data_pdfs_2m1_320_1m1[_stride_pdfs_0*ctr_0 + _stride_pdfs_0];
            const double streamed_21 = _data_pdfs_2m1_321_11[_stride_pdfs_0*ctr_0 - _stride_pdfs_0];
            const double streamed_22 = _data_pdfs_2m1_322_11[_stride_pdfs_0*ctr_0 + _stride_pdfs_0];
            const double streamed_23 = _data_pdfs_21_323_1m1[_stride_pdfs_0*ctr_0 - _stride_pdfs_0];
            const double streamed_24 = _data_pdfs_21_324_1m1[_stride_pdfs_0*ctr_0 + _stride_pdfs_0];
            const double streamed_25 = _data_pdfs_21_325_11[_stride_pdfs_0*ctr_0 - _stride_pdfs_0];
            const double streamed_26 = _data_pdfs_21_326_11[_stride_pdfs_0*ctr_0 + _stride_pdfs_0];
            _data_pdfs_tmp_20_30_10[_stride_pdfs_tmp_0*ctr_0] = streamed_0;
            _data_pdfs_tmp_20_31_10[_stride_pdfs_tmp_0*ctr_0] = streamed_1;
            _data_pdfs_tmp_20_32_10[_stride_pdfs_tmp_0*ctr_0] = streamed_2;
            _data_pdfs_tmp_20_33_10[_stride_pdfs_tmp_0*ctr_0] = streamed_3;
            _data_pdfs_tmp_20_34_10[_stride_pdfs_tmp_0*ctr_0] = streamed_4;
            _data_pdfs_tmp_20_35_10[_stride_pdfs_tmp_0*ctr_0] = streamed_5;
            _data_pdfs_tmp_20_36_10[_stride_pdfs_tmp_0*ctr_0] = streamed_6;
            _data_pdfs_tmp_20_37_10[_stride_pdfs_tmp_0*ctr_0] = streamed_7;
            _data_pdfs_tmp_20_38_10[_stride_pdfs_tmp_0*ctr_0] = streamed_8;
            _data_pdfs_tmp_20_39_10[_stride_pdfs_tmp_0*ctr_0] = streamed_9;
            _data_pdfs_tmp_20_310_10[_stride_pdfs_tmp_0*ctr_0] = streamed_10;
            _data_pdfs_tmp_20_311_10[_stride_pdfs_tmp_0*ctr_0] = streamed_11;
            _data_pdfs_tmp_20_312_10[_stride_pdfs_tmp_0*ctr_0] = streamed_12;
            _data_pdfs_tmp_20_313_10[_stride_pdfs_tmp_0*ctr_0] = streamed_13;
            _data_pdfs_tmp_20_314_10[_stride_pdfs_tmp_0*ctr_0] = streamed_14;
            _data_pdfs_tmp_20_315_10[_stride_pdfs_tmp_0*ctr_0] = streamed_15;
            _data_pdfs_tmp_20_316_10[_stride_pdfs_tmp_0*ctr_0] = streamed_16;
            _data_pdfs_tmp_20_317_10[_stride_pdfs_tmp_0*ctr_0] = streamed_17;
            _data_pdfs_tmp_20_318_10[_stride_pdfs_tmp_0*ctr_0] = streamed_18;
            _data_pdfs_tmp_20_319_10[_stride_pdfs_tmp_0*ctr_0] = streamed_19;
            _data_pdfs_tmp_20_320_10[_stride_pdfs_tmp_0*ctr_0] = streamed_20;
            _data_pdfs_tmp_20_321_10[_stride_pdfs_tmp_0*ctr_0] = streamed_21;
            _data_pdfs_tmp_20_322_10[_stride_pdfs_tmp_0*ctr_0] = streamed_22;
            _data_pdfs_tmp_20_323_10[_stride_pdfs_tmp_0*ctr_0] = streamed_23;
            _data_pdfs_tmp_20_324_10[_stride_pdfs_tmp_0*ctr_0] = streamed_24;
            _data_pdfs_tmp_20_325_10[_stride_pdfs_tmp_0*ctr_0] = streamed_25;
            _data_pdfs_tmp_20_326_10[_stride_pdfs_tmp_0*ctr_0] = streamed_26;
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
      double * RESTRICT _data_pdfs_20_30 = _data_pdfs + _stride_pdfs_2*ctr_2;
      double * RESTRICT _data_pdfs_20_31 = _data_pdfs + _stride_pdfs_2*ctr_2 + _stride_pdfs_3;
      double * RESTRICT _data_pdfs_20_32 = _data_pdfs + _stride_pdfs_2*ctr_2 + 2*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_20_33 = _data_pdfs + _stride_pdfs_2*ctr_2 + 3*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_20_34 = _data_pdfs + _stride_pdfs_2*ctr_2 + 4*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_2m1_35 = _data_pdfs + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 5*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_21_36 = _data_pdfs + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 6*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_20_37 = _data_pdfs + _stride_pdfs_2*ctr_2 + 7*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_20_38 = _data_pdfs + _stride_pdfs_2*ctr_2 + 8*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_20_39 = _data_pdfs + _stride_pdfs_2*ctr_2 + 9*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_20_310 = _data_pdfs + _stride_pdfs_2*ctr_2 + 10*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_2m1_311 = _data_pdfs + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 11*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_2m1_312 = _data_pdfs + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 12*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_2m1_313 = _data_pdfs + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 13*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_2m1_314 = _data_pdfs + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 14*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_21_315 = _data_pdfs + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 15*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_21_316 = _data_pdfs + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 16*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_21_317 = _data_pdfs + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 17*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_21_318 = _data_pdfs + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 18*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_2m1_319 = _data_pdfs + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 19*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_2m1_320 = _data_pdfs + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 20*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_2m1_321 = _data_pdfs + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 21*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_2m1_322 = _data_pdfs + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 22*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_21_323 = _data_pdfs + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 23*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_21_324 = _data_pdfs + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 24*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_21_325 = _data_pdfs + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 25*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_21_326 = _data_pdfs + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 26*_stride_pdfs_3;
      double * RESTRICT  _data_pdfs_tmp_20_30 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2;
      double * RESTRICT  _data_pdfs_tmp_20_31 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + _stride_pdfs_tmp_3;
      double * RESTRICT  _data_pdfs_tmp_20_32 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 2*_stride_pdfs_tmp_3;
      double * RESTRICT  _data_pdfs_tmp_20_33 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 3*_stride_pdfs_tmp_3;
      double * RESTRICT  _data_pdfs_tmp_20_34 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 4*_stride_pdfs_tmp_3;
      double * RESTRICT  _data_pdfs_tmp_20_35 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 5*_stride_pdfs_tmp_3;
      double * RESTRICT  _data_pdfs_tmp_20_36 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 6*_stride_pdfs_tmp_3;
      double * RESTRICT  _data_pdfs_tmp_20_37 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 7*_stride_pdfs_tmp_3;
      double * RESTRICT  _data_pdfs_tmp_20_38 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 8*_stride_pdfs_tmp_3;
      double * RESTRICT  _data_pdfs_tmp_20_39 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 9*_stride_pdfs_tmp_3;
      double * RESTRICT  _data_pdfs_tmp_20_310 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 10*_stride_pdfs_tmp_3;
      double * RESTRICT  _data_pdfs_tmp_20_311 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 11*_stride_pdfs_tmp_3;
      double * RESTRICT  _data_pdfs_tmp_20_312 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 12*_stride_pdfs_tmp_3;
      double * RESTRICT  _data_pdfs_tmp_20_313 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 13*_stride_pdfs_tmp_3;
      double * RESTRICT  _data_pdfs_tmp_20_314 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 14*_stride_pdfs_tmp_3;
      double * RESTRICT  _data_pdfs_tmp_20_315 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 15*_stride_pdfs_tmp_3;
      double * RESTRICT  _data_pdfs_tmp_20_316 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 16*_stride_pdfs_tmp_3;
      double * RESTRICT  _data_pdfs_tmp_20_317 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 17*_stride_pdfs_tmp_3;
      double * RESTRICT  _data_pdfs_tmp_20_318 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 18*_stride_pdfs_tmp_3;
      double * RESTRICT  _data_pdfs_tmp_20_319 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 19*_stride_pdfs_tmp_3;
      double * RESTRICT  _data_pdfs_tmp_20_320 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 20*_stride_pdfs_tmp_3;
      double * RESTRICT  _data_pdfs_tmp_20_321 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 21*_stride_pdfs_tmp_3;
      double * RESTRICT  _data_pdfs_tmp_20_322 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 22*_stride_pdfs_tmp_3;
      double * RESTRICT  _data_pdfs_tmp_20_323 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 23*_stride_pdfs_tmp_3;
      double * RESTRICT  _data_pdfs_tmp_20_324 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 24*_stride_pdfs_tmp_3;
      double * RESTRICT  _data_pdfs_tmp_20_325 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 25*_stride_pdfs_tmp_3;
      double * RESTRICT  _data_pdfs_tmp_20_326 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 26*_stride_pdfs_tmp_3;
      for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_1; ctr_1 += 1)
      {
         double * RESTRICT _data_pdfs_20_30_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_30;
         double * RESTRICT _data_pdfs_20_31_1m1 = _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _data_pdfs_20_31;
         double * RESTRICT _data_pdfs_20_32_11 = _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _data_pdfs_20_32;
         double * RESTRICT _data_pdfs_20_33_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_33;
         double * RESTRICT _data_pdfs_20_34_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_34;
         double * RESTRICT _data_pdfs_2m1_35_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_2m1_35;
         double * RESTRICT _data_pdfs_21_36_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_21_36;
         double * RESTRICT _data_pdfs_20_37_1m1 = _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _data_pdfs_20_37;
         double * RESTRICT _data_pdfs_20_38_1m1 = _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _data_pdfs_20_38;
         double * RESTRICT _data_pdfs_20_39_11 = _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _data_pdfs_20_39;
         double * RESTRICT _data_pdfs_20_310_11 = _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _data_pdfs_20_310;
         double * RESTRICT _data_pdfs_2m1_311_1m1 = _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _data_pdfs_2m1_311;
         double * RESTRICT _data_pdfs_2m1_312_11 = _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _data_pdfs_2m1_312;
         double * RESTRICT _data_pdfs_2m1_313_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_2m1_313;
         double * RESTRICT _data_pdfs_2m1_314_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_2m1_314;
         double * RESTRICT _data_pdfs_21_315_1m1 = _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _data_pdfs_21_315;
         double * RESTRICT _data_pdfs_21_316_11 = _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _data_pdfs_21_316;
         double * RESTRICT _data_pdfs_21_317_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_21_317;
         double * RESTRICT _data_pdfs_21_318_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_21_318;
         double * RESTRICT _data_pdfs_2m1_319_1m1 = _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _data_pdfs_2m1_319;
         double * RESTRICT _data_pdfs_2m1_320_1m1 = _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _data_pdfs_2m1_320;
         double * RESTRICT _data_pdfs_2m1_321_11 = _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _data_pdfs_2m1_321;
         double * RESTRICT _data_pdfs_2m1_322_11 = _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _data_pdfs_2m1_322;
         double * RESTRICT _data_pdfs_21_323_1m1 = _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _data_pdfs_21_323;
         double * RESTRICT _data_pdfs_21_324_1m1 = _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _data_pdfs_21_324;
         double * RESTRICT _data_pdfs_21_325_11 = _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _data_pdfs_21_325;
         double * RESTRICT _data_pdfs_21_326_11 = _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _data_pdfs_21_326;
         double * RESTRICT  _data_pdfs_tmp_20_30_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_30;
         double * RESTRICT  _data_pdfs_tmp_20_31_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_31;
         double * RESTRICT  _data_pdfs_tmp_20_32_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_32;
         double * RESTRICT  _data_pdfs_tmp_20_33_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_33;
         double * RESTRICT  _data_pdfs_tmp_20_34_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_34;
         double * RESTRICT  _data_pdfs_tmp_20_35_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_35;
         double * RESTRICT  _data_pdfs_tmp_20_36_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_36;
         double * RESTRICT  _data_pdfs_tmp_20_37_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_37;
         double * RESTRICT  _data_pdfs_tmp_20_38_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_38;
         double * RESTRICT  _data_pdfs_tmp_20_39_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_39;
         double * RESTRICT  _data_pdfs_tmp_20_310_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_310;
         double * RESTRICT  _data_pdfs_tmp_20_311_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_311;
         double * RESTRICT  _data_pdfs_tmp_20_312_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_312;
         double * RESTRICT  _data_pdfs_tmp_20_313_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_313;
         double * RESTRICT  _data_pdfs_tmp_20_314_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_314;
         double * RESTRICT  _data_pdfs_tmp_20_315_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_315;
         double * RESTRICT  _data_pdfs_tmp_20_316_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_316;
         double * RESTRICT  _data_pdfs_tmp_20_317_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_317;
         double * RESTRICT  _data_pdfs_tmp_20_318_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_318;
         double * RESTRICT  _data_pdfs_tmp_20_319_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_319;
         double * RESTRICT  _data_pdfs_tmp_20_320_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_320;
         double * RESTRICT  _data_pdfs_tmp_20_321_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_321;
         double * RESTRICT  _data_pdfs_tmp_20_322_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_322;
         double * RESTRICT  _data_pdfs_tmp_20_323_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_323;
         double * RESTRICT  _data_pdfs_tmp_20_324_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_324;
         double * RESTRICT  _data_pdfs_tmp_20_325_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_325;
         double * RESTRICT  _data_pdfs_tmp_20_326_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_326;
         for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_0; ctr_0 += 1)
         {
            const double streamed_0 = _data_pdfs_20_30_10[_stride_pdfs_0*ctr_0];
            const double streamed_1 = _data_pdfs_20_31_1m1[_stride_pdfs_0*ctr_0];
            const double streamed_2 = _data_pdfs_20_32_11[_stride_pdfs_0*ctr_0];
            const double streamed_3 = _data_pdfs_20_33_10[_stride_pdfs_0*ctr_0 + _stride_pdfs_0];
            const double streamed_4 = _data_pdfs_20_34_10[_stride_pdfs_0*ctr_0 - _stride_pdfs_0];
            const double streamed_5 = _data_pdfs_2m1_35_10[_stride_pdfs_0*ctr_0];
            const double streamed_6 = _data_pdfs_21_36_10[_stride_pdfs_0*ctr_0];
            const double streamed_7 = _data_pdfs_20_37_1m1[_stride_pdfs_0*ctr_0 + _stride_pdfs_0];
            const double streamed_8 = _data_pdfs_20_38_1m1[_stride_pdfs_0*ctr_0 - _stride_pdfs_0];
            const double streamed_9 = _data_pdfs_20_39_11[_stride_pdfs_0*ctr_0 + _stride_pdfs_0];
            const double streamed_10 = _data_pdfs_20_310_11[_stride_pdfs_0*ctr_0 - _stride_pdfs_0];
            const double streamed_11 = _data_pdfs_2m1_311_1m1[_stride_pdfs_0*ctr_0];
            const double streamed_12 = _data_pdfs_2m1_312_11[_stride_pdfs_0*ctr_0];
            const double streamed_13 = _data_pdfs_2m1_313_10[_stride_pdfs_0*ctr_0 + _stride_pdfs_0];
            const double streamed_14 = _data_pdfs_2m1_314_10[_stride_pdfs_0*ctr_0 - _stride_pdfs_0];
            const double streamed_15 = _data_pdfs_21_315_1m1[_stride_pdfs_0*ctr_0];
            const double streamed_16 = _data_pdfs_21_316_11[_stride_pdfs_0*ctr_0];
            const double streamed_17 = _data_pdfs_21_317_10[_stride_pdfs_0*ctr_0 + _stride_pdfs_0];
            const double streamed_18 = _data_pdfs_21_318_10[_stride_pdfs_0*ctr_0 - _stride_pdfs_0];
            const double streamed_19 = _data_pdfs_2m1_319_1m1[_stride_pdfs_0*ctr_0 - _stride_pdfs_0];
            const double streamed_20 = _data_pdfs_2m1_320_1m1[_stride_pdfs_0*ctr_0 + _stride_pdfs_0];
            const double streamed_21 = _data_pdfs_2m1_321_11[_stride_pdfs_0*ctr_0 - _stride_pdfs_0];
            const double streamed_22 = _data_pdfs_2m1_322_11[_stride_pdfs_0*ctr_0 + _stride_pdfs_0];
            const double streamed_23 = _data_pdfs_21_323_1m1[_stride_pdfs_0*ctr_0 - _stride_pdfs_0];
            const double streamed_24 = _data_pdfs_21_324_1m1[_stride_pdfs_0*ctr_0 + _stride_pdfs_0];
            const double streamed_25 = _data_pdfs_21_325_11[_stride_pdfs_0*ctr_0 - _stride_pdfs_0];
            const double streamed_26 = _data_pdfs_21_326_11[_stride_pdfs_0*ctr_0 + _stride_pdfs_0];
            _data_pdfs_tmp_20_30_10[_stride_pdfs_tmp_0*ctr_0] = streamed_0;
            _data_pdfs_tmp_20_31_10[_stride_pdfs_tmp_0*ctr_0] = streamed_1;
            _data_pdfs_tmp_20_32_10[_stride_pdfs_tmp_0*ctr_0] = streamed_2;
            _data_pdfs_tmp_20_33_10[_stride_pdfs_tmp_0*ctr_0] = streamed_3;
            _data_pdfs_tmp_20_34_10[_stride_pdfs_tmp_0*ctr_0] = streamed_4;
            _data_pdfs_tmp_20_35_10[_stride_pdfs_tmp_0*ctr_0] = streamed_5;
            _data_pdfs_tmp_20_36_10[_stride_pdfs_tmp_0*ctr_0] = streamed_6;
            _data_pdfs_tmp_20_37_10[_stride_pdfs_tmp_0*ctr_0] = streamed_7;
            _data_pdfs_tmp_20_38_10[_stride_pdfs_tmp_0*ctr_0] = streamed_8;
            _data_pdfs_tmp_20_39_10[_stride_pdfs_tmp_0*ctr_0] = streamed_9;
            _data_pdfs_tmp_20_310_10[_stride_pdfs_tmp_0*ctr_0] = streamed_10;
            _data_pdfs_tmp_20_311_10[_stride_pdfs_tmp_0*ctr_0] = streamed_11;
            _data_pdfs_tmp_20_312_10[_stride_pdfs_tmp_0*ctr_0] = streamed_12;
            _data_pdfs_tmp_20_313_10[_stride_pdfs_tmp_0*ctr_0] = streamed_13;
            _data_pdfs_tmp_20_314_10[_stride_pdfs_tmp_0*ctr_0] = streamed_14;
            _data_pdfs_tmp_20_315_10[_stride_pdfs_tmp_0*ctr_0] = streamed_15;
            _data_pdfs_tmp_20_316_10[_stride_pdfs_tmp_0*ctr_0] = streamed_16;
            _data_pdfs_tmp_20_317_10[_stride_pdfs_tmp_0*ctr_0] = streamed_17;
            _data_pdfs_tmp_20_318_10[_stride_pdfs_tmp_0*ctr_0] = streamed_18;
            _data_pdfs_tmp_20_319_10[_stride_pdfs_tmp_0*ctr_0] = streamed_19;
            _data_pdfs_tmp_20_320_10[_stride_pdfs_tmp_0*ctr_0] = streamed_20;
            _data_pdfs_tmp_20_321_10[_stride_pdfs_tmp_0*ctr_0] = streamed_21;
            _data_pdfs_tmp_20_322_10[_stride_pdfs_tmp_0*ctr_0] = streamed_22;
            _data_pdfs_tmp_20_323_10[_stride_pdfs_tmp_0*ctr_0] = streamed_23;
            _data_pdfs_tmp_20_324_10[_stride_pdfs_tmp_0*ctr_0] = streamed_24;
            _data_pdfs_tmp_20_325_10[_stride_pdfs_tmp_0*ctr_0] = streamed_25;
            _data_pdfs_tmp_20_326_10[_stride_pdfs_tmp_0*ctr_0] = streamed_26;
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
      double * RESTRICT _data_density_20_30 = _data_density + _stride_density_2*ctr_2;
      double * RESTRICT _data_velocity_20_30 = _data_velocity + _stride_velocity_2*ctr_2;
      double * RESTRICT _data_velocity_20_31 = _data_velocity + _stride_velocity_2*ctr_2 + _stride_velocity_3;
      double * RESTRICT _data_velocity_20_32 = _data_velocity + _stride_velocity_2*ctr_2 + 2*_stride_velocity_3;
      double * RESTRICT  _data_pdfs_20_30 = _data_pdfs + _stride_pdfs_2*ctr_2;
      double * RESTRICT  _data_pdfs_20_31 = _data_pdfs + _stride_pdfs_2*ctr_2 + _stride_pdfs_3;
      double * RESTRICT  _data_pdfs_20_32 = _data_pdfs + _stride_pdfs_2*ctr_2 + 2*_stride_pdfs_3;
      double * RESTRICT  _data_pdfs_20_33 = _data_pdfs + _stride_pdfs_2*ctr_2 + 3*_stride_pdfs_3;
      double * RESTRICT  _data_pdfs_20_34 = _data_pdfs + _stride_pdfs_2*ctr_2 + 4*_stride_pdfs_3;
      double * RESTRICT  _data_pdfs_20_35 = _data_pdfs + _stride_pdfs_2*ctr_2 + 5*_stride_pdfs_3;
      double * RESTRICT  _data_pdfs_20_36 = _data_pdfs + _stride_pdfs_2*ctr_2 + 6*_stride_pdfs_3;
      double * RESTRICT  _data_pdfs_20_37 = _data_pdfs + _stride_pdfs_2*ctr_2 + 7*_stride_pdfs_3;
      double * RESTRICT  _data_pdfs_20_38 = _data_pdfs + _stride_pdfs_2*ctr_2 + 8*_stride_pdfs_3;
      double * RESTRICT  _data_pdfs_20_39 = _data_pdfs + _stride_pdfs_2*ctr_2 + 9*_stride_pdfs_3;
      double * RESTRICT  _data_pdfs_20_310 = _data_pdfs + _stride_pdfs_2*ctr_2 + 10*_stride_pdfs_3;
      double * RESTRICT  _data_pdfs_20_311 = _data_pdfs + _stride_pdfs_2*ctr_2 + 11*_stride_pdfs_3;
      double * RESTRICT  _data_pdfs_20_312 = _data_pdfs + _stride_pdfs_2*ctr_2 + 12*_stride_pdfs_3;
      double * RESTRICT  _data_pdfs_20_313 = _data_pdfs + _stride_pdfs_2*ctr_2 + 13*_stride_pdfs_3;
      double * RESTRICT  _data_pdfs_20_314 = _data_pdfs + _stride_pdfs_2*ctr_2 + 14*_stride_pdfs_3;
      double * RESTRICT  _data_pdfs_20_315 = _data_pdfs + _stride_pdfs_2*ctr_2 + 15*_stride_pdfs_3;
      double * RESTRICT  _data_pdfs_20_316 = _data_pdfs + _stride_pdfs_2*ctr_2 + 16*_stride_pdfs_3;
      double * RESTRICT  _data_pdfs_20_317 = _data_pdfs + _stride_pdfs_2*ctr_2 + 17*_stride_pdfs_3;
      double * RESTRICT  _data_pdfs_20_318 = _data_pdfs + _stride_pdfs_2*ctr_2 + 18*_stride_pdfs_3;
      double * RESTRICT  _data_pdfs_20_319 = _data_pdfs + _stride_pdfs_2*ctr_2 + 19*_stride_pdfs_3;
      double * RESTRICT  _data_pdfs_20_320 = _data_pdfs + _stride_pdfs_2*ctr_2 + 20*_stride_pdfs_3;
      double * RESTRICT  _data_pdfs_20_321 = _data_pdfs + _stride_pdfs_2*ctr_2 + 21*_stride_pdfs_3;
      double * RESTRICT  _data_pdfs_20_322 = _data_pdfs + _stride_pdfs_2*ctr_2 + 22*_stride_pdfs_3;
      double * RESTRICT  _data_pdfs_20_323 = _data_pdfs + _stride_pdfs_2*ctr_2 + 23*_stride_pdfs_3;
      double * RESTRICT  _data_pdfs_20_324 = _data_pdfs + _stride_pdfs_2*ctr_2 + 24*_stride_pdfs_3;
      double * RESTRICT  _data_pdfs_20_325 = _data_pdfs + _stride_pdfs_2*ctr_2 + 25*_stride_pdfs_3;
      double * RESTRICT  _data_pdfs_20_326 = _data_pdfs + _stride_pdfs_2*ctr_2 + 26*_stride_pdfs_3;
      for (int64_t ctr_1 = 0; ctr_1 < _size_density_1; ctr_1 += 1)
      {
         double * RESTRICT _data_density_20_30_10 = _stride_density_1*ctr_1 + _data_density_20_30;
         double * RESTRICT _data_velocity_20_30_10 = _stride_velocity_1*ctr_1 + _data_velocity_20_30;
         double * RESTRICT _data_velocity_20_31_10 = _stride_velocity_1*ctr_1 + _data_velocity_20_31;
         double * RESTRICT _data_velocity_20_32_10 = _stride_velocity_1*ctr_1 + _data_velocity_20_32;
         double * RESTRICT  _data_pdfs_20_30_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_30;
         double * RESTRICT  _data_pdfs_20_31_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_31;
         double * RESTRICT  _data_pdfs_20_32_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_32;
         double * RESTRICT  _data_pdfs_20_33_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_33;
         double * RESTRICT  _data_pdfs_20_34_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_34;
         double * RESTRICT  _data_pdfs_20_35_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_35;
         double * RESTRICT  _data_pdfs_20_36_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_36;
         double * RESTRICT  _data_pdfs_20_37_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_37;
         double * RESTRICT  _data_pdfs_20_38_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_38;
         double * RESTRICT  _data_pdfs_20_39_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_39;
         double * RESTRICT  _data_pdfs_20_310_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_310;
         double * RESTRICT  _data_pdfs_20_311_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_311;
         double * RESTRICT  _data_pdfs_20_312_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_312;
         double * RESTRICT  _data_pdfs_20_313_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_313;
         double * RESTRICT  _data_pdfs_20_314_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_314;
         double * RESTRICT  _data_pdfs_20_315_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_315;
         double * RESTRICT  _data_pdfs_20_316_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_316;
         double * RESTRICT  _data_pdfs_20_317_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_317;
         double * RESTRICT  _data_pdfs_20_318_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_318;
         double * RESTRICT  _data_pdfs_20_319_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_319;
         double * RESTRICT  _data_pdfs_20_320_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_320;
         double * RESTRICT  _data_pdfs_20_321_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_321;
         double * RESTRICT  _data_pdfs_20_322_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_322;
         double * RESTRICT  _data_pdfs_20_323_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_323;
         double * RESTRICT  _data_pdfs_20_324_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_324;
         double * RESTRICT  _data_pdfs_20_325_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_325;
         double * RESTRICT  _data_pdfs_20_326_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_326;
         for (int64_t ctr_0 = 0; ctr_0 < _size_density_0; ctr_0 += 1)
         {
            const double rho = _data_density_20_30_10[_stride_density_0*ctr_0];
            const double delta_rho = rho - 1.0;
            const double u_0 = _data_velocity_20_30_10[_stride_velocity_0*ctr_0];
            const double u_1 = _data_velocity_20_31_10[_stride_velocity_0*ctr_0];
            const double u_2 = _data_velocity_20_32_10[_stride_velocity_0*ctr_0];
            _data_pdfs_20_30_10[_stride_pdfs_0*ctr_0] = delta_rho*0.29629629629629628 - 0.44444444444444442*(u_0*u_0) - 0.44444444444444442*(u_1*u_1) - 0.44444444444444442*(u_2*u_2);
            _data_pdfs_20_31_10[_stride_pdfs_0*ctr_0] = delta_rho*0.07407407407407407 + u_1*0.22222222222222221 - 0.1111111111111111*(u_0*u_0) - 0.1111111111111111*(u_2*u_2) + 0.22222222222222221*(u_1*u_1);
            _data_pdfs_20_32_10[_stride_pdfs_0*ctr_0] = delta_rho*0.07407407407407407 + u_1*-0.22222222222222221 - 0.1111111111111111*(u_0*u_0) - 0.1111111111111111*(u_2*u_2) + 0.22222222222222221*(u_1*u_1);
            _data_pdfs_20_33_10[_stride_pdfs_0*ctr_0] = delta_rho*0.07407407407407407 + u_0*-0.22222222222222221 - 0.1111111111111111*(u_1*u_1) - 0.1111111111111111*(u_2*u_2) + 0.22222222222222221*(u_0*u_0);
            _data_pdfs_20_34_10[_stride_pdfs_0*ctr_0] = delta_rho*0.07407407407407407 + u_0*0.22222222222222221 - 0.1111111111111111*(u_1*u_1) - 0.1111111111111111*(u_2*u_2) + 0.22222222222222221*(u_0*u_0);
            _data_pdfs_20_35_10[_stride_pdfs_0*ctr_0] = delta_rho*0.07407407407407407 + u_2*0.22222222222222221 - 0.1111111111111111*(u_0*u_0) - 0.1111111111111111*(u_1*u_1) + 0.22222222222222221*(u_2*u_2);
            _data_pdfs_20_36_10[_stride_pdfs_0*ctr_0] = delta_rho*0.07407407407407407 + u_2*-0.22222222222222221 - 0.1111111111111111*(u_0*u_0) - 0.1111111111111111*(u_1*u_1) + 0.22222222222222221*(u_2*u_2);
            _data_pdfs_20_37_10[_stride_pdfs_0*ctr_0] = delta_rho*0.018518518518518517 + u_0*u_1*-0.16666666666666666 + u_0*-0.055555555555555552 + u_1*0.055555555555555552 - 0.027777777777777776*(u_2*u_2) + 0.055555555555555552*(u_0*u_0) + 0.055555555555555552*(u_1*u_1);
            _data_pdfs_20_38_10[_stride_pdfs_0*ctr_0] = delta_rho*0.018518518518518517 + u_0*u_1*0.16666666666666666 + u_0*0.055555555555555552 + u_1*0.055555555555555552 - 0.027777777777777776*(u_2*u_2) + 0.055555555555555552*(u_0*u_0) + 0.055555555555555552*(u_1*u_1);
            _data_pdfs_20_39_10[_stride_pdfs_0*ctr_0] = delta_rho*0.018518518518518517 + u_0*u_1*0.16666666666666666 + u_0*-0.055555555555555552 + u_1*-0.055555555555555552 - 0.027777777777777776*(u_2*u_2) + 0.055555555555555552*(u_0*u_0) + 0.055555555555555552*(u_1*u_1);
            _data_pdfs_20_310_10[_stride_pdfs_0*ctr_0] = delta_rho*0.018518518518518517 + u_0*u_1*-0.16666666666666666 + u_0*0.055555555555555552 + u_1*-0.055555555555555552 - 0.027777777777777776*(u_2*u_2) + 0.055555555555555552*(u_0*u_0) + 0.055555555555555552*(u_1*u_1);
            _data_pdfs_20_311_10[_stride_pdfs_0*ctr_0] = delta_rho*0.018518518518518517 + u_1*u_2*0.16666666666666666 + u_1*0.055555555555555552 + u_2*0.055555555555555552 - 0.027777777777777776*(u_0*u_0) + 0.055555555555555552*(u_1*u_1) + 0.055555555555555552*(u_2*u_2);
            _data_pdfs_20_312_10[_stride_pdfs_0*ctr_0] = delta_rho*0.018518518518518517 + u_1*u_2*-0.16666666666666666 + u_1*-0.055555555555555552 + u_2*0.055555555555555552 - 0.027777777777777776*(u_0*u_0) + 0.055555555555555552*(u_1*u_1) + 0.055555555555555552*(u_2*u_2);
            _data_pdfs_20_313_10[_stride_pdfs_0*ctr_0] = delta_rho*0.018518518518518517 + u_0*u_2*-0.16666666666666666 + u_0*-0.055555555555555552 + u_2*0.055555555555555552 - 0.027777777777777776*(u_1*u_1) + 0.055555555555555552*(u_0*u_0) + 0.055555555555555552*(u_2*u_2);
            _data_pdfs_20_314_10[_stride_pdfs_0*ctr_0] = delta_rho*0.018518518518518517 + u_0*u_2*0.16666666666666666 + u_0*0.055555555555555552 + u_2*0.055555555555555552 - 0.027777777777777776*(u_1*u_1) + 0.055555555555555552*(u_0*u_0) + 0.055555555555555552*(u_2*u_2);
            _data_pdfs_20_315_10[_stride_pdfs_0*ctr_0] = delta_rho*0.018518518518518517 + u_1*u_2*-0.16666666666666666 + u_1*0.055555555555555552 + u_2*-0.055555555555555552 - 0.027777777777777776*(u_0*u_0) + 0.055555555555555552*(u_1*u_1) + 0.055555555555555552*(u_2*u_2);
            _data_pdfs_20_316_10[_stride_pdfs_0*ctr_0] = delta_rho*0.018518518518518517 + u_1*u_2*0.16666666666666666 + u_1*-0.055555555555555552 + u_2*-0.055555555555555552 - 0.027777777777777776*(u_0*u_0) + 0.055555555555555552*(u_1*u_1) + 0.055555555555555552*(u_2*u_2);
            _data_pdfs_20_317_10[_stride_pdfs_0*ctr_0] = delta_rho*0.018518518518518517 + u_0*u_2*0.16666666666666666 + u_0*-0.055555555555555552 + u_2*-0.055555555555555552 - 0.027777777777777776*(u_1*u_1) + 0.055555555555555552*(u_0*u_0) + 0.055555555555555552*(u_2*u_2);
            _data_pdfs_20_318_10[_stride_pdfs_0*ctr_0] = delta_rho*0.018518518518518517 + u_0*u_2*-0.16666666666666666 + u_0*0.055555555555555552 + u_2*-0.055555555555555552 - 0.027777777777777776*(u_1*u_1) + 0.055555555555555552*(u_0*u_0) + 0.055555555555555552*(u_2*u_2);
            _data_pdfs_20_319_10[_stride_pdfs_0*ctr_0] = delta_rho*0.0046296296296296294 + u_0*u_1*0.041666666666666664 + u_0*u_2*0.041666666666666664 + u_0*0.013888888888888888 + u_1*u_2*0.041666666666666664 + u_1*0.013888888888888888 + u_2*0.013888888888888888 + 0.013888888888888888*(u_0*u_0) + 0.013888888888888888*(u_1*u_1) + 0.013888888888888888*(u_2*u_2);
            _data_pdfs_20_320_10[_stride_pdfs_0*ctr_0] = delta_rho*0.0046296296296296294 + u_0*u_1*-0.041666666666666664 + u_0*u_2*-0.041666666666666664 + u_0*-0.013888888888888888 + u_1*u_2*0.041666666666666664 + u_1*0.013888888888888888 + u_2*0.013888888888888888 + 0.013888888888888888*(u_0*u_0) + 0.013888888888888888*(u_1*u_1) + 0.013888888888888888*(u_2*u_2);
            _data_pdfs_20_321_10[_stride_pdfs_0*ctr_0] = delta_rho*0.0046296296296296294 + u_0*u_1*-0.041666666666666664 + u_0*u_2*0.041666666666666664 + u_0*0.013888888888888888 + u_1*u_2*-0.041666666666666664 + u_1*-0.013888888888888888 + u_2*0.013888888888888888 + 0.013888888888888888*(u_0*u_0) + 0.013888888888888888*(u_1*u_1) + 0.013888888888888888*(u_2*u_2);
            _data_pdfs_20_322_10[_stride_pdfs_0*ctr_0] = delta_rho*0.0046296296296296294 + u_0*u_1*0.041666666666666664 + u_0*u_2*-0.041666666666666664 + u_0*-0.013888888888888888 + u_1*u_2*-0.041666666666666664 + u_1*-0.013888888888888888 + u_2*0.013888888888888888 + 0.013888888888888888*(u_0*u_0) + 0.013888888888888888*(u_1*u_1) + 0.013888888888888888*(u_2*u_2);
            _data_pdfs_20_323_10[_stride_pdfs_0*ctr_0] = delta_rho*0.0046296296296296294 + u_0*u_1*0.041666666666666664 + u_0*u_2*-0.041666666666666664 + u_0*0.013888888888888888 + u_1*u_2*-0.041666666666666664 + u_1*0.013888888888888888 + u_2*-0.013888888888888888 + 0.013888888888888888*(u_0*u_0) + 0.013888888888888888*(u_1*u_1) + 0.013888888888888888*(u_2*u_2);
            _data_pdfs_20_324_10[_stride_pdfs_0*ctr_0] = delta_rho*0.0046296296296296294 + u_0*u_1*-0.041666666666666664 + u_0*u_2*0.041666666666666664 + u_0*-0.013888888888888888 + u_1*u_2*-0.041666666666666664 + u_1*0.013888888888888888 + u_2*-0.013888888888888888 + 0.013888888888888888*(u_0*u_0) + 0.013888888888888888*(u_1*u_1) + 0.013888888888888888*(u_2*u_2);
            _data_pdfs_20_325_10[_stride_pdfs_0*ctr_0] = delta_rho*0.0046296296296296294 + u_0*u_1*-0.041666666666666664 + u_0*u_2*-0.041666666666666664 + u_0*0.013888888888888888 + u_1*u_2*0.041666666666666664 + u_1*-0.013888888888888888 + u_2*-0.013888888888888888 + 0.013888888888888888*(u_0*u_0) + 0.013888888888888888*(u_1*u_1) + 0.013888888888888888*(u_2*u_2);
            _data_pdfs_20_326_10[_stride_pdfs_0*ctr_0] = delta_rho*0.0046296296296296294 + u_0*u_1*0.041666666666666664 + u_0*u_2*0.041666666666666664 + u_0*-0.013888888888888888 + u_1*u_2*0.041666666666666664 + u_1*-0.013888888888888888 + u_2*-0.013888888888888888 + 0.013888888888888888*(u_0*u_0) + 0.013888888888888888*(u_1*u_1) + 0.013888888888888888*(u_2*u_2);
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
      double * RESTRICT _data_pdfs_20_310 = _data_pdfs + _stride_pdfs_2*ctr_2 + 10*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_20_314 = _data_pdfs + _stride_pdfs_2*ctr_2 + 14*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_20_318 = _data_pdfs + _stride_pdfs_2*ctr_2 + 18*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_20_319 = _data_pdfs + _stride_pdfs_2*ctr_2 + 19*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_20_321 = _data_pdfs + _stride_pdfs_2*ctr_2 + 21*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_20_323 = _data_pdfs + _stride_pdfs_2*ctr_2 + 23*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_20_325 = _data_pdfs + _stride_pdfs_2*ctr_2 + 25*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_20_34 = _data_pdfs + _stride_pdfs_2*ctr_2 + 4*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_20_38 = _data_pdfs + _stride_pdfs_2*ctr_2 + 8*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_20_313 = _data_pdfs + _stride_pdfs_2*ctr_2 + 13*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_20_317 = _data_pdfs + _stride_pdfs_2*ctr_2 + 17*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_20_320 = _data_pdfs + _stride_pdfs_2*ctr_2 + 20*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_20_322 = _data_pdfs + _stride_pdfs_2*ctr_2 + 22*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_20_324 = _data_pdfs + _stride_pdfs_2*ctr_2 + 24*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_20_326 = _data_pdfs + _stride_pdfs_2*ctr_2 + 26*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_20_33 = _data_pdfs + _stride_pdfs_2*ctr_2 + 3*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_20_37 = _data_pdfs + _stride_pdfs_2*ctr_2 + 7*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_20_39 = _data_pdfs + _stride_pdfs_2*ctr_2 + 9*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_20_31 = _data_pdfs + _stride_pdfs_2*ctr_2 + _stride_pdfs_3;
      double * RESTRICT _data_pdfs_20_311 = _data_pdfs + _stride_pdfs_2*ctr_2 + 11*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_20_315 = _data_pdfs + _stride_pdfs_2*ctr_2 + 15*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_20_312 = _data_pdfs + _stride_pdfs_2*ctr_2 + 12*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_20_316 = _data_pdfs + _stride_pdfs_2*ctr_2 + 16*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_20_32 = _data_pdfs + _stride_pdfs_2*ctr_2 + 2*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_20_35 = _data_pdfs + _stride_pdfs_2*ctr_2 + 5*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_20_30 = _data_pdfs + _stride_pdfs_2*ctr_2;
      double * RESTRICT _data_pdfs_20_36 = _data_pdfs + _stride_pdfs_2*ctr_2 + 6*_stride_pdfs_3;
      double * RESTRICT  _data_density_20_30 = _data_density + _stride_density_2*ctr_2;
      double * RESTRICT  _data_velocity_20_30 = _data_velocity + _stride_velocity_2*ctr_2;
      double * RESTRICT  _data_velocity_20_31 = _data_velocity + _stride_velocity_2*ctr_2 + _stride_velocity_3;
      double * RESTRICT  _data_velocity_20_32 = _data_velocity + _stride_velocity_2*ctr_2 + 2*_stride_velocity_3;
      for (int64_t ctr_1 = 0; ctr_1 < _size_density_1; ctr_1 += 1)
      {
         double * RESTRICT _data_pdfs_20_310_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_310;
         double * RESTRICT _data_pdfs_20_314_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_314;
         double * RESTRICT _data_pdfs_20_318_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_318;
         double * RESTRICT _data_pdfs_20_319_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_319;
         double * RESTRICT _data_pdfs_20_321_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_321;
         double * RESTRICT _data_pdfs_20_323_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_323;
         double * RESTRICT _data_pdfs_20_325_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_325;
         double * RESTRICT _data_pdfs_20_34_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_34;
         double * RESTRICT _data_pdfs_20_38_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_38;
         double * RESTRICT _data_pdfs_20_313_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_313;
         double * RESTRICT _data_pdfs_20_317_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_317;
         double * RESTRICT _data_pdfs_20_320_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_320;
         double * RESTRICT _data_pdfs_20_322_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_322;
         double * RESTRICT _data_pdfs_20_324_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_324;
         double * RESTRICT _data_pdfs_20_326_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_326;
         double * RESTRICT _data_pdfs_20_33_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_33;
         double * RESTRICT _data_pdfs_20_37_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_37;
         double * RESTRICT _data_pdfs_20_39_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_39;
         double * RESTRICT _data_pdfs_20_31_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_31;
         double * RESTRICT _data_pdfs_20_311_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_311;
         double * RESTRICT _data_pdfs_20_315_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_315;
         double * RESTRICT _data_pdfs_20_312_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_312;
         double * RESTRICT _data_pdfs_20_316_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_316;
         double * RESTRICT _data_pdfs_20_32_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_32;
         double * RESTRICT _data_pdfs_20_35_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_35;
         double * RESTRICT _data_pdfs_20_30_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_30;
         double * RESTRICT _data_pdfs_20_36_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_36;
         double * RESTRICT  _data_density_20_30_10 = _stride_density_1*ctr_1 + _data_density_20_30;
         double * RESTRICT  _data_velocity_20_30_10 = _stride_velocity_1*ctr_1 + _data_velocity_20_30;
         double * RESTRICT  _data_velocity_20_31_10 = _stride_velocity_1*ctr_1 + _data_velocity_20_31;
         double * RESTRICT  _data_velocity_20_32_10 = _stride_velocity_1*ctr_1 + _data_velocity_20_32;
         for (int64_t ctr_0 = 0; ctr_0 < _size_density_0; ctr_0 += 1)
         {
            const double vel0Term = _data_pdfs_20_310_10[_stride_pdfs_0*ctr_0] + _data_pdfs_20_314_10[_stride_pdfs_0*ctr_0] + _data_pdfs_20_318_10[_stride_pdfs_0*ctr_0] + _data_pdfs_20_319_10[_stride_pdfs_0*ctr_0] + _data_pdfs_20_321_10[_stride_pdfs_0*ctr_0] + _data_pdfs_20_323_10[_stride_pdfs_0*ctr_0] + _data_pdfs_20_325_10[_stride_pdfs_0*ctr_0] + _data_pdfs_20_34_10[_stride_pdfs_0*ctr_0] + _data_pdfs_20_38_10[_stride_pdfs_0*ctr_0];
            const double momdensity_0 = vel0Term - 1.0*_data_pdfs_20_313_10[_stride_pdfs_0*ctr_0] - 1.0*_data_pdfs_20_317_10[_stride_pdfs_0*ctr_0] - 1.0*_data_pdfs_20_320_10[_stride_pdfs_0*ctr_0] - 1.0*_data_pdfs_20_322_10[_stride_pdfs_0*ctr_0] - 1.0*_data_pdfs_20_324_10[_stride_pdfs_0*ctr_0] - 1.0*_data_pdfs_20_326_10[_stride_pdfs_0*ctr_0] - 1.0*_data_pdfs_20_33_10[_stride_pdfs_0*ctr_0] - 1.0*_data_pdfs_20_37_10[_stride_pdfs_0*ctr_0] - 1.0*_data_pdfs_20_39_10[_stride_pdfs_0*ctr_0];
            const double vel1Term = _data_pdfs_20_311_10[_stride_pdfs_0*ctr_0] + _data_pdfs_20_315_10[_stride_pdfs_0*ctr_0] + _data_pdfs_20_31_10[_stride_pdfs_0*ctr_0] + _data_pdfs_20_320_10[_stride_pdfs_0*ctr_0] + _data_pdfs_20_324_10[_stride_pdfs_0*ctr_0] + _data_pdfs_20_37_10[_stride_pdfs_0*ctr_0];
            const double momdensity_1 = vel1Term - 1.0*_data_pdfs_20_310_10[_stride_pdfs_0*ctr_0] - 1.0*_data_pdfs_20_312_10[_stride_pdfs_0*ctr_0] - 1.0*_data_pdfs_20_316_10[_stride_pdfs_0*ctr_0] - 1.0*_data_pdfs_20_321_10[_stride_pdfs_0*ctr_0] - 1.0*_data_pdfs_20_322_10[_stride_pdfs_0*ctr_0] - 1.0*_data_pdfs_20_325_10[_stride_pdfs_0*ctr_0] - 1.0*_data_pdfs_20_326_10[_stride_pdfs_0*ctr_0] - 1.0*_data_pdfs_20_32_10[_stride_pdfs_0*ctr_0] - 1.0*_data_pdfs_20_39_10[_stride_pdfs_0*ctr_0] + _data_pdfs_20_319_10[_stride_pdfs_0*ctr_0] + _data_pdfs_20_323_10[_stride_pdfs_0*ctr_0] + _data_pdfs_20_38_10[_stride_pdfs_0*ctr_0];
            const double vel2Term = _data_pdfs_20_312_10[_stride_pdfs_0*ctr_0] + _data_pdfs_20_313_10[_stride_pdfs_0*ctr_0] + _data_pdfs_20_322_10[_stride_pdfs_0*ctr_0] + _data_pdfs_20_35_10[_stride_pdfs_0*ctr_0];
            const double delta_rho = vel0Term + vel1Term + vel2Term + _data_pdfs_20_30_10[_stride_pdfs_0*ctr_0] + _data_pdfs_20_316_10[_stride_pdfs_0*ctr_0] + _data_pdfs_20_317_10[_stride_pdfs_0*ctr_0] + _data_pdfs_20_326_10[_stride_pdfs_0*ctr_0] + _data_pdfs_20_32_10[_stride_pdfs_0*ctr_0] + _data_pdfs_20_33_10[_stride_pdfs_0*ctr_0] + _data_pdfs_20_36_10[_stride_pdfs_0*ctr_0] + _data_pdfs_20_39_10[_stride_pdfs_0*ctr_0];
            const double momdensity_2 = vel2Term - 1.0*_data_pdfs_20_315_10[_stride_pdfs_0*ctr_0] - 1.0*_data_pdfs_20_316_10[_stride_pdfs_0*ctr_0] - 1.0*_data_pdfs_20_317_10[_stride_pdfs_0*ctr_0] - 1.0*_data_pdfs_20_318_10[_stride_pdfs_0*ctr_0] - 1.0*_data_pdfs_20_323_10[_stride_pdfs_0*ctr_0] - 1.0*_data_pdfs_20_324_10[_stride_pdfs_0*ctr_0] - 1.0*_data_pdfs_20_325_10[_stride_pdfs_0*ctr_0] - 1.0*_data_pdfs_20_326_10[_stride_pdfs_0*ctr_0] - 1.0*_data_pdfs_20_36_10[_stride_pdfs_0*ctr_0] + _data_pdfs_20_311_10[_stride_pdfs_0*ctr_0] + _data_pdfs_20_314_10[_stride_pdfs_0*ctr_0] + _data_pdfs_20_319_10[_stride_pdfs_0*ctr_0] + _data_pdfs_20_320_10[_stride_pdfs_0*ctr_0] + _data_pdfs_20_321_10[_stride_pdfs_0*ctr_0];
            const double rho = delta_rho + 1.0;
            const double u_0 = momdensity_0;
            const double u_1 = momdensity_1;
            const double u_2 = momdensity_2;
            _data_density_20_30_10[_stride_density_0*ctr_0] = rho;
            _data_velocity_20_30_10[_stride_velocity_0*ctr_0] = u_0;
            _data_velocity_20_31_10[_stride_velocity_0*ctr_0] = u_1;
            _data_velocity_20_32_10[_stride_velocity_0*ctr_0] = u_2;
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