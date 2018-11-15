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
//! \\author Martin Bauer <martin.bauer@fau.de>
//======================================================================================================================

#include <cmath>

#include "core/DataTypes.h"
#include "core/Macros.h"
#include "lbm/field/PdfField.h"
#include "lbm/sweeps/Streaming.h"
#include "UniformGridGPU_LatticeModel.h"

#ifdef _MSC_VER
#  pragma warning( disable : 4458 )
#endif

#define FUNC_PREFIX

using namespace std;

namespace walberla {
namespace lbm {

namespace internal_kernel_streamCollide {
static FUNC_PREFIX void kernel_streamCollide(double * const _data_pdfs, double * _data_pdfs_tmp, int64_t const _size_pdfs_0, int64_t const _size_pdfs_1, int64_t const _size_pdfs_2, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2, int64_t const _stride_pdfs_3, double omega)
{
   const double xi_1 = omega*0.166666666666667;
   const double xi_5 = omega*0.0416666666666667;
   for (int ctr_2 = 1; ctr_2 < _size_pdfs_2 - 1; ctr_2 += 1)
   {
      double * const _data_pdfs_21_317 = _data_pdfs + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 17*_stride_pdfs_3;
      double * const _data_pdfs_20_39 = _data_pdfs + _stride_pdfs_2*ctr_2 + 9*_stride_pdfs_3;
      double * const _data_pdfs_21_316 = _data_pdfs + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 16*_stride_pdfs_3;
      double * const _data_pdfs_2m1_314 = _data_pdfs + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 14*_stride_pdfs_3;
      double * const _data_pdfs_20_310 = _data_pdfs + _stride_pdfs_2*ctr_2 + 10*_stride_pdfs_3;
      double * const _data_pdfs_20_38 = _data_pdfs + _stride_pdfs_2*ctr_2 + 8*_stride_pdfs_3;
      double * const _data_pdfs_21_318 = _data_pdfs + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 18*_stride_pdfs_3;
      double * const _data_pdfs_20_34 = _data_pdfs + _stride_pdfs_2*ctr_2 + 4*_stride_pdfs_3;
      double * const _data_pdfs_2m1_311 = _data_pdfs + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 11*_stride_pdfs_3;
      double * const _data_pdfs_20_37 = _data_pdfs + _stride_pdfs_2*ctr_2 + 7*_stride_pdfs_3;
      double * const _data_pdfs_20_31 = _data_pdfs + _stride_pdfs_2*ctr_2 + _stride_pdfs_3;
      double * const _data_pdfs_21_315 = _data_pdfs + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 15*_stride_pdfs_3;
      double * const _data_pdfs_2m1_313 = _data_pdfs + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 13*_stride_pdfs_3;
      double * const _data_pdfs_2m1_312 = _data_pdfs + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 12*_stride_pdfs_3;
      double * const _data_pdfs_2m1_35 = _data_pdfs + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 5*_stride_pdfs_3;
      double * const _data_pdfs_20_30 = _data_pdfs + _stride_pdfs_2*ctr_2;
      double * const _data_pdfs_20_33 = _data_pdfs + _stride_pdfs_2*ctr_2 + 3*_stride_pdfs_3;
      double * const _data_pdfs_20_32 = _data_pdfs + _stride_pdfs_2*ctr_2 + 2*_stride_pdfs_3;
      double * const _data_pdfs_21_36 = _data_pdfs + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 6*_stride_pdfs_3;
      double * _data_pdfs_tmp_20_30 = _data_pdfs_tmp + _stride_pdfs_2*ctr_2;
      double * _data_pdfs_tmp_20_31 = _data_pdfs_tmp + _stride_pdfs_2*ctr_2 + _stride_pdfs_3;
      double * _data_pdfs_tmp_20_32 = _data_pdfs_tmp + _stride_pdfs_2*ctr_2 + 2*_stride_pdfs_3;
      double * _data_pdfs_tmp_20_33 = _data_pdfs_tmp + _stride_pdfs_2*ctr_2 + 3*_stride_pdfs_3;
      double * _data_pdfs_tmp_20_34 = _data_pdfs_tmp + _stride_pdfs_2*ctr_2 + 4*_stride_pdfs_3;
      double * _data_pdfs_tmp_20_35 = _data_pdfs_tmp + _stride_pdfs_2*ctr_2 + 5*_stride_pdfs_3;
      double * _data_pdfs_tmp_20_36 = _data_pdfs_tmp + _stride_pdfs_2*ctr_2 + 6*_stride_pdfs_3;
      double * _data_pdfs_tmp_20_37 = _data_pdfs_tmp + _stride_pdfs_2*ctr_2 + 7*_stride_pdfs_3;
      double * _data_pdfs_tmp_20_38 = _data_pdfs_tmp + _stride_pdfs_2*ctr_2 + 8*_stride_pdfs_3;
      double * _data_pdfs_tmp_20_39 = _data_pdfs_tmp + _stride_pdfs_2*ctr_2 + 9*_stride_pdfs_3;
      double * _data_pdfs_tmp_20_310 = _data_pdfs_tmp + _stride_pdfs_2*ctr_2 + 10*_stride_pdfs_3;
      double * _data_pdfs_tmp_20_311 = _data_pdfs_tmp + _stride_pdfs_2*ctr_2 + 11*_stride_pdfs_3;
      double * _data_pdfs_tmp_20_312 = _data_pdfs_tmp + _stride_pdfs_2*ctr_2 + 12*_stride_pdfs_3;
      double * _data_pdfs_tmp_20_313 = _data_pdfs_tmp + _stride_pdfs_2*ctr_2 + 13*_stride_pdfs_3;
      double * _data_pdfs_tmp_20_314 = _data_pdfs_tmp + _stride_pdfs_2*ctr_2 + 14*_stride_pdfs_3;
      double * _data_pdfs_tmp_20_315 = _data_pdfs_tmp + _stride_pdfs_2*ctr_2 + 15*_stride_pdfs_3;
      double * _data_pdfs_tmp_20_316 = _data_pdfs_tmp + _stride_pdfs_2*ctr_2 + 16*_stride_pdfs_3;
      double * _data_pdfs_tmp_20_317 = _data_pdfs_tmp + _stride_pdfs_2*ctr_2 + 17*_stride_pdfs_3;
      double * _data_pdfs_tmp_20_318 = _data_pdfs_tmp + _stride_pdfs_2*ctr_2 + 18*_stride_pdfs_3;
      for (int ctr_1 = 1; ctr_1 < _size_pdfs_1 - 1; ctr_1 += 1)
      {
         double * const _data_pdfs_21_317_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_21_317;
         double * const _data_pdfs_20_39_11 = _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _data_pdfs_20_39;
         double * const _data_pdfs_21_316_11 = _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _data_pdfs_21_316;
         double * const _data_pdfs_2m1_314_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_2m1_314;
         double * const _data_pdfs_20_310_11 = _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _data_pdfs_20_310;
         double * const _data_pdfs_20_38_1m1 = _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _data_pdfs_20_38;
         double * const _data_pdfs_21_318_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_21_318;
         double * const _data_pdfs_20_34_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_34;
         double * const _data_pdfs_2m1_311_1m1 = _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _data_pdfs_2m1_311;
         double * const _data_pdfs_20_37_1m1 = _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _data_pdfs_20_37;
         double * const _data_pdfs_20_31_1m1 = _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _data_pdfs_20_31;
         double * const _data_pdfs_21_315_1m1 = _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _data_pdfs_21_315;
         double * const _data_pdfs_2m1_313_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_2m1_313;
         double * const _data_pdfs_2m1_312_11 = _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _data_pdfs_2m1_312;
         double * const _data_pdfs_2m1_35_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_2m1_35;
         double * const _data_pdfs_20_30_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_30;
         double * const _data_pdfs_20_33_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_33;
         double * const _data_pdfs_20_32_11 = _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _data_pdfs_20_32;
         double * const _data_pdfs_21_36_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_21_36;
         double * _data_pdfs_tmp_20_30_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_tmp_20_30;
         double * _data_pdfs_tmp_20_31_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_tmp_20_31;
         double * _data_pdfs_tmp_20_32_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_tmp_20_32;
         double * _data_pdfs_tmp_20_33_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_tmp_20_33;
         double * _data_pdfs_tmp_20_34_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_tmp_20_34;
         double * _data_pdfs_tmp_20_35_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_tmp_20_35;
         double * _data_pdfs_tmp_20_36_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_tmp_20_36;
         double * _data_pdfs_tmp_20_37_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_tmp_20_37;
         double * _data_pdfs_tmp_20_38_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_tmp_20_38;
         double * _data_pdfs_tmp_20_39_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_tmp_20_39;
         double * _data_pdfs_tmp_20_310_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_tmp_20_310;
         double * _data_pdfs_tmp_20_311_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_tmp_20_311;
         double * _data_pdfs_tmp_20_312_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_tmp_20_312;
         double * _data_pdfs_tmp_20_313_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_tmp_20_313;
         double * _data_pdfs_tmp_20_314_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_tmp_20_314;
         double * _data_pdfs_tmp_20_315_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_tmp_20_315;
         double * _data_pdfs_tmp_20_316_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_tmp_20_316;
         double * _data_pdfs_tmp_20_317_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_tmp_20_317;
         double * _data_pdfs_tmp_20_318_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_tmp_20_318;
         for (int ctr_0 = 1; ctr_0 < _size_pdfs_0 - 1; ctr_0 += 1)
         {
            const double xi_18 = -_data_pdfs_21_317_10[_stride_pdfs_0*ctr_0 + _stride_pdfs_0];
            const double xi_19 = -_data_pdfs_20_39_11[_stride_pdfs_0*ctr_0 + _stride_pdfs_0];
            const double xi_20 = -_data_pdfs_21_316_11[_stride_pdfs_0*ctr_0];
            const double vel0Term = _data_pdfs_20_310_11[_stride_pdfs_0*ctr_0 - _stride_pdfs_0] + _data_pdfs_20_34_10[_stride_pdfs_0*ctr_0 - _stride_pdfs_0] + _data_pdfs_20_38_1m1[_stride_pdfs_0*ctr_0 - _stride_pdfs_0] + _data_pdfs_21_318_10[_stride_pdfs_0*ctr_0 - _stride_pdfs_0] + _data_pdfs_2m1_314_10[_stride_pdfs_0*ctr_0 - _stride_pdfs_0];
            const double vel1Term = _data_pdfs_20_31_1m1[_stride_pdfs_0*ctr_0] + _data_pdfs_20_37_1m1[_stride_pdfs_0*ctr_0 + _stride_pdfs_0] + _data_pdfs_21_315_1m1[_stride_pdfs_0*ctr_0] + _data_pdfs_2m1_311_1m1[_stride_pdfs_0*ctr_0];
            const double vel2Term = _data_pdfs_2m1_312_11[_stride_pdfs_0*ctr_0] + _data_pdfs_2m1_313_10[_stride_pdfs_0*ctr_0 + _stride_pdfs_0] + _data_pdfs_2m1_35_10[_stride_pdfs_0*ctr_0];
            const double rho = vel0Term + vel1Term + vel2Term + _data_pdfs_20_30_10[_stride_pdfs_0*ctr_0] + _data_pdfs_20_32_11[_stride_pdfs_0*ctr_0] + _data_pdfs_20_33_10[_stride_pdfs_0*ctr_0 + _stride_pdfs_0] + _data_pdfs_20_39_11[_stride_pdfs_0*ctr_0 + _stride_pdfs_0] + _data_pdfs_21_316_11[_stride_pdfs_0*ctr_0] + _data_pdfs_21_317_10[_stride_pdfs_0*ctr_0 + _stride_pdfs_0] + _data_pdfs_21_36_10[_stride_pdfs_0*ctr_0];
            const double xi_27 = rho*-0.333333333333333;
            const double u_0 = vel0Term + xi_18 + xi_19 - _data_pdfs_20_33_10[_stride_pdfs_0*ctr_0 + _stride_pdfs_0] - _data_pdfs_20_37_1m1[_stride_pdfs_0*ctr_0 + _stride_pdfs_0] - _data_pdfs_2m1_313_10[_stride_pdfs_0*ctr_0 + _stride_pdfs_0];
            const double xi_23 = (u_0*u_0);
            const double u_1 = vel1Term + xi_19 + xi_20 - _data_pdfs_20_310_11[_stride_pdfs_0*ctr_0 - _stride_pdfs_0] - _data_pdfs_20_32_11[_stride_pdfs_0*ctr_0] + _data_pdfs_20_38_1m1[_stride_pdfs_0*ctr_0 - _stride_pdfs_0] - _data_pdfs_2m1_312_11[_stride_pdfs_0*ctr_0];
            const double xi_21 = -u_1;
            const double xi_24 = (u_1*u_1);
            const double u_2 = vel2Term + xi_18 + xi_20 - _data_pdfs_21_315_1m1[_stride_pdfs_0*ctr_0] - _data_pdfs_21_318_10[_stride_pdfs_0*ctr_0 - _stride_pdfs_0] - _data_pdfs_21_36_10[_stride_pdfs_0*ctr_0] + _data_pdfs_2m1_311_1m1[_stride_pdfs_0*ctr_0] + _data_pdfs_2m1_314_10[_stride_pdfs_0*ctr_0 - _stride_pdfs_0];
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
            const double xi_28 = f_eq_common + xi_25 + xi_27;
            const double xi_29 = f_eq_common + xi_23 + xi_27;
            const double xi_30 = f_eq_common + xi_24 + xi_27;
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
            _data_pdfs_tmp_20_30_10[_stride_pdfs_0*ctr_0] = omega*(f_eq_common*0.333333333333333 - _data_pdfs_20_30_10[_stride_pdfs_0*ctr_0]) + _data_pdfs_20_30_10[_stride_pdfs_0*ctr_0];
            _data_pdfs_tmp_20_31_10[_stride_pdfs_0*ctr_0] = xi_1*(u_1 + xi_2 - 6*_data_pdfs_20_31_1m1[_stride_pdfs_0*ctr_0]) + _data_pdfs_20_31_1m1[_stride_pdfs_0*ctr_0];
            _data_pdfs_tmp_20_32_10[_stride_pdfs_0*ctr_0] = xi_1*(xi_2 + xi_21 - 6*_data_pdfs_20_32_11[_stride_pdfs_0*ctr_0]) + _data_pdfs_20_32_11[_stride_pdfs_0*ctr_0];
            _data_pdfs_tmp_20_33_10[_stride_pdfs_0*ctr_0] = xi_1*(-u_0 + xi_3 - 6*_data_pdfs_20_33_10[_stride_pdfs_0*ctr_0 + _stride_pdfs_0]) + _data_pdfs_20_33_10[_stride_pdfs_0*ctr_0 + _stride_pdfs_0];
            _data_pdfs_tmp_20_34_10[_stride_pdfs_0*ctr_0] = xi_1*(u_0 + xi_3 - 6*_data_pdfs_20_34_10[_stride_pdfs_0*ctr_0 - _stride_pdfs_0]) + _data_pdfs_20_34_10[_stride_pdfs_0*ctr_0 - _stride_pdfs_0];
            _data_pdfs_tmp_20_35_10[_stride_pdfs_0*ctr_0] = xi_1*(u_2 + xi_4 - 6*_data_pdfs_2m1_35_10[_stride_pdfs_0*ctr_0]) + _data_pdfs_2m1_35_10[_stride_pdfs_0*ctr_0];
            _data_pdfs_tmp_20_36_10[_stride_pdfs_0*ctr_0] = xi_1*(xi_22 + xi_4 - 6*_data_pdfs_21_36_10[_stride_pdfs_0*ctr_0]) + _data_pdfs_21_36_10[_stride_pdfs_0*ctr_0];
            _data_pdfs_tmp_20_37_10[_stride_pdfs_0*ctr_0] = xi_5*(-xi_6 + xi_7 - 24*_data_pdfs_20_37_1m1[_stride_pdfs_0*ctr_0 + _stride_pdfs_0]) + _data_pdfs_20_37_1m1[_stride_pdfs_0*ctr_0 + _stride_pdfs_0];
            _data_pdfs_tmp_20_38_10[_stride_pdfs_0*ctr_0] = xi_5*(xi_8 + xi_9 - 24*_data_pdfs_20_38_1m1[_stride_pdfs_0*ctr_0 - _stride_pdfs_0]) + _data_pdfs_20_38_1m1[_stride_pdfs_0*ctr_0 - _stride_pdfs_0];
            _data_pdfs_tmp_20_39_10[_stride_pdfs_0*ctr_0] = xi_5*(-xi_8 + xi_9 - 24*_data_pdfs_20_39_11[_stride_pdfs_0*ctr_0 + _stride_pdfs_0]) + _data_pdfs_20_39_11[_stride_pdfs_0*ctr_0 + _stride_pdfs_0];
            _data_pdfs_tmp_20_310_10[_stride_pdfs_0*ctr_0] = xi_5*(xi_6 + xi_7 - 24*_data_pdfs_20_310_11[_stride_pdfs_0*ctr_0 - _stride_pdfs_0]) + _data_pdfs_20_310_11[_stride_pdfs_0*ctr_0 - _stride_pdfs_0];
            _data_pdfs_tmp_20_311_10[_stride_pdfs_0*ctr_0] = xi_5*(xi_10 + xi_11 - 24*_data_pdfs_2m1_311_1m1[_stride_pdfs_0*ctr_0]) + _data_pdfs_2m1_311_1m1[_stride_pdfs_0*ctr_0];
            _data_pdfs_tmp_20_312_10[_stride_pdfs_0*ctr_0] = xi_5*(-xi_12 + xi_13 - 24*_data_pdfs_2m1_312_11[_stride_pdfs_0*ctr_0]) + _data_pdfs_2m1_312_11[_stride_pdfs_0*ctr_0];
            _data_pdfs_tmp_20_313_10[_stride_pdfs_0*ctr_0] = xi_5*(-xi_14 + xi_15 - 24*_data_pdfs_2m1_313_10[_stride_pdfs_0*ctr_0 + _stride_pdfs_0]) + _data_pdfs_2m1_313_10[_stride_pdfs_0*ctr_0 + _stride_pdfs_0];
            _data_pdfs_tmp_20_314_10[_stride_pdfs_0*ctr_0] = xi_5*(xi_16 + xi_17 - 24*_data_pdfs_2m1_314_10[_stride_pdfs_0*ctr_0 - _stride_pdfs_0]) + _data_pdfs_2m1_314_10[_stride_pdfs_0*ctr_0 - _stride_pdfs_0];
            _data_pdfs_tmp_20_315_10[_stride_pdfs_0*ctr_0] = xi_5*(xi_12 + xi_13 - 24*_data_pdfs_21_315_1m1[_stride_pdfs_0*ctr_0]) + _data_pdfs_21_315_1m1[_stride_pdfs_0*ctr_0];
            _data_pdfs_tmp_20_316_10[_stride_pdfs_0*ctr_0] = xi_5*(-xi_10 + xi_11 - 24*_data_pdfs_21_316_11[_stride_pdfs_0*ctr_0]) + _data_pdfs_21_316_11[_stride_pdfs_0*ctr_0];
            _data_pdfs_tmp_20_317_10[_stride_pdfs_0*ctr_0] = xi_5*(-xi_16 + xi_17 - 24*_data_pdfs_21_317_10[_stride_pdfs_0*ctr_0 + _stride_pdfs_0]) + _data_pdfs_21_317_10[_stride_pdfs_0*ctr_0 + _stride_pdfs_0];
            _data_pdfs_tmp_20_318_10[_stride_pdfs_0*ctr_0] = xi_5*(xi_14 + xi_15 - 24*_data_pdfs_21_318_10[_stride_pdfs_0*ctr_0 - _stride_pdfs_0]) + _data_pdfs_21_318_10[_stride_pdfs_0*ctr_0 - _stride_pdfs_0];
         }
      }
   }
}
}
namespace internal_kernel_collide {
static FUNC_PREFIX void kernel_collide(double * _data_pdfs, int64_t const _size_pdfs_0, int64_t const _size_pdfs_1, int64_t const _size_pdfs_2, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2, int64_t const _stride_pdfs_3, double omega)
{
   const double xi_1 = omega*0.166666666666667;
   const double xi_5 = omega*0.0416666666666667;
   for (int ctr_2 = 1; ctr_2 < _size_pdfs_2 - 1; ctr_2 += 1)
   {
      double * _data_pdfs_20_317 = _data_pdfs + _stride_pdfs_2*ctr_2 + 17*_stride_pdfs_3;
      double * _data_pdfs_20_39 = _data_pdfs + _stride_pdfs_2*ctr_2 + 9*_stride_pdfs_3;
      double * _data_pdfs_20_316 = _data_pdfs + _stride_pdfs_2*ctr_2 + 16*_stride_pdfs_3;
      double * _data_pdfs_20_310 = _data_pdfs + _stride_pdfs_2*ctr_2 + 10*_stride_pdfs_3;
      double * _data_pdfs_20_314 = _data_pdfs + _stride_pdfs_2*ctr_2 + 14*_stride_pdfs_3;
      double * _data_pdfs_20_318 = _data_pdfs + _stride_pdfs_2*ctr_2 + 18*_stride_pdfs_3;
      double * _data_pdfs_20_34 = _data_pdfs + _stride_pdfs_2*ctr_2 + 4*_stride_pdfs_3;
      double * _data_pdfs_20_38 = _data_pdfs + _stride_pdfs_2*ctr_2 + 8*_stride_pdfs_3;
      double * _data_pdfs_20_31 = _data_pdfs + _stride_pdfs_2*ctr_2 + _stride_pdfs_3;
      double * _data_pdfs_20_311 = _data_pdfs + _stride_pdfs_2*ctr_2 + 11*_stride_pdfs_3;
      double * _data_pdfs_20_315 = _data_pdfs + _stride_pdfs_2*ctr_2 + 15*_stride_pdfs_3;
      double * _data_pdfs_20_37 = _data_pdfs + _stride_pdfs_2*ctr_2 + 7*_stride_pdfs_3;
      double * _data_pdfs_20_312 = _data_pdfs + _stride_pdfs_2*ctr_2 + 12*_stride_pdfs_3;
      double * _data_pdfs_20_313 = _data_pdfs + _stride_pdfs_2*ctr_2 + 13*_stride_pdfs_3;
      double * _data_pdfs_20_35 = _data_pdfs + _stride_pdfs_2*ctr_2 + 5*_stride_pdfs_3;
      double * _data_pdfs_20_30 = _data_pdfs + _stride_pdfs_2*ctr_2;
      double * _data_pdfs_20_32 = _data_pdfs + _stride_pdfs_2*ctr_2 + 2*_stride_pdfs_3;
      double * _data_pdfs_20_33 = _data_pdfs + _stride_pdfs_2*ctr_2 + 3*_stride_pdfs_3;
      double * _data_pdfs_20_36 = _data_pdfs + _stride_pdfs_2*ctr_2 + 6*_stride_pdfs_3;
      for (int ctr_1 = 1; ctr_1 < _size_pdfs_1 - 1; ctr_1 += 1)
      {
         double * _data_pdfs_20_317_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_317;
         double * _data_pdfs_20_39_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_39;
         double * _data_pdfs_20_316_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_316;
         double * _data_pdfs_20_310_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_310;
         double * _data_pdfs_20_314_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_314;
         double * _data_pdfs_20_318_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_318;
         double * _data_pdfs_20_34_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_34;
         double * _data_pdfs_20_38_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_38;
         double * _data_pdfs_20_31_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_31;
         double * _data_pdfs_20_311_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_311;
         double * _data_pdfs_20_315_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_315;
         double * _data_pdfs_20_37_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_37;
         double * _data_pdfs_20_312_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_312;
         double * _data_pdfs_20_313_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_313;
         double * _data_pdfs_20_35_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_35;
         double * _data_pdfs_20_30_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_30;
         double * _data_pdfs_20_32_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_32;
         double * _data_pdfs_20_33_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_33;
         double * _data_pdfs_20_36_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_36;
         for (int ctr_0 = 1; ctr_0 < _size_pdfs_0 - 1; ctr_0 += 1)
         {
            const double xi_18 = -_data_pdfs_20_317_10[_stride_pdfs_0*ctr_0];
            const double xi_19 = -_data_pdfs_20_39_10[_stride_pdfs_0*ctr_0];
            const double xi_20 = -_data_pdfs_20_316_10[_stride_pdfs_0*ctr_0];
            const double vel0Term = _data_pdfs_20_310_10[_stride_pdfs_0*ctr_0] + _data_pdfs_20_314_10[_stride_pdfs_0*ctr_0] + _data_pdfs_20_318_10[_stride_pdfs_0*ctr_0] + _data_pdfs_20_34_10[_stride_pdfs_0*ctr_0] + _data_pdfs_20_38_10[_stride_pdfs_0*ctr_0];
            const double vel1Term = _data_pdfs_20_311_10[_stride_pdfs_0*ctr_0] + _data_pdfs_20_315_10[_stride_pdfs_0*ctr_0] + _data_pdfs_20_31_10[_stride_pdfs_0*ctr_0] + _data_pdfs_20_37_10[_stride_pdfs_0*ctr_0];
            const double vel2Term = _data_pdfs_20_312_10[_stride_pdfs_0*ctr_0] + _data_pdfs_20_313_10[_stride_pdfs_0*ctr_0] + _data_pdfs_20_35_10[_stride_pdfs_0*ctr_0];
            const double rho = vel0Term + vel1Term + vel2Term + _data_pdfs_20_30_10[_stride_pdfs_0*ctr_0] + _data_pdfs_20_316_10[_stride_pdfs_0*ctr_0] + _data_pdfs_20_317_10[_stride_pdfs_0*ctr_0] + _data_pdfs_20_32_10[_stride_pdfs_0*ctr_0] + _data_pdfs_20_33_10[_stride_pdfs_0*ctr_0] + _data_pdfs_20_36_10[_stride_pdfs_0*ctr_0] + _data_pdfs_20_39_10[_stride_pdfs_0*ctr_0];
            const double xi_27 = rho*-0.333333333333333;
            const double u_0 = vel0Term + xi_18 + xi_19 - _data_pdfs_20_313_10[_stride_pdfs_0*ctr_0] - _data_pdfs_20_33_10[_stride_pdfs_0*ctr_0] - _data_pdfs_20_37_10[_stride_pdfs_0*ctr_0];
            const double xi_23 = (u_0*u_0);
            const double u_1 = vel1Term + xi_19 + xi_20 - _data_pdfs_20_310_10[_stride_pdfs_0*ctr_0] - _data_pdfs_20_312_10[_stride_pdfs_0*ctr_0] - _data_pdfs_20_32_10[_stride_pdfs_0*ctr_0] + _data_pdfs_20_38_10[_stride_pdfs_0*ctr_0];
            const double xi_21 = -u_1;
            const double xi_24 = (u_1*u_1);
            const double u_2 = vel2Term + xi_18 + xi_20 + _data_pdfs_20_311_10[_stride_pdfs_0*ctr_0] + _data_pdfs_20_314_10[_stride_pdfs_0*ctr_0] - _data_pdfs_20_315_10[_stride_pdfs_0*ctr_0] - _data_pdfs_20_318_10[_stride_pdfs_0*ctr_0] - _data_pdfs_20_36_10[_stride_pdfs_0*ctr_0];
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
            const double xi_28 = f_eq_common + xi_25 + xi_27;
            const double xi_29 = f_eq_common + xi_23 + xi_27;
            const double xi_30 = f_eq_common + xi_24 + xi_27;
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
            _data_pdfs_20_30_10[_stride_pdfs_0*ctr_0] = omega*(f_eq_common*0.333333333333333 - _data_pdfs_20_30_10[_stride_pdfs_0*ctr_0]) + _data_pdfs_20_30_10[_stride_pdfs_0*ctr_0];
            _data_pdfs_20_31_10[_stride_pdfs_0*ctr_0] = xi_1*(u_1 + xi_2 - 6*_data_pdfs_20_31_10[_stride_pdfs_0*ctr_0]) + _data_pdfs_20_31_10[_stride_pdfs_0*ctr_0];
            _data_pdfs_20_32_10[_stride_pdfs_0*ctr_0] = xi_1*(xi_2 + xi_21 - 6*_data_pdfs_20_32_10[_stride_pdfs_0*ctr_0]) + _data_pdfs_20_32_10[_stride_pdfs_0*ctr_0];
            _data_pdfs_20_33_10[_stride_pdfs_0*ctr_0] = xi_1*(-u_0 + xi_3 - 6*_data_pdfs_20_33_10[_stride_pdfs_0*ctr_0]) + _data_pdfs_20_33_10[_stride_pdfs_0*ctr_0];
            _data_pdfs_20_34_10[_stride_pdfs_0*ctr_0] = xi_1*(u_0 + xi_3 - 6*_data_pdfs_20_34_10[_stride_pdfs_0*ctr_0]) + _data_pdfs_20_34_10[_stride_pdfs_0*ctr_0];
            _data_pdfs_20_35_10[_stride_pdfs_0*ctr_0] = xi_1*(u_2 + xi_4 - 6*_data_pdfs_20_35_10[_stride_pdfs_0*ctr_0]) + _data_pdfs_20_35_10[_stride_pdfs_0*ctr_0];
            _data_pdfs_20_36_10[_stride_pdfs_0*ctr_0] = xi_1*(xi_22 + xi_4 - 6*_data_pdfs_20_36_10[_stride_pdfs_0*ctr_0]) + _data_pdfs_20_36_10[_stride_pdfs_0*ctr_0];
            _data_pdfs_20_37_10[_stride_pdfs_0*ctr_0] = xi_5*(-xi_6 + xi_7 - 24*_data_pdfs_20_37_10[_stride_pdfs_0*ctr_0]) + _data_pdfs_20_37_10[_stride_pdfs_0*ctr_0];
            _data_pdfs_20_38_10[_stride_pdfs_0*ctr_0] = xi_5*(xi_8 + xi_9 - 24*_data_pdfs_20_38_10[_stride_pdfs_0*ctr_0]) + _data_pdfs_20_38_10[_stride_pdfs_0*ctr_0];
            _data_pdfs_20_39_10[_stride_pdfs_0*ctr_0] = xi_5*(-xi_8 + xi_9 - 24*_data_pdfs_20_39_10[_stride_pdfs_0*ctr_0]) + _data_pdfs_20_39_10[_stride_pdfs_0*ctr_0];
            _data_pdfs_20_310_10[_stride_pdfs_0*ctr_0] = xi_5*(xi_6 + xi_7 - 24*_data_pdfs_20_310_10[_stride_pdfs_0*ctr_0]) + _data_pdfs_20_310_10[_stride_pdfs_0*ctr_0];
            _data_pdfs_20_311_10[_stride_pdfs_0*ctr_0] = xi_5*(xi_10 + xi_11 - 24*_data_pdfs_20_311_10[_stride_pdfs_0*ctr_0]) + _data_pdfs_20_311_10[_stride_pdfs_0*ctr_0];
            _data_pdfs_20_312_10[_stride_pdfs_0*ctr_0] = xi_5*(-xi_12 + xi_13 - 24*_data_pdfs_20_312_10[_stride_pdfs_0*ctr_0]) + _data_pdfs_20_312_10[_stride_pdfs_0*ctr_0];
            _data_pdfs_20_313_10[_stride_pdfs_0*ctr_0] = xi_5*(-xi_14 + xi_15 - 24*_data_pdfs_20_313_10[_stride_pdfs_0*ctr_0]) + _data_pdfs_20_313_10[_stride_pdfs_0*ctr_0];
            _data_pdfs_20_314_10[_stride_pdfs_0*ctr_0] = xi_5*(xi_16 + xi_17 - 24*_data_pdfs_20_314_10[_stride_pdfs_0*ctr_0]) + _data_pdfs_20_314_10[_stride_pdfs_0*ctr_0];
            _data_pdfs_20_315_10[_stride_pdfs_0*ctr_0] = xi_5*(xi_12 + xi_13 - 24*_data_pdfs_20_315_10[_stride_pdfs_0*ctr_0]) + _data_pdfs_20_315_10[_stride_pdfs_0*ctr_0];
            _data_pdfs_20_316_10[_stride_pdfs_0*ctr_0] = xi_5*(-xi_10 + xi_11 - 24*_data_pdfs_20_316_10[_stride_pdfs_0*ctr_0]) + _data_pdfs_20_316_10[_stride_pdfs_0*ctr_0];
            _data_pdfs_20_317_10[_stride_pdfs_0*ctr_0] = xi_5*(-xi_16 + xi_17 - 24*_data_pdfs_20_317_10[_stride_pdfs_0*ctr_0]) + _data_pdfs_20_317_10[_stride_pdfs_0*ctr_0];
            _data_pdfs_20_318_10[_stride_pdfs_0*ctr_0] = xi_5*(xi_14 + xi_15 - 24*_data_pdfs_20_318_10[_stride_pdfs_0*ctr_0]) + _data_pdfs_20_318_10[_stride_pdfs_0*ctr_0];
         }
      }
   }
}
}
namespace internal_kernel_stream {
static FUNC_PREFIX void kernel_stream(double * const _data_pdfs, double * _data_pdfs_tmp, int64_t const _size_pdfs_0, int64_t const _size_pdfs_1, int64_t const _size_pdfs_2, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2, int64_t const _stride_pdfs_3, int64_t const _stride_pdfs_tmp_0, int64_t const _stride_pdfs_tmp_1, int64_t const _stride_pdfs_tmp_2, int64_t const _stride_pdfs_tmp_3)
{
   for (int ctr_2 = 1; ctr_2 < _size_pdfs_2 - 1; ctr_2 += 1)
   {
      double * _data_pdfs_tmp_20_30 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2;
      double * const _data_pdfs_20_30 = _data_pdfs + _stride_pdfs_2*ctr_2;
      double * _data_pdfs_tmp_20_31 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + _stride_pdfs_tmp_3;
      double * const _data_pdfs_20_31 = _data_pdfs + _stride_pdfs_2*ctr_2 + _stride_pdfs_3;
      double * _data_pdfs_tmp_20_32 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 2*_stride_pdfs_tmp_3;
      double * const _data_pdfs_20_32 = _data_pdfs + _stride_pdfs_2*ctr_2 + 2*_stride_pdfs_3;
      double * _data_pdfs_tmp_20_33 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 3*_stride_pdfs_tmp_3;
      double * const _data_pdfs_20_33 = _data_pdfs + _stride_pdfs_2*ctr_2 + 3*_stride_pdfs_3;
      double * _data_pdfs_tmp_20_34 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 4*_stride_pdfs_tmp_3;
      double * const _data_pdfs_20_34 = _data_pdfs + _stride_pdfs_2*ctr_2 + 4*_stride_pdfs_3;
      double * _data_pdfs_tmp_20_35 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 5*_stride_pdfs_tmp_3;
      double * const _data_pdfs_2m1_35 = _data_pdfs + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 5*_stride_pdfs_3;
      double * _data_pdfs_tmp_20_36 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 6*_stride_pdfs_tmp_3;
      double * const _data_pdfs_21_36 = _data_pdfs + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 6*_stride_pdfs_3;
      double * _data_pdfs_tmp_20_37 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 7*_stride_pdfs_tmp_3;
      double * const _data_pdfs_20_37 = _data_pdfs + _stride_pdfs_2*ctr_2 + 7*_stride_pdfs_3;
      double * _data_pdfs_tmp_20_38 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 8*_stride_pdfs_tmp_3;
      double * const _data_pdfs_20_38 = _data_pdfs + _stride_pdfs_2*ctr_2 + 8*_stride_pdfs_3;
      double * _data_pdfs_tmp_20_39 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 9*_stride_pdfs_tmp_3;
      double * const _data_pdfs_20_39 = _data_pdfs + _stride_pdfs_2*ctr_2 + 9*_stride_pdfs_3;
      double * _data_pdfs_tmp_20_310 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 10*_stride_pdfs_tmp_3;
      double * const _data_pdfs_20_310 = _data_pdfs + _stride_pdfs_2*ctr_2 + 10*_stride_pdfs_3;
      double * _data_pdfs_tmp_20_311 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 11*_stride_pdfs_tmp_3;
      double * const _data_pdfs_2m1_311 = _data_pdfs + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 11*_stride_pdfs_3;
      double * _data_pdfs_tmp_20_312 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 12*_stride_pdfs_tmp_3;
      double * const _data_pdfs_2m1_312 = _data_pdfs + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 12*_stride_pdfs_3;
      double * _data_pdfs_tmp_20_313 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 13*_stride_pdfs_tmp_3;
      double * const _data_pdfs_2m1_313 = _data_pdfs + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 13*_stride_pdfs_3;
      double * _data_pdfs_tmp_20_314 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 14*_stride_pdfs_tmp_3;
      double * const _data_pdfs_2m1_314 = _data_pdfs + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 14*_stride_pdfs_3;
      double * _data_pdfs_tmp_20_315 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 15*_stride_pdfs_tmp_3;
      double * const _data_pdfs_21_315 = _data_pdfs + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 15*_stride_pdfs_3;
      double * _data_pdfs_tmp_20_316 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 16*_stride_pdfs_tmp_3;
      double * const _data_pdfs_21_316 = _data_pdfs + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 16*_stride_pdfs_3;
      double * _data_pdfs_tmp_20_317 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 17*_stride_pdfs_tmp_3;
      double * const _data_pdfs_21_317 = _data_pdfs + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 17*_stride_pdfs_3;
      double * _data_pdfs_tmp_20_318 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 18*_stride_pdfs_tmp_3;
      double * const _data_pdfs_21_318 = _data_pdfs + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 18*_stride_pdfs_3;
      for (int ctr_1 = 1; ctr_1 < _size_pdfs_1 - 1; ctr_1 += 1)
      {
         double * _data_pdfs_tmp_20_30_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_30;
         double * const _data_pdfs_20_30_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_30;
         double * _data_pdfs_tmp_20_31_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_31;
         double * const _data_pdfs_20_31_1m1 = _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _data_pdfs_20_31;
         double * _data_pdfs_tmp_20_32_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_32;
         double * const _data_pdfs_20_32_11 = _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _data_pdfs_20_32;
         double * _data_pdfs_tmp_20_33_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_33;
         double * const _data_pdfs_20_33_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_33;
         double * _data_pdfs_tmp_20_34_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_34;
         double * const _data_pdfs_20_34_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_34;
         double * _data_pdfs_tmp_20_35_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_35;
         double * const _data_pdfs_2m1_35_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_2m1_35;
         double * _data_pdfs_tmp_20_36_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_36;
         double * const _data_pdfs_21_36_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_21_36;
         double * _data_pdfs_tmp_20_37_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_37;
         double * const _data_pdfs_20_37_1m1 = _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _data_pdfs_20_37;
         double * _data_pdfs_tmp_20_38_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_38;
         double * const _data_pdfs_20_38_1m1 = _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _data_pdfs_20_38;
         double * _data_pdfs_tmp_20_39_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_39;
         double * const _data_pdfs_20_39_11 = _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _data_pdfs_20_39;
         double * _data_pdfs_tmp_20_310_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_310;
         double * const _data_pdfs_20_310_11 = _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _data_pdfs_20_310;
         double * _data_pdfs_tmp_20_311_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_311;
         double * const _data_pdfs_2m1_311_1m1 = _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _data_pdfs_2m1_311;
         double * _data_pdfs_tmp_20_312_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_312;
         double * const _data_pdfs_2m1_312_11 = _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _data_pdfs_2m1_312;
         double * _data_pdfs_tmp_20_313_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_313;
         double * const _data_pdfs_2m1_313_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_2m1_313;
         double * _data_pdfs_tmp_20_314_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_314;
         double * const _data_pdfs_2m1_314_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_2m1_314;
         double * _data_pdfs_tmp_20_315_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_315;
         double * const _data_pdfs_21_315_1m1 = _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _data_pdfs_21_315;
         double * _data_pdfs_tmp_20_316_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_316;
         double * const _data_pdfs_21_316_11 = _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _data_pdfs_21_316;
         double * _data_pdfs_tmp_20_317_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_317;
         double * const _data_pdfs_21_317_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_21_317;
         double * _data_pdfs_tmp_20_318_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_318;
         double * const _data_pdfs_21_318_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_21_318;
         for (int ctr_0 = 1; ctr_0 < _size_pdfs_0 - 1; ctr_0 += 1)
         {
            _data_pdfs_tmp_20_30_10[_stride_pdfs_tmp_0*ctr_0] = _data_pdfs_20_30_10[_stride_pdfs_0*ctr_0];
            _data_pdfs_tmp_20_31_10[_stride_pdfs_tmp_0*ctr_0] = _data_pdfs_20_31_1m1[_stride_pdfs_0*ctr_0];
            _data_pdfs_tmp_20_32_10[_stride_pdfs_tmp_0*ctr_0] = _data_pdfs_20_32_11[_stride_pdfs_0*ctr_0];
            _data_pdfs_tmp_20_33_10[_stride_pdfs_tmp_0*ctr_0] = _data_pdfs_20_33_10[_stride_pdfs_0*ctr_0 + _stride_pdfs_0];
            _data_pdfs_tmp_20_34_10[_stride_pdfs_tmp_0*ctr_0] = _data_pdfs_20_34_10[_stride_pdfs_0*ctr_0 - _stride_pdfs_0];
            _data_pdfs_tmp_20_35_10[_stride_pdfs_tmp_0*ctr_0] = _data_pdfs_2m1_35_10[_stride_pdfs_0*ctr_0];
            _data_pdfs_tmp_20_36_10[_stride_pdfs_tmp_0*ctr_0] = _data_pdfs_21_36_10[_stride_pdfs_0*ctr_0];
            _data_pdfs_tmp_20_37_10[_stride_pdfs_tmp_0*ctr_0] = _data_pdfs_20_37_1m1[_stride_pdfs_0*ctr_0 + _stride_pdfs_0];
            _data_pdfs_tmp_20_38_10[_stride_pdfs_tmp_0*ctr_0] = _data_pdfs_20_38_1m1[_stride_pdfs_0*ctr_0 - _stride_pdfs_0];
            _data_pdfs_tmp_20_39_10[_stride_pdfs_tmp_0*ctr_0] = _data_pdfs_20_39_11[_stride_pdfs_0*ctr_0 + _stride_pdfs_0];
            _data_pdfs_tmp_20_310_10[_stride_pdfs_tmp_0*ctr_0] = _data_pdfs_20_310_11[_stride_pdfs_0*ctr_0 - _stride_pdfs_0];
            _data_pdfs_tmp_20_311_10[_stride_pdfs_tmp_0*ctr_0] = _data_pdfs_2m1_311_1m1[_stride_pdfs_0*ctr_0];
            _data_pdfs_tmp_20_312_10[_stride_pdfs_tmp_0*ctr_0] = _data_pdfs_2m1_312_11[_stride_pdfs_0*ctr_0];
            _data_pdfs_tmp_20_313_10[_stride_pdfs_tmp_0*ctr_0] = _data_pdfs_2m1_313_10[_stride_pdfs_0*ctr_0 + _stride_pdfs_0];
            _data_pdfs_tmp_20_314_10[_stride_pdfs_tmp_0*ctr_0] = _data_pdfs_2m1_314_10[_stride_pdfs_0*ctr_0 - _stride_pdfs_0];
            _data_pdfs_tmp_20_315_10[_stride_pdfs_tmp_0*ctr_0] = _data_pdfs_21_315_1m1[_stride_pdfs_0*ctr_0];
            _data_pdfs_tmp_20_316_10[_stride_pdfs_tmp_0*ctr_0] = _data_pdfs_21_316_11[_stride_pdfs_0*ctr_0];
            _data_pdfs_tmp_20_317_10[_stride_pdfs_tmp_0*ctr_0] = _data_pdfs_21_317_10[_stride_pdfs_0*ctr_0 + _stride_pdfs_0];
            _data_pdfs_tmp_20_318_10[_stride_pdfs_tmp_0*ctr_0] = _data_pdfs_21_318_10[_stride_pdfs_0*ctr_0 - _stride_pdfs_0];
         }
      }
   }
}
}


const real_t UniformGridGPU_LatticeModel::w[19] = { 0.333333333333333,0.0555555555555556,0.0555555555555556,0.0555555555555556,0.0555555555555556,0.0555555555555556,0.0555555555555556,0.0277777777777778,0.0277777777777778,0.0277777777777778,0.0277777777777778,0.0277777777777778,0.0277777777777778,0.0277777777777778,0.0277777777777778,0.0277777777777778,0.0277777777777778,0.0277777777777778,0.0277777777777778 };
const real_t UniformGridGPU_LatticeModel::wInv[19] = { 3.00000000000000,18.0000000000000,18.0000000000000,18.0000000000000,18.0000000000000,18.0000000000000,18.0000000000000,36.0000000000000,36.0000000000000,36.0000000000000,36.0000000000000,36.0000000000000,36.0000000000000,36.0000000000000,36.0000000000000,36.0000000000000,36.0000000000000,36.0000000000000,36.0000000000000 };

void UniformGridGPU_LatticeModel::Sweep::streamCollide( IBlock * block, const uint_t numberOfGhostLayersToInclude )
{
    auto pdfs = block->getData< GhostLayerField<double, 19> >(pdfsID);
    GhostLayerField<double, 19> * pdfs_tmp;
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


    auto & lm = dynamic_cast< lbm::PdfField<UniformGridGPU_LatticeModel> * > (pdfs)->latticeModel();
    lm.configureBlock(block);

    auto & omega = lm.omega;
    WALBERLA_ASSERT_GREATER_EQUAL(-cell_idx_c(numberOfGhostLayersToInclude) - 1, -int_c(pdfs->nrOfGhostLayers()));
    double * const _data_pdfs = pdfs->dataAt(-cell_idx_c(numberOfGhostLayersToInclude) - 1, -cell_idx_c(numberOfGhostLayersToInclude) - 1, -cell_idx_c(numberOfGhostLayersToInclude) - 1, 0);
    WALBERLA_ASSERT_GREATER_EQUAL(-cell_idx_c(numberOfGhostLayersToInclude) - 1, -int_c(pdfs_tmp->nrOfGhostLayers()));
    double * _data_pdfs_tmp = pdfs_tmp->dataAt(-cell_idx_c(numberOfGhostLayersToInclude) - 1, -cell_idx_c(numberOfGhostLayersToInclude) - 1, -cell_idx_c(numberOfGhostLayersToInclude) - 1, 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(), int64_t(pdfs->xSize() + 2*cell_idx_c(numberOfGhostLayersToInclude) + 2));
    const int64_t _size_pdfs_0 = int64_t(pdfs->xSize() + 2*cell_idx_c(numberOfGhostLayersToInclude) + 2);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(), int64_t(pdfs->ySize() + 2*cell_idx_c(numberOfGhostLayersToInclude) + 2));
    const int64_t _size_pdfs_1 = int64_t(pdfs->ySize() + 2*cell_idx_c(numberOfGhostLayersToInclude) + 2);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(), int64_t(pdfs->zSize() + 2*cell_idx_c(numberOfGhostLayersToInclude) + 2));
    const int64_t _size_pdfs_2 = int64_t(pdfs->zSize() + 2*cell_idx_c(numberOfGhostLayersToInclude) + 2);
    const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
    const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
    const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
    const int64_t _stride_pdfs_3 = int64_t(pdfs->fStride());
    internal_kernel_streamCollide::kernel_streamCollide(_data_pdfs, _data_pdfs_tmp, _size_pdfs_0, _size_pdfs_1, _size_pdfs_2, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3, omega);
    pdfs->swapDataPointers(pdfs_tmp);

}

void UniformGridGPU_LatticeModel::Sweep::collide( IBlock * block, const uint_t numberOfGhostLayersToInclude )
{
   auto pdfs = block->getData< GhostLayerField<double, 19> >(pdfsID);


    auto & lm = dynamic_cast< lbm::PdfField<UniformGridGPU_LatticeModel> * > (pdfs)->latticeModel();
    lm.configureBlock(block);

    auto & omega = lm.omega;
    WALBERLA_ASSERT_GREATER_EQUAL(-cell_idx_c(numberOfGhostLayersToInclude) - 1, -int_c(pdfs->nrOfGhostLayers()));
    double * _data_pdfs = pdfs->dataAt(-cell_idx_c(numberOfGhostLayersToInclude) - 1, -cell_idx_c(numberOfGhostLayersToInclude) - 1, -cell_idx_c(numberOfGhostLayersToInclude) - 1, 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(), int64_t(pdfs->xSize() + 2*cell_idx_c(numberOfGhostLayersToInclude) + 2));
    const int64_t _size_pdfs_0 = int64_t(pdfs->xSize() + 2*cell_idx_c(numberOfGhostLayersToInclude) + 2);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(), int64_t(pdfs->ySize() + 2*cell_idx_c(numberOfGhostLayersToInclude) + 2));
    const int64_t _size_pdfs_1 = int64_t(pdfs->ySize() + 2*cell_idx_c(numberOfGhostLayersToInclude) + 2);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(), int64_t(pdfs->zSize() + 2*cell_idx_c(numberOfGhostLayersToInclude) + 2));
    const int64_t _size_pdfs_2 = int64_t(pdfs->zSize() + 2*cell_idx_c(numberOfGhostLayersToInclude) + 2);
    const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
    const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
    const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
    const int64_t _stride_pdfs_3 = int64_t(pdfs->fStride());
    internal_kernel_collide::kernel_collide(_data_pdfs, _size_pdfs_0, _size_pdfs_1, _size_pdfs_2, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3, omega);
}


void UniformGridGPU_LatticeModel::Sweep::stream( IBlock * block, const uint_t numberOfGhostLayersToInclude )
{
    auto pdfs = block->getData< GhostLayerField<double, 19> >(pdfsID);
    GhostLayerField<double, 19> * pdfs_tmp;
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


    WALBERLA_ASSERT_GREATER_EQUAL(-cell_idx_c(numberOfGhostLayersToInclude) - 1, -int_c(pdfs->nrOfGhostLayers()));
    double * const _data_pdfs = pdfs->dataAt(-cell_idx_c(numberOfGhostLayersToInclude) - 1, -cell_idx_c(numberOfGhostLayersToInclude) - 1, -cell_idx_c(numberOfGhostLayersToInclude) - 1, 0);
    WALBERLA_ASSERT_GREATER_EQUAL(-cell_idx_c(numberOfGhostLayersToInclude) - 1, -int_c(pdfs_tmp->nrOfGhostLayers()));
    double * _data_pdfs_tmp = pdfs_tmp->dataAt(-cell_idx_c(numberOfGhostLayersToInclude) - 1, -cell_idx_c(numberOfGhostLayersToInclude) - 1, -cell_idx_c(numberOfGhostLayersToInclude) - 1, 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(), int64_t(pdfs->xSize() + 2*cell_idx_c(numberOfGhostLayersToInclude) + 2));
    const int64_t _size_pdfs_0 = int64_t(pdfs->xSize() + 2*cell_idx_c(numberOfGhostLayersToInclude) + 2);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(), int64_t(pdfs->ySize() + 2*cell_idx_c(numberOfGhostLayersToInclude) + 2));
    const int64_t _size_pdfs_1 = int64_t(pdfs->ySize() + 2*cell_idx_c(numberOfGhostLayersToInclude) + 2);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(), int64_t(pdfs->zSize() + 2*cell_idx_c(numberOfGhostLayersToInclude) + 2));
    const int64_t _size_pdfs_2 = int64_t(pdfs->zSize() + 2*cell_idx_c(numberOfGhostLayersToInclude) + 2);
    const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
    const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
    const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
    const int64_t _stride_pdfs_3 = int64_t(pdfs->fStride());
    const int64_t _stride_pdfs_tmp_0 = int64_t(pdfs_tmp->xStride());
    const int64_t _stride_pdfs_tmp_1 = int64_t(pdfs_tmp->yStride());
    const int64_t _stride_pdfs_tmp_2 = int64_t(pdfs_tmp->zStride());
    const int64_t _stride_pdfs_tmp_3 = int64_t(pdfs_tmp->fStride());
    internal_kernel_stream::kernel_stream(_data_pdfs, _data_pdfs_tmp, _size_pdfs_0, _size_pdfs_1, _size_pdfs_2, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3, _stride_pdfs_tmp_0, _stride_pdfs_tmp_1, _stride_pdfs_tmp_2, _stride_pdfs_tmp_3);
    
    pdfs->swapDataPointers(pdfs_tmp);

}


} // namespace lbm
} // namespace walberla




// Buffer Packing

namespace walberla {
namespace mpi {

mpi::SendBuffer & operator<< (mpi::SendBuffer & buf, const ::walberla::lbm::UniformGridGPU_LatticeModel & lm)
{
    buf << lm.currentLevel;
    return buf;
}

mpi::RecvBuffer & operator>> (mpi::RecvBuffer & buf, ::walberla::lbm::UniformGridGPU_LatticeModel & lm)
{
    buf >> lm.currentLevel;
    return buf;
}


} // namespace mpi
} // namespace walberla
