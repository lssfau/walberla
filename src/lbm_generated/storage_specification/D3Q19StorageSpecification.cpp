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
//! \\file D3Q19StorageSpecification.cpp
//! \\author lbmpy
//======================================================================================================================

#include "D3Q19StorageSpecification.h"

#if ( defined WALBERLA_CXX_COMPILER_IS_GNU ) || ( defined WALBERLA_CXX_COMPILER_IS_CLANG )
#   pragma GCC diagnostic push
#   pragma GCC diagnostic ignored "-Wfloat-equal"
#   pragma GCC diagnostic ignored "-Wshadow"
#   pragma GCC diagnostic ignored "-Wconversion"
#   pragma GCC diagnostic ignored "-Wunused-variable"
#endif

/*************************************************************************************
 *                                Kernel Definitions
*************************************************************************************/
namespace internal_d3q19storagespecification_pack_ALL {
static FUNC_PREFIX void d3q19storagespecification_pack_ALL(double * RESTRICT  _data_buffer, double * RESTRICT const _data_pdfs_src, int64_t const _size_pdfs_src_0, int64_t const _size_pdfs_src_1, int64_t const _size_pdfs_src_2, int64_t const _stride_pdfs_src_0, int64_t const _stride_pdfs_src_1, int64_t const _stride_pdfs_src_2, int64_t const _stride_pdfs_src_3)
{
   for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_src_0; ctr_0 += 1)
   {
      double * RESTRICT _data_pdfs_src_00_30 = _data_pdfs_src + _stride_pdfs_src_0*ctr_0;
      double * RESTRICT _data_pdfs_src_00_31 = _data_pdfs_src + _stride_pdfs_src_0*ctr_0 + _stride_pdfs_src_3;
      double * RESTRICT _data_pdfs_src_00_32 = _data_pdfs_src + _stride_pdfs_src_0*ctr_0 + 2*_stride_pdfs_src_3;
      double * RESTRICT _data_pdfs_src_00_33 = _data_pdfs_src + _stride_pdfs_src_0*ctr_0 + 3*_stride_pdfs_src_3;
      double * RESTRICT _data_pdfs_src_00_34 = _data_pdfs_src + _stride_pdfs_src_0*ctr_0 + 4*_stride_pdfs_src_3;
      double * RESTRICT _data_pdfs_src_00_35 = _data_pdfs_src + _stride_pdfs_src_0*ctr_0 + 5*_stride_pdfs_src_3;
      double * RESTRICT _data_pdfs_src_00_36 = _data_pdfs_src + _stride_pdfs_src_0*ctr_0 + 6*_stride_pdfs_src_3;
      double * RESTRICT _data_pdfs_src_00_37 = _data_pdfs_src + _stride_pdfs_src_0*ctr_0 + 7*_stride_pdfs_src_3;
      double * RESTRICT _data_pdfs_src_00_38 = _data_pdfs_src + _stride_pdfs_src_0*ctr_0 + 8*_stride_pdfs_src_3;
      double * RESTRICT _data_pdfs_src_00_39 = _data_pdfs_src + _stride_pdfs_src_0*ctr_0 + 9*_stride_pdfs_src_3;
      double * RESTRICT _data_pdfs_src_00_310 = _data_pdfs_src + _stride_pdfs_src_0*ctr_0 + 10*_stride_pdfs_src_3;
      double * RESTRICT _data_pdfs_src_00_311 = _data_pdfs_src + _stride_pdfs_src_0*ctr_0 + 11*_stride_pdfs_src_3;
      double * RESTRICT _data_pdfs_src_00_312 = _data_pdfs_src + _stride_pdfs_src_0*ctr_0 + 12*_stride_pdfs_src_3;
      double * RESTRICT _data_pdfs_src_00_313 = _data_pdfs_src + _stride_pdfs_src_0*ctr_0 + 13*_stride_pdfs_src_3;
      double * RESTRICT _data_pdfs_src_00_314 = _data_pdfs_src + _stride_pdfs_src_0*ctr_0 + 14*_stride_pdfs_src_3;
      double * RESTRICT _data_pdfs_src_00_315 = _data_pdfs_src + _stride_pdfs_src_0*ctr_0 + 15*_stride_pdfs_src_3;
      double * RESTRICT _data_pdfs_src_00_316 = _data_pdfs_src + _stride_pdfs_src_0*ctr_0 + 16*_stride_pdfs_src_3;
      double * RESTRICT _data_pdfs_src_00_317 = _data_pdfs_src + _stride_pdfs_src_0*ctr_0 + 17*_stride_pdfs_src_3;
      double * RESTRICT _data_pdfs_src_00_318 = _data_pdfs_src + _stride_pdfs_src_0*ctr_0 + 18*_stride_pdfs_src_3;
      for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_src_1; ctr_1 += 1)
      {
         double * RESTRICT _data_pdfs_src_00_30_10 = _stride_pdfs_src_1*ctr_1 + _data_pdfs_src_00_30;
         double * RESTRICT _data_pdfs_src_00_31_10 = _stride_pdfs_src_1*ctr_1 + _data_pdfs_src_00_31;
         double * RESTRICT _data_pdfs_src_00_32_10 = _stride_pdfs_src_1*ctr_1 + _data_pdfs_src_00_32;
         double * RESTRICT _data_pdfs_src_00_33_10 = _stride_pdfs_src_1*ctr_1 + _data_pdfs_src_00_33;
         double * RESTRICT _data_pdfs_src_00_34_10 = _stride_pdfs_src_1*ctr_1 + _data_pdfs_src_00_34;
         double * RESTRICT _data_pdfs_src_00_35_10 = _stride_pdfs_src_1*ctr_1 + _data_pdfs_src_00_35;
         double * RESTRICT _data_pdfs_src_00_36_10 = _stride_pdfs_src_1*ctr_1 + _data_pdfs_src_00_36;
         double * RESTRICT _data_pdfs_src_00_37_10 = _stride_pdfs_src_1*ctr_1 + _data_pdfs_src_00_37;
         double * RESTRICT _data_pdfs_src_00_38_10 = _stride_pdfs_src_1*ctr_1 + _data_pdfs_src_00_38;
         double * RESTRICT _data_pdfs_src_00_39_10 = _stride_pdfs_src_1*ctr_1 + _data_pdfs_src_00_39;
         double * RESTRICT _data_pdfs_src_00_310_10 = _stride_pdfs_src_1*ctr_1 + _data_pdfs_src_00_310;
         double * RESTRICT _data_pdfs_src_00_311_10 = _stride_pdfs_src_1*ctr_1 + _data_pdfs_src_00_311;
         double * RESTRICT _data_pdfs_src_00_312_10 = _stride_pdfs_src_1*ctr_1 + _data_pdfs_src_00_312;
         double * RESTRICT _data_pdfs_src_00_313_10 = _stride_pdfs_src_1*ctr_1 + _data_pdfs_src_00_313;
         double * RESTRICT _data_pdfs_src_00_314_10 = _stride_pdfs_src_1*ctr_1 + _data_pdfs_src_00_314;
         double * RESTRICT _data_pdfs_src_00_315_10 = _stride_pdfs_src_1*ctr_1 + _data_pdfs_src_00_315;
         double * RESTRICT _data_pdfs_src_00_316_10 = _stride_pdfs_src_1*ctr_1 + _data_pdfs_src_00_316;
         double * RESTRICT _data_pdfs_src_00_317_10 = _stride_pdfs_src_1*ctr_1 + _data_pdfs_src_00_317;
         double * RESTRICT _data_pdfs_src_00_318_10 = _stride_pdfs_src_1*ctr_1 + _data_pdfs_src_00_318;
         for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_src_2; ctr_2 += 1)
         {
            _data_buffer[19*_size_pdfs_src_1*_size_pdfs_src_2*ctr_0 + 19*_size_pdfs_src_2*ctr_1 + 19*ctr_2] = _data_pdfs_src_00_30_10[_stride_pdfs_src_2*ctr_2];
            _data_buffer[19*_size_pdfs_src_1*_size_pdfs_src_2*ctr_0 + 19*_size_pdfs_src_2*ctr_1 + 19*ctr_2 + 1] = _data_pdfs_src_00_31_10[_stride_pdfs_src_2*ctr_2];
            _data_buffer[19*_size_pdfs_src_1*_size_pdfs_src_2*ctr_0 + 19*_size_pdfs_src_2*ctr_1 + 19*ctr_2 + 2] = _data_pdfs_src_00_32_10[_stride_pdfs_src_2*ctr_2];
            _data_buffer[19*_size_pdfs_src_1*_size_pdfs_src_2*ctr_0 + 19*_size_pdfs_src_2*ctr_1 + 19*ctr_2 + 3] = _data_pdfs_src_00_33_10[_stride_pdfs_src_2*ctr_2];
            _data_buffer[19*_size_pdfs_src_1*_size_pdfs_src_2*ctr_0 + 19*_size_pdfs_src_2*ctr_1 + 19*ctr_2 + 4] = _data_pdfs_src_00_34_10[_stride_pdfs_src_2*ctr_2];
            _data_buffer[19*_size_pdfs_src_1*_size_pdfs_src_2*ctr_0 + 19*_size_pdfs_src_2*ctr_1 + 19*ctr_2 + 5] = _data_pdfs_src_00_35_10[_stride_pdfs_src_2*ctr_2];
            _data_buffer[19*_size_pdfs_src_1*_size_pdfs_src_2*ctr_0 + 19*_size_pdfs_src_2*ctr_1 + 19*ctr_2 + 6] = _data_pdfs_src_00_36_10[_stride_pdfs_src_2*ctr_2];
            _data_buffer[19*_size_pdfs_src_1*_size_pdfs_src_2*ctr_0 + 19*_size_pdfs_src_2*ctr_1 + 19*ctr_2 + 7] = _data_pdfs_src_00_37_10[_stride_pdfs_src_2*ctr_2];
            _data_buffer[19*_size_pdfs_src_1*_size_pdfs_src_2*ctr_0 + 19*_size_pdfs_src_2*ctr_1 + 19*ctr_2 + 8] = _data_pdfs_src_00_38_10[_stride_pdfs_src_2*ctr_2];
            _data_buffer[19*_size_pdfs_src_1*_size_pdfs_src_2*ctr_0 + 19*_size_pdfs_src_2*ctr_1 + 19*ctr_2 + 9] = _data_pdfs_src_00_39_10[_stride_pdfs_src_2*ctr_2];
            _data_buffer[19*_size_pdfs_src_1*_size_pdfs_src_2*ctr_0 + 19*_size_pdfs_src_2*ctr_1 + 19*ctr_2 + 10] = _data_pdfs_src_00_310_10[_stride_pdfs_src_2*ctr_2];
            _data_buffer[19*_size_pdfs_src_1*_size_pdfs_src_2*ctr_0 + 19*_size_pdfs_src_2*ctr_1 + 19*ctr_2 + 11] = _data_pdfs_src_00_311_10[_stride_pdfs_src_2*ctr_2];
            _data_buffer[19*_size_pdfs_src_1*_size_pdfs_src_2*ctr_0 + 19*_size_pdfs_src_2*ctr_1 + 19*ctr_2 + 12] = _data_pdfs_src_00_312_10[_stride_pdfs_src_2*ctr_2];
            _data_buffer[19*_size_pdfs_src_1*_size_pdfs_src_2*ctr_0 + 19*_size_pdfs_src_2*ctr_1 + 19*ctr_2 + 13] = _data_pdfs_src_00_313_10[_stride_pdfs_src_2*ctr_2];
            _data_buffer[19*_size_pdfs_src_1*_size_pdfs_src_2*ctr_0 + 19*_size_pdfs_src_2*ctr_1 + 19*ctr_2 + 14] = _data_pdfs_src_00_314_10[_stride_pdfs_src_2*ctr_2];
            _data_buffer[19*_size_pdfs_src_1*_size_pdfs_src_2*ctr_0 + 19*_size_pdfs_src_2*ctr_1 + 19*ctr_2 + 15] = _data_pdfs_src_00_315_10[_stride_pdfs_src_2*ctr_2];
            _data_buffer[19*_size_pdfs_src_1*_size_pdfs_src_2*ctr_0 + 19*_size_pdfs_src_2*ctr_1 + 19*ctr_2 + 16] = _data_pdfs_src_00_316_10[_stride_pdfs_src_2*ctr_2];
            _data_buffer[19*_size_pdfs_src_1*_size_pdfs_src_2*ctr_0 + 19*_size_pdfs_src_2*ctr_1 + 19*ctr_2 + 17] = _data_pdfs_src_00_317_10[_stride_pdfs_src_2*ctr_2];
            _data_buffer[19*_size_pdfs_src_1*_size_pdfs_src_2*ctr_0 + 19*_size_pdfs_src_2*ctr_1 + 19*ctr_2 + 18] = _data_pdfs_src_00_318_10[_stride_pdfs_src_2*ctr_2];
         }
      }
   }
}
}

namespace internal_d3q19storagespecification_unpack_ALL {
static FUNC_PREFIX void d3q19storagespecification_unpack_ALL(const double * RESTRICT const _data_buffer, double * RESTRICT  _data_pdfs_dst, int64_t const _size_pdfs_dst_0, int64_t const _size_pdfs_dst_1, int64_t const _size_pdfs_dst_2, int64_t const _stride_pdfs_dst_0, int64_t const _stride_pdfs_dst_1, int64_t const _stride_pdfs_dst_2, int64_t const _stride_pdfs_dst_3)
{
   for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_dst_0; ctr_0 += 1)
   {
      double * RESTRICT  _data_pdfs_dst_00_30 = _data_pdfs_dst + _stride_pdfs_dst_0*ctr_0;
      double * RESTRICT  _data_pdfs_dst_00_31 = _data_pdfs_dst + _stride_pdfs_dst_0*ctr_0 + _stride_pdfs_dst_3;
      double * RESTRICT  _data_pdfs_dst_00_32 = _data_pdfs_dst + _stride_pdfs_dst_0*ctr_0 + 2*_stride_pdfs_dst_3;
      double * RESTRICT  _data_pdfs_dst_00_33 = _data_pdfs_dst + _stride_pdfs_dst_0*ctr_0 + 3*_stride_pdfs_dst_3;
      double * RESTRICT  _data_pdfs_dst_00_34 = _data_pdfs_dst + _stride_pdfs_dst_0*ctr_0 + 4*_stride_pdfs_dst_3;
      double * RESTRICT  _data_pdfs_dst_00_35 = _data_pdfs_dst + _stride_pdfs_dst_0*ctr_0 + 5*_stride_pdfs_dst_3;
      double * RESTRICT  _data_pdfs_dst_00_36 = _data_pdfs_dst + _stride_pdfs_dst_0*ctr_0 + 6*_stride_pdfs_dst_3;
      double * RESTRICT  _data_pdfs_dst_00_37 = _data_pdfs_dst + _stride_pdfs_dst_0*ctr_0 + 7*_stride_pdfs_dst_3;
      double * RESTRICT  _data_pdfs_dst_00_38 = _data_pdfs_dst + _stride_pdfs_dst_0*ctr_0 + 8*_stride_pdfs_dst_3;
      double * RESTRICT  _data_pdfs_dst_00_39 = _data_pdfs_dst + _stride_pdfs_dst_0*ctr_0 + 9*_stride_pdfs_dst_3;
      double * RESTRICT  _data_pdfs_dst_00_310 = _data_pdfs_dst + _stride_pdfs_dst_0*ctr_0 + 10*_stride_pdfs_dst_3;
      double * RESTRICT  _data_pdfs_dst_00_311 = _data_pdfs_dst + _stride_pdfs_dst_0*ctr_0 + 11*_stride_pdfs_dst_3;
      double * RESTRICT  _data_pdfs_dst_00_312 = _data_pdfs_dst + _stride_pdfs_dst_0*ctr_0 + 12*_stride_pdfs_dst_3;
      double * RESTRICT  _data_pdfs_dst_00_313 = _data_pdfs_dst + _stride_pdfs_dst_0*ctr_0 + 13*_stride_pdfs_dst_3;
      double * RESTRICT  _data_pdfs_dst_00_314 = _data_pdfs_dst + _stride_pdfs_dst_0*ctr_0 + 14*_stride_pdfs_dst_3;
      double * RESTRICT  _data_pdfs_dst_00_315 = _data_pdfs_dst + _stride_pdfs_dst_0*ctr_0 + 15*_stride_pdfs_dst_3;
      double * RESTRICT  _data_pdfs_dst_00_316 = _data_pdfs_dst + _stride_pdfs_dst_0*ctr_0 + 16*_stride_pdfs_dst_3;
      double * RESTRICT  _data_pdfs_dst_00_317 = _data_pdfs_dst + _stride_pdfs_dst_0*ctr_0 + 17*_stride_pdfs_dst_3;
      double * RESTRICT  _data_pdfs_dst_00_318 = _data_pdfs_dst + _stride_pdfs_dst_0*ctr_0 + 18*_stride_pdfs_dst_3;
      for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_dst_1; ctr_1 += 1)
      {
         double * RESTRICT  _data_pdfs_dst_00_30_10 = _stride_pdfs_dst_1*ctr_1 + _data_pdfs_dst_00_30;
         double * RESTRICT  _data_pdfs_dst_00_31_10 = _stride_pdfs_dst_1*ctr_1 + _data_pdfs_dst_00_31;
         double * RESTRICT  _data_pdfs_dst_00_32_10 = _stride_pdfs_dst_1*ctr_1 + _data_pdfs_dst_00_32;
         double * RESTRICT  _data_pdfs_dst_00_33_10 = _stride_pdfs_dst_1*ctr_1 + _data_pdfs_dst_00_33;
         double * RESTRICT  _data_pdfs_dst_00_34_10 = _stride_pdfs_dst_1*ctr_1 + _data_pdfs_dst_00_34;
         double * RESTRICT  _data_pdfs_dst_00_35_10 = _stride_pdfs_dst_1*ctr_1 + _data_pdfs_dst_00_35;
         double * RESTRICT  _data_pdfs_dst_00_36_10 = _stride_pdfs_dst_1*ctr_1 + _data_pdfs_dst_00_36;
         double * RESTRICT  _data_pdfs_dst_00_37_10 = _stride_pdfs_dst_1*ctr_1 + _data_pdfs_dst_00_37;
         double * RESTRICT  _data_pdfs_dst_00_38_10 = _stride_pdfs_dst_1*ctr_1 + _data_pdfs_dst_00_38;
         double * RESTRICT  _data_pdfs_dst_00_39_10 = _stride_pdfs_dst_1*ctr_1 + _data_pdfs_dst_00_39;
         double * RESTRICT  _data_pdfs_dst_00_310_10 = _stride_pdfs_dst_1*ctr_1 + _data_pdfs_dst_00_310;
         double * RESTRICT  _data_pdfs_dst_00_311_10 = _stride_pdfs_dst_1*ctr_1 + _data_pdfs_dst_00_311;
         double * RESTRICT  _data_pdfs_dst_00_312_10 = _stride_pdfs_dst_1*ctr_1 + _data_pdfs_dst_00_312;
         double * RESTRICT  _data_pdfs_dst_00_313_10 = _stride_pdfs_dst_1*ctr_1 + _data_pdfs_dst_00_313;
         double * RESTRICT  _data_pdfs_dst_00_314_10 = _stride_pdfs_dst_1*ctr_1 + _data_pdfs_dst_00_314;
         double * RESTRICT  _data_pdfs_dst_00_315_10 = _stride_pdfs_dst_1*ctr_1 + _data_pdfs_dst_00_315;
         double * RESTRICT  _data_pdfs_dst_00_316_10 = _stride_pdfs_dst_1*ctr_1 + _data_pdfs_dst_00_316;
         double * RESTRICT  _data_pdfs_dst_00_317_10 = _stride_pdfs_dst_1*ctr_1 + _data_pdfs_dst_00_317;
         double * RESTRICT  _data_pdfs_dst_00_318_10 = _stride_pdfs_dst_1*ctr_1 + _data_pdfs_dst_00_318;
         for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_dst_2; ctr_2 += 1)
         {
            _data_pdfs_dst_00_30_10[_stride_pdfs_dst_2*ctr_2] = _data_buffer[19*_size_pdfs_dst_1*_size_pdfs_dst_2*ctr_0 + 19*_size_pdfs_dst_2*ctr_1 + 19*ctr_2];
            _data_pdfs_dst_00_31_10[_stride_pdfs_dst_2*ctr_2] = _data_buffer[19*_size_pdfs_dst_1*_size_pdfs_dst_2*ctr_0 + 19*_size_pdfs_dst_2*ctr_1 + 19*ctr_2 + 1];
            _data_pdfs_dst_00_32_10[_stride_pdfs_dst_2*ctr_2] = _data_buffer[19*_size_pdfs_dst_1*_size_pdfs_dst_2*ctr_0 + 19*_size_pdfs_dst_2*ctr_1 + 19*ctr_2 + 2];
            _data_pdfs_dst_00_33_10[_stride_pdfs_dst_2*ctr_2] = _data_buffer[19*_size_pdfs_dst_1*_size_pdfs_dst_2*ctr_0 + 19*_size_pdfs_dst_2*ctr_1 + 19*ctr_2 + 3];
            _data_pdfs_dst_00_34_10[_stride_pdfs_dst_2*ctr_2] = _data_buffer[19*_size_pdfs_dst_1*_size_pdfs_dst_2*ctr_0 + 19*_size_pdfs_dst_2*ctr_1 + 19*ctr_2 + 4];
            _data_pdfs_dst_00_35_10[_stride_pdfs_dst_2*ctr_2] = _data_buffer[19*_size_pdfs_dst_1*_size_pdfs_dst_2*ctr_0 + 19*_size_pdfs_dst_2*ctr_1 + 19*ctr_2 + 5];
            _data_pdfs_dst_00_36_10[_stride_pdfs_dst_2*ctr_2] = _data_buffer[19*_size_pdfs_dst_1*_size_pdfs_dst_2*ctr_0 + 19*_size_pdfs_dst_2*ctr_1 + 19*ctr_2 + 6];
            _data_pdfs_dst_00_37_10[_stride_pdfs_dst_2*ctr_2] = _data_buffer[19*_size_pdfs_dst_1*_size_pdfs_dst_2*ctr_0 + 19*_size_pdfs_dst_2*ctr_1 + 19*ctr_2 + 7];
            _data_pdfs_dst_00_38_10[_stride_pdfs_dst_2*ctr_2] = _data_buffer[19*_size_pdfs_dst_1*_size_pdfs_dst_2*ctr_0 + 19*_size_pdfs_dst_2*ctr_1 + 19*ctr_2 + 8];
            _data_pdfs_dst_00_39_10[_stride_pdfs_dst_2*ctr_2] = _data_buffer[19*_size_pdfs_dst_1*_size_pdfs_dst_2*ctr_0 + 19*_size_pdfs_dst_2*ctr_1 + 19*ctr_2 + 9];
            _data_pdfs_dst_00_310_10[_stride_pdfs_dst_2*ctr_2] = _data_buffer[19*_size_pdfs_dst_1*_size_pdfs_dst_2*ctr_0 + 19*_size_pdfs_dst_2*ctr_1 + 19*ctr_2 + 10];
            _data_pdfs_dst_00_311_10[_stride_pdfs_dst_2*ctr_2] = _data_buffer[19*_size_pdfs_dst_1*_size_pdfs_dst_2*ctr_0 + 19*_size_pdfs_dst_2*ctr_1 + 19*ctr_2 + 11];
            _data_pdfs_dst_00_312_10[_stride_pdfs_dst_2*ctr_2] = _data_buffer[19*_size_pdfs_dst_1*_size_pdfs_dst_2*ctr_0 + 19*_size_pdfs_dst_2*ctr_1 + 19*ctr_2 + 12];
            _data_pdfs_dst_00_313_10[_stride_pdfs_dst_2*ctr_2] = _data_buffer[19*_size_pdfs_dst_1*_size_pdfs_dst_2*ctr_0 + 19*_size_pdfs_dst_2*ctr_1 + 19*ctr_2 + 13];
            _data_pdfs_dst_00_314_10[_stride_pdfs_dst_2*ctr_2] = _data_buffer[19*_size_pdfs_dst_1*_size_pdfs_dst_2*ctr_0 + 19*_size_pdfs_dst_2*ctr_1 + 19*ctr_2 + 14];
            _data_pdfs_dst_00_315_10[_stride_pdfs_dst_2*ctr_2] = _data_buffer[19*_size_pdfs_dst_1*_size_pdfs_dst_2*ctr_0 + 19*_size_pdfs_dst_2*ctr_1 + 19*ctr_2 + 15];
            _data_pdfs_dst_00_316_10[_stride_pdfs_dst_2*ctr_2] = _data_buffer[19*_size_pdfs_dst_1*_size_pdfs_dst_2*ctr_0 + 19*_size_pdfs_dst_2*ctr_1 + 19*ctr_2 + 16];
            _data_pdfs_dst_00_317_10[_stride_pdfs_dst_2*ctr_2] = _data_buffer[19*_size_pdfs_dst_1*_size_pdfs_dst_2*ctr_0 + 19*_size_pdfs_dst_2*ctr_1 + 19*ctr_2 + 17];
            _data_pdfs_dst_00_318_10[_stride_pdfs_dst_2*ctr_2] = _data_buffer[19*_size_pdfs_dst_1*_size_pdfs_dst_2*ctr_0 + 19*_size_pdfs_dst_2*ctr_1 + 19*ctr_2 + 18];
         }
      }
   }
}
}

namespace internal_d3q19storagespecification_localCopy_ALL {
static FUNC_PREFIX void d3q19storagespecification_localCopy_ALL(double * RESTRICT  _data_pdfs_dst, double * RESTRICT const _data_pdfs_src, int64_t const _size_pdfs_dst_0, int64_t const _size_pdfs_dst_1, int64_t const _size_pdfs_dst_2, int64_t const _stride_pdfs_dst_0, int64_t const _stride_pdfs_dst_1, int64_t const _stride_pdfs_dst_2, int64_t const _stride_pdfs_dst_3, int64_t const _stride_pdfs_src_0, int64_t const _stride_pdfs_src_1, int64_t const _stride_pdfs_src_2, int64_t const _stride_pdfs_src_3)
{
   for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_dst_0; ctr_0 += 1)
   {
      double * RESTRICT  _data_pdfs_dst_00_30 = _data_pdfs_dst + _stride_pdfs_dst_0*ctr_0;
      double * RESTRICT _data_pdfs_src_00_30 = _data_pdfs_src + _stride_pdfs_src_0*ctr_0;
      double * RESTRICT  _data_pdfs_dst_00_31 = _data_pdfs_dst + _stride_pdfs_dst_0*ctr_0 + _stride_pdfs_dst_3;
      double * RESTRICT _data_pdfs_src_00_31 = _data_pdfs_src + _stride_pdfs_src_0*ctr_0 + _stride_pdfs_src_3;
      double * RESTRICT  _data_pdfs_dst_00_32 = _data_pdfs_dst + _stride_pdfs_dst_0*ctr_0 + 2*_stride_pdfs_dst_3;
      double * RESTRICT _data_pdfs_src_00_32 = _data_pdfs_src + _stride_pdfs_src_0*ctr_0 + 2*_stride_pdfs_src_3;
      double * RESTRICT  _data_pdfs_dst_00_33 = _data_pdfs_dst + _stride_pdfs_dst_0*ctr_0 + 3*_stride_pdfs_dst_3;
      double * RESTRICT _data_pdfs_src_00_33 = _data_pdfs_src + _stride_pdfs_src_0*ctr_0 + 3*_stride_pdfs_src_3;
      double * RESTRICT  _data_pdfs_dst_00_34 = _data_pdfs_dst + _stride_pdfs_dst_0*ctr_0 + 4*_stride_pdfs_dst_3;
      double * RESTRICT _data_pdfs_src_00_34 = _data_pdfs_src + _stride_pdfs_src_0*ctr_0 + 4*_stride_pdfs_src_3;
      double * RESTRICT  _data_pdfs_dst_00_35 = _data_pdfs_dst + _stride_pdfs_dst_0*ctr_0 + 5*_stride_pdfs_dst_3;
      double * RESTRICT _data_pdfs_src_00_35 = _data_pdfs_src + _stride_pdfs_src_0*ctr_0 + 5*_stride_pdfs_src_3;
      double * RESTRICT  _data_pdfs_dst_00_36 = _data_pdfs_dst + _stride_pdfs_dst_0*ctr_0 + 6*_stride_pdfs_dst_3;
      double * RESTRICT _data_pdfs_src_00_36 = _data_pdfs_src + _stride_pdfs_src_0*ctr_0 + 6*_stride_pdfs_src_3;
      double * RESTRICT  _data_pdfs_dst_00_37 = _data_pdfs_dst + _stride_pdfs_dst_0*ctr_0 + 7*_stride_pdfs_dst_3;
      double * RESTRICT _data_pdfs_src_00_37 = _data_pdfs_src + _stride_pdfs_src_0*ctr_0 + 7*_stride_pdfs_src_3;
      double * RESTRICT  _data_pdfs_dst_00_38 = _data_pdfs_dst + _stride_pdfs_dst_0*ctr_0 + 8*_stride_pdfs_dst_3;
      double * RESTRICT _data_pdfs_src_00_38 = _data_pdfs_src + _stride_pdfs_src_0*ctr_0 + 8*_stride_pdfs_src_3;
      double * RESTRICT  _data_pdfs_dst_00_39 = _data_pdfs_dst + _stride_pdfs_dst_0*ctr_0 + 9*_stride_pdfs_dst_3;
      double * RESTRICT _data_pdfs_src_00_39 = _data_pdfs_src + _stride_pdfs_src_0*ctr_0 + 9*_stride_pdfs_src_3;
      double * RESTRICT  _data_pdfs_dst_00_310 = _data_pdfs_dst + _stride_pdfs_dst_0*ctr_0 + 10*_stride_pdfs_dst_3;
      double * RESTRICT _data_pdfs_src_00_310 = _data_pdfs_src + _stride_pdfs_src_0*ctr_0 + 10*_stride_pdfs_src_3;
      double * RESTRICT  _data_pdfs_dst_00_311 = _data_pdfs_dst + _stride_pdfs_dst_0*ctr_0 + 11*_stride_pdfs_dst_3;
      double * RESTRICT _data_pdfs_src_00_311 = _data_pdfs_src + _stride_pdfs_src_0*ctr_0 + 11*_stride_pdfs_src_3;
      double * RESTRICT  _data_pdfs_dst_00_312 = _data_pdfs_dst + _stride_pdfs_dst_0*ctr_0 + 12*_stride_pdfs_dst_3;
      double * RESTRICT _data_pdfs_src_00_312 = _data_pdfs_src + _stride_pdfs_src_0*ctr_0 + 12*_stride_pdfs_src_3;
      double * RESTRICT  _data_pdfs_dst_00_313 = _data_pdfs_dst + _stride_pdfs_dst_0*ctr_0 + 13*_stride_pdfs_dst_3;
      double * RESTRICT _data_pdfs_src_00_313 = _data_pdfs_src + _stride_pdfs_src_0*ctr_0 + 13*_stride_pdfs_src_3;
      double * RESTRICT  _data_pdfs_dst_00_314 = _data_pdfs_dst + _stride_pdfs_dst_0*ctr_0 + 14*_stride_pdfs_dst_3;
      double * RESTRICT _data_pdfs_src_00_314 = _data_pdfs_src + _stride_pdfs_src_0*ctr_0 + 14*_stride_pdfs_src_3;
      double * RESTRICT  _data_pdfs_dst_00_315 = _data_pdfs_dst + _stride_pdfs_dst_0*ctr_0 + 15*_stride_pdfs_dst_3;
      double * RESTRICT _data_pdfs_src_00_315 = _data_pdfs_src + _stride_pdfs_src_0*ctr_0 + 15*_stride_pdfs_src_3;
      double * RESTRICT  _data_pdfs_dst_00_316 = _data_pdfs_dst + _stride_pdfs_dst_0*ctr_0 + 16*_stride_pdfs_dst_3;
      double * RESTRICT _data_pdfs_src_00_316 = _data_pdfs_src + _stride_pdfs_src_0*ctr_0 + 16*_stride_pdfs_src_3;
      double * RESTRICT  _data_pdfs_dst_00_317 = _data_pdfs_dst + _stride_pdfs_dst_0*ctr_0 + 17*_stride_pdfs_dst_3;
      double * RESTRICT _data_pdfs_src_00_317 = _data_pdfs_src + _stride_pdfs_src_0*ctr_0 + 17*_stride_pdfs_src_3;
      double * RESTRICT  _data_pdfs_dst_00_318 = _data_pdfs_dst + _stride_pdfs_dst_0*ctr_0 + 18*_stride_pdfs_dst_3;
      double * RESTRICT _data_pdfs_src_00_318 = _data_pdfs_src + _stride_pdfs_src_0*ctr_0 + 18*_stride_pdfs_src_3;
      for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_dst_1; ctr_1 += 1)
      {
         double * RESTRICT  _data_pdfs_dst_00_30_10 = _stride_pdfs_dst_1*ctr_1 + _data_pdfs_dst_00_30;
         double * RESTRICT _data_pdfs_src_00_30_10 = _stride_pdfs_src_1*ctr_1 + _data_pdfs_src_00_30;
         double * RESTRICT  _data_pdfs_dst_00_31_10 = _stride_pdfs_dst_1*ctr_1 + _data_pdfs_dst_00_31;
         double * RESTRICT _data_pdfs_src_00_31_10 = _stride_pdfs_src_1*ctr_1 + _data_pdfs_src_00_31;
         double * RESTRICT  _data_pdfs_dst_00_32_10 = _stride_pdfs_dst_1*ctr_1 + _data_pdfs_dst_00_32;
         double * RESTRICT _data_pdfs_src_00_32_10 = _stride_pdfs_src_1*ctr_1 + _data_pdfs_src_00_32;
         double * RESTRICT  _data_pdfs_dst_00_33_10 = _stride_pdfs_dst_1*ctr_1 + _data_pdfs_dst_00_33;
         double * RESTRICT _data_pdfs_src_00_33_10 = _stride_pdfs_src_1*ctr_1 + _data_pdfs_src_00_33;
         double * RESTRICT  _data_pdfs_dst_00_34_10 = _stride_pdfs_dst_1*ctr_1 + _data_pdfs_dst_00_34;
         double * RESTRICT _data_pdfs_src_00_34_10 = _stride_pdfs_src_1*ctr_1 + _data_pdfs_src_00_34;
         double * RESTRICT  _data_pdfs_dst_00_35_10 = _stride_pdfs_dst_1*ctr_1 + _data_pdfs_dst_00_35;
         double * RESTRICT _data_pdfs_src_00_35_10 = _stride_pdfs_src_1*ctr_1 + _data_pdfs_src_00_35;
         double * RESTRICT  _data_pdfs_dst_00_36_10 = _stride_pdfs_dst_1*ctr_1 + _data_pdfs_dst_00_36;
         double * RESTRICT _data_pdfs_src_00_36_10 = _stride_pdfs_src_1*ctr_1 + _data_pdfs_src_00_36;
         double * RESTRICT  _data_pdfs_dst_00_37_10 = _stride_pdfs_dst_1*ctr_1 + _data_pdfs_dst_00_37;
         double * RESTRICT _data_pdfs_src_00_37_10 = _stride_pdfs_src_1*ctr_1 + _data_pdfs_src_00_37;
         double * RESTRICT  _data_pdfs_dst_00_38_10 = _stride_pdfs_dst_1*ctr_1 + _data_pdfs_dst_00_38;
         double * RESTRICT _data_pdfs_src_00_38_10 = _stride_pdfs_src_1*ctr_1 + _data_pdfs_src_00_38;
         double * RESTRICT  _data_pdfs_dst_00_39_10 = _stride_pdfs_dst_1*ctr_1 + _data_pdfs_dst_00_39;
         double * RESTRICT _data_pdfs_src_00_39_10 = _stride_pdfs_src_1*ctr_1 + _data_pdfs_src_00_39;
         double * RESTRICT  _data_pdfs_dst_00_310_10 = _stride_pdfs_dst_1*ctr_1 + _data_pdfs_dst_00_310;
         double * RESTRICT _data_pdfs_src_00_310_10 = _stride_pdfs_src_1*ctr_1 + _data_pdfs_src_00_310;
         double * RESTRICT  _data_pdfs_dst_00_311_10 = _stride_pdfs_dst_1*ctr_1 + _data_pdfs_dst_00_311;
         double * RESTRICT _data_pdfs_src_00_311_10 = _stride_pdfs_src_1*ctr_1 + _data_pdfs_src_00_311;
         double * RESTRICT  _data_pdfs_dst_00_312_10 = _stride_pdfs_dst_1*ctr_1 + _data_pdfs_dst_00_312;
         double * RESTRICT _data_pdfs_src_00_312_10 = _stride_pdfs_src_1*ctr_1 + _data_pdfs_src_00_312;
         double * RESTRICT  _data_pdfs_dst_00_313_10 = _stride_pdfs_dst_1*ctr_1 + _data_pdfs_dst_00_313;
         double * RESTRICT _data_pdfs_src_00_313_10 = _stride_pdfs_src_1*ctr_1 + _data_pdfs_src_00_313;
         double * RESTRICT  _data_pdfs_dst_00_314_10 = _stride_pdfs_dst_1*ctr_1 + _data_pdfs_dst_00_314;
         double * RESTRICT _data_pdfs_src_00_314_10 = _stride_pdfs_src_1*ctr_1 + _data_pdfs_src_00_314;
         double * RESTRICT  _data_pdfs_dst_00_315_10 = _stride_pdfs_dst_1*ctr_1 + _data_pdfs_dst_00_315;
         double * RESTRICT _data_pdfs_src_00_315_10 = _stride_pdfs_src_1*ctr_1 + _data_pdfs_src_00_315;
         double * RESTRICT  _data_pdfs_dst_00_316_10 = _stride_pdfs_dst_1*ctr_1 + _data_pdfs_dst_00_316;
         double * RESTRICT _data_pdfs_src_00_316_10 = _stride_pdfs_src_1*ctr_1 + _data_pdfs_src_00_316;
         double * RESTRICT  _data_pdfs_dst_00_317_10 = _stride_pdfs_dst_1*ctr_1 + _data_pdfs_dst_00_317;
         double * RESTRICT _data_pdfs_src_00_317_10 = _stride_pdfs_src_1*ctr_1 + _data_pdfs_src_00_317;
         double * RESTRICT  _data_pdfs_dst_00_318_10 = _stride_pdfs_dst_1*ctr_1 + _data_pdfs_dst_00_318;
         double * RESTRICT _data_pdfs_src_00_318_10 = _stride_pdfs_src_1*ctr_1 + _data_pdfs_src_00_318;
         for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_dst_2; ctr_2 += 1)
         {
            _data_pdfs_dst_00_30_10[_stride_pdfs_dst_2*ctr_2] = _data_pdfs_src_00_30_10[_stride_pdfs_src_2*ctr_2];
            _data_pdfs_dst_00_31_10[_stride_pdfs_dst_2*ctr_2] = _data_pdfs_src_00_31_10[_stride_pdfs_src_2*ctr_2];
            _data_pdfs_dst_00_32_10[_stride_pdfs_dst_2*ctr_2] = _data_pdfs_src_00_32_10[_stride_pdfs_src_2*ctr_2];
            _data_pdfs_dst_00_33_10[_stride_pdfs_dst_2*ctr_2] = _data_pdfs_src_00_33_10[_stride_pdfs_src_2*ctr_2];
            _data_pdfs_dst_00_34_10[_stride_pdfs_dst_2*ctr_2] = _data_pdfs_src_00_34_10[_stride_pdfs_src_2*ctr_2];
            _data_pdfs_dst_00_35_10[_stride_pdfs_dst_2*ctr_2] = _data_pdfs_src_00_35_10[_stride_pdfs_src_2*ctr_2];
            _data_pdfs_dst_00_36_10[_stride_pdfs_dst_2*ctr_2] = _data_pdfs_src_00_36_10[_stride_pdfs_src_2*ctr_2];
            _data_pdfs_dst_00_37_10[_stride_pdfs_dst_2*ctr_2] = _data_pdfs_src_00_37_10[_stride_pdfs_src_2*ctr_2];
            _data_pdfs_dst_00_38_10[_stride_pdfs_dst_2*ctr_2] = _data_pdfs_src_00_38_10[_stride_pdfs_src_2*ctr_2];
            _data_pdfs_dst_00_39_10[_stride_pdfs_dst_2*ctr_2] = _data_pdfs_src_00_39_10[_stride_pdfs_src_2*ctr_2];
            _data_pdfs_dst_00_310_10[_stride_pdfs_dst_2*ctr_2] = _data_pdfs_src_00_310_10[_stride_pdfs_src_2*ctr_2];
            _data_pdfs_dst_00_311_10[_stride_pdfs_dst_2*ctr_2] = _data_pdfs_src_00_311_10[_stride_pdfs_src_2*ctr_2];
            _data_pdfs_dst_00_312_10[_stride_pdfs_dst_2*ctr_2] = _data_pdfs_src_00_312_10[_stride_pdfs_src_2*ctr_2];
            _data_pdfs_dst_00_313_10[_stride_pdfs_dst_2*ctr_2] = _data_pdfs_src_00_313_10[_stride_pdfs_src_2*ctr_2];
            _data_pdfs_dst_00_314_10[_stride_pdfs_dst_2*ctr_2] = _data_pdfs_src_00_314_10[_stride_pdfs_src_2*ctr_2];
            _data_pdfs_dst_00_315_10[_stride_pdfs_dst_2*ctr_2] = _data_pdfs_src_00_315_10[_stride_pdfs_src_2*ctr_2];
            _data_pdfs_dst_00_316_10[_stride_pdfs_dst_2*ctr_2] = _data_pdfs_src_00_316_10[_stride_pdfs_src_2*ctr_2];
            _data_pdfs_dst_00_317_10[_stride_pdfs_dst_2*ctr_2] = _data_pdfs_src_00_317_10[_stride_pdfs_src_2*ctr_2];
            _data_pdfs_dst_00_318_10[_stride_pdfs_dst_2*ctr_2] = _data_pdfs_src_00_318_10[_stride_pdfs_src_2*ctr_2];
         }
      }
   }
}
}


namespace internal_d3q19storagespecification_pack_TE {
static FUNC_PREFIX void d3q19storagespecification_pack_TE(double * RESTRICT  _data_buffer, double * RESTRICT const _data_pdfs_src, int64_t const _size_pdfs_src_0, int64_t const _size_pdfs_src_1, int64_t const _size_pdfs_src_2, int64_t const _stride_pdfs_src_0, int64_t const _stride_pdfs_src_1, int64_t const _stride_pdfs_src_2, int64_t const _stride_pdfs_src_3)
{
   for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_src_0; ctr_0 += 1)
   {
      double * RESTRICT _data_pdfs_src_00_314 = _data_pdfs_src + _stride_pdfs_src_0*ctr_0 + 14*_stride_pdfs_src_3;
      for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_src_1; ctr_1 += 1)
      {
         double * RESTRICT _data_pdfs_src_00_314_10 = _stride_pdfs_src_1*ctr_1 + _data_pdfs_src_00_314;
         for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_src_2; ctr_2 += 1)
         {
            _data_buffer[_size_pdfs_src_1*_size_pdfs_src_2*ctr_0 + _size_pdfs_src_2*ctr_1 + ctr_2] = _data_pdfs_src_00_314_10[_stride_pdfs_src_2*ctr_2];
         }
      }
   }
}
}

namespace internal_d3q19storagespecification_pack_SW {
static FUNC_PREFIX void d3q19storagespecification_pack_SW(double * RESTRICT  _data_buffer, double * RESTRICT const _data_pdfs_src, int64_t const _size_pdfs_src_0, int64_t const _size_pdfs_src_1, int64_t const _size_pdfs_src_2, int64_t const _stride_pdfs_src_0, int64_t const _stride_pdfs_src_1, int64_t const _stride_pdfs_src_2, int64_t const _stride_pdfs_src_3)
{
   for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_src_0; ctr_0 += 1)
   {
      double * RESTRICT _data_pdfs_src_00_39 = _data_pdfs_src + _stride_pdfs_src_0*ctr_0 + 9*_stride_pdfs_src_3;
      for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_src_1; ctr_1 += 1)
      {
         double * RESTRICT _data_pdfs_src_00_39_10 = _stride_pdfs_src_1*ctr_1 + _data_pdfs_src_00_39;
         for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_src_2; ctr_2 += 1)
         {
            _data_buffer[_size_pdfs_src_1*_size_pdfs_src_2*ctr_0 + _size_pdfs_src_2*ctr_1 + ctr_2] = _data_pdfs_src_00_39_10[_stride_pdfs_src_2*ctr_2];
         }
      }
   }
}
}

namespace internal_d3q19storagespecification_pack_T {
static FUNC_PREFIX void d3q19storagespecification_pack_T(double * RESTRICT  _data_buffer, double * RESTRICT const _data_pdfs_src, int64_t const _size_pdfs_src_0, int64_t const _size_pdfs_src_1, int64_t const _size_pdfs_src_2, int64_t const _stride_pdfs_src_0, int64_t const _stride_pdfs_src_1, int64_t const _stride_pdfs_src_2, int64_t const _stride_pdfs_src_3)
{
   for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_src_0; ctr_0 += 1)
   {
      double * RESTRICT _data_pdfs_src_00_35 = _data_pdfs_src + _stride_pdfs_src_0*ctr_0 + 5*_stride_pdfs_src_3;
      double * RESTRICT _data_pdfs_src_00_311 = _data_pdfs_src + _stride_pdfs_src_0*ctr_0 + 11*_stride_pdfs_src_3;
      double * RESTRICT _data_pdfs_src_00_312 = _data_pdfs_src + _stride_pdfs_src_0*ctr_0 + 12*_stride_pdfs_src_3;
      double * RESTRICT _data_pdfs_src_00_313 = _data_pdfs_src + _stride_pdfs_src_0*ctr_0 + 13*_stride_pdfs_src_3;
      double * RESTRICT _data_pdfs_src_00_314 = _data_pdfs_src + _stride_pdfs_src_0*ctr_0 + 14*_stride_pdfs_src_3;
      for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_src_1; ctr_1 += 1)
      {
         double * RESTRICT _data_pdfs_src_00_35_10 = _stride_pdfs_src_1*ctr_1 + _data_pdfs_src_00_35;
         double * RESTRICT _data_pdfs_src_00_311_10 = _stride_pdfs_src_1*ctr_1 + _data_pdfs_src_00_311;
         double * RESTRICT _data_pdfs_src_00_312_10 = _stride_pdfs_src_1*ctr_1 + _data_pdfs_src_00_312;
         double * RESTRICT _data_pdfs_src_00_313_10 = _stride_pdfs_src_1*ctr_1 + _data_pdfs_src_00_313;
         double * RESTRICT _data_pdfs_src_00_314_10 = _stride_pdfs_src_1*ctr_1 + _data_pdfs_src_00_314;
         for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_src_2; ctr_2 += 1)
         {
            _data_buffer[5*_size_pdfs_src_1*_size_pdfs_src_2*ctr_0 + 5*_size_pdfs_src_2*ctr_1 + 5*ctr_2] = _data_pdfs_src_00_35_10[_stride_pdfs_src_2*ctr_2];
            _data_buffer[5*_size_pdfs_src_1*_size_pdfs_src_2*ctr_0 + 5*_size_pdfs_src_2*ctr_1 + 5*ctr_2 + 1] = _data_pdfs_src_00_311_10[_stride_pdfs_src_2*ctr_2];
            _data_buffer[5*_size_pdfs_src_1*_size_pdfs_src_2*ctr_0 + 5*_size_pdfs_src_2*ctr_1 + 5*ctr_2 + 2] = _data_pdfs_src_00_312_10[_stride_pdfs_src_2*ctr_2];
            _data_buffer[5*_size_pdfs_src_1*_size_pdfs_src_2*ctr_0 + 5*_size_pdfs_src_2*ctr_1 + 5*ctr_2 + 3] = _data_pdfs_src_00_313_10[_stride_pdfs_src_2*ctr_2];
            _data_buffer[5*_size_pdfs_src_1*_size_pdfs_src_2*ctr_0 + 5*_size_pdfs_src_2*ctr_1 + 5*ctr_2 + 4] = _data_pdfs_src_00_314_10[_stride_pdfs_src_2*ctr_2];
         }
      }
   }
}
}

namespace internal_d3q19storagespecification_pack_BS {
static FUNC_PREFIX void d3q19storagespecification_pack_BS(double * RESTRICT  _data_buffer, double * RESTRICT const _data_pdfs_src, int64_t const _size_pdfs_src_0, int64_t const _size_pdfs_src_1, int64_t const _size_pdfs_src_2, int64_t const _stride_pdfs_src_0, int64_t const _stride_pdfs_src_1, int64_t const _stride_pdfs_src_2, int64_t const _stride_pdfs_src_3)
{
   for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_src_0; ctr_0 += 1)
   {
      double * RESTRICT _data_pdfs_src_00_316 = _data_pdfs_src + _stride_pdfs_src_0*ctr_0 + 16*_stride_pdfs_src_3;
      for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_src_1; ctr_1 += 1)
      {
         double * RESTRICT _data_pdfs_src_00_316_10 = _stride_pdfs_src_1*ctr_1 + _data_pdfs_src_00_316;
         for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_src_2; ctr_2 += 1)
         {
            _data_buffer[_size_pdfs_src_1*_size_pdfs_src_2*ctr_0 + _size_pdfs_src_2*ctr_1 + ctr_2] = _data_pdfs_src_00_316_10[_stride_pdfs_src_2*ctr_2];
         }
      }
   }
}
}

namespace internal_d3q19storagespecification_pack_TN {
static FUNC_PREFIX void d3q19storagespecification_pack_TN(double * RESTRICT  _data_buffer, double * RESTRICT const _data_pdfs_src, int64_t const _size_pdfs_src_0, int64_t const _size_pdfs_src_1, int64_t const _size_pdfs_src_2, int64_t const _stride_pdfs_src_0, int64_t const _stride_pdfs_src_1, int64_t const _stride_pdfs_src_2, int64_t const _stride_pdfs_src_3)
{
   for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_src_0; ctr_0 += 1)
   {
      double * RESTRICT _data_pdfs_src_00_311 = _data_pdfs_src + _stride_pdfs_src_0*ctr_0 + 11*_stride_pdfs_src_3;
      for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_src_1; ctr_1 += 1)
      {
         double * RESTRICT _data_pdfs_src_00_311_10 = _stride_pdfs_src_1*ctr_1 + _data_pdfs_src_00_311;
         for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_src_2; ctr_2 += 1)
         {
            _data_buffer[_size_pdfs_src_1*_size_pdfs_src_2*ctr_0 + _size_pdfs_src_2*ctr_1 + ctr_2] = _data_pdfs_src_00_311_10[_stride_pdfs_src_2*ctr_2];
         }
      }
   }
}
}

namespace internal_d3q19storagespecification_pack_BW {
static FUNC_PREFIX void d3q19storagespecification_pack_BW(double * RESTRICT  _data_buffer, double * RESTRICT const _data_pdfs_src, int64_t const _size_pdfs_src_0, int64_t const _size_pdfs_src_1, int64_t const _size_pdfs_src_2, int64_t const _stride_pdfs_src_0, int64_t const _stride_pdfs_src_1, int64_t const _stride_pdfs_src_2, int64_t const _stride_pdfs_src_3)
{
   for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_src_0; ctr_0 += 1)
   {
      double * RESTRICT _data_pdfs_src_00_317 = _data_pdfs_src + _stride_pdfs_src_0*ctr_0 + 17*_stride_pdfs_src_3;
      for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_src_1; ctr_1 += 1)
      {
         double * RESTRICT _data_pdfs_src_00_317_10 = _stride_pdfs_src_1*ctr_1 + _data_pdfs_src_00_317;
         for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_src_2; ctr_2 += 1)
         {
            _data_buffer[_size_pdfs_src_1*_size_pdfs_src_2*ctr_0 + _size_pdfs_src_2*ctr_1 + ctr_2] = _data_pdfs_src_00_317_10[_stride_pdfs_src_2*ctr_2];
         }
      }
   }
}
}

namespace internal_d3q19storagespecification_pack_N {
static FUNC_PREFIX void d3q19storagespecification_pack_N(double * RESTRICT  _data_buffer, double * RESTRICT const _data_pdfs_src, int64_t const _size_pdfs_src_0, int64_t const _size_pdfs_src_1, int64_t const _size_pdfs_src_2, int64_t const _stride_pdfs_src_0, int64_t const _stride_pdfs_src_1, int64_t const _stride_pdfs_src_2, int64_t const _stride_pdfs_src_3)
{
   for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_src_0; ctr_0 += 1)
   {
      double * RESTRICT _data_pdfs_src_00_31 = _data_pdfs_src + _stride_pdfs_src_0*ctr_0 + _stride_pdfs_src_3;
      double * RESTRICT _data_pdfs_src_00_37 = _data_pdfs_src + _stride_pdfs_src_0*ctr_0 + 7*_stride_pdfs_src_3;
      double * RESTRICT _data_pdfs_src_00_38 = _data_pdfs_src + _stride_pdfs_src_0*ctr_0 + 8*_stride_pdfs_src_3;
      double * RESTRICT _data_pdfs_src_00_311 = _data_pdfs_src + _stride_pdfs_src_0*ctr_0 + 11*_stride_pdfs_src_3;
      double * RESTRICT _data_pdfs_src_00_315 = _data_pdfs_src + _stride_pdfs_src_0*ctr_0 + 15*_stride_pdfs_src_3;
      for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_src_1; ctr_1 += 1)
      {
         double * RESTRICT _data_pdfs_src_00_31_10 = _stride_pdfs_src_1*ctr_1 + _data_pdfs_src_00_31;
         double * RESTRICT _data_pdfs_src_00_37_10 = _stride_pdfs_src_1*ctr_1 + _data_pdfs_src_00_37;
         double * RESTRICT _data_pdfs_src_00_38_10 = _stride_pdfs_src_1*ctr_1 + _data_pdfs_src_00_38;
         double * RESTRICT _data_pdfs_src_00_311_10 = _stride_pdfs_src_1*ctr_1 + _data_pdfs_src_00_311;
         double * RESTRICT _data_pdfs_src_00_315_10 = _stride_pdfs_src_1*ctr_1 + _data_pdfs_src_00_315;
         for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_src_2; ctr_2 += 1)
         {
            _data_buffer[5*_size_pdfs_src_1*_size_pdfs_src_2*ctr_0 + 5*_size_pdfs_src_2*ctr_1 + 5*ctr_2] = _data_pdfs_src_00_31_10[_stride_pdfs_src_2*ctr_2];
            _data_buffer[5*_size_pdfs_src_1*_size_pdfs_src_2*ctr_0 + 5*_size_pdfs_src_2*ctr_1 + 5*ctr_2 + 1] = _data_pdfs_src_00_37_10[_stride_pdfs_src_2*ctr_2];
            _data_buffer[5*_size_pdfs_src_1*_size_pdfs_src_2*ctr_0 + 5*_size_pdfs_src_2*ctr_1 + 5*ctr_2 + 2] = _data_pdfs_src_00_38_10[_stride_pdfs_src_2*ctr_2];
            _data_buffer[5*_size_pdfs_src_1*_size_pdfs_src_2*ctr_0 + 5*_size_pdfs_src_2*ctr_1 + 5*ctr_2 + 3] = _data_pdfs_src_00_311_10[_stride_pdfs_src_2*ctr_2];
            _data_buffer[5*_size_pdfs_src_1*_size_pdfs_src_2*ctr_0 + 5*_size_pdfs_src_2*ctr_1 + 5*ctr_2 + 4] = _data_pdfs_src_00_315_10[_stride_pdfs_src_2*ctr_2];
         }
      }
   }
}
}

namespace internal_d3q19storagespecification_pack_E {
static FUNC_PREFIX void d3q19storagespecification_pack_E(double * RESTRICT  _data_buffer, double * RESTRICT const _data_pdfs_src, int64_t const _size_pdfs_src_0, int64_t const _size_pdfs_src_1, int64_t const _size_pdfs_src_2, int64_t const _stride_pdfs_src_0, int64_t const _stride_pdfs_src_1, int64_t const _stride_pdfs_src_2, int64_t const _stride_pdfs_src_3)
{
   for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_src_0; ctr_0 += 1)
   {
      double * RESTRICT _data_pdfs_src_00_34 = _data_pdfs_src + _stride_pdfs_src_0*ctr_0 + 4*_stride_pdfs_src_3;
      double * RESTRICT _data_pdfs_src_00_38 = _data_pdfs_src + _stride_pdfs_src_0*ctr_0 + 8*_stride_pdfs_src_3;
      double * RESTRICT _data_pdfs_src_00_310 = _data_pdfs_src + _stride_pdfs_src_0*ctr_0 + 10*_stride_pdfs_src_3;
      double * RESTRICT _data_pdfs_src_00_314 = _data_pdfs_src + _stride_pdfs_src_0*ctr_0 + 14*_stride_pdfs_src_3;
      double * RESTRICT _data_pdfs_src_00_318 = _data_pdfs_src + _stride_pdfs_src_0*ctr_0 + 18*_stride_pdfs_src_3;
      for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_src_1; ctr_1 += 1)
      {
         double * RESTRICT _data_pdfs_src_00_34_10 = _stride_pdfs_src_1*ctr_1 + _data_pdfs_src_00_34;
         double * RESTRICT _data_pdfs_src_00_38_10 = _stride_pdfs_src_1*ctr_1 + _data_pdfs_src_00_38;
         double * RESTRICT _data_pdfs_src_00_310_10 = _stride_pdfs_src_1*ctr_1 + _data_pdfs_src_00_310;
         double * RESTRICT _data_pdfs_src_00_314_10 = _stride_pdfs_src_1*ctr_1 + _data_pdfs_src_00_314;
         double * RESTRICT _data_pdfs_src_00_318_10 = _stride_pdfs_src_1*ctr_1 + _data_pdfs_src_00_318;
         for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_src_2; ctr_2 += 1)
         {
            _data_buffer[5*_size_pdfs_src_1*_size_pdfs_src_2*ctr_0 + 5*_size_pdfs_src_2*ctr_1 + 5*ctr_2] = _data_pdfs_src_00_34_10[_stride_pdfs_src_2*ctr_2];
            _data_buffer[5*_size_pdfs_src_1*_size_pdfs_src_2*ctr_0 + 5*_size_pdfs_src_2*ctr_1 + 5*ctr_2 + 1] = _data_pdfs_src_00_38_10[_stride_pdfs_src_2*ctr_2];
            _data_buffer[5*_size_pdfs_src_1*_size_pdfs_src_2*ctr_0 + 5*_size_pdfs_src_2*ctr_1 + 5*ctr_2 + 2] = _data_pdfs_src_00_310_10[_stride_pdfs_src_2*ctr_2];
            _data_buffer[5*_size_pdfs_src_1*_size_pdfs_src_2*ctr_0 + 5*_size_pdfs_src_2*ctr_1 + 5*ctr_2 + 3] = _data_pdfs_src_00_314_10[_stride_pdfs_src_2*ctr_2];
            _data_buffer[5*_size_pdfs_src_1*_size_pdfs_src_2*ctr_0 + 5*_size_pdfs_src_2*ctr_1 + 5*ctr_2 + 4] = _data_pdfs_src_00_318_10[_stride_pdfs_src_2*ctr_2];
         }
      }
   }
}
}

namespace internal_d3q19storagespecification_pack_NW {
static FUNC_PREFIX void d3q19storagespecification_pack_NW(double * RESTRICT  _data_buffer, double * RESTRICT const _data_pdfs_src, int64_t const _size_pdfs_src_0, int64_t const _size_pdfs_src_1, int64_t const _size_pdfs_src_2, int64_t const _stride_pdfs_src_0, int64_t const _stride_pdfs_src_1, int64_t const _stride_pdfs_src_2, int64_t const _stride_pdfs_src_3)
{
   for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_src_0; ctr_0 += 1)
   {
      double * RESTRICT _data_pdfs_src_00_37 = _data_pdfs_src + _stride_pdfs_src_0*ctr_0 + 7*_stride_pdfs_src_3;
      for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_src_1; ctr_1 += 1)
      {
         double * RESTRICT _data_pdfs_src_00_37_10 = _stride_pdfs_src_1*ctr_1 + _data_pdfs_src_00_37;
         for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_src_2; ctr_2 += 1)
         {
            _data_buffer[_size_pdfs_src_1*_size_pdfs_src_2*ctr_0 + _size_pdfs_src_2*ctr_1 + ctr_2] = _data_pdfs_src_00_37_10[_stride_pdfs_src_2*ctr_2];
         }
      }
   }
}
}

namespace internal_d3q19storagespecification_pack_NE {
static FUNC_PREFIX void d3q19storagespecification_pack_NE(double * RESTRICT  _data_buffer, double * RESTRICT const _data_pdfs_src, int64_t const _size_pdfs_src_0, int64_t const _size_pdfs_src_1, int64_t const _size_pdfs_src_2, int64_t const _stride_pdfs_src_0, int64_t const _stride_pdfs_src_1, int64_t const _stride_pdfs_src_2, int64_t const _stride_pdfs_src_3)
{
   for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_src_0; ctr_0 += 1)
   {
      double * RESTRICT _data_pdfs_src_00_38 = _data_pdfs_src + _stride_pdfs_src_0*ctr_0 + 8*_stride_pdfs_src_3;
      for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_src_1; ctr_1 += 1)
      {
         double * RESTRICT _data_pdfs_src_00_38_10 = _stride_pdfs_src_1*ctr_1 + _data_pdfs_src_00_38;
         for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_src_2; ctr_2 += 1)
         {
            _data_buffer[_size_pdfs_src_1*_size_pdfs_src_2*ctr_0 + _size_pdfs_src_2*ctr_1 + ctr_2] = _data_pdfs_src_00_38_10[_stride_pdfs_src_2*ctr_2];
         }
      }
   }
}
}

namespace internal_d3q19storagespecification_pack_TW {
static FUNC_PREFIX void d3q19storagespecification_pack_TW(double * RESTRICT  _data_buffer, double * RESTRICT const _data_pdfs_src, int64_t const _size_pdfs_src_0, int64_t const _size_pdfs_src_1, int64_t const _size_pdfs_src_2, int64_t const _stride_pdfs_src_0, int64_t const _stride_pdfs_src_1, int64_t const _stride_pdfs_src_2, int64_t const _stride_pdfs_src_3)
{
   for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_src_0; ctr_0 += 1)
   {
      double * RESTRICT _data_pdfs_src_00_313 = _data_pdfs_src + _stride_pdfs_src_0*ctr_0 + 13*_stride_pdfs_src_3;
      for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_src_1; ctr_1 += 1)
      {
         double * RESTRICT _data_pdfs_src_00_313_10 = _stride_pdfs_src_1*ctr_1 + _data_pdfs_src_00_313;
         for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_src_2; ctr_2 += 1)
         {
            _data_buffer[_size_pdfs_src_1*_size_pdfs_src_2*ctr_0 + _size_pdfs_src_2*ctr_1 + ctr_2] = _data_pdfs_src_00_313_10[_stride_pdfs_src_2*ctr_2];
         }
      }
   }
}
}

namespace internal_d3q19storagespecification_pack_BE {
static FUNC_PREFIX void d3q19storagespecification_pack_BE(double * RESTRICT  _data_buffer, double * RESTRICT const _data_pdfs_src, int64_t const _size_pdfs_src_0, int64_t const _size_pdfs_src_1, int64_t const _size_pdfs_src_2, int64_t const _stride_pdfs_src_0, int64_t const _stride_pdfs_src_1, int64_t const _stride_pdfs_src_2, int64_t const _stride_pdfs_src_3)
{
   for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_src_0; ctr_0 += 1)
   {
      double * RESTRICT _data_pdfs_src_00_318 = _data_pdfs_src + _stride_pdfs_src_0*ctr_0 + 18*_stride_pdfs_src_3;
      for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_src_1; ctr_1 += 1)
      {
         double * RESTRICT _data_pdfs_src_00_318_10 = _stride_pdfs_src_1*ctr_1 + _data_pdfs_src_00_318;
         for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_src_2; ctr_2 += 1)
         {
            _data_buffer[_size_pdfs_src_1*_size_pdfs_src_2*ctr_0 + _size_pdfs_src_2*ctr_1 + ctr_2] = _data_pdfs_src_00_318_10[_stride_pdfs_src_2*ctr_2];
         }
      }
   }
}
}

namespace internal_d3q19storagespecification_pack_W {
static FUNC_PREFIX void d3q19storagespecification_pack_W(double * RESTRICT  _data_buffer, double * RESTRICT const _data_pdfs_src, int64_t const _size_pdfs_src_0, int64_t const _size_pdfs_src_1, int64_t const _size_pdfs_src_2, int64_t const _stride_pdfs_src_0, int64_t const _stride_pdfs_src_1, int64_t const _stride_pdfs_src_2, int64_t const _stride_pdfs_src_3)
{
   for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_src_0; ctr_0 += 1)
   {
      double * RESTRICT _data_pdfs_src_00_33 = _data_pdfs_src + _stride_pdfs_src_0*ctr_0 + 3*_stride_pdfs_src_3;
      double * RESTRICT _data_pdfs_src_00_37 = _data_pdfs_src + _stride_pdfs_src_0*ctr_0 + 7*_stride_pdfs_src_3;
      double * RESTRICT _data_pdfs_src_00_39 = _data_pdfs_src + _stride_pdfs_src_0*ctr_0 + 9*_stride_pdfs_src_3;
      double * RESTRICT _data_pdfs_src_00_313 = _data_pdfs_src + _stride_pdfs_src_0*ctr_0 + 13*_stride_pdfs_src_3;
      double * RESTRICT _data_pdfs_src_00_317 = _data_pdfs_src + _stride_pdfs_src_0*ctr_0 + 17*_stride_pdfs_src_3;
      for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_src_1; ctr_1 += 1)
      {
         double * RESTRICT _data_pdfs_src_00_33_10 = _stride_pdfs_src_1*ctr_1 + _data_pdfs_src_00_33;
         double * RESTRICT _data_pdfs_src_00_37_10 = _stride_pdfs_src_1*ctr_1 + _data_pdfs_src_00_37;
         double * RESTRICT _data_pdfs_src_00_39_10 = _stride_pdfs_src_1*ctr_1 + _data_pdfs_src_00_39;
         double * RESTRICT _data_pdfs_src_00_313_10 = _stride_pdfs_src_1*ctr_1 + _data_pdfs_src_00_313;
         double * RESTRICT _data_pdfs_src_00_317_10 = _stride_pdfs_src_1*ctr_1 + _data_pdfs_src_00_317;
         for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_src_2; ctr_2 += 1)
         {
            _data_buffer[5*_size_pdfs_src_1*_size_pdfs_src_2*ctr_0 + 5*_size_pdfs_src_2*ctr_1 + 5*ctr_2] = _data_pdfs_src_00_33_10[_stride_pdfs_src_2*ctr_2];
            _data_buffer[5*_size_pdfs_src_1*_size_pdfs_src_2*ctr_0 + 5*_size_pdfs_src_2*ctr_1 + 5*ctr_2 + 1] = _data_pdfs_src_00_37_10[_stride_pdfs_src_2*ctr_2];
            _data_buffer[5*_size_pdfs_src_1*_size_pdfs_src_2*ctr_0 + 5*_size_pdfs_src_2*ctr_1 + 5*ctr_2 + 2] = _data_pdfs_src_00_39_10[_stride_pdfs_src_2*ctr_2];
            _data_buffer[5*_size_pdfs_src_1*_size_pdfs_src_2*ctr_0 + 5*_size_pdfs_src_2*ctr_1 + 5*ctr_2 + 3] = _data_pdfs_src_00_313_10[_stride_pdfs_src_2*ctr_2];
            _data_buffer[5*_size_pdfs_src_1*_size_pdfs_src_2*ctr_0 + 5*_size_pdfs_src_2*ctr_1 + 5*ctr_2 + 4] = _data_pdfs_src_00_317_10[_stride_pdfs_src_2*ctr_2];
         }
      }
   }
}
}

namespace internal_d3q19storagespecification_pack_S {
static FUNC_PREFIX void d3q19storagespecification_pack_S(double * RESTRICT  _data_buffer, double * RESTRICT const _data_pdfs_src, int64_t const _size_pdfs_src_0, int64_t const _size_pdfs_src_1, int64_t const _size_pdfs_src_2, int64_t const _stride_pdfs_src_0, int64_t const _stride_pdfs_src_1, int64_t const _stride_pdfs_src_2, int64_t const _stride_pdfs_src_3)
{
   for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_src_0; ctr_0 += 1)
   {
      double * RESTRICT _data_pdfs_src_00_32 = _data_pdfs_src + _stride_pdfs_src_0*ctr_0 + 2*_stride_pdfs_src_3;
      double * RESTRICT _data_pdfs_src_00_39 = _data_pdfs_src + _stride_pdfs_src_0*ctr_0 + 9*_stride_pdfs_src_3;
      double * RESTRICT _data_pdfs_src_00_310 = _data_pdfs_src + _stride_pdfs_src_0*ctr_0 + 10*_stride_pdfs_src_3;
      double * RESTRICT _data_pdfs_src_00_312 = _data_pdfs_src + _stride_pdfs_src_0*ctr_0 + 12*_stride_pdfs_src_3;
      double * RESTRICT _data_pdfs_src_00_316 = _data_pdfs_src + _stride_pdfs_src_0*ctr_0 + 16*_stride_pdfs_src_3;
      for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_src_1; ctr_1 += 1)
      {
         double * RESTRICT _data_pdfs_src_00_32_10 = _stride_pdfs_src_1*ctr_1 + _data_pdfs_src_00_32;
         double * RESTRICT _data_pdfs_src_00_39_10 = _stride_pdfs_src_1*ctr_1 + _data_pdfs_src_00_39;
         double * RESTRICT _data_pdfs_src_00_310_10 = _stride_pdfs_src_1*ctr_1 + _data_pdfs_src_00_310;
         double * RESTRICT _data_pdfs_src_00_312_10 = _stride_pdfs_src_1*ctr_1 + _data_pdfs_src_00_312;
         double * RESTRICT _data_pdfs_src_00_316_10 = _stride_pdfs_src_1*ctr_1 + _data_pdfs_src_00_316;
         for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_src_2; ctr_2 += 1)
         {
            _data_buffer[5*_size_pdfs_src_1*_size_pdfs_src_2*ctr_0 + 5*_size_pdfs_src_2*ctr_1 + 5*ctr_2] = _data_pdfs_src_00_32_10[_stride_pdfs_src_2*ctr_2];
            _data_buffer[5*_size_pdfs_src_1*_size_pdfs_src_2*ctr_0 + 5*_size_pdfs_src_2*ctr_1 + 5*ctr_2 + 1] = _data_pdfs_src_00_39_10[_stride_pdfs_src_2*ctr_2];
            _data_buffer[5*_size_pdfs_src_1*_size_pdfs_src_2*ctr_0 + 5*_size_pdfs_src_2*ctr_1 + 5*ctr_2 + 2] = _data_pdfs_src_00_310_10[_stride_pdfs_src_2*ctr_2];
            _data_buffer[5*_size_pdfs_src_1*_size_pdfs_src_2*ctr_0 + 5*_size_pdfs_src_2*ctr_1 + 5*ctr_2 + 3] = _data_pdfs_src_00_312_10[_stride_pdfs_src_2*ctr_2];
            _data_buffer[5*_size_pdfs_src_1*_size_pdfs_src_2*ctr_0 + 5*_size_pdfs_src_2*ctr_1 + 5*ctr_2 + 4] = _data_pdfs_src_00_316_10[_stride_pdfs_src_2*ctr_2];
         }
      }
   }
}
}

namespace internal_d3q19storagespecification_pack_SE {
static FUNC_PREFIX void d3q19storagespecification_pack_SE(double * RESTRICT  _data_buffer, double * RESTRICT const _data_pdfs_src, int64_t const _size_pdfs_src_0, int64_t const _size_pdfs_src_1, int64_t const _size_pdfs_src_2, int64_t const _stride_pdfs_src_0, int64_t const _stride_pdfs_src_1, int64_t const _stride_pdfs_src_2, int64_t const _stride_pdfs_src_3)
{
   for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_src_0; ctr_0 += 1)
   {
      double * RESTRICT _data_pdfs_src_00_310 = _data_pdfs_src + _stride_pdfs_src_0*ctr_0 + 10*_stride_pdfs_src_3;
      for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_src_1; ctr_1 += 1)
      {
         double * RESTRICT _data_pdfs_src_00_310_10 = _stride_pdfs_src_1*ctr_1 + _data_pdfs_src_00_310;
         for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_src_2; ctr_2 += 1)
         {
            _data_buffer[_size_pdfs_src_1*_size_pdfs_src_2*ctr_0 + _size_pdfs_src_2*ctr_1 + ctr_2] = _data_pdfs_src_00_310_10[_stride_pdfs_src_2*ctr_2];
         }
      }
   }
}
}

namespace internal_d3q19storagespecification_pack_TS {
static FUNC_PREFIX void d3q19storagespecification_pack_TS(double * RESTRICT  _data_buffer, double * RESTRICT const _data_pdfs_src, int64_t const _size_pdfs_src_0, int64_t const _size_pdfs_src_1, int64_t const _size_pdfs_src_2, int64_t const _stride_pdfs_src_0, int64_t const _stride_pdfs_src_1, int64_t const _stride_pdfs_src_2, int64_t const _stride_pdfs_src_3)
{
   for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_src_0; ctr_0 += 1)
   {
      double * RESTRICT _data_pdfs_src_00_312 = _data_pdfs_src + _stride_pdfs_src_0*ctr_0 + 12*_stride_pdfs_src_3;
      for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_src_1; ctr_1 += 1)
      {
         double * RESTRICT _data_pdfs_src_00_312_10 = _stride_pdfs_src_1*ctr_1 + _data_pdfs_src_00_312;
         for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_src_2; ctr_2 += 1)
         {
            _data_buffer[_size_pdfs_src_1*_size_pdfs_src_2*ctr_0 + _size_pdfs_src_2*ctr_1 + ctr_2] = _data_pdfs_src_00_312_10[_stride_pdfs_src_2*ctr_2];
         }
      }
   }
}
}

namespace internal_d3q19storagespecification_pack_BN {
static FUNC_PREFIX void d3q19storagespecification_pack_BN(double * RESTRICT  _data_buffer, double * RESTRICT const _data_pdfs_src, int64_t const _size_pdfs_src_0, int64_t const _size_pdfs_src_1, int64_t const _size_pdfs_src_2, int64_t const _stride_pdfs_src_0, int64_t const _stride_pdfs_src_1, int64_t const _stride_pdfs_src_2, int64_t const _stride_pdfs_src_3)
{
   for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_src_0; ctr_0 += 1)
   {
      double * RESTRICT _data_pdfs_src_00_315 = _data_pdfs_src + _stride_pdfs_src_0*ctr_0 + 15*_stride_pdfs_src_3;
      for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_src_1; ctr_1 += 1)
      {
         double * RESTRICT _data_pdfs_src_00_315_10 = _stride_pdfs_src_1*ctr_1 + _data_pdfs_src_00_315;
         for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_src_2; ctr_2 += 1)
         {
            _data_buffer[_size_pdfs_src_1*_size_pdfs_src_2*ctr_0 + _size_pdfs_src_2*ctr_1 + ctr_2] = _data_pdfs_src_00_315_10[_stride_pdfs_src_2*ctr_2];
         }
      }
   }
}
}

namespace internal_d3q19storagespecification_pack_B {
static FUNC_PREFIX void d3q19storagespecification_pack_B(double * RESTRICT  _data_buffer, double * RESTRICT const _data_pdfs_src, int64_t const _size_pdfs_src_0, int64_t const _size_pdfs_src_1, int64_t const _size_pdfs_src_2, int64_t const _stride_pdfs_src_0, int64_t const _stride_pdfs_src_1, int64_t const _stride_pdfs_src_2, int64_t const _stride_pdfs_src_3)
{
   for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_src_0; ctr_0 += 1)
   {
      double * RESTRICT _data_pdfs_src_00_36 = _data_pdfs_src + _stride_pdfs_src_0*ctr_0 + 6*_stride_pdfs_src_3;
      double * RESTRICT _data_pdfs_src_00_315 = _data_pdfs_src + _stride_pdfs_src_0*ctr_0 + 15*_stride_pdfs_src_3;
      double * RESTRICT _data_pdfs_src_00_316 = _data_pdfs_src + _stride_pdfs_src_0*ctr_0 + 16*_stride_pdfs_src_3;
      double * RESTRICT _data_pdfs_src_00_317 = _data_pdfs_src + _stride_pdfs_src_0*ctr_0 + 17*_stride_pdfs_src_3;
      double * RESTRICT _data_pdfs_src_00_318 = _data_pdfs_src + _stride_pdfs_src_0*ctr_0 + 18*_stride_pdfs_src_3;
      for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_src_1; ctr_1 += 1)
      {
         double * RESTRICT _data_pdfs_src_00_36_10 = _stride_pdfs_src_1*ctr_1 + _data_pdfs_src_00_36;
         double * RESTRICT _data_pdfs_src_00_315_10 = _stride_pdfs_src_1*ctr_1 + _data_pdfs_src_00_315;
         double * RESTRICT _data_pdfs_src_00_316_10 = _stride_pdfs_src_1*ctr_1 + _data_pdfs_src_00_316;
         double * RESTRICT _data_pdfs_src_00_317_10 = _stride_pdfs_src_1*ctr_1 + _data_pdfs_src_00_317;
         double * RESTRICT _data_pdfs_src_00_318_10 = _stride_pdfs_src_1*ctr_1 + _data_pdfs_src_00_318;
         for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_src_2; ctr_2 += 1)
         {
            _data_buffer[5*_size_pdfs_src_1*_size_pdfs_src_2*ctr_0 + 5*_size_pdfs_src_2*ctr_1 + 5*ctr_2] = _data_pdfs_src_00_36_10[_stride_pdfs_src_2*ctr_2];
            _data_buffer[5*_size_pdfs_src_1*_size_pdfs_src_2*ctr_0 + 5*_size_pdfs_src_2*ctr_1 + 5*ctr_2 + 1] = _data_pdfs_src_00_315_10[_stride_pdfs_src_2*ctr_2];
            _data_buffer[5*_size_pdfs_src_1*_size_pdfs_src_2*ctr_0 + 5*_size_pdfs_src_2*ctr_1 + 5*ctr_2 + 2] = _data_pdfs_src_00_316_10[_stride_pdfs_src_2*ctr_2];
            _data_buffer[5*_size_pdfs_src_1*_size_pdfs_src_2*ctr_0 + 5*_size_pdfs_src_2*ctr_1 + 5*ctr_2 + 3] = _data_pdfs_src_00_317_10[_stride_pdfs_src_2*ctr_2];
            _data_buffer[5*_size_pdfs_src_1*_size_pdfs_src_2*ctr_0 + 5*_size_pdfs_src_2*ctr_1 + 5*ctr_2 + 4] = _data_pdfs_src_00_318_10[_stride_pdfs_src_2*ctr_2];
         }
      }
   }
}
}

namespace internal_d3q19storagespecification_unpack_BW {
static FUNC_PREFIX void d3q19storagespecification_unpack_BW(const double * RESTRICT const _data_buffer, double * RESTRICT  _data_pdfs_dst, int64_t const _size_pdfs_dst_0, int64_t const _size_pdfs_dst_1, int64_t const _size_pdfs_dst_2, int64_t const _stride_pdfs_dst_0, int64_t const _stride_pdfs_dst_1, int64_t const _stride_pdfs_dst_2, int64_t const _stride_pdfs_dst_3)
{
   for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_dst_0; ctr_0 += 1)
   {
      double * RESTRICT  _data_pdfs_dst_00_314 = _data_pdfs_dst + _stride_pdfs_dst_0*ctr_0 + 14*_stride_pdfs_dst_3;
      for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_dst_1; ctr_1 += 1)
      {
         double * RESTRICT  _data_pdfs_dst_00_314_10 = _stride_pdfs_dst_1*ctr_1 + _data_pdfs_dst_00_314;
         for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_dst_2; ctr_2 += 1)
         {
            _data_pdfs_dst_00_314_10[_stride_pdfs_dst_2*ctr_2] = _data_buffer[_size_pdfs_dst_1*_size_pdfs_dst_2*ctr_0 + _size_pdfs_dst_2*ctr_1 + ctr_2];
         }
      }
   }
}
}

namespace internal_d3q19storagespecification_unpack_N {
static FUNC_PREFIX void d3q19storagespecification_unpack_N(const double * RESTRICT const _data_buffer, double * RESTRICT  _data_pdfs_dst, int64_t const _size_pdfs_dst_0, int64_t const _size_pdfs_dst_1, int64_t const _size_pdfs_dst_2, int64_t const _stride_pdfs_dst_0, int64_t const _stride_pdfs_dst_1, int64_t const _stride_pdfs_dst_2, int64_t const _stride_pdfs_dst_3)
{
   for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_dst_0; ctr_0 += 1)
   {
      double * RESTRICT  _data_pdfs_dst_00_32 = _data_pdfs_dst + _stride_pdfs_dst_0*ctr_0 + 2*_stride_pdfs_dst_3;
      double * RESTRICT  _data_pdfs_dst_00_39 = _data_pdfs_dst + _stride_pdfs_dst_0*ctr_0 + 9*_stride_pdfs_dst_3;
      double * RESTRICT  _data_pdfs_dst_00_310 = _data_pdfs_dst + _stride_pdfs_dst_0*ctr_0 + 10*_stride_pdfs_dst_3;
      double * RESTRICT  _data_pdfs_dst_00_312 = _data_pdfs_dst + _stride_pdfs_dst_0*ctr_0 + 12*_stride_pdfs_dst_3;
      double * RESTRICT  _data_pdfs_dst_00_316 = _data_pdfs_dst + _stride_pdfs_dst_0*ctr_0 + 16*_stride_pdfs_dst_3;
      for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_dst_1; ctr_1 += 1)
      {
         double * RESTRICT  _data_pdfs_dst_00_32_10 = _stride_pdfs_dst_1*ctr_1 + _data_pdfs_dst_00_32;
         double * RESTRICT  _data_pdfs_dst_00_39_10 = _stride_pdfs_dst_1*ctr_1 + _data_pdfs_dst_00_39;
         double * RESTRICT  _data_pdfs_dst_00_310_10 = _stride_pdfs_dst_1*ctr_1 + _data_pdfs_dst_00_310;
         double * RESTRICT  _data_pdfs_dst_00_312_10 = _stride_pdfs_dst_1*ctr_1 + _data_pdfs_dst_00_312;
         double * RESTRICT  _data_pdfs_dst_00_316_10 = _stride_pdfs_dst_1*ctr_1 + _data_pdfs_dst_00_316;
         for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_dst_2; ctr_2 += 1)
         {
            _data_pdfs_dst_00_32_10[_stride_pdfs_dst_2*ctr_2] = _data_buffer[5*_size_pdfs_dst_1*_size_pdfs_dst_2*ctr_0 + 5*_size_pdfs_dst_2*ctr_1 + 5*ctr_2];
            _data_pdfs_dst_00_39_10[_stride_pdfs_dst_2*ctr_2] = _data_buffer[5*_size_pdfs_dst_1*_size_pdfs_dst_2*ctr_0 + 5*_size_pdfs_dst_2*ctr_1 + 5*ctr_2 + 1];
            _data_pdfs_dst_00_310_10[_stride_pdfs_dst_2*ctr_2] = _data_buffer[5*_size_pdfs_dst_1*_size_pdfs_dst_2*ctr_0 + 5*_size_pdfs_dst_2*ctr_1 + 5*ctr_2 + 2];
            _data_pdfs_dst_00_312_10[_stride_pdfs_dst_2*ctr_2] = _data_buffer[5*_size_pdfs_dst_1*_size_pdfs_dst_2*ctr_0 + 5*_size_pdfs_dst_2*ctr_1 + 5*ctr_2 + 3];
            _data_pdfs_dst_00_316_10[_stride_pdfs_dst_2*ctr_2] = _data_buffer[5*_size_pdfs_dst_1*_size_pdfs_dst_2*ctr_0 + 5*_size_pdfs_dst_2*ctr_1 + 5*ctr_2 + 4];
         }
      }
   }
}
}

namespace internal_d3q19storagespecification_unpack_SE {
static FUNC_PREFIX void d3q19storagespecification_unpack_SE(const double * RESTRICT const _data_buffer, double * RESTRICT  _data_pdfs_dst, int64_t const _size_pdfs_dst_0, int64_t const _size_pdfs_dst_1, int64_t const _size_pdfs_dst_2, int64_t const _stride_pdfs_dst_0, int64_t const _stride_pdfs_dst_1, int64_t const _stride_pdfs_dst_2, int64_t const _stride_pdfs_dst_3)
{
   for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_dst_0; ctr_0 += 1)
   {
      double * RESTRICT  _data_pdfs_dst_00_37 = _data_pdfs_dst + _stride_pdfs_dst_0*ctr_0 + 7*_stride_pdfs_dst_3;
      for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_dst_1; ctr_1 += 1)
      {
         double * RESTRICT  _data_pdfs_dst_00_37_10 = _stride_pdfs_dst_1*ctr_1 + _data_pdfs_dst_00_37;
         for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_dst_2; ctr_2 += 1)
         {
            _data_pdfs_dst_00_37_10[_stride_pdfs_dst_2*ctr_2] = _data_buffer[_size_pdfs_dst_1*_size_pdfs_dst_2*ctr_0 + _size_pdfs_dst_2*ctr_1 + ctr_2];
         }
      }
   }
}
}

namespace internal_d3q19storagespecification_unpack_TE {
static FUNC_PREFIX void d3q19storagespecification_unpack_TE(const double * RESTRICT const _data_buffer, double * RESTRICT  _data_pdfs_dst, int64_t const _size_pdfs_dst_0, int64_t const _size_pdfs_dst_1, int64_t const _size_pdfs_dst_2, int64_t const _stride_pdfs_dst_0, int64_t const _stride_pdfs_dst_1, int64_t const _stride_pdfs_dst_2, int64_t const _stride_pdfs_dst_3)
{
   for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_dst_0; ctr_0 += 1)
   {
      double * RESTRICT  _data_pdfs_dst_00_317 = _data_pdfs_dst + _stride_pdfs_dst_0*ctr_0 + 17*_stride_pdfs_dst_3;
      for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_dst_1; ctr_1 += 1)
      {
         double * RESTRICT  _data_pdfs_dst_00_317_10 = _stride_pdfs_dst_1*ctr_1 + _data_pdfs_dst_00_317;
         for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_dst_2; ctr_2 += 1)
         {
            _data_pdfs_dst_00_317_10[_stride_pdfs_dst_2*ctr_2] = _data_buffer[_size_pdfs_dst_1*_size_pdfs_dst_2*ctr_0 + _size_pdfs_dst_2*ctr_1 + ctr_2];
         }
      }
   }
}
}

namespace internal_d3q19storagespecification_unpack_T {
static FUNC_PREFIX void d3q19storagespecification_unpack_T(const double * RESTRICT const _data_buffer, double * RESTRICT  _data_pdfs_dst, int64_t const _size_pdfs_dst_0, int64_t const _size_pdfs_dst_1, int64_t const _size_pdfs_dst_2, int64_t const _stride_pdfs_dst_0, int64_t const _stride_pdfs_dst_1, int64_t const _stride_pdfs_dst_2, int64_t const _stride_pdfs_dst_3)
{
   for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_dst_0; ctr_0 += 1)
   {
      double * RESTRICT  _data_pdfs_dst_00_36 = _data_pdfs_dst + _stride_pdfs_dst_0*ctr_0 + 6*_stride_pdfs_dst_3;
      double * RESTRICT  _data_pdfs_dst_00_315 = _data_pdfs_dst + _stride_pdfs_dst_0*ctr_0 + 15*_stride_pdfs_dst_3;
      double * RESTRICT  _data_pdfs_dst_00_316 = _data_pdfs_dst + _stride_pdfs_dst_0*ctr_0 + 16*_stride_pdfs_dst_3;
      double * RESTRICT  _data_pdfs_dst_00_317 = _data_pdfs_dst + _stride_pdfs_dst_0*ctr_0 + 17*_stride_pdfs_dst_3;
      double * RESTRICT  _data_pdfs_dst_00_318 = _data_pdfs_dst + _stride_pdfs_dst_0*ctr_0 + 18*_stride_pdfs_dst_3;
      for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_dst_1; ctr_1 += 1)
      {
         double * RESTRICT  _data_pdfs_dst_00_36_10 = _stride_pdfs_dst_1*ctr_1 + _data_pdfs_dst_00_36;
         double * RESTRICT  _data_pdfs_dst_00_315_10 = _stride_pdfs_dst_1*ctr_1 + _data_pdfs_dst_00_315;
         double * RESTRICT  _data_pdfs_dst_00_316_10 = _stride_pdfs_dst_1*ctr_1 + _data_pdfs_dst_00_316;
         double * RESTRICT  _data_pdfs_dst_00_317_10 = _stride_pdfs_dst_1*ctr_1 + _data_pdfs_dst_00_317;
         double * RESTRICT  _data_pdfs_dst_00_318_10 = _stride_pdfs_dst_1*ctr_1 + _data_pdfs_dst_00_318;
         for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_dst_2; ctr_2 += 1)
         {
            _data_pdfs_dst_00_36_10[_stride_pdfs_dst_2*ctr_2] = _data_buffer[5*_size_pdfs_dst_1*_size_pdfs_dst_2*ctr_0 + 5*_size_pdfs_dst_2*ctr_1 + 5*ctr_2];
            _data_pdfs_dst_00_315_10[_stride_pdfs_dst_2*ctr_2] = _data_buffer[5*_size_pdfs_dst_1*_size_pdfs_dst_2*ctr_0 + 5*_size_pdfs_dst_2*ctr_1 + 5*ctr_2 + 1];
            _data_pdfs_dst_00_316_10[_stride_pdfs_dst_2*ctr_2] = _data_buffer[5*_size_pdfs_dst_1*_size_pdfs_dst_2*ctr_0 + 5*_size_pdfs_dst_2*ctr_1 + 5*ctr_2 + 2];
            _data_pdfs_dst_00_317_10[_stride_pdfs_dst_2*ctr_2] = _data_buffer[5*_size_pdfs_dst_1*_size_pdfs_dst_2*ctr_0 + 5*_size_pdfs_dst_2*ctr_1 + 5*ctr_2 + 3];
            _data_pdfs_dst_00_318_10[_stride_pdfs_dst_2*ctr_2] = _data_buffer[5*_size_pdfs_dst_1*_size_pdfs_dst_2*ctr_0 + 5*_size_pdfs_dst_2*ctr_1 + 5*ctr_2 + 4];
         }
      }
   }
}
}

namespace internal_d3q19storagespecification_unpack_TS {
static FUNC_PREFIX void d3q19storagespecification_unpack_TS(const double * RESTRICT const _data_buffer, double * RESTRICT  _data_pdfs_dst, int64_t const _size_pdfs_dst_0, int64_t const _size_pdfs_dst_1, int64_t const _size_pdfs_dst_2, int64_t const _stride_pdfs_dst_0, int64_t const _stride_pdfs_dst_1, int64_t const _stride_pdfs_dst_2, int64_t const _stride_pdfs_dst_3)
{
   for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_dst_0; ctr_0 += 1)
   {
      double * RESTRICT  _data_pdfs_dst_00_315 = _data_pdfs_dst + _stride_pdfs_dst_0*ctr_0 + 15*_stride_pdfs_dst_3;
      for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_dst_1; ctr_1 += 1)
      {
         double * RESTRICT  _data_pdfs_dst_00_315_10 = _stride_pdfs_dst_1*ctr_1 + _data_pdfs_dst_00_315;
         for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_dst_2; ctr_2 += 1)
         {
            _data_pdfs_dst_00_315_10[_stride_pdfs_dst_2*ctr_2] = _data_buffer[_size_pdfs_dst_1*_size_pdfs_dst_2*ctr_0 + _size_pdfs_dst_2*ctr_1 + ctr_2];
         }
      }
   }
}
}

namespace internal_d3q19storagespecification_unpack_BE {
static FUNC_PREFIX void d3q19storagespecification_unpack_BE(const double * RESTRICT const _data_buffer, double * RESTRICT  _data_pdfs_dst, int64_t const _size_pdfs_dst_0, int64_t const _size_pdfs_dst_1, int64_t const _size_pdfs_dst_2, int64_t const _stride_pdfs_dst_0, int64_t const _stride_pdfs_dst_1, int64_t const _stride_pdfs_dst_2, int64_t const _stride_pdfs_dst_3)
{
   for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_dst_0; ctr_0 += 1)
   {
      double * RESTRICT  _data_pdfs_dst_00_313 = _data_pdfs_dst + _stride_pdfs_dst_0*ctr_0 + 13*_stride_pdfs_dst_3;
      for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_dst_1; ctr_1 += 1)
      {
         double * RESTRICT  _data_pdfs_dst_00_313_10 = _stride_pdfs_dst_1*ctr_1 + _data_pdfs_dst_00_313;
         for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_dst_2; ctr_2 += 1)
         {
            _data_pdfs_dst_00_313_10[_stride_pdfs_dst_2*ctr_2] = _data_buffer[_size_pdfs_dst_1*_size_pdfs_dst_2*ctr_0 + _size_pdfs_dst_2*ctr_1 + ctr_2];
         }
      }
   }
}
}

namespace internal_d3q19storagespecification_unpack_NW {
static FUNC_PREFIX void d3q19storagespecification_unpack_NW(const double * RESTRICT const _data_buffer, double * RESTRICT  _data_pdfs_dst, int64_t const _size_pdfs_dst_0, int64_t const _size_pdfs_dst_1, int64_t const _size_pdfs_dst_2, int64_t const _stride_pdfs_dst_0, int64_t const _stride_pdfs_dst_1, int64_t const _stride_pdfs_dst_2, int64_t const _stride_pdfs_dst_3)
{
   for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_dst_0; ctr_0 += 1)
   {
      double * RESTRICT  _data_pdfs_dst_00_310 = _data_pdfs_dst + _stride_pdfs_dst_0*ctr_0 + 10*_stride_pdfs_dst_3;
      for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_dst_1; ctr_1 += 1)
      {
         double * RESTRICT  _data_pdfs_dst_00_310_10 = _stride_pdfs_dst_1*ctr_1 + _data_pdfs_dst_00_310;
         for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_dst_2; ctr_2 += 1)
         {
            _data_pdfs_dst_00_310_10[_stride_pdfs_dst_2*ctr_2] = _data_buffer[_size_pdfs_dst_1*_size_pdfs_dst_2*ctr_0 + _size_pdfs_dst_2*ctr_1 + ctr_2];
         }
      }
   }
}
}

namespace internal_d3q19storagespecification_unpack_NE {
static FUNC_PREFIX void d3q19storagespecification_unpack_NE(const double * RESTRICT const _data_buffer, double * RESTRICT  _data_pdfs_dst, int64_t const _size_pdfs_dst_0, int64_t const _size_pdfs_dst_1, int64_t const _size_pdfs_dst_2, int64_t const _stride_pdfs_dst_0, int64_t const _stride_pdfs_dst_1, int64_t const _stride_pdfs_dst_2, int64_t const _stride_pdfs_dst_3)
{
   for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_dst_0; ctr_0 += 1)
   {
      double * RESTRICT  _data_pdfs_dst_00_39 = _data_pdfs_dst + _stride_pdfs_dst_0*ctr_0 + 9*_stride_pdfs_dst_3;
      for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_dst_1; ctr_1 += 1)
      {
         double * RESTRICT  _data_pdfs_dst_00_39_10 = _stride_pdfs_dst_1*ctr_1 + _data_pdfs_dst_00_39;
         for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_dst_2; ctr_2 += 1)
         {
            _data_pdfs_dst_00_39_10[_stride_pdfs_dst_2*ctr_2] = _data_buffer[_size_pdfs_dst_1*_size_pdfs_dst_2*ctr_0 + _size_pdfs_dst_2*ctr_1 + ctr_2];
         }
      }
   }
}
}

namespace internal_d3q19storagespecification_unpack_BS {
static FUNC_PREFIX void d3q19storagespecification_unpack_BS(const double * RESTRICT const _data_buffer, double * RESTRICT  _data_pdfs_dst, int64_t const _size_pdfs_dst_0, int64_t const _size_pdfs_dst_1, int64_t const _size_pdfs_dst_2, int64_t const _stride_pdfs_dst_0, int64_t const _stride_pdfs_dst_1, int64_t const _stride_pdfs_dst_2, int64_t const _stride_pdfs_dst_3)
{
   for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_dst_0; ctr_0 += 1)
   {
      double * RESTRICT  _data_pdfs_dst_00_311 = _data_pdfs_dst + _stride_pdfs_dst_0*ctr_0 + 11*_stride_pdfs_dst_3;
      for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_dst_1; ctr_1 += 1)
      {
         double * RESTRICT  _data_pdfs_dst_00_311_10 = _stride_pdfs_dst_1*ctr_1 + _data_pdfs_dst_00_311;
         for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_dst_2; ctr_2 += 1)
         {
            _data_pdfs_dst_00_311_10[_stride_pdfs_dst_2*ctr_2] = _data_buffer[_size_pdfs_dst_1*_size_pdfs_dst_2*ctr_0 + _size_pdfs_dst_2*ctr_1 + ctr_2];
         }
      }
   }
}
}

namespace internal_d3q19storagespecification_unpack_E {
static FUNC_PREFIX void d3q19storagespecification_unpack_E(const double * RESTRICT const _data_buffer, double * RESTRICT  _data_pdfs_dst, int64_t const _size_pdfs_dst_0, int64_t const _size_pdfs_dst_1, int64_t const _size_pdfs_dst_2, int64_t const _stride_pdfs_dst_0, int64_t const _stride_pdfs_dst_1, int64_t const _stride_pdfs_dst_2, int64_t const _stride_pdfs_dst_3)
{
   for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_dst_0; ctr_0 += 1)
   {
      double * RESTRICT  _data_pdfs_dst_00_33 = _data_pdfs_dst + _stride_pdfs_dst_0*ctr_0 + 3*_stride_pdfs_dst_3;
      double * RESTRICT  _data_pdfs_dst_00_37 = _data_pdfs_dst + _stride_pdfs_dst_0*ctr_0 + 7*_stride_pdfs_dst_3;
      double * RESTRICT  _data_pdfs_dst_00_39 = _data_pdfs_dst + _stride_pdfs_dst_0*ctr_0 + 9*_stride_pdfs_dst_3;
      double * RESTRICT  _data_pdfs_dst_00_313 = _data_pdfs_dst + _stride_pdfs_dst_0*ctr_0 + 13*_stride_pdfs_dst_3;
      double * RESTRICT  _data_pdfs_dst_00_317 = _data_pdfs_dst + _stride_pdfs_dst_0*ctr_0 + 17*_stride_pdfs_dst_3;
      for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_dst_1; ctr_1 += 1)
      {
         double * RESTRICT  _data_pdfs_dst_00_33_10 = _stride_pdfs_dst_1*ctr_1 + _data_pdfs_dst_00_33;
         double * RESTRICT  _data_pdfs_dst_00_37_10 = _stride_pdfs_dst_1*ctr_1 + _data_pdfs_dst_00_37;
         double * RESTRICT  _data_pdfs_dst_00_39_10 = _stride_pdfs_dst_1*ctr_1 + _data_pdfs_dst_00_39;
         double * RESTRICT  _data_pdfs_dst_00_313_10 = _stride_pdfs_dst_1*ctr_1 + _data_pdfs_dst_00_313;
         double * RESTRICT  _data_pdfs_dst_00_317_10 = _stride_pdfs_dst_1*ctr_1 + _data_pdfs_dst_00_317;
         for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_dst_2; ctr_2 += 1)
         {
            _data_pdfs_dst_00_33_10[_stride_pdfs_dst_2*ctr_2] = _data_buffer[5*_size_pdfs_dst_1*_size_pdfs_dst_2*ctr_0 + 5*_size_pdfs_dst_2*ctr_1 + 5*ctr_2];
            _data_pdfs_dst_00_37_10[_stride_pdfs_dst_2*ctr_2] = _data_buffer[5*_size_pdfs_dst_1*_size_pdfs_dst_2*ctr_0 + 5*_size_pdfs_dst_2*ctr_1 + 5*ctr_2 + 1];
            _data_pdfs_dst_00_39_10[_stride_pdfs_dst_2*ctr_2] = _data_buffer[5*_size_pdfs_dst_1*_size_pdfs_dst_2*ctr_0 + 5*_size_pdfs_dst_2*ctr_1 + 5*ctr_2 + 2];
            _data_pdfs_dst_00_313_10[_stride_pdfs_dst_2*ctr_2] = _data_buffer[5*_size_pdfs_dst_1*_size_pdfs_dst_2*ctr_0 + 5*_size_pdfs_dst_2*ctr_1 + 5*ctr_2 + 3];
            _data_pdfs_dst_00_317_10[_stride_pdfs_dst_2*ctr_2] = _data_buffer[5*_size_pdfs_dst_1*_size_pdfs_dst_2*ctr_0 + 5*_size_pdfs_dst_2*ctr_1 + 5*ctr_2 + 4];
         }
      }
   }
}
}

namespace internal_d3q19storagespecification_unpack_S {
static FUNC_PREFIX void d3q19storagespecification_unpack_S(const double * RESTRICT const _data_buffer, double * RESTRICT  _data_pdfs_dst, int64_t const _size_pdfs_dst_0, int64_t const _size_pdfs_dst_1, int64_t const _size_pdfs_dst_2, int64_t const _stride_pdfs_dst_0, int64_t const _stride_pdfs_dst_1, int64_t const _stride_pdfs_dst_2, int64_t const _stride_pdfs_dst_3)
{
   for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_dst_0; ctr_0 += 1)
   {
      double * RESTRICT  _data_pdfs_dst_00_31 = _data_pdfs_dst + _stride_pdfs_dst_0*ctr_0 + _stride_pdfs_dst_3;
      double * RESTRICT  _data_pdfs_dst_00_37 = _data_pdfs_dst + _stride_pdfs_dst_0*ctr_0 + 7*_stride_pdfs_dst_3;
      double * RESTRICT  _data_pdfs_dst_00_38 = _data_pdfs_dst + _stride_pdfs_dst_0*ctr_0 + 8*_stride_pdfs_dst_3;
      double * RESTRICT  _data_pdfs_dst_00_311 = _data_pdfs_dst + _stride_pdfs_dst_0*ctr_0 + 11*_stride_pdfs_dst_3;
      double * RESTRICT  _data_pdfs_dst_00_315 = _data_pdfs_dst + _stride_pdfs_dst_0*ctr_0 + 15*_stride_pdfs_dst_3;
      for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_dst_1; ctr_1 += 1)
      {
         double * RESTRICT  _data_pdfs_dst_00_31_10 = _stride_pdfs_dst_1*ctr_1 + _data_pdfs_dst_00_31;
         double * RESTRICT  _data_pdfs_dst_00_37_10 = _stride_pdfs_dst_1*ctr_1 + _data_pdfs_dst_00_37;
         double * RESTRICT  _data_pdfs_dst_00_38_10 = _stride_pdfs_dst_1*ctr_1 + _data_pdfs_dst_00_38;
         double * RESTRICT  _data_pdfs_dst_00_311_10 = _stride_pdfs_dst_1*ctr_1 + _data_pdfs_dst_00_311;
         double * RESTRICT  _data_pdfs_dst_00_315_10 = _stride_pdfs_dst_1*ctr_1 + _data_pdfs_dst_00_315;
         for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_dst_2; ctr_2 += 1)
         {
            _data_pdfs_dst_00_31_10[_stride_pdfs_dst_2*ctr_2] = _data_buffer[5*_size_pdfs_dst_1*_size_pdfs_dst_2*ctr_0 + 5*_size_pdfs_dst_2*ctr_1 + 5*ctr_2];
            _data_pdfs_dst_00_37_10[_stride_pdfs_dst_2*ctr_2] = _data_buffer[5*_size_pdfs_dst_1*_size_pdfs_dst_2*ctr_0 + 5*_size_pdfs_dst_2*ctr_1 + 5*ctr_2 + 1];
            _data_pdfs_dst_00_38_10[_stride_pdfs_dst_2*ctr_2] = _data_buffer[5*_size_pdfs_dst_1*_size_pdfs_dst_2*ctr_0 + 5*_size_pdfs_dst_2*ctr_1 + 5*ctr_2 + 2];
            _data_pdfs_dst_00_311_10[_stride_pdfs_dst_2*ctr_2] = _data_buffer[5*_size_pdfs_dst_1*_size_pdfs_dst_2*ctr_0 + 5*_size_pdfs_dst_2*ctr_1 + 5*ctr_2 + 3];
            _data_pdfs_dst_00_315_10[_stride_pdfs_dst_2*ctr_2] = _data_buffer[5*_size_pdfs_dst_1*_size_pdfs_dst_2*ctr_0 + 5*_size_pdfs_dst_2*ctr_1 + 5*ctr_2 + 4];
         }
      }
   }
}
}

namespace internal_d3q19storagespecification_unpack_W {
static FUNC_PREFIX void d3q19storagespecification_unpack_W(const double * RESTRICT const _data_buffer, double * RESTRICT  _data_pdfs_dst, int64_t const _size_pdfs_dst_0, int64_t const _size_pdfs_dst_1, int64_t const _size_pdfs_dst_2, int64_t const _stride_pdfs_dst_0, int64_t const _stride_pdfs_dst_1, int64_t const _stride_pdfs_dst_2, int64_t const _stride_pdfs_dst_3)
{
   for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_dst_0; ctr_0 += 1)
   {
      double * RESTRICT  _data_pdfs_dst_00_34 = _data_pdfs_dst + _stride_pdfs_dst_0*ctr_0 + 4*_stride_pdfs_dst_3;
      double * RESTRICT  _data_pdfs_dst_00_38 = _data_pdfs_dst + _stride_pdfs_dst_0*ctr_0 + 8*_stride_pdfs_dst_3;
      double * RESTRICT  _data_pdfs_dst_00_310 = _data_pdfs_dst + _stride_pdfs_dst_0*ctr_0 + 10*_stride_pdfs_dst_3;
      double * RESTRICT  _data_pdfs_dst_00_314 = _data_pdfs_dst + _stride_pdfs_dst_0*ctr_0 + 14*_stride_pdfs_dst_3;
      double * RESTRICT  _data_pdfs_dst_00_318 = _data_pdfs_dst + _stride_pdfs_dst_0*ctr_0 + 18*_stride_pdfs_dst_3;
      for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_dst_1; ctr_1 += 1)
      {
         double * RESTRICT  _data_pdfs_dst_00_34_10 = _stride_pdfs_dst_1*ctr_1 + _data_pdfs_dst_00_34;
         double * RESTRICT  _data_pdfs_dst_00_38_10 = _stride_pdfs_dst_1*ctr_1 + _data_pdfs_dst_00_38;
         double * RESTRICT  _data_pdfs_dst_00_310_10 = _stride_pdfs_dst_1*ctr_1 + _data_pdfs_dst_00_310;
         double * RESTRICT  _data_pdfs_dst_00_314_10 = _stride_pdfs_dst_1*ctr_1 + _data_pdfs_dst_00_314;
         double * RESTRICT  _data_pdfs_dst_00_318_10 = _stride_pdfs_dst_1*ctr_1 + _data_pdfs_dst_00_318;
         for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_dst_2; ctr_2 += 1)
         {
            _data_pdfs_dst_00_34_10[_stride_pdfs_dst_2*ctr_2] = _data_buffer[5*_size_pdfs_dst_1*_size_pdfs_dst_2*ctr_0 + 5*_size_pdfs_dst_2*ctr_1 + 5*ctr_2];
            _data_pdfs_dst_00_38_10[_stride_pdfs_dst_2*ctr_2] = _data_buffer[5*_size_pdfs_dst_1*_size_pdfs_dst_2*ctr_0 + 5*_size_pdfs_dst_2*ctr_1 + 5*ctr_2 + 1];
            _data_pdfs_dst_00_310_10[_stride_pdfs_dst_2*ctr_2] = _data_buffer[5*_size_pdfs_dst_1*_size_pdfs_dst_2*ctr_0 + 5*_size_pdfs_dst_2*ctr_1 + 5*ctr_2 + 2];
            _data_pdfs_dst_00_314_10[_stride_pdfs_dst_2*ctr_2] = _data_buffer[5*_size_pdfs_dst_1*_size_pdfs_dst_2*ctr_0 + 5*_size_pdfs_dst_2*ctr_1 + 5*ctr_2 + 3];
            _data_pdfs_dst_00_318_10[_stride_pdfs_dst_2*ctr_2] = _data_buffer[5*_size_pdfs_dst_1*_size_pdfs_dst_2*ctr_0 + 5*_size_pdfs_dst_2*ctr_1 + 5*ctr_2 + 4];
         }
      }
   }
}
}

namespace internal_d3q19storagespecification_unpack_SW {
static FUNC_PREFIX void d3q19storagespecification_unpack_SW(const double * RESTRICT const _data_buffer, double * RESTRICT  _data_pdfs_dst, int64_t const _size_pdfs_dst_0, int64_t const _size_pdfs_dst_1, int64_t const _size_pdfs_dst_2, int64_t const _stride_pdfs_dst_0, int64_t const _stride_pdfs_dst_1, int64_t const _stride_pdfs_dst_2, int64_t const _stride_pdfs_dst_3)
{
   for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_dst_0; ctr_0 += 1)
   {
      double * RESTRICT  _data_pdfs_dst_00_38 = _data_pdfs_dst + _stride_pdfs_dst_0*ctr_0 + 8*_stride_pdfs_dst_3;
      for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_dst_1; ctr_1 += 1)
      {
         double * RESTRICT  _data_pdfs_dst_00_38_10 = _stride_pdfs_dst_1*ctr_1 + _data_pdfs_dst_00_38;
         for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_dst_2; ctr_2 += 1)
         {
            _data_pdfs_dst_00_38_10[_stride_pdfs_dst_2*ctr_2] = _data_buffer[_size_pdfs_dst_1*_size_pdfs_dst_2*ctr_0 + _size_pdfs_dst_2*ctr_1 + ctr_2];
         }
      }
   }
}
}

namespace internal_d3q19storagespecification_unpack_B {
static FUNC_PREFIX void d3q19storagespecification_unpack_B(const double * RESTRICT const _data_buffer, double * RESTRICT  _data_pdfs_dst, int64_t const _size_pdfs_dst_0, int64_t const _size_pdfs_dst_1, int64_t const _size_pdfs_dst_2, int64_t const _stride_pdfs_dst_0, int64_t const _stride_pdfs_dst_1, int64_t const _stride_pdfs_dst_2, int64_t const _stride_pdfs_dst_3)
{
   for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_dst_0; ctr_0 += 1)
   {
      double * RESTRICT  _data_pdfs_dst_00_35 = _data_pdfs_dst + _stride_pdfs_dst_0*ctr_0 + 5*_stride_pdfs_dst_3;
      double * RESTRICT  _data_pdfs_dst_00_311 = _data_pdfs_dst + _stride_pdfs_dst_0*ctr_0 + 11*_stride_pdfs_dst_3;
      double * RESTRICT  _data_pdfs_dst_00_312 = _data_pdfs_dst + _stride_pdfs_dst_0*ctr_0 + 12*_stride_pdfs_dst_3;
      double * RESTRICT  _data_pdfs_dst_00_313 = _data_pdfs_dst + _stride_pdfs_dst_0*ctr_0 + 13*_stride_pdfs_dst_3;
      double * RESTRICT  _data_pdfs_dst_00_314 = _data_pdfs_dst + _stride_pdfs_dst_0*ctr_0 + 14*_stride_pdfs_dst_3;
      for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_dst_1; ctr_1 += 1)
      {
         double * RESTRICT  _data_pdfs_dst_00_35_10 = _stride_pdfs_dst_1*ctr_1 + _data_pdfs_dst_00_35;
         double * RESTRICT  _data_pdfs_dst_00_311_10 = _stride_pdfs_dst_1*ctr_1 + _data_pdfs_dst_00_311;
         double * RESTRICT  _data_pdfs_dst_00_312_10 = _stride_pdfs_dst_1*ctr_1 + _data_pdfs_dst_00_312;
         double * RESTRICT  _data_pdfs_dst_00_313_10 = _stride_pdfs_dst_1*ctr_1 + _data_pdfs_dst_00_313;
         double * RESTRICT  _data_pdfs_dst_00_314_10 = _stride_pdfs_dst_1*ctr_1 + _data_pdfs_dst_00_314;
         for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_dst_2; ctr_2 += 1)
         {
            _data_pdfs_dst_00_35_10[_stride_pdfs_dst_2*ctr_2] = _data_buffer[5*_size_pdfs_dst_1*_size_pdfs_dst_2*ctr_0 + 5*_size_pdfs_dst_2*ctr_1 + 5*ctr_2];
            _data_pdfs_dst_00_311_10[_stride_pdfs_dst_2*ctr_2] = _data_buffer[5*_size_pdfs_dst_1*_size_pdfs_dst_2*ctr_0 + 5*_size_pdfs_dst_2*ctr_1 + 5*ctr_2 + 1];
            _data_pdfs_dst_00_312_10[_stride_pdfs_dst_2*ctr_2] = _data_buffer[5*_size_pdfs_dst_1*_size_pdfs_dst_2*ctr_0 + 5*_size_pdfs_dst_2*ctr_1 + 5*ctr_2 + 2];
            _data_pdfs_dst_00_313_10[_stride_pdfs_dst_2*ctr_2] = _data_buffer[5*_size_pdfs_dst_1*_size_pdfs_dst_2*ctr_0 + 5*_size_pdfs_dst_2*ctr_1 + 5*ctr_2 + 3];
            _data_pdfs_dst_00_314_10[_stride_pdfs_dst_2*ctr_2] = _data_buffer[5*_size_pdfs_dst_1*_size_pdfs_dst_2*ctr_0 + 5*_size_pdfs_dst_2*ctr_1 + 5*ctr_2 + 4];
         }
      }
   }
}
}

namespace internal_d3q19storagespecification_unpack_TN {
static FUNC_PREFIX void d3q19storagespecification_unpack_TN(const double * RESTRICT const _data_buffer, double * RESTRICT  _data_pdfs_dst, int64_t const _size_pdfs_dst_0, int64_t const _size_pdfs_dst_1, int64_t const _size_pdfs_dst_2, int64_t const _stride_pdfs_dst_0, int64_t const _stride_pdfs_dst_1, int64_t const _stride_pdfs_dst_2, int64_t const _stride_pdfs_dst_3)
{
   for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_dst_0; ctr_0 += 1)
   {
      double * RESTRICT  _data_pdfs_dst_00_316 = _data_pdfs_dst + _stride_pdfs_dst_0*ctr_0 + 16*_stride_pdfs_dst_3;
      for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_dst_1; ctr_1 += 1)
      {
         double * RESTRICT  _data_pdfs_dst_00_316_10 = _stride_pdfs_dst_1*ctr_1 + _data_pdfs_dst_00_316;
         for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_dst_2; ctr_2 += 1)
         {
            _data_pdfs_dst_00_316_10[_stride_pdfs_dst_2*ctr_2] = _data_buffer[_size_pdfs_dst_1*_size_pdfs_dst_2*ctr_0 + _size_pdfs_dst_2*ctr_1 + ctr_2];
         }
      }
   }
}
}

namespace internal_d3q19storagespecification_unpack_BN {
static FUNC_PREFIX void d3q19storagespecification_unpack_BN(const double * RESTRICT const _data_buffer, double * RESTRICT  _data_pdfs_dst, int64_t const _size_pdfs_dst_0, int64_t const _size_pdfs_dst_1, int64_t const _size_pdfs_dst_2, int64_t const _stride_pdfs_dst_0, int64_t const _stride_pdfs_dst_1, int64_t const _stride_pdfs_dst_2, int64_t const _stride_pdfs_dst_3)
{
   for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_dst_0; ctr_0 += 1)
   {
      double * RESTRICT  _data_pdfs_dst_00_312 = _data_pdfs_dst + _stride_pdfs_dst_0*ctr_0 + 12*_stride_pdfs_dst_3;
      for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_dst_1; ctr_1 += 1)
      {
         double * RESTRICT  _data_pdfs_dst_00_312_10 = _stride_pdfs_dst_1*ctr_1 + _data_pdfs_dst_00_312;
         for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_dst_2; ctr_2 += 1)
         {
            _data_pdfs_dst_00_312_10[_stride_pdfs_dst_2*ctr_2] = _data_buffer[_size_pdfs_dst_1*_size_pdfs_dst_2*ctr_0 + _size_pdfs_dst_2*ctr_1 + ctr_2];
         }
      }
   }
}
}

namespace internal_d3q19storagespecification_unpack_TW {
static FUNC_PREFIX void d3q19storagespecification_unpack_TW(const double * RESTRICT const _data_buffer, double * RESTRICT  _data_pdfs_dst, int64_t const _size_pdfs_dst_0, int64_t const _size_pdfs_dst_1, int64_t const _size_pdfs_dst_2, int64_t const _stride_pdfs_dst_0, int64_t const _stride_pdfs_dst_1, int64_t const _stride_pdfs_dst_2, int64_t const _stride_pdfs_dst_3)
{
   for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_dst_0; ctr_0 += 1)
   {
      double * RESTRICT  _data_pdfs_dst_00_318 = _data_pdfs_dst + _stride_pdfs_dst_0*ctr_0 + 18*_stride_pdfs_dst_3;
      for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_dst_1; ctr_1 += 1)
      {
         double * RESTRICT  _data_pdfs_dst_00_318_10 = _stride_pdfs_dst_1*ctr_1 + _data_pdfs_dst_00_318;
         for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_dst_2; ctr_2 += 1)
         {
            _data_pdfs_dst_00_318_10[_stride_pdfs_dst_2*ctr_2] = _data_buffer[_size_pdfs_dst_1*_size_pdfs_dst_2*ctr_0 + _size_pdfs_dst_2*ctr_1 + ctr_2];
         }
      }
   }
}
}

namespace internal_d3q19storagespecification_localCopy_NE {
static FUNC_PREFIX void d3q19storagespecification_localCopy_NE(double * RESTRICT  _data_pdfs_dst, double * RESTRICT const _data_pdfs_src, int64_t const _size_pdfs_dst_0, int64_t const _size_pdfs_dst_1, int64_t const _size_pdfs_dst_2, int64_t const _stride_pdfs_dst_0, int64_t const _stride_pdfs_dst_1, int64_t const _stride_pdfs_dst_2, int64_t const _stride_pdfs_dst_3, int64_t const _stride_pdfs_src_0, int64_t const _stride_pdfs_src_1, int64_t const _stride_pdfs_src_2, int64_t const _stride_pdfs_src_3)
{
   for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_dst_0; ctr_0 += 1)
   {
      double * RESTRICT  _data_pdfs_dst_00_38 = _data_pdfs_dst + _stride_pdfs_dst_0*ctr_0 + 8*_stride_pdfs_dst_3;
      double * RESTRICT _data_pdfs_src_00_38 = _data_pdfs_src + _stride_pdfs_src_0*ctr_0 + 8*_stride_pdfs_src_3;
      for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_dst_1; ctr_1 += 1)
      {
         double * RESTRICT  _data_pdfs_dst_00_38_10 = _stride_pdfs_dst_1*ctr_1 + _data_pdfs_dst_00_38;
         double * RESTRICT _data_pdfs_src_00_38_10 = _stride_pdfs_src_1*ctr_1 + _data_pdfs_src_00_38;
         for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_dst_2; ctr_2 += 1)
         {
            _data_pdfs_dst_00_38_10[_stride_pdfs_dst_2*ctr_2] = _data_pdfs_src_00_38_10[_stride_pdfs_src_2*ctr_2];
         }
      }
   }
}
}

namespace internal_d3q19storagespecification_localCopy_TS {
static FUNC_PREFIX void d3q19storagespecification_localCopy_TS(double * RESTRICT  _data_pdfs_dst, double * RESTRICT const _data_pdfs_src, int64_t const _size_pdfs_dst_0, int64_t const _size_pdfs_dst_1, int64_t const _size_pdfs_dst_2, int64_t const _stride_pdfs_dst_0, int64_t const _stride_pdfs_dst_1, int64_t const _stride_pdfs_dst_2, int64_t const _stride_pdfs_dst_3, int64_t const _stride_pdfs_src_0, int64_t const _stride_pdfs_src_1, int64_t const _stride_pdfs_src_2, int64_t const _stride_pdfs_src_3)
{
   for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_dst_0; ctr_0 += 1)
   {
      double * RESTRICT  _data_pdfs_dst_00_312 = _data_pdfs_dst + _stride_pdfs_dst_0*ctr_0 + 12*_stride_pdfs_dst_3;
      double * RESTRICT _data_pdfs_src_00_312 = _data_pdfs_src + _stride_pdfs_src_0*ctr_0 + 12*_stride_pdfs_src_3;
      for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_dst_1; ctr_1 += 1)
      {
         double * RESTRICT  _data_pdfs_dst_00_312_10 = _stride_pdfs_dst_1*ctr_1 + _data_pdfs_dst_00_312;
         double * RESTRICT _data_pdfs_src_00_312_10 = _stride_pdfs_src_1*ctr_1 + _data_pdfs_src_00_312;
         for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_dst_2; ctr_2 += 1)
         {
            _data_pdfs_dst_00_312_10[_stride_pdfs_dst_2*ctr_2] = _data_pdfs_src_00_312_10[_stride_pdfs_src_2*ctr_2];
         }
      }
   }
}
}

namespace internal_d3q19storagespecification_localCopy_BE {
static FUNC_PREFIX void d3q19storagespecification_localCopy_BE(double * RESTRICT  _data_pdfs_dst, double * RESTRICT const _data_pdfs_src, int64_t const _size_pdfs_dst_0, int64_t const _size_pdfs_dst_1, int64_t const _size_pdfs_dst_2, int64_t const _stride_pdfs_dst_0, int64_t const _stride_pdfs_dst_1, int64_t const _stride_pdfs_dst_2, int64_t const _stride_pdfs_dst_3, int64_t const _stride_pdfs_src_0, int64_t const _stride_pdfs_src_1, int64_t const _stride_pdfs_src_2, int64_t const _stride_pdfs_src_3)
{
   for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_dst_0; ctr_0 += 1)
   {
      double * RESTRICT  _data_pdfs_dst_00_318 = _data_pdfs_dst + _stride_pdfs_dst_0*ctr_0 + 18*_stride_pdfs_dst_3;
      double * RESTRICT _data_pdfs_src_00_318 = _data_pdfs_src + _stride_pdfs_src_0*ctr_0 + 18*_stride_pdfs_src_3;
      for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_dst_1; ctr_1 += 1)
      {
         double * RESTRICT  _data_pdfs_dst_00_318_10 = _stride_pdfs_dst_1*ctr_1 + _data_pdfs_dst_00_318;
         double * RESTRICT _data_pdfs_src_00_318_10 = _stride_pdfs_src_1*ctr_1 + _data_pdfs_src_00_318;
         for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_dst_2; ctr_2 += 1)
         {
            _data_pdfs_dst_00_318_10[_stride_pdfs_dst_2*ctr_2] = _data_pdfs_src_00_318_10[_stride_pdfs_src_2*ctr_2];
         }
      }
   }
}
}

namespace internal_d3q19storagespecification_localCopy_BS {
static FUNC_PREFIX void d3q19storagespecification_localCopy_BS(double * RESTRICT  _data_pdfs_dst, double * RESTRICT const _data_pdfs_src, int64_t const _size_pdfs_dst_0, int64_t const _size_pdfs_dst_1, int64_t const _size_pdfs_dst_2, int64_t const _stride_pdfs_dst_0, int64_t const _stride_pdfs_dst_1, int64_t const _stride_pdfs_dst_2, int64_t const _stride_pdfs_dst_3, int64_t const _stride_pdfs_src_0, int64_t const _stride_pdfs_src_1, int64_t const _stride_pdfs_src_2, int64_t const _stride_pdfs_src_3)
{
   for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_dst_0; ctr_0 += 1)
   {
      double * RESTRICT  _data_pdfs_dst_00_316 = _data_pdfs_dst + _stride_pdfs_dst_0*ctr_0 + 16*_stride_pdfs_dst_3;
      double * RESTRICT _data_pdfs_src_00_316 = _data_pdfs_src + _stride_pdfs_src_0*ctr_0 + 16*_stride_pdfs_src_3;
      for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_dst_1; ctr_1 += 1)
      {
         double * RESTRICT  _data_pdfs_dst_00_316_10 = _stride_pdfs_dst_1*ctr_1 + _data_pdfs_dst_00_316;
         double * RESTRICT _data_pdfs_src_00_316_10 = _stride_pdfs_src_1*ctr_1 + _data_pdfs_src_00_316;
         for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_dst_2; ctr_2 += 1)
         {
            _data_pdfs_dst_00_316_10[_stride_pdfs_dst_2*ctr_2] = _data_pdfs_src_00_316_10[_stride_pdfs_src_2*ctr_2];
         }
      }
   }
}
}

namespace internal_d3q19storagespecification_localCopy_BW {
static FUNC_PREFIX void d3q19storagespecification_localCopy_BW(double * RESTRICT  _data_pdfs_dst, double * RESTRICT const _data_pdfs_src, int64_t const _size_pdfs_dst_0, int64_t const _size_pdfs_dst_1, int64_t const _size_pdfs_dst_2, int64_t const _stride_pdfs_dst_0, int64_t const _stride_pdfs_dst_1, int64_t const _stride_pdfs_dst_2, int64_t const _stride_pdfs_dst_3, int64_t const _stride_pdfs_src_0, int64_t const _stride_pdfs_src_1, int64_t const _stride_pdfs_src_2, int64_t const _stride_pdfs_src_3)
{
   for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_dst_0; ctr_0 += 1)
   {
      double * RESTRICT  _data_pdfs_dst_00_317 = _data_pdfs_dst + _stride_pdfs_dst_0*ctr_0 + 17*_stride_pdfs_dst_3;
      double * RESTRICT _data_pdfs_src_00_317 = _data_pdfs_src + _stride_pdfs_src_0*ctr_0 + 17*_stride_pdfs_src_3;
      for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_dst_1; ctr_1 += 1)
      {
         double * RESTRICT  _data_pdfs_dst_00_317_10 = _stride_pdfs_dst_1*ctr_1 + _data_pdfs_dst_00_317;
         double * RESTRICT _data_pdfs_src_00_317_10 = _stride_pdfs_src_1*ctr_1 + _data_pdfs_src_00_317;
         for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_dst_2; ctr_2 += 1)
         {
            _data_pdfs_dst_00_317_10[_stride_pdfs_dst_2*ctr_2] = _data_pdfs_src_00_317_10[_stride_pdfs_src_2*ctr_2];
         }
      }
   }
}
}

namespace internal_d3q19storagespecification_localCopy_T {
static FUNC_PREFIX void d3q19storagespecification_localCopy_T(double * RESTRICT  _data_pdfs_dst, double * RESTRICT const _data_pdfs_src, int64_t const _size_pdfs_dst_0, int64_t const _size_pdfs_dst_1, int64_t const _size_pdfs_dst_2, int64_t const _stride_pdfs_dst_0, int64_t const _stride_pdfs_dst_1, int64_t const _stride_pdfs_dst_2, int64_t const _stride_pdfs_dst_3, int64_t const _stride_pdfs_src_0, int64_t const _stride_pdfs_src_1, int64_t const _stride_pdfs_src_2, int64_t const _stride_pdfs_src_3)
{
   for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_dst_0; ctr_0 += 1)
   {
      double * RESTRICT  _data_pdfs_dst_00_35 = _data_pdfs_dst + _stride_pdfs_dst_0*ctr_0 + 5*_stride_pdfs_dst_3;
      double * RESTRICT _data_pdfs_src_00_35 = _data_pdfs_src + _stride_pdfs_src_0*ctr_0 + 5*_stride_pdfs_src_3;
      double * RESTRICT  _data_pdfs_dst_00_311 = _data_pdfs_dst + _stride_pdfs_dst_0*ctr_0 + 11*_stride_pdfs_dst_3;
      double * RESTRICT _data_pdfs_src_00_311 = _data_pdfs_src + _stride_pdfs_src_0*ctr_0 + 11*_stride_pdfs_src_3;
      double * RESTRICT  _data_pdfs_dst_00_312 = _data_pdfs_dst + _stride_pdfs_dst_0*ctr_0 + 12*_stride_pdfs_dst_3;
      double * RESTRICT _data_pdfs_src_00_312 = _data_pdfs_src + _stride_pdfs_src_0*ctr_0 + 12*_stride_pdfs_src_3;
      double * RESTRICT  _data_pdfs_dst_00_313 = _data_pdfs_dst + _stride_pdfs_dst_0*ctr_0 + 13*_stride_pdfs_dst_3;
      double * RESTRICT _data_pdfs_src_00_313 = _data_pdfs_src + _stride_pdfs_src_0*ctr_0 + 13*_stride_pdfs_src_3;
      double * RESTRICT  _data_pdfs_dst_00_314 = _data_pdfs_dst + _stride_pdfs_dst_0*ctr_0 + 14*_stride_pdfs_dst_3;
      double * RESTRICT _data_pdfs_src_00_314 = _data_pdfs_src + _stride_pdfs_src_0*ctr_0 + 14*_stride_pdfs_src_3;
      for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_dst_1; ctr_1 += 1)
      {
         double * RESTRICT  _data_pdfs_dst_00_35_10 = _stride_pdfs_dst_1*ctr_1 + _data_pdfs_dst_00_35;
         double * RESTRICT _data_pdfs_src_00_35_10 = _stride_pdfs_src_1*ctr_1 + _data_pdfs_src_00_35;
         double * RESTRICT  _data_pdfs_dst_00_311_10 = _stride_pdfs_dst_1*ctr_1 + _data_pdfs_dst_00_311;
         double * RESTRICT _data_pdfs_src_00_311_10 = _stride_pdfs_src_1*ctr_1 + _data_pdfs_src_00_311;
         double * RESTRICT  _data_pdfs_dst_00_312_10 = _stride_pdfs_dst_1*ctr_1 + _data_pdfs_dst_00_312;
         double * RESTRICT _data_pdfs_src_00_312_10 = _stride_pdfs_src_1*ctr_1 + _data_pdfs_src_00_312;
         double * RESTRICT  _data_pdfs_dst_00_313_10 = _stride_pdfs_dst_1*ctr_1 + _data_pdfs_dst_00_313;
         double * RESTRICT _data_pdfs_src_00_313_10 = _stride_pdfs_src_1*ctr_1 + _data_pdfs_src_00_313;
         double * RESTRICT  _data_pdfs_dst_00_314_10 = _stride_pdfs_dst_1*ctr_1 + _data_pdfs_dst_00_314;
         double * RESTRICT _data_pdfs_src_00_314_10 = _stride_pdfs_src_1*ctr_1 + _data_pdfs_src_00_314;
         for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_dst_2; ctr_2 += 1)
         {
            _data_pdfs_dst_00_35_10[_stride_pdfs_dst_2*ctr_2] = _data_pdfs_src_00_35_10[_stride_pdfs_src_2*ctr_2];
            _data_pdfs_dst_00_311_10[_stride_pdfs_dst_2*ctr_2] = _data_pdfs_src_00_311_10[_stride_pdfs_src_2*ctr_2];
            _data_pdfs_dst_00_312_10[_stride_pdfs_dst_2*ctr_2] = _data_pdfs_src_00_312_10[_stride_pdfs_src_2*ctr_2];
            _data_pdfs_dst_00_313_10[_stride_pdfs_dst_2*ctr_2] = _data_pdfs_src_00_313_10[_stride_pdfs_src_2*ctr_2];
            _data_pdfs_dst_00_314_10[_stride_pdfs_dst_2*ctr_2] = _data_pdfs_src_00_314_10[_stride_pdfs_src_2*ctr_2];
         }
      }
   }
}
}

namespace internal_d3q19storagespecification_localCopy_TN {
static FUNC_PREFIX void d3q19storagespecification_localCopy_TN(double * RESTRICT  _data_pdfs_dst, double * RESTRICT const _data_pdfs_src, int64_t const _size_pdfs_dst_0, int64_t const _size_pdfs_dst_1, int64_t const _size_pdfs_dst_2, int64_t const _stride_pdfs_dst_0, int64_t const _stride_pdfs_dst_1, int64_t const _stride_pdfs_dst_2, int64_t const _stride_pdfs_dst_3, int64_t const _stride_pdfs_src_0, int64_t const _stride_pdfs_src_1, int64_t const _stride_pdfs_src_2, int64_t const _stride_pdfs_src_3)
{
   for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_dst_0; ctr_0 += 1)
   {
      double * RESTRICT  _data_pdfs_dst_00_311 = _data_pdfs_dst + _stride_pdfs_dst_0*ctr_0 + 11*_stride_pdfs_dst_3;
      double * RESTRICT _data_pdfs_src_00_311 = _data_pdfs_src + _stride_pdfs_src_0*ctr_0 + 11*_stride_pdfs_src_3;
      for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_dst_1; ctr_1 += 1)
      {
         double * RESTRICT  _data_pdfs_dst_00_311_10 = _stride_pdfs_dst_1*ctr_1 + _data_pdfs_dst_00_311;
         double * RESTRICT _data_pdfs_src_00_311_10 = _stride_pdfs_src_1*ctr_1 + _data_pdfs_src_00_311;
         for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_dst_2; ctr_2 += 1)
         {
            _data_pdfs_dst_00_311_10[_stride_pdfs_dst_2*ctr_2] = _data_pdfs_src_00_311_10[_stride_pdfs_src_2*ctr_2];
         }
      }
   }
}
}

namespace internal_d3q19storagespecification_localCopy_W {
static FUNC_PREFIX void d3q19storagespecification_localCopy_W(double * RESTRICT  _data_pdfs_dst, double * RESTRICT const _data_pdfs_src, int64_t const _size_pdfs_dst_0, int64_t const _size_pdfs_dst_1, int64_t const _size_pdfs_dst_2, int64_t const _stride_pdfs_dst_0, int64_t const _stride_pdfs_dst_1, int64_t const _stride_pdfs_dst_2, int64_t const _stride_pdfs_dst_3, int64_t const _stride_pdfs_src_0, int64_t const _stride_pdfs_src_1, int64_t const _stride_pdfs_src_2, int64_t const _stride_pdfs_src_3)
{
   for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_dst_0; ctr_0 += 1)
   {
      double * RESTRICT  _data_pdfs_dst_00_33 = _data_pdfs_dst + _stride_pdfs_dst_0*ctr_0 + 3*_stride_pdfs_dst_3;
      double * RESTRICT _data_pdfs_src_00_33 = _data_pdfs_src + _stride_pdfs_src_0*ctr_0 + 3*_stride_pdfs_src_3;
      double * RESTRICT  _data_pdfs_dst_00_37 = _data_pdfs_dst + _stride_pdfs_dst_0*ctr_0 + 7*_stride_pdfs_dst_3;
      double * RESTRICT _data_pdfs_src_00_37 = _data_pdfs_src + _stride_pdfs_src_0*ctr_0 + 7*_stride_pdfs_src_3;
      double * RESTRICT  _data_pdfs_dst_00_39 = _data_pdfs_dst + _stride_pdfs_dst_0*ctr_0 + 9*_stride_pdfs_dst_3;
      double * RESTRICT _data_pdfs_src_00_39 = _data_pdfs_src + _stride_pdfs_src_0*ctr_0 + 9*_stride_pdfs_src_3;
      double * RESTRICT  _data_pdfs_dst_00_313 = _data_pdfs_dst + _stride_pdfs_dst_0*ctr_0 + 13*_stride_pdfs_dst_3;
      double * RESTRICT _data_pdfs_src_00_313 = _data_pdfs_src + _stride_pdfs_src_0*ctr_0 + 13*_stride_pdfs_src_3;
      double * RESTRICT  _data_pdfs_dst_00_317 = _data_pdfs_dst + _stride_pdfs_dst_0*ctr_0 + 17*_stride_pdfs_dst_3;
      double * RESTRICT _data_pdfs_src_00_317 = _data_pdfs_src + _stride_pdfs_src_0*ctr_0 + 17*_stride_pdfs_src_3;
      for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_dst_1; ctr_1 += 1)
      {
         double * RESTRICT  _data_pdfs_dst_00_33_10 = _stride_pdfs_dst_1*ctr_1 + _data_pdfs_dst_00_33;
         double * RESTRICT _data_pdfs_src_00_33_10 = _stride_pdfs_src_1*ctr_1 + _data_pdfs_src_00_33;
         double * RESTRICT  _data_pdfs_dst_00_37_10 = _stride_pdfs_dst_1*ctr_1 + _data_pdfs_dst_00_37;
         double * RESTRICT _data_pdfs_src_00_37_10 = _stride_pdfs_src_1*ctr_1 + _data_pdfs_src_00_37;
         double * RESTRICT  _data_pdfs_dst_00_39_10 = _stride_pdfs_dst_1*ctr_1 + _data_pdfs_dst_00_39;
         double * RESTRICT _data_pdfs_src_00_39_10 = _stride_pdfs_src_1*ctr_1 + _data_pdfs_src_00_39;
         double * RESTRICT  _data_pdfs_dst_00_313_10 = _stride_pdfs_dst_1*ctr_1 + _data_pdfs_dst_00_313;
         double * RESTRICT _data_pdfs_src_00_313_10 = _stride_pdfs_src_1*ctr_1 + _data_pdfs_src_00_313;
         double * RESTRICT  _data_pdfs_dst_00_317_10 = _stride_pdfs_dst_1*ctr_1 + _data_pdfs_dst_00_317;
         double * RESTRICT _data_pdfs_src_00_317_10 = _stride_pdfs_src_1*ctr_1 + _data_pdfs_src_00_317;
         for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_dst_2; ctr_2 += 1)
         {
            _data_pdfs_dst_00_33_10[_stride_pdfs_dst_2*ctr_2] = _data_pdfs_src_00_33_10[_stride_pdfs_src_2*ctr_2];
            _data_pdfs_dst_00_37_10[_stride_pdfs_dst_2*ctr_2] = _data_pdfs_src_00_37_10[_stride_pdfs_src_2*ctr_2];
            _data_pdfs_dst_00_39_10[_stride_pdfs_dst_2*ctr_2] = _data_pdfs_src_00_39_10[_stride_pdfs_src_2*ctr_2];
            _data_pdfs_dst_00_313_10[_stride_pdfs_dst_2*ctr_2] = _data_pdfs_src_00_313_10[_stride_pdfs_src_2*ctr_2];
            _data_pdfs_dst_00_317_10[_stride_pdfs_dst_2*ctr_2] = _data_pdfs_src_00_317_10[_stride_pdfs_src_2*ctr_2];
         }
      }
   }
}
}

namespace internal_d3q19storagespecification_localCopy_E {
static FUNC_PREFIX void d3q19storagespecification_localCopy_E(double * RESTRICT  _data_pdfs_dst, double * RESTRICT const _data_pdfs_src, int64_t const _size_pdfs_dst_0, int64_t const _size_pdfs_dst_1, int64_t const _size_pdfs_dst_2, int64_t const _stride_pdfs_dst_0, int64_t const _stride_pdfs_dst_1, int64_t const _stride_pdfs_dst_2, int64_t const _stride_pdfs_dst_3, int64_t const _stride_pdfs_src_0, int64_t const _stride_pdfs_src_1, int64_t const _stride_pdfs_src_2, int64_t const _stride_pdfs_src_3)
{
   for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_dst_0; ctr_0 += 1)
   {
      double * RESTRICT  _data_pdfs_dst_00_34 = _data_pdfs_dst + _stride_pdfs_dst_0*ctr_0 + 4*_stride_pdfs_dst_3;
      double * RESTRICT _data_pdfs_src_00_34 = _data_pdfs_src + _stride_pdfs_src_0*ctr_0 + 4*_stride_pdfs_src_3;
      double * RESTRICT  _data_pdfs_dst_00_38 = _data_pdfs_dst + _stride_pdfs_dst_0*ctr_0 + 8*_stride_pdfs_dst_3;
      double * RESTRICT _data_pdfs_src_00_38 = _data_pdfs_src + _stride_pdfs_src_0*ctr_0 + 8*_stride_pdfs_src_3;
      double * RESTRICT  _data_pdfs_dst_00_310 = _data_pdfs_dst + _stride_pdfs_dst_0*ctr_0 + 10*_stride_pdfs_dst_3;
      double * RESTRICT _data_pdfs_src_00_310 = _data_pdfs_src + _stride_pdfs_src_0*ctr_0 + 10*_stride_pdfs_src_3;
      double * RESTRICT  _data_pdfs_dst_00_314 = _data_pdfs_dst + _stride_pdfs_dst_0*ctr_0 + 14*_stride_pdfs_dst_3;
      double * RESTRICT _data_pdfs_src_00_314 = _data_pdfs_src + _stride_pdfs_src_0*ctr_0 + 14*_stride_pdfs_src_3;
      double * RESTRICT  _data_pdfs_dst_00_318 = _data_pdfs_dst + _stride_pdfs_dst_0*ctr_0 + 18*_stride_pdfs_dst_3;
      double * RESTRICT _data_pdfs_src_00_318 = _data_pdfs_src + _stride_pdfs_src_0*ctr_0 + 18*_stride_pdfs_src_3;
      for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_dst_1; ctr_1 += 1)
      {
         double * RESTRICT  _data_pdfs_dst_00_34_10 = _stride_pdfs_dst_1*ctr_1 + _data_pdfs_dst_00_34;
         double * RESTRICT _data_pdfs_src_00_34_10 = _stride_pdfs_src_1*ctr_1 + _data_pdfs_src_00_34;
         double * RESTRICT  _data_pdfs_dst_00_38_10 = _stride_pdfs_dst_1*ctr_1 + _data_pdfs_dst_00_38;
         double * RESTRICT _data_pdfs_src_00_38_10 = _stride_pdfs_src_1*ctr_1 + _data_pdfs_src_00_38;
         double * RESTRICT  _data_pdfs_dst_00_310_10 = _stride_pdfs_dst_1*ctr_1 + _data_pdfs_dst_00_310;
         double * RESTRICT _data_pdfs_src_00_310_10 = _stride_pdfs_src_1*ctr_1 + _data_pdfs_src_00_310;
         double * RESTRICT  _data_pdfs_dst_00_314_10 = _stride_pdfs_dst_1*ctr_1 + _data_pdfs_dst_00_314;
         double * RESTRICT _data_pdfs_src_00_314_10 = _stride_pdfs_src_1*ctr_1 + _data_pdfs_src_00_314;
         double * RESTRICT  _data_pdfs_dst_00_318_10 = _stride_pdfs_dst_1*ctr_1 + _data_pdfs_dst_00_318;
         double * RESTRICT _data_pdfs_src_00_318_10 = _stride_pdfs_src_1*ctr_1 + _data_pdfs_src_00_318;
         for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_dst_2; ctr_2 += 1)
         {
            _data_pdfs_dst_00_34_10[_stride_pdfs_dst_2*ctr_2] = _data_pdfs_src_00_34_10[_stride_pdfs_src_2*ctr_2];
            _data_pdfs_dst_00_38_10[_stride_pdfs_dst_2*ctr_2] = _data_pdfs_src_00_38_10[_stride_pdfs_src_2*ctr_2];
            _data_pdfs_dst_00_310_10[_stride_pdfs_dst_2*ctr_2] = _data_pdfs_src_00_310_10[_stride_pdfs_src_2*ctr_2];
            _data_pdfs_dst_00_314_10[_stride_pdfs_dst_2*ctr_2] = _data_pdfs_src_00_314_10[_stride_pdfs_src_2*ctr_2];
            _data_pdfs_dst_00_318_10[_stride_pdfs_dst_2*ctr_2] = _data_pdfs_src_00_318_10[_stride_pdfs_src_2*ctr_2];
         }
      }
   }
}
}

namespace internal_d3q19storagespecification_localCopy_TW {
static FUNC_PREFIX void d3q19storagespecification_localCopy_TW(double * RESTRICT  _data_pdfs_dst, double * RESTRICT const _data_pdfs_src, int64_t const _size_pdfs_dst_0, int64_t const _size_pdfs_dst_1, int64_t const _size_pdfs_dst_2, int64_t const _stride_pdfs_dst_0, int64_t const _stride_pdfs_dst_1, int64_t const _stride_pdfs_dst_2, int64_t const _stride_pdfs_dst_3, int64_t const _stride_pdfs_src_0, int64_t const _stride_pdfs_src_1, int64_t const _stride_pdfs_src_2, int64_t const _stride_pdfs_src_3)
{
   for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_dst_0; ctr_0 += 1)
   {
      double * RESTRICT  _data_pdfs_dst_00_313 = _data_pdfs_dst + _stride_pdfs_dst_0*ctr_0 + 13*_stride_pdfs_dst_3;
      double * RESTRICT _data_pdfs_src_00_313 = _data_pdfs_src + _stride_pdfs_src_0*ctr_0 + 13*_stride_pdfs_src_3;
      for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_dst_1; ctr_1 += 1)
      {
         double * RESTRICT  _data_pdfs_dst_00_313_10 = _stride_pdfs_dst_1*ctr_1 + _data_pdfs_dst_00_313;
         double * RESTRICT _data_pdfs_src_00_313_10 = _stride_pdfs_src_1*ctr_1 + _data_pdfs_src_00_313;
         for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_dst_2; ctr_2 += 1)
         {
            _data_pdfs_dst_00_313_10[_stride_pdfs_dst_2*ctr_2] = _data_pdfs_src_00_313_10[_stride_pdfs_src_2*ctr_2];
         }
      }
   }
}
}

namespace internal_d3q19storagespecification_localCopy_SW {
static FUNC_PREFIX void d3q19storagespecification_localCopy_SW(double * RESTRICT  _data_pdfs_dst, double * RESTRICT const _data_pdfs_src, int64_t const _size_pdfs_dst_0, int64_t const _size_pdfs_dst_1, int64_t const _size_pdfs_dst_2, int64_t const _stride_pdfs_dst_0, int64_t const _stride_pdfs_dst_1, int64_t const _stride_pdfs_dst_2, int64_t const _stride_pdfs_dst_3, int64_t const _stride_pdfs_src_0, int64_t const _stride_pdfs_src_1, int64_t const _stride_pdfs_src_2, int64_t const _stride_pdfs_src_3)
{
   for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_dst_0; ctr_0 += 1)
   {
      double * RESTRICT  _data_pdfs_dst_00_39 = _data_pdfs_dst + _stride_pdfs_dst_0*ctr_0 + 9*_stride_pdfs_dst_3;
      double * RESTRICT _data_pdfs_src_00_39 = _data_pdfs_src + _stride_pdfs_src_0*ctr_0 + 9*_stride_pdfs_src_3;
      for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_dst_1; ctr_1 += 1)
      {
         double * RESTRICT  _data_pdfs_dst_00_39_10 = _stride_pdfs_dst_1*ctr_1 + _data_pdfs_dst_00_39;
         double * RESTRICT _data_pdfs_src_00_39_10 = _stride_pdfs_src_1*ctr_1 + _data_pdfs_src_00_39;
         for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_dst_2; ctr_2 += 1)
         {
            _data_pdfs_dst_00_39_10[_stride_pdfs_dst_2*ctr_2] = _data_pdfs_src_00_39_10[_stride_pdfs_src_2*ctr_2];
         }
      }
   }
}
}

namespace internal_d3q19storagespecification_localCopy_NW {
static FUNC_PREFIX void d3q19storagespecification_localCopy_NW(double * RESTRICT  _data_pdfs_dst, double * RESTRICT const _data_pdfs_src, int64_t const _size_pdfs_dst_0, int64_t const _size_pdfs_dst_1, int64_t const _size_pdfs_dst_2, int64_t const _stride_pdfs_dst_0, int64_t const _stride_pdfs_dst_1, int64_t const _stride_pdfs_dst_2, int64_t const _stride_pdfs_dst_3, int64_t const _stride_pdfs_src_0, int64_t const _stride_pdfs_src_1, int64_t const _stride_pdfs_src_2, int64_t const _stride_pdfs_src_3)
{
   for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_dst_0; ctr_0 += 1)
   {
      double * RESTRICT  _data_pdfs_dst_00_37 = _data_pdfs_dst + _stride_pdfs_dst_0*ctr_0 + 7*_stride_pdfs_dst_3;
      double * RESTRICT _data_pdfs_src_00_37 = _data_pdfs_src + _stride_pdfs_src_0*ctr_0 + 7*_stride_pdfs_src_3;
      for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_dst_1; ctr_1 += 1)
      {
         double * RESTRICT  _data_pdfs_dst_00_37_10 = _stride_pdfs_dst_1*ctr_1 + _data_pdfs_dst_00_37;
         double * RESTRICT _data_pdfs_src_00_37_10 = _stride_pdfs_src_1*ctr_1 + _data_pdfs_src_00_37;
         for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_dst_2; ctr_2 += 1)
         {
            _data_pdfs_dst_00_37_10[_stride_pdfs_dst_2*ctr_2] = _data_pdfs_src_00_37_10[_stride_pdfs_src_2*ctr_2];
         }
      }
   }
}
}

namespace internal_d3q19storagespecification_localCopy_BN {
static FUNC_PREFIX void d3q19storagespecification_localCopy_BN(double * RESTRICT  _data_pdfs_dst, double * RESTRICT const _data_pdfs_src, int64_t const _size_pdfs_dst_0, int64_t const _size_pdfs_dst_1, int64_t const _size_pdfs_dst_2, int64_t const _stride_pdfs_dst_0, int64_t const _stride_pdfs_dst_1, int64_t const _stride_pdfs_dst_2, int64_t const _stride_pdfs_dst_3, int64_t const _stride_pdfs_src_0, int64_t const _stride_pdfs_src_1, int64_t const _stride_pdfs_src_2, int64_t const _stride_pdfs_src_3)
{
   for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_dst_0; ctr_0 += 1)
   {
      double * RESTRICT  _data_pdfs_dst_00_315 = _data_pdfs_dst + _stride_pdfs_dst_0*ctr_0 + 15*_stride_pdfs_dst_3;
      double * RESTRICT _data_pdfs_src_00_315 = _data_pdfs_src + _stride_pdfs_src_0*ctr_0 + 15*_stride_pdfs_src_3;
      for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_dst_1; ctr_1 += 1)
      {
         double * RESTRICT  _data_pdfs_dst_00_315_10 = _stride_pdfs_dst_1*ctr_1 + _data_pdfs_dst_00_315;
         double * RESTRICT _data_pdfs_src_00_315_10 = _stride_pdfs_src_1*ctr_1 + _data_pdfs_src_00_315;
         for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_dst_2; ctr_2 += 1)
         {
            _data_pdfs_dst_00_315_10[_stride_pdfs_dst_2*ctr_2] = _data_pdfs_src_00_315_10[_stride_pdfs_src_2*ctr_2];
         }
      }
   }
}
}

namespace internal_d3q19storagespecification_localCopy_TE {
static FUNC_PREFIX void d3q19storagespecification_localCopy_TE(double * RESTRICT  _data_pdfs_dst, double * RESTRICT const _data_pdfs_src, int64_t const _size_pdfs_dst_0, int64_t const _size_pdfs_dst_1, int64_t const _size_pdfs_dst_2, int64_t const _stride_pdfs_dst_0, int64_t const _stride_pdfs_dst_1, int64_t const _stride_pdfs_dst_2, int64_t const _stride_pdfs_dst_3, int64_t const _stride_pdfs_src_0, int64_t const _stride_pdfs_src_1, int64_t const _stride_pdfs_src_2, int64_t const _stride_pdfs_src_3)
{
   for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_dst_0; ctr_0 += 1)
   {
      double * RESTRICT  _data_pdfs_dst_00_314 = _data_pdfs_dst + _stride_pdfs_dst_0*ctr_0 + 14*_stride_pdfs_dst_3;
      double * RESTRICT _data_pdfs_src_00_314 = _data_pdfs_src + _stride_pdfs_src_0*ctr_0 + 14*_stride_pdfs_src_3;
      for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_dst_1; ctr_1 += 1)
      {
         double * RESTRICT  _data_pdfs_dst_00_314_10 = _stride_pdfs_dst_1*ctr_1 + _data_pdfs_dst_00_314;
         double * RESTRICT _data_pdfs_src_00_314_10 = _stride_pdfs_src_1*ctr_1 + _data_pdfs_src_00_314;
         for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_dst_2; ctr_2 += 1)
         {
            _data_pdfs_dst_00_314_10[_stride_pdfs_dst_2*ctr_2] = _data_pdfs_src_00_314_10[_stride_pdfs_src_2*ctr_2];
         }
      }
   }
}
}

namespace internal_d3q19storagespecification_localCopy_B {
static FUNC_PREFIX void d3q19storagespecification_localCopy_B(double * RESTRICT  _data_pdfs_dst, double * RESTRICT const _data_pdfs_src, int64_t const _size_pdfs_dst_0, int64_t const _size_pdfs_dst_1, int64_t const _size_pdfs_dst_2, int64_t const _stride_pdfs_dst_0, int64_t const _stride_pdfs_dst_1, int64_t const _stride_pdfs_dst_2, int64_t const _stride_pdfs_dst_3, int64_t const _stride_pdfs_src_0, int64_t const _stride_pdfs_src_1, int64_t const _stride_pdfs_src_2, int64_t const _stride_pdfs_src_3)
{
   for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_dst_0; ctr_0 += 1)
   {
      double * RESTRICT  _data_pdfs_dst_00_36 = _data_pdfs_dst + _stride_pdfs_dst_0*ctr_0 + 6*_stride_pdfs_dst_3;
      double * RESTRICT _data_pdfs_src_00_36 = _data_pdfs_src + _stride_pdfs_src_0*ctr_0 + 6*_stride_pdfs_src_3;
      double * RESTRICT  _data_pdfs_dst_00_315 = _data_pdfs_dst + _stride_pdfs_dst_0*ctr_0 + 15*_stride_pdfs_dst_3;
      double * RESTRICT _data_pdfs_src_00_315 = _data_pdfs_src + _stride_pdfs_src_0*ctr_0 + 15*_stride_pdfs_src_3;
      double * RESTRICT  _data_pdfs_dst_00_316 = _data_pdfs_dst + _stride_pdfs_dst_0*ctr_0 + 16*_stride_pdfs_dst_3;
      double * RESTRICT _data_pdfs_src_00_316 = _data_pdfs_src + _stride_pdfs_src_0*ctr_0 + 16*_stride_pdfs_src_3;
      double * RESTRICT  _data_pdfs_dst_00_317 = _data_pdfs_dst + _stride_pdfs_dst_0*ctr_0 + 17*_stride_pdfs_dst_3;
      double * RESTRICT _data_pdfs_src_00_317 = _data_pdfs_src + _stride_pdfs_src_0*ctr_0 + 17*_stride_pdfs_src_3;
      double * RESTRICT  _data_pdfs_dst_00_318 = _data_pdfs_dst + _stride_pdfs_dst_0*ctr_0 + 18*_stride_pdfs_dst_3;
      double * RESTRICT _data_pdfs_src_00_318 = _data_pdfs_src + _stride_pdfs_src_0*ctr_0 + 18*_stride_pdfs_src_3;
      for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_dst_1; ctr_1 += 1)
      {
         double * RESTRICT  _data_pdfs_dst_00_36_10 = _stride_pdfs_dst_1*ctr_1 + _data_pdfs_dst_00_36;
         double * RESTRICT _data_pdfs_src_00_36_10 = _stride_pdfs_src_1*ctr_1 + _data_pdfs_src_00_36;
         double * RESTRICT  _data_pdfs_dst_00_315_10 = _stride_pdfs_dst_1*ctr_1 + _data_pdfs_dst_00_315;
         double * RESTRICT _data_pdfs_src_00_315_10 = _stride_pdfs_src_1*ctr_1 + _data_pdfs_src_00_315;
         double * RESTRICT  _data_pdfs_dst_00_316_10 = _stride_pdfs_dst_1*ctr_1 + _data_pdfs_dst_00_316;
         double * RESTRICT _data_pdfs_src_00_316_10 = _stride_pdfs_src_1*ctr_1 + _data_pdfs_src_00_316;
         double * RESTRICT  _data_pdfs_dst_00_317_10 = _stride_pdfs_dst_1*ctr_1 + _data_pdfs_dst_00_317;
         double * RESTRICT _data_pdfs_src_00_317_10 = _stride_pdfs_src_1*ctr_1 + _data_pdfs_src_00_317;
         double * RESTRICT  _data_pdfs_dst_00_318_10 = _stride_pdfs_dst_1*ctr_1 + _data_pdfs_dst_00_318;
         double * RESTRICT _data_pdfs_src_00_318_10 = _stride_pdfs_src_1*ctr_1 + _data_pdfs_src_00_318;
         for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_dst_2; ctr_2 += 1)
         {
            _data_pdfs_dst_00_36_10[_stride_pdfs_dst_2*ctr_2] = _data_pdfs_src_00_36_10[_stride_pdfs_src_2*ctr_2];
            _data_pdfs_dst_00_315_10[_stride_pdfs_dst_2*ctr_2] = _data_pdfs_src_00_315_10[_stride_pdfs_src_2*ctr_2];
            _data_pdfs_dst_00_316_10[_stride_pdfs_dst_2*ctr_2] = _data_pdfs_src_00_316_10[_stride_pdfs_src_2*ctr_2];
            _data_pdfs_dst_00_317_10[_stride_pdfs_dst_2*ctr_2] = _data_pdfs_src_00_317_10[_stride_pdfs_src_2*ctr_2];
            _data_pdfs_dst_00_318_10[_stride_pdfs_dst_2*ctr_2] = _data_pdfs_src_00_318_10[_stride_pdfs_src_2*ctr_2];
         }
      }
   }
}
}

namespace internal_d3q19storagespecification_localCopy_N {
static FUNC_PREFIX void d3q19storagespecification_localCopy_N(double * RESTRICT  _data_pdfs_dst, double * RESTRICT const _data_pdfs_src, int64_t const _size_pdfs_dst_0, int64_t const _size_pdfs_dst_1, int64_t const _size_pdfs_dst_2, int64_t const _stride_pdfs_dst_0, int64_t const _stride_pdfs_dst_1, int64_t const _stride_pdfs_dst_2, int64_t const _stride_pdfs_dst_3, int64_t const _stride_pdfs_src_0, int64_t const _stride_pdfs_src_1, int64_t const _stride_pdfs_src_2, int64_t const _stride_pdfs_src_3)
{
   for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_dst_0; ctr_0 += 1)
   {
      double * RESTRICT  _data_pdfs_dst_00_31 = _data_pdfs_dst + _stride_pdfs_dst_0*ctr_0 + _stride_pdfs_dst_3;
      double * RESTRICT _data_pdfs_src_00_31 = _data_pdfs_src + _stride_pdfs_src_0*ctr_0 + _stride_pdfs_src_3;
      double * RESTRICT  _data_pdfs_dst_00_37 = _data_pdfs_dst + _stride_pdfs_dst_0*ctr_0 + 7*_stride_pdfs_dst_3;
      double * RESTRICT _data_pdfs_src_00_37 = _data_pdfs_src + _stride_pdfs_src_0*ctr_0 + 7*_stride_pdfs_src_3;
      double * RESTRICT  _data_pdfs_dst_00_38 = _data_pdfs_dst + _stride_pdfs_dst_0*ctr_0 + 8*_stride_pdfs_dst_3;
      double * RESTRICT _data_pdfs_src_00_38 = _data_pdfs_src + _stride_pdfs_src_0*ctr_0 + 8*_stride_pdfs_src_3;
      double * RESTRICT  _data_pdfs_dst_00_311 = _data_pdfs_dst + _stride_pdfs_dst_0*ctr_0 + 11*_stride_pdfs_dst_3;
      double * RESTRICT _data_pdfs_src_00_311 = _data_pdfs_src + _stride_pdfs_src_0*ctr_0 + 11*_stride_pdfs_src_3;
      double * RESTRICT  _data_pdfs_dst_00_315 = _data_pdfs_dst + _stride_pdfs_dst_0*ctr_0 + 15*_stride_pdfs_dst_3;
      double * RESTRICT _data_pdfs_src_00_315 = _data_pdfs_src + _stride_pdfs_src_0*ctr_0 + 15*_stride_pdfs_src_3;
      for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_dst_1; ctr_1 += 1)
      {
         double * RESTRICT  _data_pdfs_dst_00_31_10 = _stride_pdfs_dst_1*ctr_1 + _data_pdfs_dst_00_31;
         double * RESTRICT _data_pdfs_src_00_31_10 = _stride_pdfs_src_1*ctr_1 + _data_pdfs_src_00_31;
         double * RESTRICT  _data_pdfs_dst_00_37_10 = _stride_pdfs_dst_1*ctr_1 + _data_pdfs_dst_00_37;
         double * RESTRICT _data_pdfs_src_00_37_10 = _stride_pdfs_src_1*ctr_1 + _data_pdfs_src_00_37;
         double * RESTRICT  _data_pdfs_dst_00_38_10 = _stride_pdfs_dst_1*ctr_1 + _data_pdfs_dst_00_38;
         double * RESTRICT _data_pdfs_src_00_38_10 = _stride_pdfs_src_1*ctr_1 + _data_pdfs_src_00_38;
         double * RESTRICT  _data_pdfs_dst_00_311_10 = _stride_pdfs_dst_1*ctr_1 + _data_pdfs_dst_00_311;
         double * RESTRICT _data_pdfs_src_00_311_10 = _stride_pdfs_src_1*ctr_1 + _data_pdfs_src_00_311;
         double * RESTRICT  _data_pdfs_dst_00_315_10 = _stride_pdfs_dst_1*ctr_1 + _data_pdfs_dst_00_315;
         double * RESTRICT _data_pdfs_src_00_315_10 = _stride_pdfs_src_1*ctr_1 + _data_pdfs_src_00_315;
         for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_dst_2; ctr_2 += 1)
         {
            _data_pdfs_dst_00_31_10[_stride_pdfs_dst_2*ctr_2] = _data_pdfs_src_00_31_10[_stride_pdfs_src_2*ctr_2];
            _data_pdfs_dst_00_37_10[_stride_pdfs_dst_2*ctr_2] = _data_pdfs_src_00_37_10[_stride_pdfs_src_2*ctr_2];
            _data_pdfs_dst_00_38_10[_stride_pdfs_dst_2*ctr_2] = _data_pdfs_src_00_38_10[_stride_pdfs_src_2*ctr_2];
            _data_pdfs_dst_00_311_10[_stride_pdfs_dst_2*ctr_2] = _data_pdfs_src_00_311_10[_stride_pdfs_src_2*ctr_2];
            _data_pdfs_dst_00_315_10[_stride_pdfs_dst_2*ctr_2] = _data_pdfs_src_00_315_10[_stride_pdfs_src_2*ctr_2];
         }
      }
   }
}
}

namespace internal_d3q19storagespecification_localCopy_S {
static FUNC_PREFIX void d3q19storagespecification_localCopy_S(double * RESTRICT  _data_pdfs_dst, double * RESTRICT const _data_pdfs_src, int64_t const _size_pdfs_dst_0, int64_t const _size_pdfs_dst_1, int64_t const _size_pdfs_dst_2, int64_t const _stride_pdfs_dst_0, int64_t const _stride_pdfs_dst_1, int64_t const _stride_pdfs_dst_2, int64_t const _stride_pdfs_dst_3, int64_t const _stride_pdfs_src_0, int64_t const _stride_pdfs_src_1, int64_t const _stride_pdfs_src_2, int64_t const _stride_pdfs_src_3)
{
   for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_dst_0; ctr_0 += 1)
   {
      double * RESTRICT  _data_pdfs_dst_00_32 = _data_pdfs_dst + _stride_pdfs_dst_0*ctr_0 + 2*_stride_pdfs_dst_3;
      double * RESTRICT _data_pdfs_src_00_32 = _data_pdfs_src + _stride_pdfs_src_0*ctr_0 + 2*_stride_pdfs_src_3;
      double * RESTRICT  _data_pdfs_dst_00_39 = _data_pdfs_dst + _stride_pdfs_dst_0*ctr_0 + 9*_stride_pdfs_dst_3;
      double * RESTRICT _data_pdfs_src_00_39 = _data_pdfs_src + _stride_pdfs_src_0*ctr_0 + 9*_stride_pdfs_src_3;
      double * RESTRICT  _data_pdfs_dst_00_310 = _data_pdfs_dst + _stride_pdfs_dst_0*ctr_0 + 10*_stride_pdfs_dst_3;
      double * RESTRICT _data_pdfs_src_00_310 = _data_pdfs_src + _stride_pdfs_src_0*ctr_0 + 10*_stride_pdfs_src_3;
      double * RESTRICT  _data_pdfs_dst_00_312 = _data_pdfs_dst + _stride_pdfs_dst_0*ctr_0 + 12*_stride_pdfs_dst_3;
      double * RESTRICT _data_pdfs_src_00_312 = _data_pdfs_src + _stride_pdfs_src_0*ctr_0 + 12*_stride_pdfs_src_3;
      double * RESTRICT  _data_pdfs_dst_00_316 = _data_pdfs_dst + _stride_pdfs_dst_0*ctr_0 + 16*_stride_pdfs_dst_3;
      double * RESTRICT _data_pdfs_src_00_316 = _data_pdfs_src + _stride_pdfs_src_0*ctr_0 + 16*_stride_pdfs_src_3;
      for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_dst_1; ctr_1 += 1)
      {
         double * RESTRICT  _data_pdfs_dst_00_32_10 = _stride_pdfs_dst_1*ctr_1 + _data_pdfs_dst_00_32;
         double * RESTRICT _data_pdfs_src_00_32_10 = _stride_pdfs_src_1*ctr_1 + _data_pdfs_src_00_32;
         double * RESTRICT  _data_pdfs_dst_00_39_10 = _stride_pdfs_dst_1*ctr_1 + _data_pdfs_dst_00_39;
         double * RESTRICT _data_pdfs_src_00_39_10 = _stride_pdfs_src_1*ctr_1 + _data_pdfs_src_00_39;
         double * RESTRICT  _data_pdfs_dst_00_310_10 = _stride_pdfs_dst_1*ctr_1 + _data_pdfs_dst_00_310;
         double * RESTRICT _data_pdfs_src_00_310_10 = _stride_pdfs_src_1*ctr_1 + _data_pdfs_src_00_310;
         double * RESTRICT  _data_pdfs_dst_00_312_10 = _stride_pdfs_dst_1*ctr_1 + _data_pdfs_dst_00_312;
         double * RESTRICT _data_pdfs_src_00_312_10 = _stride_pdfs_src_1*ctr_1 + _data_pdfs_src_00_312;
         double * RESTRICT  _data_pdfs_dst_00_316_10 = _stride_pdfs_dst_1*ctr_1 + _data_pdfs_dst_00_316;
         double * RESTRICT _data_pdfs_src_00_316_10 = _stride_pdfs_src_1*ctr_1 + _data_pdfs_src_00_316;
         for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_dst_2; ctr_2 += 1)
         {
            _data_pdfs_dst_00_32_10[_stride_pdfs_dst_2*ctr_2] = _data_pdfs_src_00_32_10[_stride_pdfs_src_2*ctr_2];
            _data_pdfs_dst_00_39_10[_stride_pdfs_dst_2*ctr_2] = _data_pdfs_src_00_39_10[_stride_pdfs_src_2*ctr_2];
            _data_pdfs_dst_00_310_10[_stride_pdfs_dst_2*ctr_2] = _data_pdfs_src_00_310_10[_stride_pdfs_src_2*ctr_2];
            _data_pdfs_dst_00_312_10[_stride_pdfs_dst_2*ctr_2] = _data_pdfs_src_00_312_10[_stride_pdfs_src_2*ctr_2];
            _data_pdfs_dst_00_316_10[_stride_pdfs_dst_2*ctr_2] = _data_pdfs_src_00_316_10[_stride_pdfs_src_2*ctr_2];
         }
      }
   }
}
}

namespace internal_d3q19storagespecification_localCopy_SE {
static FUNC_PREFIX void d3q19storagespecification_localCopy_SE(double * RESTRICT  _data_pdfs_dst, double * RESTRICT const _data_pdfs_src, int64_t const _size_pdfs_dst_0, int64_t const _size_pdfs_dst_1, int64_t const _size_pdfs_dst_2, int64_t const _stride_pdfs_dst_0, int64_t const _stride_pdfs_dst_1, int64_t const _stride_pdfs_dst_2, int64_t const _stride_pdfs_dst_3, int64_t const _stride_pdfs_src_0, int64_t const _stride_pdfs_src_1, int64_t const _stride_pdfs_src_2, int64_t const _stride_pdfs_src_3)
{
   for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_dst_0; ctr_0 += 1)
   {
      double * RESTRICT  _data_pdfs_dst_00_310 = _data_pdfs_dst + _stride_pdfs_dst_0*ctr_0 + 10*_stride_pdfs_dst_3;
      double * RESTRICT _data_pdfs_src_00_310 = _data_pdfs_src + _stride_pdfs_src_0*ctr_0 + 10*_stride_pdfs_src_3;
      for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_dst_1; ctr_1 += 1)
      {
         double * RESTRICT  _data_pdfs_dst_00_310_10 = _stride_pdfs_dst_1*ctr_1 + _data_pdfs_dst_00_310;
         double * RESTRICT _data_pdfs_src_00_310_10 = _stride_pdfs_src_1*ctr_1 + _data_pdfs_src_00_310;
         for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_dst_2; ctr_2 += 1)
         {
            _data_pdfs_dst_00_310_10[_stride_pdfs_dst_2*ctr_2] = _data_pdfs_src_00_310_10[_stride_pdfs_src_2*ctr_2];
         }
      }
   }
}
}




/*************************************************************************************
 *                                 Kernel Wrappers
*************************************************************************************/

namespace walberla {
namespace lbm {

   void D3Q19StorageSpecification::PackKernels::packAll(PdfField_T * pdfs_src, CellInterval & ci, unsigned char * outBuffer) const
   {
      double * buffer = reinterpret_cast<double*>(outBuffer);
      double * RESTRICT  _data_buffer = buffer;
      WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs_src->nrOfGhostLayers()))
      WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs_src->nrOfGhostLayers()))
      WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs_src->nrOfGhostLayers()))
      double * RESTRICT const _data_pdfs_src = pdfs_src->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
      WALBERLA_ASSERT_GREATER_EQUAL(pdfs_src->xSizeWithGhostLayer(), int64_t(int64_c(ci.xSize()) + 0))
      const int64_t _size_pdfs_src_0 = int64_t(int64_c(ci.xSize()) + 0);
      WALBERLA_ASSERT_GREATER_EQUAL(pdfs_src->ySizeWithGhostLayer(), int64_t(int64_c(ci.ySize()) + 0))
      const int64_t _size_pdfs_src_1 = int64_t(int64_c(ci.ySize()) + 0);
      WALBERLA_ASSERT_GREATER_EQUAL(pdfs_src->zSizeWithGhostLayer(), int64_t(int64_c(ci.zSize()) + 0))
      const int64_t _size_pdfs_src_2 = int64_t(int64_c(ci.zSize()) + 0);
      const int64_t _stride_pdfs_src_0 = int64_t(pdfs_src->xStride());
      const int64_t _stride_pdfs_src_1 = int64_t(pdfs_src->yStride());
      const int64_t _stride_pdfs_src_2 = int64_t(pdfs_src->zStride());
      const int64_t _stride_pdfs_src_3 = int64_t(1 * int64_t(pdfs_src->fStride()));
      internal_d3q19storagespecification_pack_ALL::d3q19storagespecification_pack_ALL(_data_buffer, _data_pdfs_src, _size_pdfs_src_0, _size_pdfs_src_1, _size_pdfs_src_2, _stride_pdfs_src_0, _stride_pdfs_src_1, _stride_pdfs_src_2, _stride_pdfs_src_3);
   }


   void D3Q19StorageSpecification::PackKernels::unpackAll(PdfField_T * pdfs_dst, CellInterval & ci, unsigned char * inBuffer) const
   {
      double * buffer = reinterpret_cast<double*>(inBuffer);
      double * RESTRICT const _data_buffer = buffer;
      WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs_dst->nrOfGhostLayers()))
      WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs_dst->nrOfGhostLayers()))
      WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs_dst->nrOfGhostLayers()))
      double * RESTRICT  _data_pdfs_dst = pdfs_dst->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
      WALBERLA_ASSERT_GREATER_EQUAL(pdfs_dst->xSizeWithGhostLayer(), int64_t(int64_c(ci.xSize()) + 0))
      const int64_t _size_pdfs_dst_0 = int64_t(int64_c(ci.xSize()) + 0);
      WALBERLA_ASSERT_GREATER_EQUAL(pdfs_dst->ySizeWithGhostLayer(), int64_t(int64_c(ci.ySize()) + 0))
      const int64_t _size_pdfs_dst_1 = int64_t(int64_c(ci.ySize()) + 0);
      WALBERLA_ASSERT_GREATER_EQUAL(pdfs_dst->zSizeWithGhostLayer(), int64_t(int64_c(ci.zSize()) + 0))
      const int64_t _size_pdfs_dst_2 = int64_t(int64_c(ci.zSize()) + 0);
      const int64_t _stride_pdfs_dst_0 = int64_t(pdfs_dst->xStride());
      const int64_t _stride_pdfs_dst_1 = int64_t(pdfs_dst->yStride());
      const int64_t _stride_pdfs_dst_2 = int64_t(pdfs_dst->zStride());
      const int64_t _stride_pdfs_dst_3 = int64_t(1 * int64_t(pdfs_dst->fStride()));
      internal_d3q19storagespecification_unpack_ALL::d3q19storagespecification_unpack_ALL(_data_buffer, _data_pdfs_dst, _size_pdfs_dst_0, _size_pdfs_dst_1, _size_pdfs_dst_2, _stride_pdfs_dst_0, _stride_pdfs_dst_1, _stride_pdfs_dst_2, _stride_pdfs_dst_3);
   }


   void D3Q19StorageSpecification::PackKernels::localCopyAll(PdfField_T * pdfs_src, CellInterval & srcInterval, PdfField_T * pdfs_dst, CellInterval & dstInterval) const
   {
      WALBERLA_ASSERT_EQUAL(srcInterval.xSize(), dstInterval.xSize())
      WALBERLA_ASSERT_EQUAL(srcInterval.ySize(), dstInterval.ySize())
      WALBERLA_ASSERT_EQUAL(srcInterval.zSize(), dstInterval.zSize())

      WALBERLA_ASSERT_GREATER_EQUAL(dstInterval.xMin(), -int_c(pdfs_dst->nrOfGhostLayers()))
      WALBERLA_ASSERT_GREATER_EQUAL(dstInterval.yMin(), -int_c(pdfs_dst->nrOfGhostLayers()))
      WALBERLA_ASSERT_GREATER_EQUAL(dstInterval.zMin(), -int_c(pdfs_dst->nrOfGhostLayers()))
      double * RESTRICT  _data_pdfs_dst = pdfs_dst->dataAt(dstInterval.xMin(), dstInterval.yMin(), dstInterval.zMin(), 0);
      WALBERLA_ASSERT_GREATER_EQUAL(srcInterval.xMin(), -int_c(pdfs_src->nrOfGhostLayers()))
      WALBERLA_ASSERT_GREATER_EQUAL(srcInterval.yMin(), -int_c(pdfs_src->nrOfGhostLayers()))
      WALBERLA_ASSERT_GREATER_EQUAL(srcInterval.zMin(), -int_c(pdfs_src->nrOfGhostLayers()))
      double * RESTRICT const _data_pdfs_src = pdfs_src->dataAt(srcInterval.xMin(), srcInterval.yMin(), srcInterval.zMin(), 0);
      WALBERLA_ASSERT_GREATER_EQUAL(pdfs_dst->xSizeWithGhostLayer(), int64_t(int64_c(dstInterval.xSize()) + 0))
      const int64_t _size_pdfs_dst_0 = int64_t(int64_c(dstInterval.xSize()) + 0);
      WALBERLA_ASSERT_GREATER_EQUAL(pdfs_dst->ySizeWithGhostLayer(), int64_t(int64_c(dstInterval.ySize()) + 0))
      const int64_t _size_pdfs_dst_1 = int64_t(int64_c(dstInterval.ySize()) + 0);
      WALBERLA_ASSERT_GREATER_EQUAL(pdfs_dst->zSizeWithGhostLayer(), int64_t(int64_c(dstInterval.zSize()) + 0))
      const int64_t _size_pdfs_dst_2 = int64_t(int64_c(dstInterval.zSize()) + 0);
      const int64_t _stride_pdfs_dst_0 = int64_t(pdfs_dst->xStride());
      const int64_t _stride_pdfs_dst_1 = int64_t(pdfs_dst->yStride());
      const int64_t _stride_pdfs_dst_2 = int64_t(pdfs_dst->zStride());
      const int64_t _stride_pdfs_dst_3 = int64_t(1 * int64_t(pdfs_dst->fStride()));
      const int64_t _stride_pdfs_src_0 = int64_t(pdfs_src->xStride());
      const int64_t _stride_pdfs_src_1 = int64_t(pdfs_src->yStride());
      const int64_t _stride_pdfs_src_2 = int64_t(pdfs_src->zStride());
      const int64_t _stride_pdfs_src_3 = int64_t(1 * int64_t(pdfs_src->fStride()));
      internal_d3q19storagespecification_localCopy_ALL::d3q19storagespecification_localCopy_ALL(_data_pdfs_dst, _data_pdfs_src, _size_pdfs_dst_0, _size_pdfs_dst_1, _size_pdfs_dst_2, _stride_pdfs_dst_0, _stride_pdfs_dst_1, _stride_pdfs_dst_2, _stride_pdfs_dst_3, _stride_pdfs_src_0, _stride_pdfs_src_1, _stride_pdfs_src_2, _stride_pdfs_src_3);
   }

   void D3Q19StorageSpecification::PackKernels::packDirection(PdfField_T * pdfs_src, CellInterval & ci, unsigned char * outBuffer, stencil::Direction dir) const
   {
      double * buffer = reinterpret_cast<double*>(outBuffer);
      double * RESTRICT  _data_buffer = buffer;
      WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs_src->nrOfGhostLayers()))
      WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs_src->nrOfGhostLayers()))
      WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs_src->nrOfGhostLayers()))
      double * RESTRICT const _data_pdfs_src = pdfs_src->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
      WALBERLA_ASSERT_GREATER_EQUAL(pdfs_src->xSizeWithGhostLayer(), int64_t(int64_c(ci.xSize()) + 0))
      const int64_t _size_pdfs_src_0 = int64_t(int64_c(ci.xSize()) + 0);
      WALBERLA_ASSERT_GREATER_EQUAL(pdfs_src->ySizeWithGhostLayer(), int64_t(int64_c(ci.ySize()) + 0))
      const int64_t _size_pdfs_src_1 = int64_t(int64_c(ci.ySize()) + 0);
      WALBERLA_ASSERT_GREATER_EQUAL(pdfs_src->zSizeWithGhostLayer(), int64_t(int64_c(ci.zSize()) + 0))
      const int64_t _size_pdfs_src_2 = int64_t(int64_c(ci.zSize()) + 0);
      const int64_t _stride_pdfs_src_0 = int64_t(pdfs_src->xStride());
      const int64_t _stride_pdfs_src_1 = int64_t(pdfs_src->yStride());
      const int64_t _stride_pdfs_src_2 = int64_t(pdfs_src->zStride());
      const int64_t _stride_pdfs_src_3 = int64_t(1 * int64_t(pdfs_src->fStride()));
      switch (dir) {
          case stencil::N : {
              internal_d3q19storagespecification_pack_N::d3q19storagespecification_pack_N(_data_buffer, _data_pdfs_src, _size_pdfs_src_0, _size_pdfs_src_1, _size_pdfs_src_2, _stride_pdfs_src_0, _stride_pdfs_src_1, _stride_pdfs_src_2, _stride_pdfs_src_3);
              break;
          }
          case stencil::S : {
              internal_d3q19storagespecification_pack_S::d3q19storagespecification_pack_S(_data_buffer, _data_pdfs_src, _size_pdfs_src_0, _size_pdfs_src_1, _size_pdfs_src_2, _stride_pdfs_src_0, _stride_pdfs_src_1, _stride_pdfs_src_2, _stride_pdfs_src_3);
              break;
          }
          case stencil::W : {
              internal_d3q19storagespecification_pack_W::d3q19storagespecification_pack_W(_data_buffer, _data_pdfs_src, _size_pdfs_src_0, _size_pdfs_src_1, _size_pdfs_src_2, _stride_pdfs_src_0, _stride_pdfs_src_1, _stride_pdfs_src_2, _stride_pdfs_src_3);
              break;
          }
          case stencil::E : {
              internal_d3q19storagespecification_pack_E::d3q19storagespecification_pack_E(_data_buffer, _data_pdfs_src, _size_pdfs_src_0, _size_pdfs_src_1, _size_pdfs_src_2, _stride_pdfs_src_0, _stride_pdfs_src_1, _stride_pdfs_src_2, _stride_pdfs_src_3);
              break;
          }
          case stencil::T : {
              internal_d3q19storagespecification_pack_T::d3q19storagespecification_pack_T(_data_buffer, _data_pdfs_src, _size_pdfs_src_0, _size_pdfs_src_1, _size_pdfs_src_2, _stride_pdfs_src_0, _stride_pdfs_src_1, _stride_pdfs_src_2, _stride_pdfs_src_3);
              break;
          }
          case stencil::B : {
              internal_d3q19storagespecification_pack_B::d3q19storagespecification_pack_B(_data_buffer, _data_pdfs_src, _size_pdfs_src_0, _size_pdfs_src_1, _size_pdfs_src_2, _stride_pdfs_src_0, _stride_pdfs_src_1, _stride_pdfs_src_2, _stride_pdfs_src_3);
              break;
          }
          case stencil::NW : {
              internal_d3q19storagespecification_pack_NW::d3q19storagespecification_pack_NW(_data_buffer, _data_pdfs_src, _size_pdfs_src_0, _size_pdfs_src_1, _size_pdfs_src_2, _stride_pdfs_src_0, _stride_pdfs_src_1, _stride_pdfs_src_2, _stride_pdfs_src_3);
              break;
          }
          case stencil::NE : {
              internal_d3q19storagespecification_pack_NE::d3q19storagespecification_pack_NE(_data_buffer, _data_pdfs_src, _size_pdfs_src_0, _size_pdfs_src_1, _size_pdfs_src_2, _stride_pdfs_src_0, _stride_pdfs_src_1, _stride_pdfs_src_2, _stride_pdfs_src_3);
              break;
          }
          case stencil::SW : {
              internal_d3q19storagespecification_pack_SW::d3q19storagespecification_pack_SW(_data_buffer, _data_pdfs_src, _size_pdfs_src_0, _size_pdfs_src_1, _size_pdfs_src_2, _stride_pdfs_src_0, _stride_pdfs_src_1, _stride_pdfs_src_2, _stride_pdfs_src_3);
              break;
          }
          case stencil::SE : {
              internal_d3q19storagespecification_pack_SE::d3q19storagespecification_pack_SE(_data_buffer, _data_pdfs_src, _size_pdfs_src_0, _size_pdfs_src_1, _size_pdfs_src_2, _stride_pdfs_src_0, _stride_pdfs_src_1, _stride_pdfs_src_2, _stride_pdfs_src_3);
              break;
          }
          case stencil::TN : {
              internal_d3q19storagespecification_pack_TN::d3q19storagespecification_pack_TN(_data_buffer, _data_pdfs_src, _size_pdfs_src_0, _size_pdfs_src_1, _size_pdfs_src_2, _stride_pdfs_src_0, _stride_pdfs_src_1, _stride_pdfs_src_2, _stride_pdfs_src_3);
              break;
          }
          case stencil::TS : {
              internal_d3q19storagespecification_pack_TS::d3q19storagespecification_pack_TS(_data_buffer, _data_pdfs_src, _size_pdfs_src_0, _size_pdfs_src_1, _size_pdfs_src_2, _stride_pdfs_src_0, _stride_pdfs_src_1, _stride_pdfs_src_2, _stride_pdfs_src_3);
              break;
          }
          case stencil::TW : {
              internal_d3q19storagespecification_pack_TW::d3q19storagespecification_pack_TW(_data_buffer, _data_pdfs_src, _size_pdfs_src_0, _size_pdfs_src_1, _size_pdfs_src_2, _stride_pdfs_src_0, _stride_pdfs_src_1, _stride_pdfs_src_2, _stride_pdfs_src_3);
              break;
          }
          case stencil::TE : {
              internal_d3q19storagespecification_pack_TE::d3q19storagespecification_pack_TE(_data_buffer, _data_pdfs_src, _size_pdfs_src_0, _size_pdfs_src_1, _size_pdfs_src_2, _stride_pdfs_src_0, _stride_pdfs_src_1, _stride_pdfs_src_2, _stride_pdfs_src_3);
              break;
          }
          case stencil::BN : {
              internal_d3q19storagespecification_pack_BN::d3q19storagespecification_pack_BN(_data_buffer, _data_pdfs_src, _size_pdfs_src_0, _size_pdfs_src_1, _size_pdfs_src_2, _stride_pdfs_src_0, _stride_pdfs_src_1, _stride_pdfs_src_2, _stride_pdfs_src_3);
              break;
          }
          case stencil::BS : {
              internal_d3q19storagespecification_pack_BS::d3q19storagespecification_pack_BS(_data_buffer, _data_pdfs_src, _size_pdfs_src_0, _size_pdfs_src_1, _size_pdfs_src_2, _stride_pdfs_src_0, _stride_pdfs_src_1, _stride_pdfs_src_2, _stride_pdfs_src_3);
              break;
          }
          case stencil::BW : {
              internal_d3q19storagespecification_pack_BW::d3q19storagespecification_pack_BW(_data_buffer, _data_pdfs_src, _size_pdfs_src_0, _size_pdfs_src_1, _size_pdfs_src_2, _stride_pdfs_src_0, _stride_pdfs_src_1, _stride_pdfs_src_2, _stride_pdfs_src_3);
              break;
          }
          case stencil::BE : {
              internal_d3q19storagespecification_pack_BE::d3q19storagespecification_pack_BE(_data_buffer, _data_pdfs_src, _size_pdfs_src_0, _size_pdfs_src_1, _size_pdfs_src_2, _stride_pdfs_src_0, _stride_pdfs_src_1, _stride_pdfs_src_2, _stride_pdfs_src_3);
              break;
          }default: break; 
      }
   }

   void D3Q19StorageSpecification::PackKernels::unpackDirection(PdfField_T * pdfs_dst, CellInterval & ci, unsigned char * inBuffer, stencil::Direction dir) const
   {
      double * buffer = reinterpret_cast<double*>(inBuffer);
      double * RESTRICT const _data_buffer = buffer;
      WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs_dst->nrOfGhostLayers()))
      WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs_dst->nrOfGhostLayers()))
      WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs_dst->nrOfGhostLayers()))
      double * RESTRICT  _data_pdfs_dst = pdfs_dst->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
      WALBERLA_ASSERT_GREATER_EQUAL(pdfs_dst->xSizeWithGhostLayer(), int64_t(int64_c(ci.xSize()) + 0))
      const int64_t _size_pdfs_dst_0 = int64_t(int64_c(ci.xSize()) + 0);
      WALBERLA_ASSERT_GREATER_EQUAL(pdfs_dst->ySizeWithGhostLayer(), int64_t(int64_c(ci.ySize()) + 0))
      const int64_t _size_pdfs_dst_1 = int64_t(int64_c(ci.ySize()) + 0);
      WALBERLA_ASSERT_GREATER_EQUAL(pdfs_dst->zSizeWithGhostLayer(), int64_t(int64_c(ci.zSize()) + 0))
      const int64_t _size_pdfs_dst_2 = int64_t(int64_c(ci.zSize()) + 0);
      const int64_t _stride_pdfs_dst_0 = int64_t(pdfs_dst->xStride());
      const int64_t _stride_pdfs_dst_1 = int64_t(pdfs_dst->yStride());
      const int64_t _stride_pdfs_dst_2 = int64_t(pdfs_dst->zStride());
      const int64_t _stride_pdfs_dst_3 = int64_t(1 * int64_t(pdfs_dst->fStride()));
      switch (dir) {
          case stencil::N : {
              internal_d3q19storagespecification_unpack_N::d3q19storagespecification_unpack_N(_data_buffer, _data_pdfs_dst, _size_pdfs_dst_0, _size_pdfs_dst_1, _size_pdfs_dst_2, _stride_pdfs_dst_0, _stride_pdfs_dst_1, _stride_pdfs_dst_2, _stride_pdfs_dst_3);
              break;
          }
          case stencil::S : {
              internal_d3q19storagespecification_unpack_S::d3q19storagespecification_unpack_S(_data_buffer, _data_pdfs_dst, _size_pdfs_dst_0, _size_pdfs_dst_1, _size_pdfs_dst_2, _stride_pdfs_dst_0, _stride_pdfs_dst_1, _stride_pdfs_dst_2, _stride_pdfs_dst_3);
              break;
          }
          case stencil::W : {
              internal_d3q19storagespecification_unpack_W::d3q19storagespecification_unpack_W(_data_buffer, _data_pdfs_dst, _size_pdfs_dst_0, _size_pdfs_dst_1, _size_pdfs_dst_2, _stride_pdfs_dst_0, _stride_pdfs_dst_1, _stride_pdfs_dst_2, _stride_pdfs_dst_3);
              break;
          }
          case stencil::E : {
              internal_d3q19storagespecification_unpack_E::d3q19storagespecification_unpack_E(_data_buffer, _data_pdfs_dst, _size_pdfs_dst_0, _size_pdfs_dst_1, _size_pdfs_dst_2, _stride_pdfs_dst_0, _stride_pdfs_dst_1, _stride_pdfs_dst_2, _stride_pdfs_dst_3);
              break;
          }
          case stencil::T : {
              internal_d3q19storagespecification_unpack_T::d3q19storagespecification_unpack_T(_data_buffer, _data_pdfs_dst, _size_pdfs_dst_0, _size_pdfs_dst_1, _size_pdfs_dst_2, _stride_pdfs_dst_0, _stride_pdfs_dst_1, _stride_pdfs_dst_2, _stride_pdfs_dst_3);
              break;
          }
          case stencil::B : {
              internal_d3q19storagespecification_unpack_B::d3q19storagespecification_unpack_B(_data_buffer, _data_pdfs_dst, _size_pdfs_dst_0, _size_pdfs_dst_1, _size_pdfs_dst_2, _stride_pdfs_dst_0, _stride_pdfs_dst_1, _stride_pdfs_dst_2, _stride_pdfs_dst_3);
              break;
          }
          case stencil::NW : {
              internal_d3q19storagespecification_unpack_NW::d3q19storagespecification_unpack_NW(_data_buffer, _data_pdfs_dst, _size_pdfs_dst_0, _size_pdfs_dst_1, _size_pdfs_dst_2, _stride_pdfs_dst_0, _stride_pdfs_dst_1, _stride_pdfs_dst_2, _stride_pdfs_dst_3);
              break;
          }
          case stencil::NE : {
              internal_d3q19storagespecification_unpack_NE::d3q19storagespecification_unpack_NE(_data_buffer, _data_pdfs_dst, _size_pdfs_dst_0, _size_pdfs_dst_1, _size_pdfs_dst_2, _stride_pdfs_dst_0, _stride_pdfs_dst_1, _stride_pdfs_dst_2, _stride_pdfs_dst_3);
              break;
          }
          case stencil::SW : {
              internal_d3q19storagespecification_unpack_SW::d3q19storagespecification_unpack_SW(_data_buffer, _data_pdfs_dst, _size_pdfs_dst_0, _size_pdfs_dst_1, _size_pdfs_dst_2, _stride_pdfs_dst_0, _stride_pdfs_dst_1, _stride_pdfs_dst_2, _stride_pdfs_dst_3);
              break;
          }
          case stencil::SE : {
              internal_d3q19storagespecification_unpack_SE::d3q19storagespecification_unpack_SE(_data_buffer, _data_pdfs_dst, _size_pdfs_dst_0, _size_pdfs_dst_1, _size_pdfs_dst_2, _stride_pdfs_dst_0, _stride_pdfs_dst_1, _stride_pdfs_dst_2, _stride_pdfs_dst_3);
              break;
          }
          case stencil::TN : {
              internal_d3q19storagespecification_unpack_TN::d3q19storagespecification_unpack_TN(_data_buffer, _data_pdfs_dst, _size_pdfs_dst_0, _size_pdfs_dst_1, _size_pdfs_dst_2, _stride_pdfs_dst_0, _stride_pdfs_dst_1, _stride_pdfs_dst_2, _stride_pdfs_dst_3);
              break;
          }
          case stencil::TS : {
              internal_d3q19storagespecification_unpack_TS::d3q19storagespecification_unpack_TS(_data_buffer, _data_pdfs_dst, _size_pdfs_dst_0, _size_pdfs_dst_1, _size_pdfs_dst_2, _stride_pdfs_dst_0, _stride_pdfs_dst_1, _stride_pdfs_dst_2, _stride_pdfs_dst_3);
              break;
          }
          case stencil::TW : {
              internal_d3q19storagespecification_unpack_TW::d3q19storagespecification_unpack_TW(_data_buffer, _data_pdfs_dst, _size_pdfs_dst_0, _size_pdfs_dst_1, _size_pdfs_dst_2, _stride_pdfs_dst_0, _stride_pdfs_dst_1, _stride_pdfs_dst_2, _stride_pdfs_dst_3);
              break;
          }
          case stencil::TE : {
              internal_d3q19storagespecification_unpack_TE::d3q19storagespecification_unpack_TE(_data_buffer, _data_pdfs_dst, _size_pdfs_dst_0, _size_pdfs_dst_1, _size_pdfs_dst_2, _stride_pdfs_dst_0, _stride_pdfs_dst_1, _stride_pdfs_dst_2, _stride_pdfs_dst_3);
              break;
          }
          case stencil::BN : {
              internal_d3q19storagespecification_unpack_BN::d3q19storagespecification_unpack_BN(_data_buffer, _data_pdfs_dst, _size_pdfs_dst_0, _size_pdfs_dst_1, _size_pdfs_dst_2, _stride_pdfs_dst_0, _stride_pdfs_dst_1, _stride_pdfs_dst_2, _stride_pdfs_dst_3);
              break;
          }
          case stencil::BS : {
              internal_d3q19storagespecification_unpack_BS::d3q19storagespecification_unpack_BS(_data_buffer, _data_pdfs_dst, _size_pdfs_dst_0, _size_pdfs_dst_1, _size_pdfs_dst_2, _stride_pdfs_dst_0, _stride_pdfs_dst_1, _stride_pdfs_dst_2, _stride_pdfs_dst_3);
              break;
          }
          case stencil::BW : {
              internal_d3q19storagespecification_unpack_BW::d3q19storagespecification_unpack_BW(_data_buffer, _data_pdfs_dst, _size_pdfs_dst_0, _size_pdfs_dst_1, _size_pdfs_dst_2, _stride_pdfs_dst_0, _stride_pdfs_dst_1, _stride_pdfs_dst_2, _stride_pdfs_dst_3);
              break;
          }
          case stencil::BE : {
              internal_d3q19storagespecification_unpack_BE::d3q19storagespecification_unpack_BE(_data_buffer, _data_pdfs_dst, _size_pdfs_dst_0, _size_pdfs_dst_1, _size_pdfs_dst_2, _stride_pdfs_dst_0, _stride_pdfs_dst_1, _stride_pdfs_dst_2, _stride_pdfs_dst_3);
              break;
          }default: break; 
      }
   }

   void D3Q19StorageSpecification::PackKernels::localCopyDirection(PdfField_T * pdfs_src, CellInterval & srcInterval, PdfField_T * pdfs_dst, CellInterval & dstInterval, stencil::Direction dir) const
   {
      WALBERLA_ASSERT_EQUAL(srcInterval.xSize(), dstInterval.xSize())
      WALBERLA_ASSERT_EQUAL(srcInterval.ySize(), dstInterval.ySize())
      WALBERLA_ASSERT_EQUAL(srcInterval.zSize(), dstInterval.zSize())

      WALBERLA_ASSERT_GREATER_EQUAL(dstInterval.xMin(), -int_c(pdfs_dst->nrOfGhostLayers()))
      WALBERLA_ASSERT_GREATER_EQUAL(dstInterval.yMin(), -int_c(pdfs_dst->nrOfGhostLayers()))
      WALBERLA_ASSERT_GREATER_EQUAL(dstInterval.zMin(), -int_c(pdfs_dst->nrOfGhostLayers()))
      double * RESTRICT  _data_pdfs_dst = pdfs_dst->dataAt(dstInterval.xMin(), dstInterval.yMin(), dstInterval.zMin(), 0);
      WALBERLA_ASSERT_GREATER_EQUAL(srcInterval.xMin(), -int_c(pdfs_src->nrOfGhostLayers()))
      WALBERLA_ASSERT_GREATER_EQUAL(srcInterval.yMin(), -int_c(pdfs_src->nrOfGhostLayers()))
      WALBERLA_ASSERT_GREATER_EQUAL(srcInterval.zMin(), -int_c(pdfs_src->nrOfGhostLayers()))
      double * RESTRICT const _data_pdfs_src = pdfs_src->dataAt(srcInterval.xMin(), srcInterval.yMin(), srcInterval.zMin(), 0);
      WALBERLA_ASSERT_GREATER_EQUAL(pdfs_dst->xSizeWithGhostLayer(), int64_t(int64_c(dstInterval.xSize()) + 0))
      const int64_t _size_pdfs_dst_0 = int64_t(int64_c(dstInterval.xSize()) + 0);
      WALBERLA_ASSERT_GREATER_EQUAL(pdfs_dst->ySizeWithGhostLayer(), int64_t(int64_c(dstInterval.ySize()) + 0))
      const int64_t _size_pdfs_dst_1 = int64_t(int64_c(dstInterval.ySize()) + 0);
      WALBERLA_ASSERT_GREATER_EQUAL(pdfs_dst->zSizeWithGhostLayer(), int64_t(int64_c(dstInterval.zSize()) + 0))
      const int64_t _size_pdfs_dst_2 = int64_t(int64_c(dstInterval.zSize()) + 0);
      const int64_t _stride_pdfs_dst_0 = int64_t(pdfs_dst->xStride());
      const int64_t _stride_pdfs_dst_1 = int64_t(pdfs_dst->yStride());
      const int64_t _stride_pdfs_dst_2 = int64_t(pdfs_dst->zStride());
      const int64_t _stride_pdfs_dst_3 = int64_t(1 * int64_t(pdfs_dst->fStride()));
      const int64_t _stride_pdfs_src_0 = int64_t(pdfs_src->xStride());
      const int64_t _stride_pdfs_src_1 = int64_t(pdfs_src->yStride());
      const int64_t _stride_pdfs_src_2 = int64_t(pdfs_src->zStride());
      const int64_t _stride_pdfs_src_3 = int64_t(1 * int64_t(pdfs_src->fStride()));
      switch (dir) {
          case stencil::N : {
              internal_d3q19storagespecification_localCopy_N::d3q19storagespecification_localCopy_N(_data_pdfs_dst, _data_pdfs_src, _size_pdfs_dst_0, _size_pdfs_dst_1, _size_pdfs_dst_2, _stride_pdfs_dst_0, _stride_pdfs_dst_1, _stride_pdfs_dst_2, _stride_pdfs_dst_3, _stride_pdfs_src_0, _stride_pdfs_src_1, _stride_pdfs_src_2, _stride_pdfs_src_3);
              break;
          }
          case stencil::S : {
              internal_d3q19storagespecification_localCopy_S::d3q19storagespecification_localCopy_S(_data_pdfs_dst, _data_pdfs_src, _size_pdfs_dst_0, _size_pdfs_dst_1, _size_pdfs_dst_2, _stride_pdfs_dst_0, _stride_pdfs_dst_1, _stride_pdfs_dst_2, _stride_pdfs_dst_3, _stride_pdfs_src_0, _stride_pdfs_src_1, _stride_pdfs_src_2, _stride_pdfs_src_3);
              break;
          }
          case stencil::W : {
              internal_d3q19storagespecification_localCopy_W::d3q19storagespecification_localCopy_W(_data_pdfs_dst, _data_pdfs_src, _size_pdfs_dst_0, _size_pdfs_dst_1, _size_pdfs_dst_2, _stride_pdfs_dst_0, _stride_pdfs_dst_1, _stride_pdfs_dst_2, _stride_pdfs_dst_3, _stride_pdfs_src_0, _stride_pdfs_src_1, _stride_pdfs_src_2, _stride_pdfs_src_3);
              break;
          }
          case stencil::E : {
              internal_d3q19storagespecification_localCopy_E::d3q19storagespecification_localCopy_E(_data_pdfs_dst, _data_pdfs_src, _size_pdfs_dst_0, _size_pdfs_dst_1, _size_pdfs_dst_2, _stride_pdfs_dst_0, _stride_pdfs_dst_1, _stride_pdfs_dst_2, _stride_pdfs_dst_3, _stride_pdfs_src_0, _stride_pdfs_src_1, _stride_pdfs_src_2, _stride_pdfs_src_3);
              break;
          }
          case stencil::T : {
              internal_d3q19storagespecification_localCopy_T::d3q19storagespecification_localCopy_T(_data_pdfs_dst, _data_pdfs_src, _size_pdfs_dst_0, _size_pdfs_dst_1, _size_pdfs_dst_2, _stride_pdfs_dst_0, _stride_pdfs_dst_1, _stride_pdfs_dst_2, _stride_pdfs_dst_3, _stride_pdfs_src_0, _stride_pdfs_src_1, _stride_pdfs_src_2, _stride_pdfs_src_3);
              break;
          }
          case stencil::B : {
              internal_d3q19storagespecification_localCopy_B::d3q19storagespecification_localCopy_B(_data_pdfs_dst, _data_pdfs_src, _size_pdfs_dst_0, _size_pdfs_dst_1, _size_pdfs_dst_2, _stride_pdfs_dst_0, _stride_pdfs_dst_1, _stride_pdfs_dst_2, _stride_pdfs_dst_3, _stride_pdfs_src_0, _stride_pdfs_src_1, _stride_pdfs_src_2, _stride_pdfs_src_3);
              break;
          }
          case stencil::NW : {
              internal_d3q19storagespecification_localCopy_NW::d3q19storagespecification_localCopy_NW(_data_pdfs_dst, _data_pdfs_src, _size_pdfs_dst_0, _size_pdfs_dst_1, _size_pdfs_dst_2, _stride_pdfs_dst_0, _stride_pdfs_dst_1, _stride_pdfs_dst_2, _stride_pdfs_dst_3, _stride_pdfs_src_0, _stride_pdfs_src_1, _stride_pdfs_src_2, _stride_pdfs_src_3);
              break;
          }
          case stencil::NE : {
              internal_d3q19storagespecification_localCopy_NE::d3q19storagespecification_localCopy_NE(_data_pdfs_dst, _data_pdfs_src, _size_pdfs_dst_0, _size_pdfs_dst_1, _size_pdfs_dst_2, _stride_pdfs_dst_0, _stride_pdfs_dst_1, _stride_pdfs_dst_2, _stride_pdfs_dst_3, _stride_pdfs_src_0, _stride_pdfs_src_1, _stride_pdfs_src_2, _stride_pdfs_src_3);
              break;
          }
          case stencil::SW : {
              internal_d3q19storagespecification_localCopy_SW::d3q19storagespecification_localCopy_SW(_data_pdfs_dst, _data_pdfs_src, _size_pdfs_dst_0, _size_pdfs_dst_1, _size_pdfs_dst_2, _stride_pdfs_dst_0, _stride_pdfs_dst_1, _stride_pdfs_dst_2, _stride_pdfs_dst_3, _stride_pdfs_src_0, _stride_pdfs_src_1, _stride_pdfs_src_2, _stride_pdfs_src_3);
              break;
          }
          case stencil::SE : {
              internal_d3q19storagespecification_localCopy_SE::d3q19storagespecification_localCopy_SE(_data_pdfs_dst, _data_pdfs_src, _size_pdfs_dst_0, _size_pdfs_dst_1, _size_pdfs_dst_2, _stride_pdfs_dst_0, _stride_pdfs_dst_1, _stride_pdfs_dst_2, _stride_pdfs_dst_3, _stride_pdfs_src_0, _stride_pdfs_src_1, _stride_pdfs_src_2, _stride_pdfs_src_3);
              break;
          }
          case stencil::TN : {
              internal_d3q19storagespecification_localCopy_TN::d3q19storagespecification_localCopy_TN(_data_pdfs_dst, _data_pdfs_src, _size_pdfs_dst_0, _size_pdfs_dst_1, _size_pdfs_dst_2, _stride_pdfs_dst_0, _stride_pdfs_dst_1, _stride_pdfs_dst_2, _stride_pdfs_dst_3, _stride_pdfs_src_0, _stride_pdfs_src_1, _stride_pdfs_src_2, _stride_pdfs_src_3);
              break;
          }
          case stencil::TS : {
              internal_d3q19storagespecification_localCopy_TS::d3q19storagespecification_localCopy_TS(_data_pdfs_dst, _data_pdfs_src, _size_pdfs_dst_0, _size_pdfs_dst_1, _size_pdfs_dst_2, _stride_pdfs_dst_0, _stride_pdfs_dst_1, _stride_pdfs_dst_2, _stride_pdfs_dst_3, _stride_pdfs_src_0, _stride_pdfs_src_1, _stride_pdfs_src_2, _stride_pdfs_src_3);
              break;
          }
          case stencil::TW : {
              internal_d3q19storagespecification_localCopy_TW::d3q19storagespecification_localCopy_TW(_data_pdfs_dst, _data_pdfs_src, _size_pdfs_dst_0, _size_pdfs_dst_1, _size_pdfs_dst_2, _stride_pdfs_dst_0, _stride_pdfs_dst_1, _stride_pdfs_dst_2, _stride_pdfs_dst_3, _stride_pdfs_src_0, _stride_pdfs_src_1, _stride_pdfs_src_2, _stride_pdfs_src_3);
              break;
          }
          case stencil::TE : {
              internal_d3q19storagespecification_localCopy_TE::d3q19storagespecification_localCopy_TE(_data_pdfs_dst, _data_pdfs_src, _size_pdfs_dst_0, _size_pdfs_dst_1, _size_pdfs_dst_2, _stride_pdfs_dst_0, _stride_pdfs_dst_1, _stride_pdfs_dst_2, _stride_pdfs_dst_3, _stride_pdfs_src_0, _stride_pdfs_src_1, _stride_pdfs_src_2, _stride_pdfs_src_3);
              break;
          }
          case stencil::BN : {
              internal_d3q19storagespecification_localCopy_BN::d3q19storagespecification_localCopy_BN(_data_pdfs_dst, _data_pdfs_src, _size_pdfs_dst_0, _size_pdfs_dst_1, _size_pdfs_dst_2, _stride_pdfs_dst_0, _stride_pdfs_dst_1, _stride_pdfs_dst_2, _stride_pdfs_dst_3, _stride_pdfs_src_0, _stride_pdfs_src_1, _stride_pdfs_src_2, _stride_pdfs_src_3);
              break;
          }
          case stencil::BS : {
              internal_d3q19storagespecification_localCopy_BS::d3q19storagespecification_localCopy_BS(_data_pdfs_dst, _data_pdfs_src, _size_pdfs_dst_0, _size_pdfs_dst_1, _size_pdfs_dst_2, _stride_pdfs_dst_0, _stride_pdfs_dst_1, _stride_pdfs_dst_2, _stride_pdfs_dst_3, _stride_pdfs_src_0, _stride_pdfs_src_1, _stride_pdfs_src_2, _stride_pdfs_src_3);
              break;
          }
          case stencil::BW : {
              internal_d3q19storagespecification_localCopy_BW::d3q19storagespecification_localCopy_BW(_data_pdfs_dst, _data_pdfs_src, _size_pdfs_dst_0, _size_pdfs_dst_1, _size_pdfs_dst_2, _stride_pdfs_dst_0, _stride_pdfs_dst_1, _stride_pdfs_dst_2, _stride_pdfs_dst_3, _stride_pdfs_src_0, _stride_pdfs_src_1, _stride_pdfs_src_2, _stride_pdfs_src_3);
              break;
          }
          case stencil::BE : {
              internal_d3q19storagespecification_localCopy_BE::d3q19storagespecification_localCopy_BE(_data_pdfs_dst, _data_pdfs_src, _size_pdfs_dst_0, _size_pdfs_dst_1, _size_pdfs_dst_2, _stride_pdfs_dst_0, _stride_pdfs_dst_1, _stride_pdfs_dst_2, _stride_pdfs_dst_3, _stride_pdfs_src_0, _stride_pdfs_src_1, _stride_pdfs_src_2, _stride_pdfs_src_3);
              break;
          }default: break; 
      }
   }

   
}  // namespace lbm
}  // namespace walberla