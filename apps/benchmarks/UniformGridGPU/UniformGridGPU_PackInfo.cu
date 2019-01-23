#include "stencil/Directions.h"
#include "core/cell/CellInterval.h"
#include "cuda/GPUField.h"
#include "core/DataTypes.h"
#include "UniformGridGPU_PackInfo.h"


#define FUNC_PREFIX __global__


namespace walberla {
namespace pystencils {

using walberla::cell::CellInterval;
using walberla::stencil::Direction;



namespace internal_pack_SW {
static FUNC_PREFIX void pack_SW(double * _data_buffer, double * const _data_pdfs, int64_t const _size_pdfs_0, int64_t const _size_pdfs_1, int64_t const _size_pdfs_2, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2, int64_t const _stride_pdfs_3)
{
   if (blockDim.x*blockIdx.x + threadIdx.x < _size_pdfs_0 && blockDim.y*blockIdx.y + threadIdx.y < _size_pdfs_1 && blockDim.z*blockIdx.z + threadIdx.z < _size_pdfs_2)
   {
      const int64_t ctr_0 = blockDim.x*blockIdx.x + threadIdx.x;
      const int64_t ctr_1 = blockDim.y*blockIdx.y + threadIdx.y;
      const int64_t ctr_2 = blockDim.z*blockIdx.z + threadIdx.z;
      double * const _data_pdfs_10_20_39 = _data_pdfs + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 9*_stride_pdfs_3;
      _data_buffer[_size_pdfs_0*_size_pdfs_1*(blockDim.z*blockIdx.z + threadIdx.z) + _size_pdfs_0*(blockDim.y*blockIdx.y + threadIdx.y) + blockDim.x*blockIdx.x + threadIdx.x] = _data_pdfs_10_20_39[_stride_pdfs_0*ctr_0];
   } 
}
}

namespace internal_pack_BW {
static FUNC_PREFIX void pack_BW(double * _data_buffer, double * const _data_pdfs, int64_t const _size_pdfs_0, int64_t const _size_pdfs_1, int64_t const _size_pdfs_2, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2, int64_t const _stride_pdfs_3)
{
   if (blockDim.x*blockIdx.x + threadIdx.x < _size_pdfs_0 && blockDim.y*blockIdx.y + threadIdx.y < _size_pdfs_1 && blockDim.z*blockIdx.z + threadIdx.z < _size_pdfs_2)
   {
      const int64_t ctr_0 = blockDim.x*blockIdx.x + threadIdx.x;
      const int64_t ctr_1 = blockDim.y*blockIdx.y + threadIdx.y;
      const int64_t ctr_2 = blockDim.z*blockIdx.z + threadIdx.z;
      double * const _data_pdfs_10_20_317 = _data_pdfs + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 17*_stride_pdfs_3;
      _data_buffer[_size_pdfs_0*_size_pdfs_1*(blockDim.z*blockIdx.z + threadIdx.z) + _size_pdfs_0*(blockDim.y*blockIdx.y + threadIdx.y) + blockDim.x*blockIdx.x + threadIdx.x] = _data_pdfs_10_20_317[_stride_pdfs_0*ctr_0];
   } 
}
}

namespace internal_pack_W {
static FUNC_PREFIX void pack_W(double * _data_buffer, double * const _data_pdfs, int64_t const _size_pdfs_0, int64_t const _size_pdfs_1, int64_t const _size_pdfs_2, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2, int64_t const _stride_pdfs_3)
{
   if (blockDim.x*blockIdx.x + threadIdx.x < _size_pdfs_0 && blockDim.y*blockIdx.y + threadIdx.y < _size_pdfs_1 && blockDim.z*blockIdx.z + threadIdx.z < _size_pdfs_2)
   {
      const int64_t ctr_0 = blockDim.x*blockIdx.x + threadIdx.x;
      const int64_t ctr_1 = blockDim.y*blockIdx.y + threadIdx.y;
      const int64_t ctr_2 = blockDim.z*blockIdx.z + threadIdx.z;
      double * const _data_pdfs_10_20_313 = _data_pdfs + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 13*_stride_pdfs_3;
      _data_buffer[_size_pdfs_0*_size_pdfs_1*(5*blockDim.z*blockIdx.z + 5*threadIdx.z) + _size_pdfs_0*(5*blockDim.y*blockIdx.y + 5*threadIdx.y) + 5*blockDim.x*blockIdx.x + 5*threadIdx.x] = _data_pdfs_10_20_313[_stride_pdfs_0*ctr_0];
      double * const _data_pdfs_10_20_317 = _data_pdfs + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 17*_stride_pdfs_3;
      _data_buffer[_size_pdfs_0*_size_pdfs_1*(5*blockDim.z*blockIdx.z + 5*threadIdx.z) + _size_pdfs_0*(5*blockDim.y*blockIdx.y + 5*threadIdx.y) + 5*blockDim.x*blockIdx.x + 5*threadIdx.x + 1] = _data_pdfs_10_20_317[_stride_pdfs_0*ctr_0];
      double * const _data_pdfs_10_20_33 = _data_pdfs + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 3*_stride_pdfs_3;
      _data_buffer[_size_pdfs_0*_size_pdfs_1*(5*blockDim.z*blockIdx.z + 5*threadIdx.z) + _size_pdfs_0*(5*blockDim.y*blockIdx.y + 5*threadIdx.y) + 5*blockDim.x*blockIdx.x + 5*threadIdx.x + 2] = _data_pdfs_10_20_33[_stride_pdfs_0*ctr_0];
      double * const _data_pdfs_10_20_37 = _data_pdfs + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 7*_stride_pdfs_3;
      _data_buffer[_size_pdfs_0*_size_pdfs_1*(5*blockDim.z*blockIdx.z + 5*threadIdx.z) + _size_pdfs_0*(5*blockDim.y*blockIdx.y + 5*threadIdx.y) + 5*blockDim.x*blockIdx.x + 5*threadIdx.x + 3] = _data_pdfs_10_20_37[_stride_pdfs_0*ctr_0];
      double * const _data_pdfs_10_20_39 = _data_pdfs + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 9*_stride_pdfs_3;
      _data_buffer[_size_pdfs_0*_size_pdfs_1*(5*blockDim.z*blockIdx.z + 5*threadIdx.z) + _size_pdfs_0*(5*blockDim.y*blockIdx.y + 5*threadIdx.y) + 5*blockDim.x*blockIdx.x + 5*threadIdx.x + 4] = _data_pdfs_10_20_39[_stride_pdfs_0*ctr_0];
   } 
}
}

namespace internal_pack_TW {
static FUNC_PREFIX void pack_TW(double * _data_buffer, double * const _data_pdfs, int64_t const _size_pdfs_0, int64_t const _size_pdfs_1, int64_t const _size_pdfs_2, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2, int64_t const _stride_pdfs_3)
{
   if (blockDim.x*blockIdx.x + threadIdx.x < _size_pdfs_0 && blockDim.y*blockIdx.y + threadIdx.y < _size_pdfs_1 && blockDim.z*blockIdx.z + threadIdx.z < _size_pdfs_2)
   {
      const int64_t ctr_0 = blockDim.x*blockIdx.x + threadIdx.x;
      const int64_t ctr_1 = blockDim.y*blockIdx.y + threadIdx.y;
      const int64_t ctr_2 = blockDim.z*blockIdx.z + threadIdx.z;
      double * const _data_pdfs_10_20_313 = _data_pdfs + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 13*_stride_pdfs_3;
      _data_buffer[_size_pdfs_0*_size_pdfs_1*(blockDim.z*blockIdx.z + threadIdx.z) + _size_pdfs_0*(blockDim.y*blockIdx.y + threadIdx.y) + blockDim.x*blockIdx.x + threadIdx.x] = _data_pdfs_10_20_313[_stride_pdfs_0*ctr_0];
   } 
}
}

namespace internal_pack_NW {
static FUNC_PREFIX void pack_NW(double * _data_buffer, double * const _data_pdfs, int64_t const _size_pdfs_0, int64_t const _size_pdfs_1, int64_t const _size_pdfs_2, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2, int64_t const _stride_pdfs_3)
{
   if (blockDim.x*blockIdx.x + threadIdx.x < _size_pdfs_0 && blockDim.y*blockIdx.y + threadIdx.y < _size_pdfs_1 && blockDim.z*blockIdx.z + threadIdx.z < _size_pdfs_2)
   {
      const int64_t ctr_0 = blockDim.x*blockIdx.x + threadIdx.x;
      const int64_t ctr_1 = blockDim.y*blockIdx.y + threadIdx.y;
      const int64_t ctr_2 = blockDim.z*blockIdx.z + threadIdx.z;
      double * const _data_pdfs_10_20_37 = _data_pdfs + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 7*_stride_pdfs_3;
      _data_buffer[_size_pdfs_0*_size_pdfs_1*(blockDim.z*blockIdx.z + threadIdx.z) + _size_pdfs_0*(blockDim.y*blockIdx.y + threadIdx.y) + blockDim.x*blockIdx.x + threadIdx.x] = _data_pdfs_10_20_37[_stride_pdfs_0*ctr_0];
   } 
}
}

namespace internal_pack_BS {
static FUNC_PREFIX void pack_BS(double * _data_buffer, double * const _data_pdfs, int64_t const _size_pdfs_0, int64_t const _size_pdfs_1, int64_t const _size_pdfs_2, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2, int64_t const _stride_pdfs_3)
{
   if (blockDim.x*blockIdx.x + threadIdx.x < _size_pdfs_0 && blockDim.y*blockIdx.y + threadIdx.y < _size_pdfs_1 && blockDim.z*blockIdx.z + threadIdx.z < _size_pdfs_2)
   {
      const int64_t ctr_0 = blockDim.x*blockIdx.x + threadIdx.x;
      const int64_t ctr_1 = blockDim.y*blockIdx.y + threadIdx.y;
      const int64_t ctr_2 = blockDim.z*blockIdx.z + threadIdx.z;
      double * const _data_pdfs_10_20_316 = _data_pdfs + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 16*_stride_pdfs_3;
      _data_buffer[_size_pdfs_0*_size_pdfs_1*(blockDim.z*blockIdx.z + threadIdx.z) + _size_pdfs_0*(blockDim.y*blockIdx.y + threadIdx.y) + blockDim.x*blockIdx.x + threadIdx.x] = _data_pdfs_10_20_316[_stride_pdfs_0*ctr_0];
   } 
}
}

namespace internal_pack_S {
static FUNC_PREFIX void pack_S(double * _data_buffer, double * const _data_pdfs, int64_t const _size_pdfs_0, int64_t const _size_pdfs_1, int64_t const _size_pdfs_2, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2, int64_t const _stride_pdfs_3)
{
   if (blockDim.x*blockIdx.x + threadIdx.x < _size_pdfs_0 && blockDim.y*blockIdx.y + threadIdx.y < _size_pdfs_1 && blockDim.z*blockIdx.z + threadIdx.z < _size_pdfs_2)
   {
      const int64_t ctr_0 = blockDim.x*blockIdx.x + threadIdx.x;
      const int64_t ctr_1 = blockDim.y*blockIdx.y + threadIdx.y;
      const int64_t ctr_2 = blockDim.z*blockIdx.z + threadIdx.z;
      double * const _data_pdfs_10_20_310 = _data_pdfs + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 10*_stride_pdfs_3;
      _data_buffer[_size_pdfs_0*_size_pdfs_1*(5*blockDim.z*blockIdx.z + 5*threadIdx.z) + _size_pdfs_0*(5*blockDim.y*blockIdx.y + 5*threadIdx.y) + 5*blockDim.x*blockIdx.x + 5*threadIdx.x] = _data_pdfs_10_20_310[_stride_pdfs_0*ctr_0];
      double * const _data_pdfs_10_20_312 = _data_pdfs + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 12*_stride_pdfs_3;
      _data_buffer[_size_pdfs_0*_size_pdfs_1*(5*blockDim.z*blockIdx.z + 5*threadIdx.z) + _size_pdfs_0*(5*blockDim.y*blockIdx.y + 5*threadIdx.y) + 5*blockDim.x*blockIdx.x + 5*threadIdx.x + 1] = _data_pdfs_10_20_312[_stride_pdfs_0*ctr_0];
      double * const _data_pdfs_10_20_316 = _data_pdfs + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 16*_stride_pdfs_3;
      _data_buffer[_size_pdfs_0*_size_pdfs_1*(5*blockDim.z*blockIdx.z + 5*threadIdx.z) + _size_pdfs_0*(5*blockDim.y*blockIdx.y + 5*threadIdx.y) + 5*blockDim.x*blockIdx.x + 5*threadIdx.x + 2] = _data_pdfs_10_20_316[_stride_pdfs_0*ctr_0];
      double * const _data_pdfs_10_20_32 = _data_pdfs + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 2*_stride_pdfs_3;
      _data_buffer[_size_pdfs_0*_size_pdfs_1*(5*blockDim.z*blockIdx.z + 5*threadIdx.z) + _size_pdfs_0*(5*blockDim.y*blockIdx.y + 5*threadIdx.y) + 5*blockDim.x*blockIdx.x + 5*threadIdx.x + 3] = _data_pdfs_10_20_32[_stride_pdfs_0*ctr_0];
      double * const _data_pdfs_10_20_39 = _data_pdfs + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 9*_stride_pdfs_3;
      _data_buffer[_size_pdfs_0*_size_pdfs_1*(5*blockDim.z*blockIdx.z + 5*threadIdx.z) + _size_pdfs_0*(5*blockDim.y*blockIdx.y + 5*threadIdx.y) + 5*blockDim.x*blockIdx.x + 5*threadIdx.x + 4] = _data_pdfs_10_20_39[_stride_pdfs_0*ctr_0];
   } 
}
}

namespace internal_pack_TS {
static FUNC_PREFIX void pack_TS(double * _data_buffer, double * const _data_pdfs, int64_t const _size_pdfs_0, int64_t const _size_pdfs_1, int64_t const _size_pdfs_2, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2, int64_t const _stride_pdfs_3)
{
   if (blockDim.x*blockIdx.x + threadIdx.x < _size_pdfs_0 && blockDim.y*blockIdx.y + threadIdx.y < _size_pdfs_1 && blockDim.z*blockIdx.z + threadIdx.z < _size_pdfs_2)
   {
      const int64_t ctr_0 = blockDim.x*blockIdx.x + threadIdx.x;
      const int64_t ctr_1 = blockDim.y*blockIdx.y + threadIdx.y;
      const int64_t ctr_2 = blockDim.z*blockIdx.z + threadIdx.z;
      double * const _data_pdfs_10_20_312 = _data_pdfs + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 12*_stride_pdfs_3;
      _data_buffer[_size_pdfs_0*_size_pdfs_1*(blockDim.z*blockIdx.z + threadIdx.z) + _size_pdfs_0*(blockDim.y*blockIdx.y + threadIdx.y) + blockDim.x*blockIdx.x + threadIdx.x] = _data_pdfs_10_20_312[_stride_pdfs_0*ctr_0];
   } 
}
}

namespace internal_pack_B {
static FUNC_PREFIX void pack_B(double * _data_buffer, double * const _data_pdfs, int64_t const _size_pdfs_0, int64_t const _size_pdfs_1, int64_t const _size_pdfs_2, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2, int64_t const _stride_pdfs_3)
{
   if (blockDim.x*blockIdx.x + threadIdx.x < _size_pdfs_0 && blockDim.y*blockIdx.y + threadIdx.y < _size_pdfs_1 && blockDim.z*blockIdx.z + threadIdx.z < _size_pdfs_2)
   {
      const int64_t ctr_0 = blockDim.x*blockIdx.x + threadIdx.x;
      const int64_t ctr_1 = blockDim.y*blockIdx.y + threadIdx.y;
      const int64_t ctr_2 = blockDim.z*blockIdx.z + threadIdx.z;
      double * const _data_pdfs_10_20_315 = _data_pdfs + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 15*_stride_pdfs_3;
      _data_buffer[_size_pdfs_0*_size_pdfs_1*(5*blockDim.z*blockIdx.z + 5*threadIdx.z) + _size_pdfs_0*(5*blockDim.y*blockIdx.y + 5*threadIdx.y) + 5*blockDim.x*blockIdx.x + 5*threadIdx.x] = _data_pdfs_10_20_315[_stride_pdfs_0*ctr_0];
      double * const _data_pdfs_10_20_316 = _data_pdfs + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 16*_stride_pdfs_3;
      _data_buffer[_size_pdfs_0*_size_pdfs_1*(5*blockDim.z*blockIdx.z + 5*threadIdx.z) + _size_pdfs_0*(5*blockDim.y*blockIdx.y + 5*threadIdx.y) + 5*blockDim.x*blockIdx.x + 5*threadIdx.x + 1] = _data_pdfs_10_20_316[_stride_pdfs_0*ctr_0];
      double * const _data_pdfs_10_20_317 = _data_pdfs + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 17*_stride_pdfs_3;
      _data_buffer[_size_pdfs_0*_size_pdfs_1*(5*blockDim.z*blockIdx.z + 5*threadIdx.z) + _size_pdfs_0*(5*blockDim.y*blockIdx.y + 5*threadIdx.y) + 5*blockDim.x*blockIdx.x + 5*threadIdx.x + 2] = _data_pdfs_10_20_317[_stride_pdfs_0*ctr_0];
      double * const _data_pdfs_10_20_318 = _data_pdfs + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 18*_stride_pdfs_3;
      _data_buffer[_size_pdfs_0*_size_pdfs_1*(5*blockDim.z*blockIdx.z + 5*threadIdx.z) + _size_pdfs_0*(5*blockDim.y*blockIdx.y + 5*threadIdx.y) + 5*blockDim.x*blockIdx.x + 5*threadIdx.x + 3] = _data_pdfs_10_20_318[_stride_pdfs_0*ctr_0];
      double * const _data_pdfs_10_20_36 = _data_pdfs + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 6*_stride_pdfs_3;
      _data_buffer[_size_pdfs_0*_size_pdfs_1*(5*blockDim.z*blockIdx.z + 5*threadIdx.z) + _size_pdfs_0*(5*blockDim.y*blockIdx.y + 5*threadIdx.y) + 5*blockDim.x*blockIdx.x + 5*threadIdx.x + 4] = _data_pdfs_10_20_36[_stride_pdfs_0*ctr_0];
   } 
}
}

namespace internal_pack_C {
static FUNC_PREFIX void pack_C(double * _data_buffer, double * const _data_pdfs, int64_t const _size_pdfs_0, int64_t const _size_pdfs_1, int64_t const _size_pdfs_2, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2)
{
   if (blockDim.x*blockIdx.x + threadIdx.x < _size_pdfs_0 && blockDim.y*blockIdx.y + threadIdx.y < _size_pdfs_1 && blockDim.z*blockIdx.z + threadIdx.z < _size_pdfs_2)
   {
      const int64_t ctr_0 = blockDim.x*blockIdx.x + threadIdx.x;
      const int64_t ctr_1 = blockDim.y*blockIdx.y + threadIdx.y;
      const int64_t ctr_2 = blockDim.z*blockIdx.z + threadIdx.z;
      double * const _data_pdfs_10_20_30 = _data_pdfs + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2;
      _data_buffer[_size_pdfs_0*_size_pdfs_1*(blockDim.z*blockIdx.z + threadIdx.z) + _size_pdfs_0*(blockDim.y*blockIdx.y + threadIdx.y) + blockDim.x*blockIdx.x + threadIdx.x] = _data_pdfs_10_20_30[_stride_pdfs_0*ctr_0];
   } 
}
}

namespace internal_pack_T {
static FUNC_PREFIX void pack_T(double * _data_buffer, double * const _data_pdfs, int64_t const _size_pdfs_0, int64_t const _size_pdfs_1, int64_t const _size_pdfs_2, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2, int64_t const _stride_pdfs_3)
{
   if (blockDim.x*blockIdx.x + threadIdx.x < _size_pdfs_0 && blockDim.y*blockIdx.y + threadIdx.y < _size_pdfs_1 && blockDim.z*blockIdx.z + threadIdx.z < _size_pdfs_2)
   {
      const int64_t ctr_0 = blockDim.x*blockIdx.x + threadIdx.x;
      const int64_t ctr_1 = blockDim.y*blockIdx.y + threadIdx.y;
      const int64_t ctr_2 = blockDim.z*blockIdx.z + threadIdx.z;
      double * const _data_pdfs_10_20_311 = _data_pdfs + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 11*_stride_pdfs_3;
      _data_buffer[_size_pdfs_0*_size_pdfs_1*(5*blockDim.z*blockIdx.z + 5*threadIdx.z) + _size_pdfs_0*(5*blockDim.y*blockIdx.y + 5*threadIdx.y) + 5*blockDim.x*blockIdx.x + 5*threadIdx.x] = _data_pdfs_10_20_311[_stride_pdfs_0*ctr_0];
      double * const _data_pdfs_10_20_312 = _data_pdfs + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 12*_stride_pdfs_3;
      _data_buffer[_size_pdfs_0*_size_pdfs_1*(5*blockDim.z*blockIdx.z + 5*threadIdx.z) + _size_pdfs_0*(5*blockDim.y*blockIdx.y + 5*threadIdx.y) + 5*blockDim.x*blockIdx.x + 5*threadIdx.x + 1] = _data_pdfs_10_20_312[_stride_pdfs_0*ctr_0];
      double * const _data_pdfs_10_20_313 = _data_pdfs + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 13*_stride_pdfs_3;
      _data_buffer[_size_pdfs_0*_size_pdfs_1*(5*blockDim.z*blockIdx.z + 5*threadIdx.z) + _size_pdfs_0*(5*blockDim.y*blockIdx.y + 5*threadIdx.y) + 5*blockDim.x*blockIdx.x + 5*threadIdx.x + 2] = _data_pdfs_10_20_313[_stride_pdfs_0*ctr_0];
      double * const _data_pdfs_10_20_314 = _data_pdfs + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 14*_stride_pdfs_3;
      _data_buffer[_size_pdfs_0*_size_pdfs_1*(5*blockDim.z*blockIdx.z + 5*threadIdx.z) + _size_pdfs_0*(5*blockDim.y*blockIdx.y + 5*threadIdx.y) + 5*blockDim.x*blockIdx.x + 5*threadIdx.x + 3] = _data_pdfs_10_20_314[_stride_pdfs_0*ctr_0];
      double * const _data_pdfs_10_20_35 = _data_pdfs + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 5*_stride_pdfs_3;
      _data_buffer[_size_pdfs_0*_size_pdfs_1*(5*blockDim.z*blockIdx.z + 5*threadIdx.z) + _size_pdfs_0*(5*blockDim.y*blockIdx.y + 5*threadIdx.y) + 5*blockDim.x*blockIdx.x + 5*threadIdx.x + 4] = _data_pdfs_10_20_35[_stride_pdfs_0*ctr_0];
   } 
}
}

namespace internal_pack_BN {
static FUNC_PREFIX void pack_BN(double * _data_buffer, double * const _data_pdfs, int64_t const _size_pdfs_0, int64_t const _size_pdfs_1, int64_t const _size_pdfs_2, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2, int64_t const _stride_pdfs_3)
{
   if (blockDim.x*blockIdx.x + threadIdx.x < _size_pdfs_0 && blockDim.y*blockIdx.y + threadIdx.y < _size_pdfs_1 && blockDim.z*blockIdx.z + threadIdx.z < _size_pdfs_2)
   {
      const int64_t ctr_0 = blockDim.x*blockIdx.x + threadIdx.x;
      const int64_t ctr_1 = blockDim.y*blockIdx.y + threadIdx.y;
      const int64_t ctr_2 = blockDim.z*blockIdx.z + threadIdx.z;
      double * const _data_pdfs_10_20_315 = _data_pdfs + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 15*_stride_pdfs_3;
      _data_buffer[_size_pdfs_0*_size_pdfs_1*(blockDim.z*blockIdx.z + threadIdx.z) + _size_pdfs_0*(blockDim.y*blockIdx.y + threadIdx.y) + blockDim.x*blockIdx.x + threadIdx.x] = _data_pdfs_10_20_315[_stride_pdfs_0*ctr_0];
   } 
}
}

namespace internal_pack_N {
static FUNC_PREFIX void pack_N(double * _data_buffer, double * const _data_pdfs, int64_t const _size_pdfs_0, int64_t const _size_pdfs_1, int64_t const _size_pdfs_2, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2, int64_t const _stride_pdfs_3)
{
   if (blockDim.x*blockIdx.x + threadIdx.x < _size_pdfs_0 && blockDim.y*blockIdx.y + threadIdx.y < _size_pdfs_1 && blockDim.z*blockIdx.z + threadIdx.z < _size_pdfs_2)
   {
      const int64_t ctr_0 = blockDim.x*blockIdx.x + threadIdx.x;
      const int64_t ctr_1 = blockDim.y*blockIdx.y + threadIdx.y;
      const int64_t ctr_2 = blockDim.z*blockIdx.z + threadIdx.z;
      double * const _data_pdfs_10_20_31 = _data_pdfs + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + _stride_pdfs_3;
      _data_buffer[_size_pdfs_0*_size_pdfs_1*(5*blockDim.z*blockIdx.z + 5*threadIdx.z) + _size_pdfs_0*(5*blockDim.y*blockIdx.y + 5*threadIdx.y) + 5*blockDim.x*blockIdx.x + 5*threadIdx.x] = _data_pdfs_10_20_31[_stride_pdfs_0*ctr_0];
      double * const _data_pdfs_10_20_311 = _data_pdfs + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 11*_stride_pdfs_3;
      _data_buffer[_size_pdfs_0*_size_pdfs_1*(5*blockDim.z*blockIdx.z + 5*threadIdx.z) + _size_pdfs_0*(5*blockDim.y*blockIdx.y + 5*threadIdx.y) + 5*blockDim.x*blockIdx.x + 5*threadIdx.x + 1] = _data_pdfs_10_20_311[_stride_pdfs_0*ctr_0];
      double * const _data_pdfs_10_20_315 = _data_pdfs + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 15*_stride_pdfs_3;
      _data_buffer[_size_pdfs_0*_size_pdfs_1*(5*blockDim.z*blockIdx.z + 5*threadIdx.z) + _size_pdfs_0*(5*blockDim.y*blockIdx.y + 5*threadIdx.y) + 5*blockDim.x*blockIdx.x + 5*threadIdx.x + 2] = _data_pdfs_10_20_315[_stride_pdfs_0*ctr_0];
      double * const _data_pdfs_10_20_37 = _data_pdfs + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 7*_stride_pdfs_3;
      _data_buffer[_size_pdfs_0*_size_pdfs_1*(5*blockDim.z*blockIdx.z + 5*threadIdx.z) + _size_pdfs_0*(5*blockDim.y*blockIdx.y + 5*threadIdx.y) + 5*blockDim.x*blockIdx.x + 5*threadIdx.x + 3] = _data_pdfs_10_20_37[_stride_pdfs_0*ctr_0];
      double * const _data_pdfs_10_20_38 = _data_pdfs + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 8*_stride_pdfs_3;
      _data_buffer[_size_pdfs_0*_size_pdfs_1*(5*blockDim.z*blockIdx.z + 5*threadIdx.z) + _size_pdfs_0*(5*blockDim.y*blockIdx.y + 5*threadIdx.y) + 5*blockDim.x*blockIdx.x + 5*threadIdx.x + 4] = _data_pdfs_10_20_38[_stride_pdfs_0*ctr_0];
   } 
}
}

namespace internal_pack_TN {
static FUNC_PREFIX void pack_TN(double * _data_buffer, double * const _data_pdfs, int64_t const _size_pdfs_0, int64_t const _size_pdfs_1, int64_t const _size_pdfs_2, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2, int64_t const _stride_pdfs_3)
{
   if (blockDim.x*blockIdx.x + threadIdx.x < _size_pdfs_0 && blockDim.y*blockIdx.y + threadIdx.y < _size_pdfs_1 && blockDim.z*blockIdx.z + threadIdx.z < _size_pdfs_2)
   {
      const int64_t ctr_0 = blockDim.x*blockIdx.x + threadIdx.x;
      const int64_t ctr_1 = blockDim.y*blockIdx.y + threadIdx.y;
      const int64_t ctr_2 = blockDim.z*blockIdx.z + threadIdx.z;
      double * const _data_pdfs_10_20_311 = _data_pdfs + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 11*_stride_pdfs_3;
      _data_buffer[_size_pdfs_0*_size_pdfs_1*(blockDim.z*blockIdx.z + threadIdx.z) + _size_pdfs_0*(blockDim.y*blockIdx.y + threadIdx.y) + blockDim.x*blockIdx.x + threadIdx.x] = _data_pdfs_10_20_311[_stride_pdfs_0*ctr_0];
   } 
}
}

namespace internal_pack_SE {
static FUNC_PREFIX void pack_SE(double * _data_buffer, double * const _data_pdfs, int64_t const _size_pdfs_0, int64_t const _size_pdfs_1, int64_t const _size_pdfs_2, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2, int64_t const _stride_pdfs_3)
{
   if (blockDim.x*blockIdx.x + threadIdx.x < _size_pdfs_0 && blockDim.y*blockIdx.y + threadIdx.y < _size_pdfs_1 && blockDim.z*blockIdx.z + threadIdx.z < _size_pdfs_2)
   {
      const int64_t ctr_0 = blockDim.x*blockIdx.x + threadIdx.x;
      const int64_t ctr_1 = blockDim.y*blockIdx.y + threadIdx.y;
      const int64_t ctr_2 = blockDim.z*blockIdx.z + threadIdx.z;
      double * const _data_pdfs_10_20_310 = _data_pdfs + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 10*_stride_pdfs_3;
      _data_buffer[_size_pdfs_0*_size_pdfs_1*(blockDim.z*blockIdx.z + threadIdx.z) + _size_pdfs_0*(blockDim.y*blockIdx.y + threadIdx.y) + blockDim.x*blockIdx.x + threadIdx.x] = _data_pdfs_10_20_310[_stride_pdfs_0*ctr_0];
   } 
}
}

namespace internal_pack_BE {
static FUNC_PREFIX void pack_BE(double * _data_buffer, double * const _data_pdfs, int64_t const _size_pdfs_0, int64_t const _size_pdfs_1, int64_t const _size_pdfs_2, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2, int64_t const _stride_pdfs_3)
{
   if (blockDim.x*blockIdx.x + threadIdx.x < _size_pdfs_0 && blockDim.y*blockIdx.y + threadIdx.y < _size_pdfs_1 && blockDim.z*blockIdx.z + threadIdx.z < _size_pdfs_2)
   {
      const int64_t ctr_0 = blockDim.x*blockIdx.x + threadIdx.x;
      const int64_t ctr_1 = blockDim.y*blockIdx.y + threadIdx.y;
      const int64_t ctr_2 = blockDim.z*blockIdx.z + threadIdx.z;
      double * const _data_pdfs_10_20_318 = _data_pdfs + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 18*_stride_pdfs_3;
      _data_buffer[_size_pdfs_0*_size_pdfs_1*(blockDim.z*blockIdx.z + threadIdx.z) + _size_pdfs_0*(blockDim.y*blockIdx.y + threadIdx.y) + blockDim.x*blockIdx.x + threadIdx.x] = _data_pdfs_10_20_318[_stride_pdfs_0*ctr_0];
   } 
}
}

namespace internal_pack_E {
static FUNC_PREFIX void pack_E(double * _data_buffer, double * const _data_pdfs, int64_t const _size_pdfs_0, int64_t const _size_pdfs_1, int64_t const _size_pdfs_2, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2, int64_t const _stride_pdfs_3)
{
   if (blockDim.x*blockIdx.x + threadIdx.x < _size_pdfs_0 && blockDim.y*blockIdx.y + threadIdx.y < _size_pdfs_1 && blockDim.z*blockIdx.z + threadIdx.z < _size_pdfs_2)
   {
      const int64_t ctr_0 = blockDim.x*blockIdx.x + threadIdx.x;
      const int64_t ctr_1 = blockDim.y*blockIdx.y + threadIdx.y;
      const int64_t ctr_2 = blockDim.z*blockIdx.z + threadIdx.z;
      double * const _data_pdfs_10_20_310 = _data_pdfs + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 10*_stride_pdfs_3;
      _data_buffer[_size_pdfs_0*_size_pdfs_1*(5*blockDim.z*blockIdx.z + 5*threadIdx.z) + _size_pdfs_0*(5*blockDim.y*blockIdx.y + 5*threadIdx.y) + 5*blockDim.x*blockIdx.x + 5*threadIdx.x] = _data_pdfs_10_20_310[_stride_pdfs_0*ctr_0];
      double * const _data_pdfs_10_20_314 = _data_pdfs + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 14*_stride_pdfs_3;
      _data_buffer[_size_pdfs_0*_size_pdfs_1*(5*blockDim.z*blockIdx.z + 5*threadIdx.z) + _size_pdfs_0*(5*blockDim.y*blockIdx.y + 5*threadIdx.y) + 5*blockDim.x*blockIdx.x + 5*threadIdx.x + 1] = _data_pdfs_10_20_314[_stride_pdfs_0*ctr_0];
      double * const _data_pdfs_10_20_318 = _data_pdfs + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 18*_stride_pdfs_3;
      _data_buffer[_size_pdfs_0*_size_pdfs_1*(5*blockDim.z*blockIdx.z + 5*threadIdx.z) + _size_pdfs_0*(5*blockDim.y*blockIdx.y + 5*threadIdx.y) + 5*blockDim.x*blockIdx.x + 5*threadIdx.x + 2] = _data_pdfs_10_20_318[_stride_pdfs_0*ctr_0];
      double * const _data_pdfs_10_20_34 = _data_pdfs + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 4*_stride_pdfs_3;
      _data_buffer[_size_pdfs_0*_size_pdfs_1*(5*blockDim.z*blockIdx.z + 5*threadIdx.z) + _size_pdfs_0*(5*blockDim.y*blockIdx.y + 5*threadIdx.y) + 5*blockDim.x*blockIdx.x + 5*threadIdx.x + 3] = _data_pdfs_10_20_34[_stride_pdfs_0*ctr_0];
      double * const _data_pdfs_10_20_38 = _data_pdfs + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 8*_stride_pdfs_3;
      _data_buffer[_size_pdfs_0*_size_pdfs_1*(5*blockDim.z*blockIdx.z + 5*threadIdx.z) + _size_pdfs_0*(5*blockDim.y*blockIdx.y + 5*threadIdx.y) + 5*blockDim.x*blockIdx.x + 5*threadIdx.x + 4] = _data_pdfs_10_20_38[_stride_pdfs_0*ctr_0];
   } 
}
}

namespace internal_pack_TE {
static FUNC_PREFIX void pack_TE(double * _data_buffer, double * const _data_pdfs, int64_t const _size_pdfs_0, int64_t const _size_pdfs_1, int64_t const _size_pdfs_2, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2, int64_t const _stride_pdfs_3)
{
   if (blockDim.x*blockIdx.x + threadIdx.x < _size_pdfs_0 && blockDim.y*blockIdx.y + threadIdx.y < _size_pdfs_1 && blockDim.z*blockIdx.z + threadIdx.z < _size_pdfs_2)
   {
      const int64_t ctr_0 = blockDim.x*blockIdx.x + threadIdx.x;
      const int64_t ctr_1 = blockDim.y*blockIdx.y + threadIdx.y;
      const int64_t ctr_2 = blockDim.z*blockIdx.z + threadIdx.z;
      double * const _data_pdfs_10_20_314 = _data_pdfs + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 14*_stride_pdfs_3;
      _data_buffer[_size_pdfs_0*_size_pdfs_1*(blockDim.z*blockIdx.z + threadIdx.z) + _size_pdfs_0*(blockDim.y*blockIdx.y + threadIdx.y) + blockDim.x*blockIdx.x + threadIdx.x] = _data_pdfs_10_20_314[_stride_pdfs_0*ctr_0];
   } 
}
}

namespace internal_pack_NE {
static FUNC_PREFIX void pack_NE(double * _data_buffer, double * const _data_pdfs, int64_t const _size_pdfs_0, int64_t const _size_pdfs_1, int64_t const _size_pdfs_2, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2, int64_t const _stride_pdfs_3)
{
   if (blockDim.x*blockIdx.x + threadIdx.x < _size_pdfs_0 && blockDim.y*blockIdx.y + threadIdx.y < _size_pdfs_1 && blockDim.z*blockIdx.z + threadIdx.z < _size_pdfs_2)
   {
      const int64_t ctr_0 = blockDim.x*blockIdx.x + threadIdx.x;
      const int64_t ctr_1 = blockDim.y*blockIdx.y + threadIdx.y;
      const int64_t ctr_2 = blockDim.z*blockIdx.z + threadIdx.z;
      double * const _data_pdfs_10_20_38 = _data_pdfs + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 8*_stride_pdfs_3;
      _data_buffer[_size_pdfs_0*_size_pdfs_1*(blockDim.z*blockIdx.z + threadIdx.z) + _size_pdfs_0*(blockDim.y*blockIdx.y + threadIdx.y) + blockDim.x*blockIdx.x + threadIdx.x] = _data_pdfs_10_20_38[_stride_pdfs_0*ctr_0];
   } 
}
}



namespace internal_unpack_NE {
static FUNC_PREFIX void unpack_NE(double * const _data_buffer, double * _data_pdfs, int64_t const _size_pdfs_0, int64_t const _size_pdfs_1, int64_t const _size_pdfs_2, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2, int64_t const _stride_pdfs_3)
{
   if (blockDim.x*blockIdx.x + threadIdx.x < _size_pdfs_0 && blockDim.y*blockIdx.y + threadIdx.y < _size_pdfs_1 && blockDim.z*blockIdx.z + threadIdx.z < _size_pdfs_2)
   {
      const int64_t ctr_0 = blockDim.x*blockIdx.x + threadIdx.x;
      const int64_t ctr_1 = blockDim.y*blockIdx.y + threadIdx.y;
      const int64_t ctr_2 = blockDim.z*blockIdx.z + threadIdx.z;
      double * _data_pdfs_10_20_39 = _data_pdfs + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 9*_stride_pdfs_3;
      _data_pdfs_10_20_39[_stride_pdfs_0*ctr_0] = _data_buffer[_size_pdfs_0*_size_pdfs_1*(blockDim.z*blockIdx.z + threadIdx.z) + _size_pdfs_0*(blockDim.y*blockIdx.y + threadIdx.y) + blockDim.x*blockIdx.x + threadIdx.x];
   } 
}
}

namespace internal_unpack_TE {
static FUNC_PREFIX void unpack_TE(double * const _data_buffer, double * _data_pdfs, int64_t const _size_pdfs_0, int64_t const _size_pdfs_1, int64_t const _size_pdfs_2, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2, int64_t const _stride_pdfs_3)
{
   if (blockDim.x*blockIdx.x + threadIdx.x < _size_pdfs_0 && blockDim.y*blockIdx.y + threadIdx.y < _size_pdfs_1 && blockDim.z*blockIdx.z + threadIdx.z < _size_pdfs_2)
   {
      const int64_t ctr_0 = blockDim.x*blockIdx.x + threadIdx.x;
      const int64_t ctr_1 = blockDim.y*blockIdx.y + threadIdx.y;
      const int64_t ctr_2 = blockDim.z*blockIdx.z + threadIdx.z;
      double * _data_pdfs_10_20_317 = _data_pdfs + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 17*_stride_pdfs_3;
      _data_pdfs_10_20_317[_stride_pdfs_0*ctr_0] = _data_buffer[_size_pdfs_0*_size_pdfs_1*(blockDim.z*blockIdx.z + threadIdx.z) + _size_pdfs_0*(blockDim.y*blockIdx.y + threadIdx.y) + blockDim.x*blockIdx.x + threadIdx.x];
   } 
}
}

namespace internal_unpack_E {
static FUNC_PREFIX void unpack_E(double * const _data_buffer, double * _data_pdfs, int64_t const _size_pdfs_0, int64_t const _size_pdfs_1, int64_t const _size_pdfs_2, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2, int64_t const _stride_pdfs_3)
{
   if (blockDim.x*blockIdx.x + threadIdx.x < _size_pdfs_0 && blockDim.y*blockIdx.y + threadIdx.y < _size_pdfs_1 && blockDim.z*blockIdx.z + threadIdx.z < _size_pdfs_2)
   {
      const int64_t ctr_0 = blockDim.x*blockIdx.x + threadIdx.x;
      const int64_t ctr_1 = blockDim.y*blockIdx.y + threadIdx.y;
      const int64_t ctr_2 = blockDim.z*blockIdx.z + threadIdx.z;
      double * _data_pdfs_10_20_313 = _data_pdfs + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 13*_stride_pdfs_3;
      _data_pdfs_10_20_313[_stride_pdfs_0*ctr_0] = _data_buffer[_size_pdfs_0*_size_pdfs_1*(5*blockDim.z*blockIdx.z + 5*threadIdx.z) + _size_pdfs_0*(5*blockDim.y*blockIdx.y + 5*threadIdx.y) + 5*blockDim.x*blockIdx.x + 5*threadIdx.x];
      double * _data_pdfs_10_20_317 = _data_pdfs + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 17*_stride_pdfs_3;
      _data_pdfs_10_20_317[_stride_pdfs_0*ctr_0] = _data_buffer[_size_pdfs_0*_size_pdfs_1*(5*blockDim.z*blockIdx.z + 5*threadIdx.z) + _size_pdfs_0*(5*blockDim.y*blockIdx.y + 5*threadIdx.y) + 5*blockDim.x*blockIdx.x + 5*threadIdx.x + 1];
      double * _data_pdfs_10_20_33 = _data_pdfs + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 3*_stride_pdfs_3;
      _data_pdfs_10_20_33[_stride_pdfs_0*ctr_0] = _data_buffer[_size_pdfs_0*_size_pdfs_1*(5*blockDim.z*blockIdx.z + 5*threadIdx.z) + _size_pdfs_0*(5*blockDim.y*blockIdx.y + 5*threadIdx.y) + 5*blockDim.x*blockIdx.x + 5*threadIdx.x + 2];
      double * _data_pdfs_10_20_37 = _data_pdfs + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 7*_stride_pdfs_3;
      _data_pdfs_10_20_37[_stride_pdfs_0*ctr_0] = _data_buffer[_size_pdfs_0*_size_pdfs_1*(5*blockDim.z*blockIdx.z + 5*threadIdx.z) + _size_pdfs_0*(5*blockDim.y*blockIdx.y + 5*threadIdx.y) + 5*blockDim.x*blockIdx.x + 5*threadIdx.x + 3];
      double * _data_pdfs_10_20_39 = _data_pdfs + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 9*_stride_pdfs_3;
      _data_pdfs_10_20_39[_stride_pdfs_0*ctr_0] = _data_buffer[_size_pdfs_0*_size_pdfs_1*(5*blockDim.z*blockIdx.z + 5*threadIdx.z) + _size_pdfs_0*(5*blockDim.y*blockIdx.y + 5*threadIdx.y) + 5*blockDim.x*blockIdx.x + 5*threadIdx.x + 4];
   } 
}
}

namespace internal_unpack_BE {
static FUNC_PREFIX void unpack_BE(double * const _data_buffer, double * _data_pdfs, int64_t const _size_pdfs_0, int64_t const _size_pdfs_1, int64_t const _size_pdfs_2, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2, int64_t const _stride_pdfs_3)
{
   if (blockDim.x*blockIdx.x + threadIdx.x < _size_pdfs_0 && blockDim.y*blockIdx.y + threadIdx.y < _size_pdfs_1 && blockDim.z*blockIdx.z + threadIdx.z < _size_pdfs_2)
   {
      const int64_t ctr_0 = blockDim.x*blockIdx.x + threadIdx.x;
      const int64_t ctr_1 = blockDim.y*blockIdx.y + threadIdx.y;
      const int64_t ctr_2 = blockDim.z*blockIdx.z + threadIdx.z;
      double * _data_pdfs_10_20_313 = _data_pdfs + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 13*_stride_pdfs_3;
      _data_pdfs_10_20_313[_stride_pdfs_0*ctr_0] = _data_buffer[_size_pdfs_0*_size_pdfs_1*(blockDim.z*blockIdx.z + threadIdx.z) + _size_pdfs_0*(blockDim.y*blockIdx.y + threadIdx.y) + blockDim.x*blockIdx.x + threadIdx.x];
   } 
}
}

namespace internal_unpack_SE {
static FUNC_PREFIX void unpack_SE(double * const _data_buffer, double * _data_pdfs, int64_t const _size_pdfs_0, int64_t const _size_pdfs_1, int64_t const _size_pdfs_2, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2, int64_t const _stride_pdfs_3)
{
   if (blockDim.x*blockIdx.x + threadIdx.x < _size_pdfs_0 && blockDim.y*blockIdx.y + threadIdx.y < _size_pdfs_1 && blockDim.z*blockIdx.z + threadIdx.z < _size_pdfs_2)
   {
      const int64_t ctr_0 = blockDim.x*blockIdx.x + threadIdx.x;
      const int64_t ctr_1 = blockDim.y*blockIdx.y + threadIdx.y;
      const int64_t ctr_2 = blockDim.z*blockIdx.z + threadIdx.z;
      double * _data_pdfs_10_20_37 = _data_pdfs + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 7*_stride_pdfs_3;
      _data_pdfs_10_20_37[_stride_pdfs_0*ctr_0] = _data_buffer[_size_pdfs_0*_size_pdfs_1*(blockDim.z*blockIdx.z + threadIdx.z) + _size_pdfs_0*(blockDim.y*blockIdx.y + threadIdx.y) + blockDim.x*blockIdx.x + threadIdx.x];
   } 
}
}

namespace internal_unpack_TN {
static FUNC_PREFIX void unpack_TN(double * const _data_buffer, double * _data_pdfs, int64_t const _size_pdfs_0, int64_t const _size_pdfs_1, int64_t const _size_pdfs_2, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2, int64_t const _stride_pdfs_3)
{
   if (blockDim.x*blockIdx.x + threadIdx.x < _size_pdfs_0 && blockDim.y*blockIdx.y + threadIdx.y < _size_pdfs_1 && blockDim.z*blockIdx.z + threadIdx.z < _size_pdfs_2)
   {
      const int64_t ctr_0 = blockDim.x*blockIdx.x + threadIdx.x;
      const int64_t ctr_1 = blockDim.y*blockIdx.y + threadIdx.y;
      const int64_t ctr_2 = blockDim.z*blockIdx.z + threadIdx.z;
      double * _data_pdfs_10_20_316 = _data_pdfs + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 16*_stride_pdfs_3;
      _data_pdfs_10_20_316[_stride_pdfs_0*ctr_0] = _data_buffer[_size_pdfs_0*_size_pdfs_1*(blockDim.z*blockIdx.z + threadIdx.z) + _size_pdfs_0*(blockDim.y*blockIdx.y + threadIdx.y) + blockDim.x*blockIdx.x + threadIdx.x];
   } 
}
}

namespace internal_unpack_N {
static FUNC_PREFIX void unpack_N(double * const _data_buffer, double * _data_pdfs, int64_t const _size_pdfs_0, int64_t const _size_pdfs_1, int64_t const _size_pdfs_2, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2, int64_t const _stride_pdfs_3)
{
   if (blockDim.x*blockIdx.x + threadIdx.x < _size_pdfs_0 && blockDim.y*blockIdx.y + threadIdx.y < _size_pdfs_1 && blockDim.z*blockIdx.z + threadIdx.z < _size_pdfs_2)
   {
      const int64_t ctr_0 = blockDim.x*blockIdx.x + threadIdx.x;
      const int64_t ctr_1 = blockDim.y*blockIdx.y + threadIdx.y;
      const int64_t ctr_2 = blockDim.z*blockIdx.z + threadIdx.z;
      double * _data_pdfs_10_20_310 = _data_pdfs + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 10*_stride_pdfs_3;
      _data_pdfs_10_20_310[_stride_pdfs_0*ctr_0] = _data_buffer[_size_pdfs_0*_size_pdfs_1*(5*blockDim.z*blockIdx.z + 5*threadIdx.z) + _size_pdfs_0*(5*blockDim.y*blockIdx.y + 5*threadIdx.y) + 5*blockDim.x*blockIdx.x + 5*threadIdx.x];
      double * _data_pdfs_10_20_312 = _data_pdfs + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 12*_stride_pdfs_3;
      _data_pdfs_10_20_312[_stride_pdfs_0*ctr_0] = _data_buffer[_size_pdfs_0*_size_pdfs_1*(5*blockDim.z*blockIdx.z + 5*threadIdx.z) + _size_pdfs_0*(5*blockDim.y*blockIdx.y + 5*threadIdx.y) + 5*blockDim.x*blockIdx.x + 5*threadIdx.x + 1];
      double * _data_pdfs_10_20_316 = _data_pdfs + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 16*_stride_pdfs_3;
      _data_pdfs_10_20_316[_stride_pdfs_0*ctr_0] = _data_buffer[_size_pdfs_0*_size_pdfs_1*(5*blockDim.z*blockIdx.z + 5*threadIdx.z) + _size_pdfs_0*(5*blockDim.y*blockIdx.y + 5*threadIdx.y) + 5*blockDim.x*blockIdx.x + 5*threadIdx.x + 2];
      double * _data_pdfs_10_20_32 = _data_pdfs + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 2*_stride_pdfs_3;
      _data_pdfs_10_20_32[_stride_pdfs_0*ctr_0] = _data_buffer[_size_pdfs_0*_size_pdfs_1*(5*blockDim.z*blockIdx.z + 5*threadIdx.z) + _size_pdfs_0*(5*blockDim.y*blockIdx.y + 5*threadIdx.y) + 5*blockDim.x*blockIdx.x + 5*threadIdx.x + 3];
      double * _data_pdfs_10_20_39 = _data_pdfs + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 9*_stride_pdfs_3;
      _data_pdfs_10_20_39[_stride_pdfs_0*ctr_0] = _data_buffer[_size_pdfs_0*_size_pdfs_1*(5*blockDim.z*blockIdx.z + 5*threadIdx.z) + _size_pdfs_0*(5*blockDim.y*blockIdx.y + 5*threadIdx.y) + 5*blockDim.x*blockIdx.x + 5*threadIdx.x + 4];
   } 
}
}

namespace internal_unpack_BN {
static FUNC_PREFIX void unpack_BN(double * const _data_buffer, double * _data_pdfs, int64_t const _size_pdfs_0, int64_t const _size_pdfs_1, int64_t const _size_pdfs_2, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2, int64_t const _stride_pdfs_3)
{
   if (blockDim.x*blockIdx.x + threadIdx.x < _size_pdfs_0 && blockDim.y*blockIdx.y + threadIdx.y < _size_pdfs_1 && blockDim.z*blockIdx.z + threadIdx.z < _size_pdfs_2)
   {
      const int64_t ctr_0 = blockDim.x*blockIdx.x + threadIdx.x;
      const int64_t ctr_1 = blockDim.y*blockIdx.y + threadIdx.y;
      const int64_t ctr_2 = blockDim.z*blockIdx.z + threadIdx.z;
      double * _data_pdfs_10_20_312 = _data_pdfs + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 12*_stride_pdfs_3;
      _data_pdfs_10_20_312[_stride_pdfs_0*ctr_0] = _data_buffer[_size_pdfs_0*_size_pdfs_1*(blockDim.z*blockIdx.z + threadIdx.z) + _size_pdfs_0*(blockDim.y*blockIdx.y + threadIdx.y) + blockDim.x*blockIdx.x + threadIdx.x];
   } 
}
}

namespace internal_unpack_T {
static FUNC_PREFIX void unpack_T(double * const _data_buffer, double * _data_pdfs, int64_t const _size_pdfs_0, int64_t const _size_pdfs_1, int64_t const _size_pdfs_2, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2, int64_t const _stride_pdfs_3)
{
   if (blockDim.x*blockIdx.x + threadIdx.x < _size_pdfs_0 && blockDim.y*blockIdx.y + threadIdx.y < _size_pdfs_1 && blockDim.z*blockIdx.z + threadIdx.z < _size_pdfs_2)
   {
      const int64_t ctr_0 = blockDim.x*blockIdx.x + threadIdx.x;
      const int64_t ctr_1 = blockDim.y*blockIdx.y + threadIdx.y;
      const int64_t ctr_2 = blockDim.z*blockIdx.z + threadIdx.z;
      double * _data_pdfs_10_20_315 = _data_pdfs + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 15*_stride_pdfs_3;
      _data_pdfs_10_20_315[_stride_pdfs_0*ctr_0] = _data_buffer[_size_pdfs_0*_size_pdfs_1*(5*blockDim.z*blockIdx.z + 5*threadIdx.z) + _size_pdfs_0*(5*blockDim.y*blockIdx.y + 5*threadIdx.y) + 5*blockDim.x*blockIdx.x + 5*threadIdx.x];
      double * _data_pdfs_10_20_316 = _data_pdfs + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 16*_stride_pdfs_3;
      _data_pdfs_10_20_316[_stride_pdfs_0*ctr_0] = _data_buffer[_size_pdfs_0*_size_pdfs_1*(5*blockDim.z*blockIdx.z + 5*threadIdx.z) + _size_pdfs_0*(5*blockDim.y*blockIdx.y + 5*threadIdx.y) + 5*blockDim.x*blockIdx.x + 5*threadIdx.x + 1];
      double * _data_pdfs_10_20_317 = _data_pdfs + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 17*_stride_pdfs_3;
      _data_pdfs_10_20_317[_stride_pdfs_0*ctr_0] = _data_buffer[_size_pdfs_0*_size_pdfs_1*(5*blockDim.z*blockIdx.z + 5*threadIdx.z) + _size_pdfs_0*(5*blockDim.y*blockIdx.y + 5*threadIdx.y) + 5*blockDim.x*blockIdx.x + 5*threadIdx.x + 2];
      double * _data_pdfs_10_20_318 = _data_pdfs + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 18*_stride_pdfs_3;
      _data_pdfs_10_20_318[_stride_pdfs_0*ctr_0] = _data_buffer[_size_pdfs_0*_size_pdfs_1*(5*blockDim.z*blockIdx.z + 5*threadIdx.z) + _size_pdfs_0*(5*blockDim.y*blockIdx.y + 5*threadIdx.y) + 5*blockDim.x*blockIdx.x + 5*threadIdx.x + 3];
      double * _data_pdfs_10_20_36 = _data_pdfs + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 6*_stride_pdfs_3;
      _data_pdfs_10_20_36[_stride_pdfs_0*ctr_0] = _data_buffer[_size_pdfs_0*_size_pdfs_1*(5*blockDim.z*blockIdx.z + 5*threadIdx.z) + _size_pdfs_0*(5*blockDim.y*blockIdx.y + 5*threadIdx.y) + 5*blockDim.x*blockIdx.x + 5*threadIdx.x + 4];
   } 
}
}

namespace internal_unpack_C {
static FUNC_PREFIX void unpack_C(double * const _data_buffer, double * _data_pdfs, int64_t const _size_pdfs_0, int64_t const _size_pdfs_1, int64_t const _size_pdfs_2, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2)
{
   if (blockDim.x*blockIdx.x + threadIdx.x < _size_pdfs_0 && blockDim.y*blockIdx.y + threadIdx.y < _size_pdfs_1 && blockDim.z*blockIdx.z + threadIdx.z < _size_pdfs_2)
   {
      const int64_t ctr_0 = blockDim.x*blockIdx.x + threadIdx.x;
      const int64_t ctr_1 = blockDim.y*blockIdx.y + threadIdx.y;
      const int64_t ctr_2 = blockDim.z*blockIdx.z + threadIdx.z;
      double * _data_pdfs_10_20_30 = _data_pdfs + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2;
      _data_pdfs_10_20_30[_stride_pdfs_0*ctr_0] = _data_buffer[_size_pdfs_0*_size_pdfs_1*(blockDim.z*blockIdx.z + threadIdx.z) + _size_pdfs_0*(blockDim.y*blockIdx.y + threadIdx.y) + blockDim.x*blockIdx.x + threadIdx.x];
   } 
}
}

namespace internal_unpack_B {
static FUNC_PREFIX void unpack_B(double * const _data_buffer, double * _data_pdfs, int64_t const _size_pdfs_0, int64_t const _size_pdfs_1, int64_t const _size_pdfs_2, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2, int64_t const _stride_pdfs_3)
{
   if (blockDim.x*blockIdx.x + threadIdx.x < _size_pdfs_0 && blockDim.y*blockIdx.y + threadIdx.y < _size_pdfs_1 && blockDim.z*blockIdx.z + threadIdx.z < _size_pdfs_2)
   {
      const int64_t ctr_0 = blockDim.x*blockIdx.x + threadIdx.x;
      const int64_t ctr_1 = blockDim.y*blockIdx.y + threadIdx.y;
      const int64_t ctr_2 = blockDim.z*blockIdx.z + threadIdx.z;
      double * _data_pdfs_10_20_311 = _data_pdfs + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 11*_stride_pdfs_3;
      _data_pdfs_10_20_311[_stride_pdfs_0*ctr_0] = _data_buffer[_size_pdfs_0*_size_pdfs_1*(5*blockDim.z*blockIdx.z + 5*threadIdx.z) + _size_pdfs_0*(5*blockDim.y*blockIdx.y + 5*threadIdx.y) + 5*blockDim.x*blockIdx.x + 5*threadIdx.x];
      double * _data_pdfs_10_20_312 = _data_pdfs + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 12*_stride_pdfs_3;
      _data_pdfs_10_20_312[_stride_pdfs_0*ctr_0] = _data_buffer[_size_pdfs_0*_size_pdfs_1*(5*blockDim.z*blockIdx.z + 5*threadIdx.z) + _size_pdfs_0*(5*blockDim.y*blockIdx.y + 5*threadIdx.y) + 5*blockDim.x*blockIdx.x + 5*threadIdx.x + 1];
      double * _data_pdfs_10_20_313 = _data_pdfs + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 13*_stride_pdfs_3;
      _data_pdfs_10_20_313[_stride_pdfs_0*ctr_0] = _data_buffer[_size_pdfs_0*_size_pdfs_1*(5*blockDim.z*blockIdx.z + 5*threadIdx.z) + _size_pdfs_0*(5*blockDim.y*blockIdx.y + 5*threadIdx.y) + 5*blockDim.x*blockIdx.x + 5*threadIdx.x + 2];
      double * _data_pdfs_10_20_314 = _data_pdfs + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 14*_stride_pdfs_3;
      _data_pdfs_10_20_314[_stride_pdfs_0*ctr_0] = _data_buffer[_size_pdfs_0*_size_pdfs_1*(5*blockDim.z*blockIdx.z + 5*threadIdx.z) + _size_pdfs_0*(5*blockDim.y*blockIdx.y + 5*threadIdx.y) + 5*blockDim.x*blockIdx.x + 5*threadIdx.x + 3];
      double * _data_pdfs_10_20_35 = _data_pdfs + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 5*_stride_pdfs_3;
      _data_pdfs_10_20_35[_stride_pdfs_0*ctr_0] = _data_buffer[_size_pdfs_0*_size_pdfs_1*(5*blockDim.z*blockIdx.z + 5*threadIdx.z) + _size_pdfs_0*(5*blockDim.y*blockIdx.y + 5*threadIdx.y) + 5*blockDim.x*blockIdx.x + 5*threadIdx.x + 4];
   } 
}
}

namespace internal_unpack_TS {
static FUNC_PREFIX void unpack_TS(double * const _data_buffer, double * _data_pdfs, int64_t const _size_pdfs_0, int64_t const _size_pdfs_1, int64_t const _size_pdfs_2, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2, int64_t const _stride_pdfs_3)
{
   if (blockDim.x*blockIdx.x + threadIdx.x < _size_pdfs_0 && blockDim.y*blockIdx.y + threadIdx.y < _size_pdfs_1 && blockDim.z*blockIdx.z + threadIdx.z < _size_pdfs_2)
   {
      const int64_t ctr_0 = blockDim.x*blockIdx.x + threadIdx.x;
      const int64_t ctr_1 = blockDim.y*blockIdx.y + threadIdx.y;
      const int64_t ctr_2 = blockDim.z*blockIdx.z + threadIdx.z;
      double * _data_pdfs_10_20_315 = _data_pdfs + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 15*_stride_pdfs_3;
      _data_pdfs_10_20_315[_stride_pdfs_0*ctr_0] = _data_buffer[_size_pdfs_0*_size_pdfs_1*(blockDim.z*blockIdx.z + threadIdx.z) + _size_pdfs_0*(blockDim.y*blockIdx.y + threadIdx.y) + blockDim.x*blockIdx.x + threadIdx.x];
   } 
}
}

namespace internal_unpack_S {
static FUNC_PREFIX void unpack_S(double * const _data_buffer, double * _data_pdfs, int64_t const _size_pdfs_0, int64_t const _size_pdfs_1, int64_t const _size_pdfs_2, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2, int64_t const _stride_pdfs_3)
{
   if (blockDim.x*blockIdx.x + threadIdx.x < _size_pdfs_0 && blockDim.y*blockIdx.y + threadIdx.y < _size_pdfs_1 && blockDim.z*blockIdx.z + threadIdx.z < _size_pdfs_2)
   {
      const int64_t ctr_0 = blockDim.x*blockIdx.x + threadIdx.x;
      const int64_t ctr_1 = blockDim.y*blockIdx.y + threadIdx.y;
      const int64_t ctr_2 = blockDim.z*blockIdx.z + threadIdx.z;
      double * _data_pdfs_10_20_31 = _data_pdfs + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + _stride_pdfs_3;
      _data_pdfs_10_20_31[_stride_pdfs_0*ctr_0] = _data_buffer[_size_pdfs_0*_size_pdfs_1*(5*blockDim.z*blockIdx.z + 5*threadIdx.z) + _size_pdfs_0*(5*blockDim.y*blockIdx.y + 5*threadIdx.y) + 5*blockDim.x*blockIdx.x + 5*threadIdx.x];
      double * _data_pdfs_10_20_311 = _data_pdfs + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 11*_stride_pdfs_3;
      _data_pdfs_10_20_311[_stride_pdfs_0*ctr_0] = _data_buffer[_size_pdfs_0*_size_pdfs_1*(5*blockDim.z*blockIdx.z + 5*threadIdx.z) + _size_pdfs_0*(5*blockDim.y*blockIdx.y + 5*threadIdx.y) + 5*blockDim.x*blockIdx.x + 5*threadIdx.x + 1];
      double * _data_pdfs_10_20_315 = _data_pdfs + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 15*_stride_pdfs_3;
      _data_pdfs_10_20_315[_stride_pdfs_0*ctr_0] = _data_buffer[_size_pdfs_0*_size_pdfs_1*(5*blockDim.z*blockIdx.z + 5*threadIdx.z) + _size_pdfs_0*(5*blockDim.y*blockIdx.y + 5*threadIdx.y) + 5*blockDim.x*blockIdx.x + 5*threadIdx.x + 2];
      double * _data_pdfs_10_20_37 = _data_pdfs + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 7*_stride_pdfs_3;
      _data_pdfs_10_20_37[_stride_pdfs_0*ctr_0] = _data_buffer[_size_pdfs_0*_size_pdfs_1*(5*blockDim.z*blockIdx.z + 5*threadIdx.z) + _size_pdfs_0*(5*blockDim.y*blockIdx.y + 5*threadIdx.y) + 5*blockDim.x*blockIdx.x + 5*threadIdx.x + 3];
      double * _data_pdfs_10_20_38 = _data_pdfs + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 8*_stride_pdfs_3;
      _data_pdfs_10_20_38[_stride_pdfs_0*ctr_0] = _data_buffer[_size_pdfs_0*_size_pdfs_1*(5*blockDim.z*blockIdx.z + 5*threadIdx.z) + _size_pdfs_0*(5*blockDim.y*blockIdx.y + 5*threadIdx.y) + 5*blockDim.x*blockIdx.x + 5*threadIdx.x + 4];
   } 
}
}

namespace internal_unpack_BS {
static FUNC_PREFIX void unpack_BS(double * const _data_buffer, double * _data_pdfs, int64_t const _size_pdfs_0, int64_t const _size_pdfs_1, int64_t const _size_pdfs_2, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2, int64_t const _stride_pdfs_3)
{
   if (blockDim.x*blockIdx.x + threadIdx.x < _size_pdfs_0 && blockDim.y*blockIdx.y + threadIdx.y < _size_pdfs_1 && blockDim.z*blockIdx.z + threadIdx.z < _size_pdfs_2)
   {
      const int64_t ctr_0 = blockDim.x*blockIdx.x + threadIdx.x;
      const int64_t ctr_1 = blockDim.y*blockIdx.y + threadIdx.y;
      const int64_t ctr_2 = blockDim.z*blockIdx.z + threadIdx.z;
      double * _data_pdfs_10_20_311 = _data_pdfs + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 11*_stride_pdfs_3;
      _data_pdfs_10_20_311[_stride_pdfs_0*ctr_0] = _data_buffer[_size_pdfs_0*_size_pdfs_1*(blockDim.z*blockIdx.z + threadIdx.z) + _size_pdfs_0*(blockDim.y*blockIdx.y + threadIdx.y) + blockDim.x*blockIdx.x + threadIdx.x];
   } 
}
}

namespace internal_unpack_NW {
static FUNC_PREFIX void unpack_NW(double * const _data_buffer, double * _data_pdfs, int64_t const _size_pdfs_0, int64_t const _size_pdfs_1, int64_t const _size_pdfs_2, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2, int64_t const _stride_pdfs_3)
{
   if (blockDim.x*blockIdx.x + threadIdx.x < _size_pdfs_0 && blockDim.y*blockIdx.y + threadIdx.y < _size_pdfs_1 && blockDim.z*blockIdx.z + threadIdx.z < _size_pdfs_2)
   {
      const int64_t ctr_0 = blockDim.x*blockIdx.x + threadIdx.x;
      const int64_t ctr_1 = blockDim.y*blockIdx.y + threadIdx.y;
      const int64_t ctr_2 = blockDim.z*blockIdx.z + threadIdx.z;
      double * _data_pdfs_10_20_310 = _data_pdfs + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 10*_stride_pdfs_3;
      _data_pdfs_10_20_310[_stride_pdfs_0*ctr_0] = _data_buffer[_size_pdfs_0*_size_pdfs_1*(blockDim.z*blockIdx.z + threadIdx.z) + _size_pdfs_0*(blockDim.y*blockIdx.y + threadIdx.y) + blockDim.x*blockIdx.x + threadIdx.x];
   } 
}
}

namespace internal_unpack_TW {
static FUNC_PREFIX void unpack_TW(double * const _data_buffer, double * _data_pdfs, int64_t const _size_pdfs_0, int64_t const _size_pdfs_1, int64_t const _size_pdfs_2, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2, int64_t const _stride_pdfs_3)
{
   if (blockDim.x*blockIdx.x + threadIdx.x < _size_pdfs_0 && blockDim.y*blockIdx.y + threadIdx.y < _size_pdfs_1 && blockDim.z*blockIdx.z + threadIdx.z < _size_pdfs_2)
   {
      const int64_t ctr_0 = blockDim.x*blockIdx.x + threadIdx.x;
      const int64_t ctr_1 = blockDim.y*blockIdx.y + threadIdx.y;
      const int64_t ctr_2 = blockDim.z*blockIdx.z + threadIdx.z;
      double * _data_pdfs_10_20_318 = _data_pdfs + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 18*_stride_pdfs_3;
      _data_pdfs_10_20_318[_stride_pdfs_0*ctr_0] = _data_buffer[_size_pdfs_0*_size_pdfs_1*(blockDim.z*blockIdx.z + threadIdx.z) + _size_pdfs_0*(blockDim.y*blockIdx.y + threadIdx.y) + blockDim.x*blockIdx.x + threadIdx.x];
   } 
}
}

namespace internal_unpack_W {
static FUNC_PREFIX void unpack_W(double * const _data_buffer, double * _data_pdfs, int64_t const _size_pdfs_0, int64_t const _size_pdfs_1, int64_t const _size_pdfs_2, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2, int64_t const _stride_pdfs_3)
{
   if (blockDim.x*blockIdx.x + threadIdx.x < _size_pdfs_0 && blockDim.y*blockIdx.y + threadIdx.y < _size_pdfs_1 && blockDim.z*blockIdx.z + threadIdx.z < _size_pdfs_2)
   {
      const int64_t ctr_0 = blockDim.x*blockIdx.x + threadIdx.x;
      const int64_t ctr_1 = blockDim.y*blockIdx.y + threadIdx.y;
      const int64_t ctr_2 = blockDim.z*blockIdx.z + threadIdx.z;
      double * _data_pdfs_10_20_310 = _data_pdfs + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 10*_stride_pdfs_3;
      _data_pdfs_10_20_310[_stride_pdfs_0*ctr_0] = _data_buffer[_size_pdfs_0*_size_pdfs_1*(5*blockDim.z*blockIdx.z + 5*threadIdx.z) + _size_pdfs_0*(5*blockDim.y*blockIdx.y + 5*threadIdx.y) + 5*blockDim.x*blockIdx.x + 5*threadIdx.x];
      double * _data_pdfs_10_20_314 = _data_pdfs + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 14*_stride_pdfs_3;
      _data_pdfs_10_20_314[_stride_pdfs_0*ctr_0] = _data_buffer[_size_pdfs_0*_size_pdfs_1*(5*blockDim.z*blockIdx.z + 5*threadIdx.z) + _size_pdfs_0*(5*blockDim.y*blockIdx.y + 5*threadIdx.y) + 5*blockDim.x*blockIdx.x + 5*threadIdx.x + 1];
      double * _data_pdfs_10_20_318 = _data_pdfs + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 18*_stride_pdfs_3;
      _data_pdfs_10_20_318[_stride_pdfs_0*ctr_0] = _data_buffer[_size_pdfs_0*_size_pdfs_1*(5*blockDim.z*blockIdx.z + 5*threadIdx.z) + _size_pdfs_0*(5*blockDim.y*blockIdx.y + 5*threadIdx.y) + 5*blockDim.x*blockIdx.x + 5*threadIdx.x + 2];
      double * _data_pdfs_10_20_34 = _data_pdfs + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 4*_stride_pdfs_3;
      _data_pdfs_10_20_34[_stride_pdfs_0*ctr_0] = _data_buffer[_size_pdfs_0*_size_pdfs_1*(5*blockDim.z*blockIdx.z + 5*threadIdx.z) + _size_pdfs_0*(5*blockDim.y*blockIdx.y + 5*threadIdx.y) + 5*blockDim.x*blockIdx.x + 5*threadIdx.x + 3];
      double * _data_pdfs_10_20_38 = _data_pdfs + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 8*_stride_pdfs_3;
      _data_pdfs_10_20_38[_stride_pdfs_0*ctr_0] = _data_buffer[_size_pdfs_0*_size_pdfs_1*(5*blockDim.z*blockIdx.z + 5*threadIdx.z) + _size_pdfs_0*(5*blockDim.y*blockIdx.y + 5*threadIdx.y) + 5*blockDim.x*blockIdx.x + 5*threadIdx.x + 4];
   } 
}
}

namespace internal_unpack_BW {
static FUNC_PREFIX void unpack_BW(double * const _data_buffer, double * _data_pdfs, int64_t const _size_pdfs_0, int64_t const _size_pdfs_1, int64_t const _size_pdfs_2, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2, int64_t const _stride_pdfs_3)
{
   if (blockDim.x*blockIdx.x + threadIdx.x < _size_pdfs_0 && blockDim.y*blockIdx.y + threadIdx.y < _size_pdfs_1 && blockDim.z*blockIdx.z + threadIdx.z < _size_pdfs_2)
   {
      const int64_t ctr_0 = blockDim.x*blockIdx.x + threadIdx.x;
      const int64_t ctr_1 = blockDim.y*blockIdx.y + threadIdx.y;
      const int64_t ctr_2 = blockDim.z*blockIdx.z + threadIdx.z;
      double * _data_pdfs_10_20_314 = _data_pdfs + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 14*_stride_pdfs_3;
      _data_pdfs_10_20_314[_stride_pdfs_0*ctr_0] = _data_buffer[_size_pdfs_0*_size_pdfs_1*(blockDim.z*blockIdx.z + threadIdx.z) + _size_pdfs_0*(blockDim.y*blockIdx.y + threadIdx.y) + blockDim.x*blockIdx.x + threadIdx.x];
   } 
}
}

namespace internal_unpack_SW {
static FUNC_PREFIX void unpack_SW(double * const _data_buffer, double * _data_pdfs, int64_t const _size_pdfs_0, int64_t const _size_pdfs_1, int64_t const _size_pdfs_2, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2, int64_t const _stride_pdfs_3)
{
   if (blockDim.x*blockIdx.x + threadIdx.x < _size_pdfs_0 && blockDim.y*blockIdx.y + threadIdx.y < _size_pdfs_1 && blockDim.z*blockIdx.z + threadIdx.z < _size_pdfs_2)
   {
      const int64_t ctr_0 = blockDim.x*blockIdx.x + threadIdx.x;
      const int64_t ctr_1 = blockDim.y*blockIdx.y + threadIdx.y;
      const int64_t ctr_2 = blockDim.z*blockIdx.z + threadIdx.z;
      double * _data_pdfs_10_20_38 = _data_pdfs + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 8*_stride_pdfs_3;
      _data_pdfs_10_20_38[_stride_pdfs_0*ctr_0] = _data_buffer[_size_pdfs_0*_size_pdfs_1*(blockDim.z*blockIdx.z + threadIdx.z) + _size_pdfs_0*(blockDim.y*blockIdx.y + threadIdx.y) + blockDim.x*blockIdx.x + threadIdx.x];
   } 
}
}




void UniformGridGPU_PackInfo::pack(Direction dir, unsigned char * byte_buffer, IBlock * block, cudaStream_t stream)
{
    double * buffer = reinterpret_cast<double*>(byte_buffer);

    auto pdfs = block->getData< cuda::GPUField<double> >(pdfsID);

    CellInterval ci;
    pdfs->getSliceBeforeGhostLayer(dir, ci, 1, false);

    switch( dir )
    {
        case stencil::SW:
        {
            double * _data_buffer = buffer;
            WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()));
            WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()));
            WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()));
            double * const _data_pdfs = pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(), int64_t(cell_idx_c(ci.xSize()) + 0));
            const int64_t _size_pdfs_0 = int64_t(cell_idx_c(ci.xSize()) + 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(), int64_t(cell_idx_c(ci.ySize()) + 0));
            const int64_t _size_pdfs_1 = int64_t(cell_idx_c(ci.ySize()) + 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(), int64_t(cell_idx_c(ci.zSize()) + 0));
            const int64_t _size_pdfs_2 = int64_t(cell_idx_c(ci.zSize()) + 0);
            const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
            const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
            const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
            const int64_t _stride_pdfs_3 = int64_t(pdfs->fStride());
            dim3 _block(int(((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)), int(((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)), int(((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)));
            dim3 _grid(int(( (_size_pdfs_0) % (((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)) == 0 ? (int64_t)(_size_pdfs_0) / (int64_t)(((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)) : ( (int64_t)(_size_pdfs_0) / (int64_t)(((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)) ) +1 )), int(( (_size_pdfs_1) % (((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)) == 0 ? (int64_t)(_size_pdfs_1) / (int64_t)(((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)) : ( (int64_t)(_size_pdfs_1) / (int64_t)(((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)) ) +1 )), int(( (_size_pdfs_2) % (((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)) == 0 ? (int64_t)(_size_pdfs_2) / (int64_t)(((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)) : ( (int64_t)(_size_pdfs_2) / (int64_t)(((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)) ) +1 )));
            internal_pack_SW::pack_SW<<<_grid, _block, 0, stream>>>(_data_buffer, _data_pdfs, _size_pdfs_0, _size_pdfs_1, _size_pdfs_2, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3);
            break;
        }
        
        case stencil::BW:
        {
            double * _data_buffer = buffer;
            WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()));
            WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()));
            WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()));
            double * const _data_pdfs = pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(), int64_t(cell_idx_c(ci.xSize()) + 0));
            const int64_t _size_pdfs_0 = int64_t(cell_idx_c(ci.xSize()) + 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(), int64_t(cell_idx_c(ci.ySize()) + 0));
            const int64_t _size_pdfs_1 = int64_t(cell_idx_c(ci.ySize()) + 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(), int64_t(cell_idx_c(ci.zSize()) + 0));
            const int64_t _size_pdfs_2 = int64_t(cell_idx_c(ci.zSize()) + 0);
            const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
            const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
            const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
            const int64_t _stride_pdfs_3 = int64_t(pdfs->fStride());
            dim3 _block(int(((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)), int(((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)), int(((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)));
            dim3 _grid(int(( (_size_pdfs_0) % (((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)) == 0 ? (int64_t)(_size_pdfs_0) / (int64_t)(((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)) : ( (int64_t)(_size_pdfs_0) / (int64_t)(((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)) ) +1 )), int(( (_size_pdfs_1) % (((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)) == 0 ? (int64_t)(_size_pdfs_1) / (int64_t)(((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)) : ( (int64_t)(_size_pdfs_1) / (int64_t)(((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)) ) +1 )), int(( (_size_pdfs_2) % (((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)) == 0 ? (int64_t)(_size_pdfs_2) / (int64_t)(((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)) : ( (int64_t)(_size_pdfs_2) / (int64_t)(((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)) ) +1 )));
            internal_pack_BW::pack_BW<<<_grid, _block, 0, stream>>>(_data_buffer, _data_pdfs, _size_pdfs_0, _size_pdfs_1, _size_pdfs_2, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3);
            break;
        }
        
        case stencil::W:
        {
            double * _data_buffer = buffer;
            WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()));
            WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()));
            WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()));
            double * const _data_pdfs = pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(), int64_t(cell_idx_c(ci.xSize()) + 0));
            const int64_t _size_pdfs_0 = int64_t(cell_idx_c(ci.xSize()) + 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(), int64_t(cell_idx_c(ci.ySize()) + 0));
            const int64_t _size_pdfs_1 = int64_t(cell_idx_c(ci.ySize()) + 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(), int64_t(cell_idx_c(ci.zSize()) + 0));
            const int64_t _size_pdfs_2 = int64_t(cell_idx_c(ci.zSize()) + 0);
            const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
            const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
            const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
            const int64_t _stride_pdfs_3 = int64_t(pdfs->fStride());
            dim3 _block(int(((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)), int(((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)), int(((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)));
            dim3 _grid(int(( (_size_pdfs_0) % (((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)) == 0 ? (int64_t)(_size_pdfs_0) / (int64_t)(((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)) : ( (int64_t)(_size_pdfs_0) / (int64_t)(((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)) ) +1 )), int(( (_size_pdfs_1) % (((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)) == 0 ? (int64_t)(_size_pdfs_1) / (int64_t)(((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)) : ( (int64_t)(_size_pdfs_1) / (int64_t)(((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)) ) +1 )), int(( (_size_pdfs_2) % (((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)) == 0 ? (int64_t)(_size_pdfs_2) / (int64_t)(((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)) : ( (int64_t)(_size_pdfs_2) / (int64_t)(((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)) ) +1 )));
            internal_pack_W::pack_W<<<_grid, _block, 0, stream>>>(_data_buffer, _data_pdfs, _size_pdfs_0, _size_pdfs_1, _size_pdfs_2, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3);
            break;
        }
        
        case stencil::TW:
        {
            double * _data_buffer = buffer;
            WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()));
            WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()));
            WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()));
            double * const _data_pdfs = pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(), int64_t(cell_idx_c(ci.xSize()) + 0));
            const int64_t _size_pdfs_0 = int64_t(cell_idx_c(ci.xSize()) + 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(), int64_t(cell_idx_c(ci.ySize()) + 0));
            const int64_t _size_pdfs_1 = int64_t(cell_idx_c(ci.ySize()) + 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(), int64_t(cell_idx_c(ci.zSize()) + 0));
            const int64_t _size_pdfs_2 = int64_t(cell_idx_c(ci.zSize()) + 0);
            const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
            const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
            const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
            const int64_t _stride_pdfs_3 = int64_t(pdfs->fStride());
            dim3 _block(int(((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)), int(((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)), int(((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)));
            dim3 _grid(int(( (_size_pdfs_0) % (((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)) == 0 ? (int64_t)(_size_pdfs_0) / (int64_t)(((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)) : ( (int64_t)(_size_pdfs_0) / (int64_t)(((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)) ) +1 )), int(( (_size_pdfs_1) % (((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)) == 0 ? (int64_t)(_size_pdfs_1) / (int64_t)(((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)) : ( (int64_t)(_size_pdfs_1) / (int64_t)(((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)) ) +1 )), int(( (_size_pdfs_2) % (((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)) == 0 ? (int64_t)(_size_pdfs_2) / (int64_t)(((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)) : ( (int64_t)(_size_pdfs_2) / (int64_t)(((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)) ) +1 )));
            internal_pack_TW::pack_TW<<<_grid, _block, 0, stream>>>(_data_buffer, _data_pdfs, _size_pdfs_0, _size_pdfs_1, _size_pdfs_2, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3);
            break;
        }
        
        case stencil::NW:
        {
            double * _data_buffer = buffer;
            WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()));
            WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()));
            WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()));
            double * const _data_pdfs = pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(), int64_t(cell_idx_c(ci.xSize()) + 0));
            const int64_t _size_pdfs_0 = int64_t(cell_idx_c(ci.xSize()) + 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(), int64_t(cell_idx_c(ci.ySize()) + 0));
            const int64_t _size_pdfs_1 = int64_t(cell_idx_c(ci.ySize()) + 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(), int64_t(cell_idx_c(ci.zSize()) + 0));
            const int64_t _size_pdfs_2 = int64_t(cell_idx_c(ci.zSize()) + 0);
            const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
            const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
            const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
            const int64_t _stride_pdfs_3 = int64_t(pdfs->fStride());
            dim3 _block(int(((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)), int(((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)), int(((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)));
            dim3 _grid(int(( (_size_pdfs_0) % (((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)) == 0 ? (int64_t)(_size_pdfs_0) / (int64_t)(((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)) : ( (int64_t)(_size_pdfs_0) / (int64_t)(((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)) ) +1 )), int(( (_size_pdfs_1) % (((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)) == 0 ? (int64_t)(_size_pdfs_1) / (int64_t)(((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)) : ( (int64_t)(_size_pdfs_1) / (int64_t)(((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)) ) +1 )), int(( (_size_pdfs_2) % (((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)) == 0 ? (int64_t)(_size_pdfs_2) / (int64_t)(((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)) : ( (int64_t)(_size_pdfs_2) / (int64_t)(((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)) ) +1 )));
            internal_pack_NW::pack_NW<<<_grid, _block, 0, stream>>>(_data_buffer, _data_pdfs, _size_pdfs_0, _size_pdfs_1, _size_pdfs_2, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3);
            break;
        }
        
        case stencil::BS:
        {
            double * _data_buffer = buffer;
            WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()));
            WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()));
            WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()));
            double * const _data_pdfs = pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(), int64_t(cell_idx_c(ci.xSize()) + 0));
            const int64_t _size_pdfs_0 = int64_t(cell_idx_c(ci.xSize()) + 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(), int64_t(cell_idx_c(ci.ySize()) + 0));
            const int64_t _size_pdfs_1 = int64_t(cell_idx_c(ci.ySize()) + 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(), int64_t(cell_idx_c(ci.zSize()) + 0));
            const int64_t _size_pdfs_2 = int64_t(cell_idx_c(ci.zSize()) + 0);
            const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
            const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
            const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
            const int64_t _stride_pdfs_3 = int64_t(pdfs->fStride());
            dim3 _block(int(((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)), int(((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)), int(((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)));
            dim3 _grid(int(( (_size_pdfs_0) % (((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)) == 0 ? (int64_t)(_size_pdfs_0) / (int64_t)(((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)) : ( (int64_t)(_size_pdfs_0) / (int64_t)(((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)) ) +1 )), int(( (_size_pdfs_1) % (((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)) == 0 ? (int64_t)(_size_pdfs_1) / (int64_t)(((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)) : ( (int64_t)(_size_pdfs_1) / (int64_t)(((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)) ) +1 )), int(( (_size_pdfs_2) % (((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)) == 0 ? (int64_t)(_size_pdfs_2) / (int64_t)(((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)) : ( (int64_t)(_size_pdfs_2) / (int64_t)(((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)) ) +1 )));
            internal_pack_BS::pack_BS<<<_grid, _block, 0, stream>>>(_data_buffer, _data_pdfs, _size_pdfs_0, _size_pdfs_1, _size_pdfs_2, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3);
            break;
        }
        
        case stencil::S:
        {
            double * _data_buffer = buffer;
            WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()));
            WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()));
            WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()));
            double * const _data_pdfs = pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(), int64_t(cell_idx_c(ci.xSize()) + 0));
            const int64_t _size_pdfs_0 = int64_t(cell_idx_c(ci.xSize()) + 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(), int64_t(cell_idx_c(ci.ySize()) + 0));
            const int64_t _size_pdfs_1 = int64_t(cell_idx_c(ci.ySize()) + 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(), int64_t(cell_idx_c(ci.zSize()) + 0));
            const int64_t _size_pdfs_2 = int64_t(cell_idx_c(ci.zSize()) + 0);
            const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
            const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
            const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
            const int64_t _stride_pdfs_3 = int64_t(pdfs->fStride());
            dim3 _block(int(((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)), int(((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)), int(((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)));
            dim3 _grid(int(( (_size_pdfs_0) % (((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)) == 0 ? (int64_t)(_size_pdfs_0) / (int64_t)(((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)) : ( (int64_t)(_size_pdfs_0) / (int64_t)(((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)) ) +1 )), int(( (_size_pdfs_1) % (((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)) == 0 ? (int64_t)(_size_pdfs_1) / (int64_t)(((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)) : ( (int64_t)(_size_pdfs_1) / (int64_t)(((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)) ) +1 )), int(( (_size_pdfs_2) % (((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)) == 0 ? (int64_t)(_size_pdfs_2) / (int64_t)(((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)) : ( (int64_t)(_size_pdfs_2) / (int64_t)(((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)) ) +1 )));
            internal_pack_S::pack_S<<<_grid, _block, 0, stream>>>(_data_buffer, _data_pdfs, _size_pdfs_0, _size_pdfs_1, _size_pdfs_2, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3);
            break;
        }
        
        case stencil::TS:
        {
            double * _data_buffer = buffer;
            WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()));
            WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()));
            WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()));
            double * const _data_pdfs = pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(), int64_t(cell_idx_c(ci.xSize()) + 0));
            const int64_t _size_pdfs_0 = int64_t(cell_idx_c(ci.xSize()) + 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(), int64_t(cell_idx_c(ci.ySize()) + 0));
            const int64_t _size_pdfs_1 = int64_t(cell_idx_c(ci.ySize()) + 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(), int64_t(cell_idx_c(ci.zSize()) + 0));
            const int64_t _size_pdfs_2 = int64_t(cell_idx_c(ci.zSize()) + 0);
            const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
            const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
            const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
            const int64_t _stride_pdfs_3 = int64_t(pdfs->fStride());
            dim3 _block(int(((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)), int(((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)), int(((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)));
            dim3 _grid(int(( (_size_pdfs_0) % (((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)) == 0 ? (int64_t)(_size_pdfs_0) / (int64_t)(((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)) : ( (int64_t)(_size_pdfs_0) / (int64_t)(((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)) ) +1 )), int(( (_size_pdfs_1) % (((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)) == 0 ? (int64_t)(_size_pdfs_1) / (int64_t)(((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)) : ( (int64_t)(_size_pdfs_1) / (int64_t)(((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)) ) +1 )), int(( (_size_pdfs_2) % (((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)) == 0 ? (int64_t)(_size_pdfs_2) / (int64_t)(((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)) : ( (int64_t)(_size_pdfs_2) / (int64_t)(((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)) ) +1 )));
            internal_pack_TS::pack_TS<<<_grid, _block, 0, stream>>>(_data_buffer, _data_pdfs, _size_pdfs_0, _size_pdfs_1, _size_pdfs_2, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3);
            break;
        }
        
        case stencil::B:
        {
            double * _data_buffer = buffer;
            WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()));
            WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()));
            WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()));
            double * const _data_pdfs = pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(), int64_t(cell_idx_c(ci.xSize()) + 0));
            const int64_t _size_pdfs_0 = int64_t(cell_idx_c(ci.xSize()) + 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(), int64_t(cell_idx_c(ci.ySize()) + 0));
            const int64_t _size_pdfs_1 = int64_t(cell_idx_c(ci.ySize()) + 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(), int64_t(cell_idx_c(ci.zSize()) + 0));
            const int64_t _size_pdfs_2 = int64_t(cell_idx_c(ci.zSize()) + 0);
            const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
            const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
            const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
            const int64_t _stride_pdfs_3 = int64_t(pdfs->fStride());
            dim3 _block(int(((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)), int(((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)), int(((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)));
            dim3 _grid(int(( (_size_pdfs_0) % (((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)) == 0 ? (int64_t)(_size_pdfs_0) / (int64_t)(((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)) : ( (int64_t)(_size_pdfs_0) / (int64_t)(((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)) ) +1 )), int(( (_size_pdfs_1) % (((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)) == 0 ? (int64_t)(_size_pdfs_1) / (int64_t)(((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)) : ( (int64_t)(_size_pdfs_1) / (int64_t)(((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)) ) +1 )), int(( (_size_pdfs_2) % (((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)) == 0 ? (int64_t)(_size_pdfs_2) / (int64_t)(((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)) : ( (int64_t)(_size_pdfs_2) / (int64_t)(((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)) ) +1 )));
            internal_pack_B::pack_B<<<_grid, _block, 0, stream>>>(_data_buffer, _data_pdfs, _size_pdfs_0, _size_pdfs_1, _size_pdfs_2, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3);
            break;
        }
        
        case stencil::C:
        {
            double * _data_buffer = buffer;
            WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()));
            WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()));
            WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()));
            double * const _data_pdfs = pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(), int64_t(cell_idx_c(ci.xSize()) + 0));
            const int64_t _size_pdfs_0 = int64_t(cell_idx_c(ci.xSize()) + 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(), int64_t(cell_idx_c(ci.ySize()) + 0));
            const int64_t _size_pdfs_1 = int64_t(cell_idx_c(ci.ySize()) + 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(), int64_t(cell_idx_c(ci.zSize()) + 0));
            const int64_t _size_pdfs_2 = int64_t(cell_idx_c(ci.zSize()) + 0);
            const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
            const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
            const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
            dim3 _block(int(((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)), int(((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)), int(((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)));
            dim3 _grid(int(( (_size_pdfs_0) % (((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)) == 0 ? (int64_t)(_size_pdfs_0) / (int64_t)(((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)) : ( (int64_t)(_size_pdfs_0) / (int64_t)(((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)) ) +1 )), int(( (_size_pdfs_1) % (((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)) == 0 ? (int64_t)(_size_pdfs_1) / (int64_t)(((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)) : ( (int64_t)(_size_pdfs_1) / (int64_t)(((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)) ) +1 )), int(( (_size_pdfs_2) % (((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)) == 0 ? (int64_t)(_size_pdfs_2) / (int64_t)(((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)) : ( (int64_t)(_size_pdfs_2) / (int64_t)(((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)) ) +1 )));
            internal_pack_C::pack_C<<<_grid, _block, 0, stream>>>(_data_buffer, _data_pdfs, _size_pdfs_0, _size_pdfs_1, _size_pdfs_2, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2);
            break;
        }
        
        case stencil::T:
        {
            double * _data_buffer = buffer;
            WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()));
            WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()));
            WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()));
            double * const _data_pdfs = pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(), int64_t(cell_idx_c(ci.xSize()) + 0));
            const int64_t _size_pdfs_0 = int64_t(cell_idx_c(ci.xSize()) + 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(), int64_t(cell_idx_c(ci.ySize()) + 0));
            const int64_t _size_pdfs_1 = int64_t(cell_idx_c(ci.ySize()) + 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(), int64_t(cell_idx_c(ci.zSize()) + 0));
            const int64_t _size_pdfs_2 = int64_t(cell_idx_c(ci.zSize()) + 0);
            const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
            const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
            const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
            const int64_t _stride_pdfs_3 = int64_t(pdfs->fStride());
            dim3 _block(int(((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)), int(((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)), int(((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)));
            dim3 _grid(int(( (_size_pdfs_0) % (((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)) == 0 ? (int64_t)(_size_pdfs_0) / (int64_t)(((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)) : ( (int64_t)(_size_pdfs_0) / (int64_t)(((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)) ) +1 )), int(( (_size_pdfs_1) % (((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)) == 0 ? (int64_t)(_size_pdfs_1) / (int64_t)(((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)) : ( (int64_t)(_size_pdfs_1) / (int64_t)(((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)) ) +1 )), int(( (_size_pdfs_2) % (((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)) == 0 ? (int64_t)(_size_pdfs_2) / (int64_t)(((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)) : ( (int64_t)(_size_pdfs_2) / (int64_t)(((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)) ) +1 )));
            internal_pack_T::pack_T<<<_grid, _block, 0, stream>>>(_data_buffer, _data_pdfs, _size_pdfs_0, _size_pdfs_1, _size_pdfs_2, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3);
            break;
        }
        
        case stencil::BN:
        {
            double * _data_buffer = buffer;
            WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()));
            WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()));
            WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()));
            double * const _data_pdfs = pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(), int64_t(cell_idx_c(ci.xSize()) + 0));
            const int64_t _size_pdfs_0 = int64_t(cell_idx_c(ci.xSize()) + 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(), int64_t(cell_idx_c(ci.ySize()) + 0));
            const int64_t _size_pdfs_1 = int64_t(cell_idx_c(ci.ySize()) + 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(), int64_t(cell_idx_c(ci.zSize()) + 0));
            const int64_t _size_pdfs_2 = int64_t(cell_idx_c(ci.zSize()) + 0);
            const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
            const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
            const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
            const int64_t _stride_pdfs_3 = int64_t(pdfs->fStride());
            dim3 _block(int(((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)), int(((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)), int(((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)));
            dim3 _grid(int(( (_size_pdfs_0) % (((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)) == 0 ? (int64_t)(_size_pdfs_0) / (int64_t)(((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)) : ( (int64_t)(_size_pdfs_0) / (int64_t)(((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)) ) +1 )), int(( (_size_pdfs_1) % (((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)) == 0 ? (int64_t)(_size_pdfs_1) / (int64_t)(((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)) : ( (int64_t)(_size_pdfs_1) / (int64_t)(((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)) ) +1 )), int(( (_size_pdfs_2) % (((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)) == 0 ? (int64_t)(_size_pdfs_2) / (int64_t)(((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)) : ( (int64_t)(_size_pdfs_2) / (int64_t)(((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)) ) +1 )));
            internal_pack_BN::pack_BN<<<_grid, _block, 0, stream>>>(_data_buffer, _data_pdfs, _size_pdfs_0, _size_pdfs_1, _size_pdfs_2, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3);
            break;
        }
        
        case stencil::N:
        {
            double * _data_buffer = buffer;
            WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()));
            WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()));
            WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()));
            double * const _data_pdfs = pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(), int64_t(cell_idx_c(ci.xSize()) + 0));
            const int64_t _size_pdfs_0 = int64_t(cell_idx_c(ci.xSize()) + 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(), int64_t(cell_idx_c(ci.ySize()) + 0));
            const int64_t _size_pdfs_1 = int64_t(cell_idx_c(ci.ySize()) + 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(), int64_t(cell_idx_c(ci.zSize()) + 0));
            const int64_t _size_pdfs_2 = int64_t(cell_idx_c(ci.zSize()) + 0);
            const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
            const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
            const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
            const int64_t _stride_pdfs_3 = int64_t(pdfs->fStride());
            dim3 _block(int(((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)), int(((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)), int(((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)));
            dim3 _grid(int(( (_size_pdfs_0) % (((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)) == 0 ? (int64_t)(_size_pdfs_0) / (int64_t)(((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)) : ( (int64_t)(_size_pdfs_0) / (int64_t)(((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)) ) +1 )), int(( (_size_pdfs_1) % (((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)) == 0 ? (int64_t)(_size_pdfs_1) / (int64_t)(((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)) : ( (int64_t)(_size_pdfs_1) / (int64_t)(((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)) ) +1 )), int(( (_size_pdfs_2) % (((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)) == 0 ? (int64_t)(_size_pdfs_2) / (int64_t)(((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)) : ( (int64_t)(_size_pdfs_2) / (int64_t)(((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)) ) +1 )));
            internal_pack_N::pack_N<<<_grid, _block, 0, stream>>>(_data_buffer, _data_pdfs, _size_pdfs_0, _size_pdfs_1, _size_pdfs_2, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3);
            break;
        }
        
        case stencil::TN:
        {
            double * _data_buffer = buffer;
            WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()));
            WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()));
            WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()));
            double * const _data_pdfs = pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(), int64_t(cell_idx_c(ci.xSize()) + 0));
            const int64_t _size_pdfs_0 = int64_t(cell_idx_c(ci.xSize()) + 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(), int64_t(cell_idx_c(ci.ySize()) + 0));
            const int64_t _size_pdfs_1 = int64_t(cell_idx_c(ci.ySize()) + 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(), int64_t(cell_idx_c(ci.zSize()) + 0));
            const int64_t _size_pdfs_2 = int64_t(cell_idx_c(ci.zSize()) + 0);
            const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
            const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
            const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
            const int64_t _stride_pdfs_3 = int64_t(pdfs->fStride());
            dim3 _block(int(((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)), int(((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)), int(((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)));
            dim3 _grid(int(( (_size_pdfs_0) % (((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)) == 0 ? (int64_t)(_size_pdfs_0) / (int64_t)(((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)) : ( (int64_t)(_size_pdfs_0) / (int64_t)(((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)) ) +1 )), int(( (_size_pdfs_1) % (((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)) == 0 ? (int64_t)(_size_pdfs_1) / (int64_t)(((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)) : ( (int64_t)(_size_pdfs_1) / (int64_t)(((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)) ) +1 )), int(( (_size_pdfs_2) % (((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)) == 0 ? (int64_t)(_size_pdfs_2) / (int64_t)(((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)) : ( (int64_t)(_size_pdfs_2) / (int64_t)(((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)) ) +1 )));
            internal_pack_TN::pack_TN<<<_grid, _block, 0, stream>>>(_data_buffer, _data_pdfs, _size_pdfs_0, _size_pdfs_1, _size_pdfs_2, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3);
            break;
        }
        
        case stencil::SE:
        {
            double * _data_buffer = buffer;
            WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()));
            WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()));
            WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()));
            double * const _data_pdfs = pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(), int64_t(cell_idx_c(ci.xSize()) + 0));
            const int64_t _size_pdfs_0 = int64_t(cell_idx_c(ci.xSize()) + 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(), int64_t(cell_idx_c(ci.ySize()) + 0));
            const int64_t _size_pdfs_1 = int64_t(cell_idx_c(ci.ySize()) + 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(), int64_t(cell_idx_c(ci.zSize()) + 0));
            const int64_t _size_pdfs_2 = int64_t(cell_idx_c(ci.zSize()) + 0);
            const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
            const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
            const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
            const int64_t _stride_pdfs_3 = int64_t(pdfs->fStride());
            dim3 _block(int(((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)), int(((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)), int(((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)));
            dim3 _grid(int(( (_size_pdfs_0) % (((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)) == 0 ? (int64_t)(_size_pdfs_0) / (int64_t)(((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)) : ( (int64_t)(_size_pdfs_0) / (int64_t)(((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)) ) +1 )), int(( (_size_pdfs_1) % (((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)) == 0 ? (int64_t)(_size_pdfs_1) / (int64_t)(((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)) : ( (int64_t)(_size_pdfs_1) / (int64_t)(((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)) ) +1 )), int(( (_size_pdfs_2) % (((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)) == 0 ? (int64_t)(_size_pdfs_2) / (int64_t)(((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)) : ( (int64_t)(_size_pdfs_2) / (int64_t)(((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)) ) +1 )));
            internal_pack_SE::pack_SE<<<_grid, _block, 0, stream>>>(_data_buffer, _data_pdfs, _size_pdfs_0, _size_pdfs_1, _size_pdfs_2, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3);
            break;
        }
        
        case stencil::BE:
        {
            double * _data_buffer = buffer;
            WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()));
            WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()));
            WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()));
            double * const _data_pdfs = pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(), int64_t(cell_idx_c(ci.xSize()) + 0));
            const int64_t _size_pdfs_0 = int64_t(cell_idx_c(ci.xSize()) + 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(), int64_t(cell_idx_c(ci.ySize()) + 0));
            const int64_t _size_pdfs_1 = int64_t(cell_idx_c(ci.ySize()) + 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(), int64_t(cell_idx_c(ci.zSize()) + 0));
            const int64_t _size_pdfs_2 = int64_t(cell_idx_c(ci.zSize()) + 0);
            const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
            const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
            const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
            const int64_t _stride_pdfs_3 = int64_t(pdfs->fStride());
            dim3 _block(int(((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)), int(((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)), int(((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)));
            dim3 _grid(int(( (_size_pdfs_0) % (((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)) == 0 ? (int64_t)(_size_pdfs_0) / (int64_t)(((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)) : ( (int64_t)(_size_pdfs_0) / (int64_t)(((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)) ) +1 )), int(( (_size_pdfs_1) % (((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)) == 0 ? (int64_t)(_size_pdfs_1) / (int64_t)(((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)) : ( (int64_t)(_size_pdfs_1) / (int64_t)(((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)) ) +1 )), int(( (_size_pdfs_2) % (((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)) == 0 ? (int64_t)(_size_pdfs_2) / (int64_t)(((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)) : ( (int64_t)(_size_pdfs_2) / (int64_t)(((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)) ) +1 )));
            internal_pack_BE::pack_BE<<<_grid, _block, 0, stream>>>(_data_buffer, _data_pdfs, _size_pdfs_0, _size_pdfs_1, _size_pdfs_2, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3);
            break;
        }
        
        case stencil::E:
        {
            double * _data_buffer = buffer;
            WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()));
            WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()));
            WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()));
            double * const _data_pdfs = pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(), int64_t(cell_idx_c(ci.xSize()) + 0));
            const int64_t _size_pdfs_0 = int64_t(cell_idx_c(ci.xSize()) + 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(), int64_t(cell_idx_c(ci.ySize()) + 0));
            const int64_t _size_pdfs_1 = int64_t(cell_idx_c(ci.ySize()) + 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(), int64_t(cell_idx_c(ci.zSize()) + 0));
            const int64_t _size_pdfs_2 = int64_t(cell_idx_c(ci.zSize()) + 0);
            const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
            const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
            const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
            const int64_t _stride_pdfs_3 = int64_t(pdfs->fStride());
            dim3 _block(int(((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)), int(((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)), int(((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)));
            dim3 _grid(int(( (_size_pdfs_0) % (((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)) == 0 ? (int64_t)(_size_pdfs_0) / (int64_t)(((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)) : ( (int64_t)(_size_pdfs_0) / (int64_t)(((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)) ) +1 )), int(( (_size_pdfs_1) % (((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)) == 0 ? (int64_t)(_size_pdfs_1) / (int64_t)(((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)) : ( (int64_t)(_size_pdfs_1) / (int64_t)(((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)) ) +1 )), int(( (_size_pdfs_2) % (((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)) == 0 ? (int64_t)(_size_pdfs_2) / (int64_t)(((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)) : ( (int64_t)(_size_pdfs_2) / (int64_t)(((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)) ) +1 )));
            internal_pack_E::pack_E<<<_grid, _block, 0, stream>>>(_data_buffer, _data_pdfs, _size_pdfs_0, _size_pdfs_1, _size_pdfs_2, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3);
            break;
        }
        
        case stencil::TE:
        {
            double * _data_buffer = buffer;
            WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()));
            WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()));
            WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()));
            double * const _data_pdfs = pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(), int64_t(cell_idx_c(ci.xSize()) + 0));
            const int64_t _size_pdfs_0 = int64_t(cell_idx_c(ci.xSize()) + 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(), int64_t(cell_idx_c(ci.ySize()) + 0));
            const int64_t _size_pdfs_1 = int64_t(cell_idx_c(ci.ySize()) + 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(), int64_t(cell_idx_c(ci.zSize()) + 0));
            const int64_t _size_pdfs_2 = int64_t(cell_idx_c(ci.zSize()) + 0);
            const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
            const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
            const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
            const int64_t _stride_pdfs_3 = int64_t(pdfs->fStride());
            dim3 _block(int(((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)), int(((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)), int(((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)));
            dim3 _grid(int(( (_size_pdfs_0) % (((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)) == 0 ? (int64_t)(_size_pdfs_0) / (int64_t)(((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)) : ( (int64_t)(_size_pdfs_0) / (int64_t)(((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)) ) +1 )), int(( (_size_pdfs_1) % (((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)) == 0 ? (int64_t)(_size_pdfs_1) / (int64_t)(((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)) : ( (int64_t)(_size_pdfs_1) / (int64_t)(((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)) ) +1 )), int(( (_size_pdfs_2) % (((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)) == 0 ? (int64_t)(_size_pdfs_2) / (int64_t)(((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)) : ( (int64_t)(_size_pdfs_2) / (int64_t)(((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)) ) +1 )));
            internal_pack_TE::pack_TE<<<_grid, _block, 0, stream>>>(_data_buffer, _data_pdfs, _size_pdfs_0, _size_pdfs_1, _size_pdfs_2, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3);
            break;
        }
        
        case stencil::NE:
        {
            double * _data_buffer = buffer;
            WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()));
            WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()));
            WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()));
            double * const _data_pdfs = pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(), int64_t(cell_idx_c(ci.xSize()) + 0));
            const int64_t _size_pdfs_0 = int64_t(cell_idx_c(ci.xSize()) + 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(), int64_t(cell_idx_c(ci.ySize()) + 0));
            const int64_t _size_pdfs_1 = int64_t(cell_idx_c(ci.ySize()) + 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(), int64_t(cell_idx_c(ci.zSize()) + 0));
            const int64_t _size_pdfs_2 = int64_t(cell_idx_c(ci.zSize()) + 0);
            const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
            const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
            const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
            const int64_t _stride_pdfs_3 = int64_t(pdfs->fStride());
            dim3 _block(int(((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)), int(((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)), int(((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)));
            dim3 _grid(int(( (_size_pdfs_0) % (((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)) == 0 ? (int64_t)(_size_pdfs_0) / (int64_t)(((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)) : ( (int64_t)(_size_pdfs_0) / (int64_t)(((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)) ) +1 )), int(( (_size_pdfs_1) % (((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)) == 0 ? (int64_t)(_size_pdfs_1) / (int64_t)(((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)) : ( (int64_t)(_size_pdfs_1) / (int64_t)(((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)) ) +1 )), int(( (_size_pdfs_2) % (((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)) == 0 ? (int64_t)(_size_pdfs_2) / (int64_t)(((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)) : ( (int64_t)(_size_pdfs_2) / (int64_t)(((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)) ) +1 )));
            internal_pack_NE::pack_NE<<<_grid, _block, 0, stream>>>(_data_buffer, _data_pdfs, _size_pdfs_0, _size_pdfs_1, _size_pdfs_2, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3);
            break;
        }
        

        default:
            WALBERLA_ASSERT(false);
    }
}


void UniformGridGPU_PackInfo::unpack(Direction dir, unsigned char * byte_buffer, IBlock * block, cudaStream_t stream)
{
    double * buffer = reinterpret_cast<double*>(byte_buffer);

    auto pdfs = block->getData< cuda::GPUField<double> >(pdfsID);

    CellInterval ci;
    pdfs->getGhostRegion(dir, ci, 1, false);

    switch( dir )
    {
        case stencil::NE:
        {
            double * const _data_buffer = buffer;
            WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()));
            WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()));
            WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()));
            double * _data_pdfs = pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(), int64_t(cell_idx_c(ci.xSize()) + 0));
            const int64_t _size_pdfs_0 = int64_t(cell_idx_c(ci.xSize()) + 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(), int64_t(cell_idx_c(ci.ySize()) + 0));
            const int64_t _size_pdfs_1 = int64_t(cell_idx_c(ci.ySize()) + 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(), int64_t(cell_idx_c(ci.zSize()) + 0));
            const int64_t _size_pdfs_2 = int64_t(cell_idx_c(ci.zSize()) + 0);
            const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
            const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
            const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
            const int64_t _stride_pdfs_3 = int64_t(pdfs->fStride());
            dim3 _block(int(((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)), int(((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)), int(((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)));
            dim3 _grid(int(( (_size_pdfs_0) % (((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)) == 0 ? (int64_t)(_size_pdfs_0) / (int64_t)(((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)) : ( (int64_t)(_size_pdfs_0) / (int64_t)(((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)) ) +1 )), int(( (_size_pdfs_1) % (((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)) == 0 ? (int64_t)(_size_pdfs_1) / (int64_t)(((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)) : ( (int64_t)(_size_pdfs_1) / (int64_t)(((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)) ) +1 )), int(( (_size_pdfs_2) % (((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)) == 0 ? (int64_t)(_size_pdfs_2) / (int64_t)(((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)) : ( (int64_t)(_size_pdfs_2) / (int64_t)(((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)) ) +1 )));
            internal_unpack_NE::unpack_NE<<<_grid, _block, 0, stream>>>(_data_buffer, _data_pdfs, _size_pdfs_0, _size_pdfs_1, _size_pdfs_2, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3);
            break;
        }
        
        case stencil::TE:
        {
            double * const _data_buffer = buffer;
            WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()));
            WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()));
            WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()));
            double * _data_pdfs = pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(), int64_t(cell_idx_c(ci.xSize()) + 0));
            const int64_t _size_pdfs_0 = int64_t(cell_idx_c(ci.xSize()) + 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(), int64_t(cell_idx_c(ci.ySize()) + 0));
            const int64_t _size_pdfs_1 = int64_t(cell_idx_c(ci.ySize()) + 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(), int64_t(cell_idx_c(ci.zSize()) + 0));
            const int64_t _size_pdfs_2 = int64_t(cell_idx_c(ci.zSize()) + 0);
            const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
            const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
            const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
            const int64_t _stride_pdfs_3 = int64_t(pdfs->fStride());
            dim3 _block(int(((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)), int(((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)), int(((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)));
            dim3 _grid(int(( (_size_pdfs_0) % (((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)) == 0 ? (int64_t)(_size_pdfs_0) / (int64_t)(((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)) : ( (int64_t)(_size_pdfs_0) / (int64_t)(((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)) ) +1 )), int(( (_size_pdfs_1) % (((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)) == 0 ? (int64_t)(_size_pdfs_1) / (int64_t)(((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)) : ( (int64_t)(_size_pdfs_1) / (int64_t)(((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)) ) +1 )), int(( (_size_pdfs_2) % (((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)) == 0 ? (int64_t)(_size_pdfs_2) / (int64_t)(((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)) : ( (int64_t)(_size_pdfs_2) / (int64_t)(((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)) ) +1 )));
            internal_unpack_TE::unpack_TE<<<_grid, _block, 0, stream>>>(_data_buffer, _data_pdfs, _size_pdfs_0, _size_pdfs_1, _size_pdfs_2, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3);
            break;
        }
        
        case stencil::E:
        {
            double * const _data_buffer = buffer;
            WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()));
            WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()));
            WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()));
            double * _data_pdfs = pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(), int64_t(cell_idx_c(ci.xSize()) + 0));
            const int64_t _size_pdfs_0 = int64_t(cell_idx_c(ci.xSize()) + 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(), int64_t(cell_idx_c(ci.ySize()) + 0));
            const int64_t _size_pdfs_1 = int64_t(cell_idx_c(ci.ySize()) + 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(), int64_t(cell_idx_c(ci.zSize()) + 0));
            const int64_t _size_pdfs_2 = int64_t(cell_idx_c(ci.zSize()) + 0);
            const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
            const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
            const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
            const int64_t _stride_pdfs_3 = int64_t(pdfs->fStride());
            dim3 _block(int(((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)), int(((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)), int(((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)));
            dim3 _grid(int(( (_size_pdfs_0) % (((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)) == 0 ? (int64_t)(_size_pdfs_0) / (int64_t)(((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)) : ( (int64_t)(_size_pdfs_0) / (int64_t)(((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)) ) +1 )), int(( (_size_pdfs_1) % (((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)) == 0 ? (int64_t)(_size_pdfs_1) / (int64_t)(((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)) : ( (int64_t)(_size_pdfs_1) / (int64_t)(((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)) ) +1 )), int(( (_size_pdfs_2) % (((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)) == 0 ? (int64_t)(_size_pdfs_2) / (int64_t)(((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)) : ( (int64_t)(_size_pdfs_2) / (int64_t)(((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)) ) +1 )));
            internal_unpack_E::unpack_E<<<_grid, _block, 0, stream>>>(_data_buffer, _data_pdfs, _size_pdfs_0, _size_pdfs_1, _size_pdfs_2, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3);
            break;
        }
        
        case stencil::BE:
        {
            double * const _data_buffer = buffer;
            WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()));
            WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()));
            WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()));
            double * _data_pdfs = pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(), int64_t(cell_idx_c(ci.xSize()) + 0));
            const int64_t _size_pdfs_0 = int64_t(cell_idx_c(ci.xSize()) + 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(), int64_t(cell_idx_c(ci.ySize()) + 0));
            const int64_t _size_pdfs_1 = int64_t(cell_idx_c(ci.ySize()) + 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(), int64_t(cell_idx_c(ci.zSize()) + 0));
            const int64_t _size_pdfs_2 = int64_t(cell_idx_c(ci.zSize()) + 0);
            const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
            const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
            const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
            const int64_t _stride_pdfs_3 = int64_t(pdfs->fStride());
            dim3 _block(int(((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)), int(((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)), int(((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)));
            dim3 _grid(int(( (_size_pdfs_0) % (((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)) == 0 ? (int64_t)(_size_pdfs_0) / (int64_t)(((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)) : ( (int64_t)(_size_pdfs_0) / (int64_t)(((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)) ) +1 )), int(( (_size_pdfs_1) % (((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)) == 0 ? (int64_t)(_size_pdfs_1) / (int64_t)(((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)) : ( (int64_t)(_size_pdfs_1) / (int64_t)(((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)) ) +1 )), int(( (_size_pdfs_2) % (((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)) == 0 ? (int64_t)(_size_pdfs_2) / (int64_t)(((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)) : ( (int64_t)(_size_pdfs_2) / (int64_t)(((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)) ) +1 )));
            internal_unpack_BE::unpack_BE<<<_grid, _block, 0, stream>>>(_data_buffer, _data_pdfs, _size_pdfs_0, _size_pdfs_1, _size_pdfs_2, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3);
            break;
        }
        
        case stencil::SE:
        {
            double * const _data_buffer = buffer;
            WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()));
            WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()));
            WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()));
            double * _data_pdfs = pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(), int64_t(cell_idx_c(ci.xSize()) + 0));
            const int64_t _size_pdfs_0 = int64_t(cell_idx_c(ci.xSize()) + 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(), int64_t(cell_idx_c(ci.ySize()) + 0));
            const int64_t _size_pdfs_1 = int64_t(cell_idx_c(ci.ySize()) + 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(), int64_t(cell_idx_c(ci.zSize()) + 0));
            const int64_t _size_pdfs_2 = int64_t(cell_idx_c(ci.zSize()) + 0);
            const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
            const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
            const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
            const int64_t _stride_pdfs_3 = int64_t(pdfs->fStride());
            dim3 _block(int(((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)), int(((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)), int(((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)));
            dim3 _grid(int(( (_size_pdfs_0) % (((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)) == 0 ? (int64_t)(_size_pdfs_0) / (int64_t)(((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)) : ( (int64_t)(_size_pdfs_0) / (int64_t)(((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)) ) +1 )), int(( (_size_pdfs_1) % (((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)) == 0 ? (int64_t)(_size_pdfs_1) / (int64_t)(((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)) : ( (int64_t)(_size_pdfs_1) / (int64_t)(((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)) ) +1 )), int(( (_size_pdfs_2) % (((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)) == 0 ? (int64_t)(_size_pdfs_2) / (int64_t)(((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)) : ( (int64_t)(_size_pdfs_2) / (int64_t)(((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)) ) +1 )));
            internal_unpack_SE::unpack_SE<<<_grid, _block, 0, stream>>>(_data_buffer, _data_pdfs, _size_pdfs_0, _size_pdfs_1, _size_pdfs_2, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3);
            break;
        }
        
        case stencil::TN:
        {
            double * const _data_buffer = buffer;
            WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()));
            WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()));
            WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()));
            double * _data_pdfs = pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(), int64_t(cell_idx_c(ci.xSize()) + 0));
            const int64_t _size_pdfs_0 = int64_t(cell_idx_c(ci.xSize()) + 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(), int64_t(cell_idx_c(ci.ySize()) + 0));
            const int64_t _size_pdfs_1 = int64_t(cell_idx_c(ci.ySize()) + 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(), int64_t(cell_idx_c(ci.zSize()) + 0));
            const int64_t _size_pdfs_2 = int64_t(cell_idx_c(ci.zSize()) + 0);
            const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
            const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
            const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
            const int64_t _stride_pdfs_3 = int64_t(pdfs->fStride());
            dim3 _block(int(((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)), int(((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)), int(((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)));
            dim3 _grid(int(( (_size_pdfs_0) % (((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)) == 0 ? (int64_t)(_size_pdfs_0) / (int64_t)(((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)) : ( (int64_t)(_size_pdfs_0) / (int64_t)(((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)) ) +1 )), int(( (_size_pdfs_1) % (((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)) == 0 ? (int64_t)(_size_pdfs_1) / (int64_t)(((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)) : ( (int64_t)(_size_pdfs_1) / (int64_t)(((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)) ) +1 )), int(( (_size_pdfs_2) % (((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)) == 0 ? (int64_t)(_size_pdfs_2) / (int64_t)(((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)) : ( (int64_t)(_size_pdfs_2) / (int64_t)(((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)) ) +1 )));
            internal_unpack_TN::unpack_TN<<<_grid, _block, 0, stream>>>(_data_buffer, _data_pdfs, _size_pdfs_0, _size_pdfs_1, _size_pdfs_2, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3);
            break;
        }
        
        case stencil::N:
        {
            double * const _data_buffer = buffer;
            WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()));
            WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()));
            WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()));
            double * _data_pdfs = pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(), int64_t(cell_idx_c(ci.xSize()) + 0));
            const int64_t _size_pdfs_0 = int64_t(cell_idx_c(ci.xSize()) + 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(), int64_t(cell_idx_c(ci.ySize()) + 0));
            const int64_t _size_pdfs_1 = int64_t(cell_idx_c(ci.ySize()) + 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(), int64_t(cell_idx_c(ci.zSize()) + 0));
            const int64_t _size_pdfs_2 = int64_t(cell_idx_c(ci.zSize()) + 0);
            const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
            const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
            const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
            const int64_t _stride_pdfs_3 = int64_t(pdfs->fStride());
            dim3 _block(int(((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)), int(((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)), int(((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)));
            dim3 _grid(int(( (_size_pdfs_0) % (((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)) == 0 ? (int64_t)(_size_pdfs_0) / (int64_t)(((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)) : ( (int64_t)(_size_pdfs_0) / (int64_t)(((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)) ) +1 )), int(( (_size_pdfs_1) % (((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)) == 0 ? (int64_t)(_size_pdfs_1) / (int64_t)(((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)) : ( (int64_t)(_size_pdfs_1) / (int64_t)(((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)) ) +1 )), int(( (_size_pdfs_2) % (((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)) == 0 ? (int64_t)(_size_pdfs_2) / (int64_t)(((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)) : ( (int64_t)(_size_pdfs_2) / (int64_t)(((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)) ) +1 )));
            internal_unpack_N::unpack_N<<<_grid, _block, 0, stream>>>(_data_buffer, _data_pdfs, _size_pdfs_0, _size_pdfs_1, _size_pdfs_2, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3);
            break;
        }
        
        case stencil::BN:
        {
            double * const _data_buffer = buffer;
            WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()));
            WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()));
            WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()));
            double * _data_pdfs = pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(), int64_t(cell_idx_c(ci.xSize()) + 0));
            const int64_t _size_pdfs_0 = int64_t(cell_idx_c(ci.xSize()) + 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(), int64_t(cell_idx_c(ci.ySize()) + 0));
            const int64_t _size_pdfs_1 = int64_t(cell_idx_c(ci.ySize()) + 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(), int64_t(cell_idx_c(ci.zSize()) + 0));
            const int64_t _size_pdfs_2 = int64_t(cell_idx_c(ci.zSize()) + 0);
            const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
            const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
            const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
            const int64_t _stride_pdfs_3 = int64_t(pdfs->fStride());
            dim3 _block(int(((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)), int(((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)), int(((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)));
            dim3 _grid(int(( (_size_pdfs_0) % (((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)) == 0 ? (int64_t)(_size_pdfs_0) / (int64_t)(((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)) : ( (int64_t)(_size_pdfs_0) / (int64_t)(((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)) ) +1 )), int(( (_size_pdfs_1) % (((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)) == 0 ? (int64_t)(_size_pdfs_1) / (int64_t)(((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)) : ( (int64_t)(_size_pdfs_1) / (int64_t)(((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)) ) +1 )), int(( (_size_pdfs_2) % (((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)) == 0 ? (int64_t)(_size_pdfs_2) / (int64_t)(((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)) : ( (int64_t)(_size_pdfs_2) / (int64_t)(((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)) ) +1 )));
            internal_unpack_BN::unpack_BN<<<_grid, _block, 0, stream>>>(_data_buffer, _data_pdfs, _size_pdfs_0, _size_pdfs_1, _size_pdfs_2, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3);
            break;
        }
        
        case stencil::T:
        {
            double * const _data_buffer = buffer;
            WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()));
            WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()));
            WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()));
            double * _data_pdfs = pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(), int64_t(cell_idx_c(ci.xSize()) + 0));
            const int64_t _size_pdfs_0 = int64_t(cell_idx_c(ci.xSize()) + 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(), int64_t(cell_idx_c(ci.ySize()) + 0));
            const int64_t _size_pdfs_1 = int64_t(cell_idx_c(ci.ySize()) + 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(), int64_t(cell_idx_c(ci.zSize()) + 0));
            const int64_t _size_pdfs_2 = int64_t(cell_idx_c(ci.zSize()) + 0);
            const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
            const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
            const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
            const int64_t _stride_pdfs_3 = int64_t(pdfs->fStride());
            dim3 _block(int(((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)), int(((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)), int(((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)));
            dim3 _grid(int(( (_size_pdfs_0) % (((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)) == 0 ? (int64_t)(_size_pdfs_0) / (int64_t)(((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)) : ( (int64_t)(_size_pdfs_0) / (int64_t)(((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)) ) +1 )), int(( (_size_pdfs_1) % (((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)) == 0 ? (int64_t)(_size_pdfs_1) / (int64_t)(((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)) : ( (int64_t)(_size_pdfs_1) / (int64_t)(((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)) ) +1 )), int(( (_size_pdfs_2) % (((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)) == 0 ? (int64_t)(_size_pdfs_2) / (int64_t)(((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)) : ( (int64_t)(_size_pdfs_2) / (int64_t)(((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)) ) +1 )));
            internal_unpack_T::unpack_T<<<_grid, _block, 0, stream>>>(_data_buffer, _data_pdfs, _size_pdfs_0, _size_pdfs_1, _size_pdfs_2, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3);
            break;
        }
        
        case stencil::C:
        {
            double * const _data_buffer = buffer;
            WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()));
            WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()));
            WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()));
            double * _data_pdfs = pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(), int64_t(cell_idx_c(ci.xSize()) + 0));
            const int64_t _size_pdfs_0 = int64_t(cell_idx_c(ci.xSize()) + 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(), int64_t(cell_idx_c(ci.ySize()) + 0));
            const int64_t _size_pdfs_1 = int64_t(cell_idx_c(ci.ySize()) + 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(), int64_t(cell_idx_c(ci.zSize()) + 0));
            const int64_t _size_pdfs_2 = int64_t(cell_idx_c(ci.zSize()) + 0);
            const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
            const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
            const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
            dim3 _block(int(((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)), int(((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)), int(((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)));
            dim3 _grid(int(( (_size_pdfs_0) % (((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)) == 0 ? (int64_t)(_size_pdfs_0) / (int64_t)(((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)) : ( (int64_t)(_size_pdfs_0) / (int64_t)(((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)) ) +1 )), int(( (_size_pdfs_1) % (((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)) == 0 ? (int64_t)(_size_pdfs_1) / (int64_t)(((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)) : ( (int64_t)(_size_pdfs_1) / (int64_t)(((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)) ) +1 )), int(( (_size_pdfs_2) % (((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)) == 0 ? (int64_t)(_size_pdfs_2) / (int64_t)(((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)) : ( (int64_t)(_size_pdfs_2) / (int64_t)(((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)) ) +1 )));
            internal_unpack_C::unpack_C<<<_grid, _block, 0, stream>>>(_data_buffer, _data_pdfs, _size_pdfs_0, _size_pdfs_1, _size_pdfs_2, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2);
            break;
        }
        
        case stencil::B:
        {
            double * const _data_buffer = buffer;
            WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()));
            WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()));
            WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()));
            double * _data_pdfs = pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(), int64_t(cell_idx_c(ci.xSize()) + 0));
            const int64_t _size_pdfs_0 = int64_t(cell_idx_c(ci.xSize()) + 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(), int64_t(cell_idx_c(ci.ySize()) + 0));
            const int64_t _size_pdfs_1 = int64_t(cell_idx_c(ci.ySize()) + 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(), int64_t(cell_idx_c(ci.zSize()) + 0));
            const int64_t _size_pdfs_2 = int64_t(cell_idx_c(ci.zSize()) + 0);
            const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
            const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
            const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
            const int64_t _stride_pdfs_3 = int64_t(pdfs->fStride());
            dim3 _block(int(((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)), int(((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)), int(((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)));
            dim3 _grid(int(( (_size_pdfs_0) % (((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)) == 0 ? (int64_t)(_size_pdfs_0) / (int64_t)(((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)) : ( (int64_t)(_size_pdfs_0) / (int64_t)(((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)) ) +1 )), int(( (_size_pdfs_1) % (((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)) == 0 ? (int64_t)(_size_pdfs_1) / (int64_t)(((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)) : ( (int64_t)(_size_pdfs_1) / (int64_t)(((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)) ) +1 )), int(( (_size_pdfs_2) % (((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)) == 0 ? (int64_t)(_size_pdfs_2) / (int64_t)(((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)) : ( (int64_t)(_size_pdfs_2) / (int64_t)(((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)) ) +1 )));
            internal_unpack_B::unpack_B<<<_grid, _block, 0, stream>>>(_data_buffer, _data_pdfs, _size_pdfs_0, _size_pdfs_1, _size_pdfs_2, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3);
            break;
        }
        
        case stencil::TS:
        {
            double * const _data_buffer = buffer;
            WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()));
            WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()));
            WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()));
            double * _data_pdfs = pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(), int64_t(cell_idx_c(ci.xSize()) + 0));
            const int64_t _size_pdfs_0 = int64_t(cell_idx_c(ci.xSize()) + 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(), int64_t(cell_idx_c(ci.ySize()) + 0));
            const int64_t _size_pdfs_1 = int64_t(cell_idx_c(ci.ySize()) + 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(), int64_t(cell_idx_c(ci.zSize()) + 0));
            const int64_t _size_pdfs_2 = int64_t(cell_idx_c(ci.zSize()) + 0);
            const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
            const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
            const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
            const int64_t _stride_pdfs_3 = int64_t(pdfs->fStride());
            dim3 _block(int(((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)), int(((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)), int(((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)));
            dim3 _grid(int(( (_size_pdfs_0) % (((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)) == 0 ? (int64_t)(_size_pdfs_0) / (int64_t)(((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)) : ( (int64_t)(_size_pdfs_0) / (int64_t)(((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)) ) +1 )), int(( (_size_pdfs_1) % (((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)) == 0 ? (int64_t)(_size_pdfs_1) / (int64_t)(((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)) : ( (int64_t)(_size_pdfs_1) / (int64_t)(((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)) ) +1 )), int(( (_size_pdfs_2) % (((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)) == 0 ? (int64_t)(_size_pdfs_2) / (int64_t)(((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)) : ( (int64_t)(_size_pdfs_2) / (int64_t)(((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)) ) +1 )));
            internal_unpack_TS::unpack_TS<<<_grid, _block, 0, stream>>>(_data_buffer, _data_pdfs, _size_pdfs_0, _size_pdfs_1, _size_pdfs_2, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3);
            break;
        }
        
        case stencil::S:
        {
            double * const _data_buffer = buffer;
            WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()));
            WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()));
            WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()));
            double * _data_pdfs = pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(), int64_t(cell_idx_c(ci.xSize()) + 0));
            const int64_t _size_pdfs_0 = int64_t(cell_idx_c(ci.xSize()) + 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(), int64_t(cell_idx_c(ci.ySize()) + 0));
            const int64_t _size_pdfs_1 = int64_t(cell_idx_c(ci.ySize()) + 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(), int64_t(cell_idx_c(ci.zSize()) + 0));
            const int64_t _size_pdfs_2 = int64_t(cell_idx_c(ci.zSize()) + 0);
            const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
            const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
            const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
            const int64_t _stride_pdfs_3 = int64_t(pdfs->fStride());
            dim3 _block(int(((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)), int(((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)), int(((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)));
            dim3 _grid(int(( (_size_pdfs_0) % (((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)) == 0 ? (int64_t)(_size_pdfs_0) / (int64_t)(((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)) : ( (int64_t)(_size_pdfs_0) / (int64_t)(((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)) ) +1 )), int(( (_size_pdfs_1) % (((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)) == 0 ? (int64_t)(_size_pdfs_1) / (int64_t)(((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)) : ( (int64_t)(_size_pdfs_1) / (int64_t)(((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)) ) +1 )), int(( (_size_pdfs_2) % (((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)) == 0 ? (int64_t)(_size_pdfs_2) / (int64_t)(((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)) : ( (int64_t)(_size_pdfs_2) / (int64_t)(((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)) ) +1 )));
            internal_unpack_S::unpack_S<<<_grid, _block, 0, stream>>>(_data_buffer, _data_pdfs, _size_pdfs_0, _size_pdfs_1, _size_pdfs_2, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3);
            break;
        }
        
        case stencil::BS:
        {
            double * const _data_buffer = buffer;
            WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()));
            WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()));
            WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()));
            double * _data_pdfs = pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(), int64_t(cell_idx_c(ci.xSize()) + 0));
            const int64_t _size_pdfs_0 = int64_t(cell_idx_c(ci.xSize()) + 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(), int64_t(cell_idx_c(ci.ySize()) + 0));
            const int64_t _size_pdfs_1 = int64_t(cell_idx_c(ci.ySize()) + 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(), int64_t(cell_idx_c(ci.zSize()) + 0));
            const int64_t _size_pdfs_2 = int64_t(cell_idx_c(ci.zSize()) + 0);
            const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
            const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
            const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
            const int64_t _stride_pdfs_3 = int64_t(pdfs->fStride());
            dim3 _block(int(((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)), int(((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)), int(((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)));
            dim3 _grid(int(( (_size_pdfs_0) % (((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)) == 0 ? (int64_t)(_size_pdfs_0) / (int64_t)(((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)) : ( (int64_t)(_size_pdfs_0) / (int64_t)(((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)) ) +1 )), int(( (_size_pdfs_1) % (((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)) == 0 ? (int64_t)(_size_pdfs_1) / (int64_t)(((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)) : ( (int64_t)(_size_pdfs_1) / (int64_t)(((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)) ) +1 )), int(( (_size_pdfs_2) % (((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)) == 0 ? (int64_t)(_size_pdfs_2) / (int64_t)(((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)) : ( (int64_t)(_size_pdfs_2) / (int64_t)(((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)) ) +1 )));
            internal_unpack_BS::unpack_BS<<<_grid, _block, 0, stream>>>(_data_buffer, _data_pdfs, _size_pdfs_0, _size_pdfs_1, _size_pdfs_2, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3);
            break;
        }
        
        case stencil::NW:
        {
            double * const _data_buffer = buffer;
            WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()));
            WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()));
            WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()));
            double * _data_pdfs = pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(), int64_t(cell_idx_c(ci.xSize()) + 0));
            const int64_t _size_pdfs_0 = int64_t(cell_idx_c(ci.xSize()) + 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(), int64_t(cell_idx_c(ci.ySize()) + 0));
            const int64_t _size_pdfs_1 = int64_t(cell_idx_c(ci.ySize()) + 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(), int64_t(cell_idx_c(ci.zSize()) + 0));
            const int64_t _size_pdfs_2 = int64_t(cell_idx_c(ci.zSize()) + 0);
            const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
            const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
            const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
            const int64_t _stride_pdfs_3 = int64_t(pdfs->fStride());
            dim3 _block(int(((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)), int(((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)), int(((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)));
            dim3 _grid(int(( (_size_pdfs_0) % (((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)) == 0 ? (int64_t)(_size_pdfs_0) / (int64_t)(((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)) : ( (int64_t)(_size_pdfs_0) / (int64_t)(((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)) ) +1 )), int(( (_size_pdfs_1) % (((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)) == 0 ? (int64_t)(_size_pdfs_1) / (int64_t)(((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)) : ( (int64_t)(_size_pdfs_1) / (int64_t)(((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)) ) +1 )), int(( (_size_pdfs_2) % (((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)) == 0 ? (int64_t)(_size_pdfs_2) / (int64_t)(((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)) : ( (int64_t)(_size_pdfs_2) / (int64_t)(((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)) ) +1 )));
            internal_unpack_NW::unpack_NW<<<_grid, _block, 0, stream>>>(_data_buffer, _data_pdfs, _size_pdfs_0, _size_pdfs_1, _size_pdfs_2, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3);
            break;
        }
        
        case stencil::TW:
        {
            double * const _data_buffer = buffer;
            WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()));
            WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()));
            WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()));
            double * _data_pdfs = pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(), int64_t(cell_idx_c(ci.xSize()) + 0));
            const int64_t _size_pdfs_0 = int64_t(cell_idx_c(ci.xSize()) + 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(), int64_t(cell_idx_c(ci.ySize()) + 0));
            const int64_t _size_pdfs_1 = int64_t(cell_idx_c(ci.ySize()) + 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(), int64_t(cell_idx_c(ci.zSize()) + 0));
            const int64_t _size_pdfs_2 = int64_t(cell_idx_c(ci.zSize()) + 0);
            const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
            const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
            const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
            const int64_t _stride_pdfs_3 = int64_t(pdfs->fStride());
            dim3 _block(int(((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)), int(((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)), int(((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)));
            dim3 _grid(int(( (_size_pdfs_0) % (((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)) == 0 ? (int64_t)(_size_pdfs_0) / (int64_t)(((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)) : ( (int64_t)(_size_pdfs_0) / (int64_t)(((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)) ) +1 )), int(( (_size_pdfs_1) % (((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)) == 0 ? (int64_t)(_size_pdfs_1) / (int64_t)(((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)) : ( (int64_t)(_size_pdfs_1) / (int64_t)(((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)) ) +1 )), int(( (_size_pdfs_2) % (((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)) == 0 ? (int64_t)(_size_pdfs_2) / (int64_t)(((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)) : ( (int64_t)(_size_pdfs_2) / (int64_t)(((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)) ) +1 )));
            internal_unpack_TW::unpack_TW<<<_grid, _block, 0, stream>>>(_data_buffer, _data_pdfs, _size_pdfs_0, _size_pdfs_1, _size_pdfs_2, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3);
            break;
        }
        
        case stencil::W:
        {
            double * const _data_buffer = buffer;
            WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()));
            WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()));
            WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()));
            double * _data_pdfs = pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(), int64_t(cell_idx_c(ci.xSize()) + 0));
            const int64_t _size_pdfs_0 = int64_t(cell_idx_c(ci.xSize()) + 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(), int64_t(cell_idx_c(ci.ySize()) + 0));
            const int64_t _size_pdfs_1 = int64_t(cell_idx_c(ci.ySize()) + 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(), int64_t(cell_idx_c(ci.zSize()) + 0));
            const int64_t _size_pdfs_2 = int64_t(cell_idx_c(ci.zSize()) + 0);
            const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
            const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
            const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
            const int64_t _stride_pdfs_3 = int64_t(pdfs->fStride());
            dim3 _block(int(((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)), int(((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)), int(((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)));
            dim3 _grid(int(( (_size_pdfs_0) % (((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)) == 0 ? (int64_t)(_size_pdfs_0) / (int64_t)(((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)) : ( (int64_t)(_size_pdfs_0) / (int64_t)(((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)) ) +1 )), int(( (_size_pdfs_1) % (((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)) == 0 ? (int64_t)(_size_pdfs_1) / (int64_t)(((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)) : ( (int64_t)(_size_pdfs_1) / (int64_t)(((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)) ) +1 )), int(( (_size_pdfs_2) % (((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)) == 0 ? (int64_t)(_size_pdfs_2) / (int64_t)(((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)) : ( (int64_t)(_size_pdfs_2) / (int64_t)(((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)) ) +1 )));
            internal_unpack_W::unpack_W<<<_grid, _block, 0, stream>>>(_data_buffer, _data_pdfs, _size_pdfs_0, _size_pdfs_1, _size_pdfs_2, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3);
            break;
        }
        
        case stencil::BW:
        {
            double * const _data_buffer = buffer;
            WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()));
            WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()));
            WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()));
            double * _data_pdfs = pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(), int64_t(cell_idx_c(ci.xSize()) + 0));
            const int64_t _size_pdfs_0 = int64_t(cell_idx_c(ci.xSize()) + 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(), int64_t(cell_idx_c(ci.ySize()) + 0));
            const int64_t _size_pdfs_1 = int64_t(cell_idx_c(ci.ySize()) + 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(), int64_t(cell_idx_c(ci.zSize()) + 0));
            const int64_t _size_pdfs_2 = int64_t(cell_idx_c(ci.zSize()) + 0);
            const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
            const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
            const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
            const int64_t _stride_pdfs_3 = int64_t(pdfs->fStride());
            dim3 _block(int(((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)), int(((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)), int(((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)));
            dim3 _grid(int(( (_size_pdfs_0) % (((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)) == 0 ? (int64_t)(_size_pdfs_0) / (int64_t)(((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)) : ( (int64_t)(_size_pdfs_0) / (int64_t)(((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)) ) +1 )), int(( (_size_pdfs_1) % (((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)) == 0 ? (int64_t)(_size_pdfs_1) / (int64_t)(((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)) : ( (int64_t)(_size_pdfs_1) / (int64_t)(((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)) ) +1 )), int(( (_size_pdfs_2) % (((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)) == 0 ? (int64_t)(_size_pdfs_2) / (int64_t)(((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)) : ( (int64_t)(_size_pdfs_2) / (int64_t)(((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)) ) +1 )));
            internal_unpack_BW::unpack_BW<<<_grid, _block, 0, stream>>>(_data_buffer, _data_pdfs, _size_pdfs_0, _size_pdfs_1, _size_pdfs_2, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3);
            break;
        }
        
        case stencil::SW:
        {
            double * const _data_buffer = buffer;
            WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()));
            WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()));
            WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()));
            double * _data_pdfs = pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(), int64_t(cell_idx_c(ci.xSize()) + 0));
            const int64_t _size_pdfs_0 = int64_t(cell_idx_c(ci.xSize()) + 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(), int64_t(cell_idx_c(ci.ySize()) + 0));
            const int64_t _size_pdfs_1 = int64_t(cell_idx_c(ci.ySize()) + 0);
            WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(), int64_t(cell_idx_c(ci.zSize()) + 0));
            const int64_t _size_pdfs_2 = int64_t(cell_idx_c(ci.zSize()) + 0);
            const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
            const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
            const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
            const int64_t _stride_pdfs_3 = int64_t(pdfs->fStride());
            dim3 _block(int(((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)), int(((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)), int(((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)));
            dim3 _grid(int(( (_size_pdfs_0) % (((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)) == 0 ? (int64_t)(_size_pdfs_0) / (int64_t)(((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)) : ( (int64_t)(_size_pdfs_0) / (int64_t)(((16 < _size_pdfs_0) ? 16 : _size_pdfs_0)) ) +1 )), int(( (_size_pdfs_1) % (((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)) == 0 ? (int64_t)(_size_pdfs_1) / (int64_t)(((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)) : ( (int64_t)(_size_pdfs_1) / (int64_t)(((16 < _size_pdfs_1) ? 16 : _size_pdfs_1)) ) +1 )), int(( (_size_pdfs_2) % (((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)) == 0 ? (int64_t)(_size_pdfs_2) / (int64_t)(((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)) : ( (int64_t)(_size_pdfs_2) / (int64_t)(((1 < _size_pdfs_2) ? 1 : _size_pdfs_2)) ) +1 )));
            internal_unpack_SW::unpack_SW<<<_grid, _block, 0, stream>>>(_data_buffer, _data_pdfs, _size_pdfs_0, _size_pdfs_1, _size_pdfs_2, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3);
            break;
        }
        

        default:
            WALBERLA_ASSERT(false);
    }
}


uint_t UniformGridGPU_PackInfo::size(stencil::Direction dir, IBlock * block)
{
    auto pdfs = block->getData< cuda::GPUField<double> >(pdfsID);

    CellInterval ci;
    pdfs->getGhostRegion(dir, ci, 1, false);

    uint_t elementsPerCell = 0;

    switch( dir )
    {
        case stencil::SW:
            elementsPerCell = 1;
            break;
        
        case stencil::BW:
            elementsPerCell = 1;
            break;
        
        case stencil::W:
            elementsPerCell = 5;
            break;
        
        case stencil::TW:
            elementsPerCell = 1;
            break;
        
        case stencil::NW:
            elementsPerCell = 1;
            break;
        
        case stencil::BS:
            elementsPerCell = 1;
            break;
        
        case stencil::S:
            elementsPerCell = 5;
            break;
        
        case stencil::TS:
            elementsPerCell = 1;
            break;
        
        case stencil::B:
            elementsPerCell = 5;
            break;
        
        case stencil::C:
            elementsPerCell = 1;
            break;
        
        case stencil::T:
            elementsPerCell = 5;
            break;
        
        case stencil::BN:
            elementsPerCell = 1;
            break;
        
        case stencil::N:
            elementsPerCell = 5;
            break;
        
        case stencil::TN:
            elementsPerCell = 1;
            break;
        
        case stencil::SE:
            elementsPerCell = 1;
            break;
        
        case stencil::BE:
            elementsPerCell = 1;
            break;
        
        case stencil::E:
            elementsPerCell = 5;
            break;
        
        case stencil::TE:
            elementsPerCell = 1;
            break;
        
        case stencil::NE:
            elementsPerCell = 1;
            break;
        
        default:
            elementsPerCell = 0;
    }
    return ci.numCells() * elementsPerCell * sizeof( double );
}



} // namespace pystencils
} // namespace walberla