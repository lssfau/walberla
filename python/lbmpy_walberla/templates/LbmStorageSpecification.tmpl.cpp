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
//! \\file {{class_name}}.cpp
//! \\author lbmpy
//======================================================================================================================

#include "{{class_name}}.h"


{% if target is equalto 'cpu' -%}
#define FUNC_PREFIX
{%- elif target is equalto 'gpu' -%}
#define FUNC_PREFIX __global__
#include "gpu/GPUWrapper.h"
#include "gpu/GPUField.h"
{%- endif %}

#ifdef __GNUC__
#define RESTRICT __restrict__
#elif _MSC_VER
#define RESTRICT __restrict
#else
#define RESTRICT
#endif


#if ( defined WALBERLA_CXX_COMPILER_IS_GNU ) || ( defined WALBERLA_CXX_COMPILER_IS_CLANG )
#   pragma GCC diagnostic push
#   pragma GCC diagnostic ignored "-Wfloat-equal"
#   pragma GCC diagnostic ignored "-Wshadow"
#   pragma GCC diagnostic ignored "-Wconversion"
#   pragma GCC diagnostic ignored "-Wunused-variable"
#   pragma GCC diagnostic ignored "-Wignored-qualifiers"
#endif

#if ( defined WALBERLA_CXX_COMPILER_IS_INTEL )
#pragma warning push
#pragma warning( disable :  1599 )
#endif

#ifdef __CUDACC__
#pragma push
#ifdef __NVCC_DIAG_PRAGMA_SUPPORT__
#pragma nv_diag_suppress 191
#else
#pragma diag_suppress 191
#endif
#endif

namespace walberla {
namespace {{namespace}} {

   /*************************************************************************************
 *                                Kernel Definitions
*************************************************************************************/
   {{ kernels['packAll']      | generate_definitions }}
   {{ kernels['unpackAll']    | generate_definitions }}
   {{ kernels['localCopyAll'] | generate_definitions }}

   {{ kernels['packDirection']      | generate_definitions }}
   {{ kernels['unpackDirection']    | generate_definitions }}
   {{ kernels['localCopyDirection'] | generate_definitions }}

   {% if nonuniform -%}
   {{ kernels['localCopyRedistribute']    | generate_definitions }}
   {{ kernels['localPartialCoalescence']    | generate_definitions }}
   {{ kernels['unpackRedistribute']    | generate_definitions }}
   {{ kernels['packPartialCoalescence']    | generate_definitions }}
   {{ kernels['zeroCoalescenceRegion']    | generate_definitions }}
   {{ kernels['unpackCoalescence']    | generate_definitions }}
   {%- endif %}

   /*************************************************************************************
 *                                 Kernel Wrappers
*************************************************************************************/

   void {{class_name}}::PackKernels::packAll(
      {{- [ "PdfField_T * " + src_field.name, "CellInterval & ci",
             "unsigned char * outBuffer", kernels['packAll'].kernel_selection_parameters,
             ["gpuStream_t stream"] if is_gpu else []]
          | type_identifier_list -}}
   ) const
   {
      {{dtype}} * buffer = reinterpret_cast<{{dtype}}*>(outBuffer);
      {{kernels['packAll'] | generate_call(cell_interval='ci', stream='stream') | indent(6) }}
   }


   void {{class_name}}::PackKernels::unpackAll(
      {{- [ "PdfField_T * " + dst_field.name, "CellInterval & ci",
             "unsigned char * inBuffer", kernels['unpackAll'].kernel_selection_parameters,
             ["gpuStream_t stream"] if is_gpu else []]
          | type_identifier_list -}}
   ) const
   {
      {{dtype}} * buffer = reinterpret_cast<{{dtype}}*>(inBuffer);
      {{kernels['unpackAll'] | generate_call(cell_interval='ci', stream='stream') | indent(6) }}
   }


   void {{class_name}}::PackKernels::localCopyAll(
      {{- [ "PdfField_T * " + src_field.name, "CellInterval & srcInterval",
             "PdfField_T * " + dst_field.name, "CellInterval & dstInterval",
             kernels['localCopyAll'].kernel_selection_parameters,
             ["gpuStream_t stream"] if is_gpu else []]
          | type_identifier_list -}}
   ) const
   {
      WALBERLA_ASSERT_EQUAL(srcInterval.xSize(), dstInterval.xSize())
      WALBERLA_ASSERT_EQUAL(srcInterval.ySize(), dstInterval.ySize())
      WALBERLA_ASSERT_EQUAL(srcInterval.zSize(), dstInterval.zSize())

      {{kernels['localCopyAll']
               | generate_call(cell_interval={src_field.name : 'srcInterval', dst_field.name : 'dstInterval'}, stream='stream')
               | indent(6) }}
   }

   void {{class_name}}::PackKernels::packDirection(
      {{- [ "PdfField_T * " + src_field.name, "CellInterval & ci",
             "unsigned char * outBuffer", kernels['packDirection'].kernel_selection_parameters,
             ["gpuStream_t stream"] if is_gpu else []]
          | type_identifier_list -}}
   ) const
   {
      {{dtype}} * buffer = reinterpret_cast<{{dtype}}*>(outBuffer);
      {{kernels['packDirection'] | generate_call(cell_interval='ci', stream='stream') | indent(6) }}
   }

   void {{class_name}}::PackKernels::unpackDirection(
      {{- [ "PdfField_T * " + dst_field.name, "CellInterval & ci",
             "unsigned char * inBuffer", kernels['unpackDirection'].kernel_selection_parameters,
             ["gpuStream_t stream"] if is_gpu else []]
          | type_identifier_list -}}
   ) const
   {
      {{dtype}} * buffer = reinterpret_cast<{{dtype}}*>(inBuffer);
      {{kernels['unpackDirection'] | generate_call(cell_interval='ci', stream='stream') | indent(6) }}
   }

   void {{class_name}}::PackKernels::localCopyDirection(
      {{- [src_field.dtype.c_name + "** _data_" + src_field.name + "_dp", dst_field.dtype.c_name + "** _data_" + dst_field.name + "_dp",
            kernels['localCopyDirection'].kernel_selection_parameters,
            ["gpuStream_t stream"] if is_gpu else [], "std::array<int64_t, 4>& _sizes", "std::array<int64_t, 4>& _strides"]
          | type_identifier_list -}}
   ) const {
      {% if block_wise -%}

      {% if target is equalto 'gpu' -%}

      const int64_t indexingX = _sizes[0] * _sizes[1];
      const int64_t indexingY = _sizes[2];
      const int64_t indexingZ = _sizes[3];

      const int64_t cudaBlockSize0 = 128;
      const int64_t cudaBlockSize1 = 1;
      const int64_t cudaBlockSize2 = 1;

      const int64_t _size_0 = _sizes[0];

      const int64_t _size_{{src_field.name}}_0 = _sizes[1];
      const int64_t _size_{{src_field.name}}_1 = _sizes[2];
      const int64_t _size_{{src_field.name}}_2 = _sizes[3];
      const int64_t _size_{{dst_field.name}}_0 = _sizes[1];
      const int64_t _size_{{dst_field.name}}_1 = _sizes[2];
      const int64_t _size_{{dst_field.name}}_2 = _sizes[3];

      const int64_t _stride_{{src_field.name}}_0 = _strides[0];
      const int64_t _stride_{{src_field.name}}_1 = _strides[1];
      const int64_t _stride_{{src_field.name}}_2 = _strides[2];
      const int64_t _stride_{{src_field.name}}_3 = _strides[3];

      const int64_t _stride_{{dst_field.name}}_0 = _strides[0];
      const int64_t _stride_{{dst_field.name}}_1 = _strides[1];
      const int64_t _stride_{{dst_field.name}}_2 = _strides[2];
      const int64_t _stride_{{dst_field.name}}_3 = _strides[3];

      const dim3 _block = dim3((unsigned int)((1024 < ((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)) ? 1024 : ((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)), (unsigned int)((1024 < ((indexingY < cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)))) ? indexingY : cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0))))) ? 1024 : ((indexingY < cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)))) ? indexingY : cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0))))), (unsigned int)((64 < ((indexingZ < cudaBlockSize2*((int64_t)(cudaBlockSize0*cudaBlockSize1) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)*((indexingY < cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)))) ? indexingY : cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0))))))) ? indexingZ : cudaBlockSize2*((int64_t)(cudaBlockSize0*cudaBlockSize1) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)*((indexingY < cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)))) ? indexingY : cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)))))))) ? 64 : ((indexingZ < cudaBlockSize2*((int64_t)(cudaBlockSize0*cudaBlockSize1) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)*((indexingY < cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)))) ? indexingY : cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0))))))) ? indexingZ : cudaBlockSize2*((int64_t)(cudaBlockSize0*cudaBlockSize1) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)*((indexingY < cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)))) ? indexingY : cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)))))))));
      const dim3 _grid  = dim3((unsigned int)(( (indexingX) % (((1024 < ((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)) ? 1024 : ((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0))) == 0 ? (int64_t)(indexingX) / (int64_t)(((1024 < ((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)) ? 1024 : ((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0))) : ( (int64_t)(indexingX) / (int64_t)(((1024 < ((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)) ? 1024 : ((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0))) ) +1 )), (unsigned int)(( (indexingY) % (((1024 < ((indexingY < cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)))) ? indexingY : cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0))))) ? 1024 : ((indexingY < cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)))) ? indexingY : cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)))))) == 0 ? (int64_t)(indexingY) / (int64_t)(((1024 < ((indexingY < cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)))) ? indexingY : cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0))))) ? 1024 : ((indexingY < cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)))) ? indexingY : cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)))))) : ( (int64_t)(indexingY) / (int64_t)(((1024 < ((indexingY < cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)))) ? indexingY : cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0))))) ? 1024 : ((indexingY < cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)))) ? indexingY : cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)))))) ) +1 )), (unsigned int)(( (indexingZ) % (((64 < ((indexingZ < cudaBlockSize2*((int64_t)(cudaBlockSize0*cudaBlockSize1) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)*((indexingY < cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)))) ? indexingY : cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0))))))) ? indexingZ : cudaBlockSize2*((int64_t)(cudaBlockSize0*cudaBlockSize1) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)*((indexingY < cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)))) ? indexingY : cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)))))))) ? 64 : ((indexingZ < cudaBlockSize2*((int64_t)(cudaBlockSize0*cudaBlockSize1) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)*((indexingY < cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)))) ? indexingY : cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0))))))) ? indexingZ : cudaBlockSize2*((int64_t)(cudaBlockSize0*cudaBlockSize1) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)*((indexingY < cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)))) ? indexingY : cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0))))))))) == 0 ? (int64_t)(indexingZ) / (int64_t)(((64 < ((indexingZ < cudaBlockSize2*((int64_t)(cudaBlockSize0*cudaBlockSize1) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)*((indexingY < cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)))) ? indexingY : cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0))))))) ? indexingZ : cudaBlockSize2*((int64_t)(cudaBlockSize0*cudaBlockSize1) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)*((indexingY < cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)))) ? indexingY : cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)))))))) ? 64 : ((indexingZ < cudaBlockSize2*((int64_t)(cudaBlockSize0*cudaBlockSize1) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)*((indexingY < cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)))) ? indexingY : cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0))))))) ? indexingZ : cudaBlockSize2*((int64_t)(cudaBlockSize0*cudaBlockSize1) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)*((indexingY < cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)))) ? indexingY : cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0))))))))) : ( (int64_t)(indexingZ) / (int64_t)(((64 < ((indexingZ < cudaBlockSize2*((int64_t)(cudaBlockSize0*cudaBlockSize1) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)*((indexingY < cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)))) ? indexingY : cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0))))))) ? indexingZ : cudaBlockSize2*((int64_t)(cudaBlockSize0*cudaBlockSize1) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)*((indexingY < cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)))) ? indexingY : cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)))))))) ? 64 : ((indexingZ < cudaBlockSize2*((int64_t)(cudaBlockSize0*cudaBlockSize1) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)*((indexingY < cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)))) ? indexingY : cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0))))))) ? indexingZ : cudaBlockSize2*((int64_t)(cudaBlockSize0*cudaBlockSize1) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)*((indexingY < cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)))) ? indexingY : cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0))))))))) ) +1 )));
      {%- endif %}

      {{kernels['localCopyDirection']
               | generate_call(plain_kernel_call=True, stream='stream')
               | indent(6) }}

      {%else%}
      WALBERLA_ABORT("Block wise local communication is not implemented")
      {%- endif %}
   }

   void {{class_name}}::PackKernels::localCopyDirection(
      {{- [ "PdfField_T * " + src_field.name, "CellInterval & srcInterval",
             "PdfField_T * " + dst_field.name, "CellInterval & dstInterval",
             kernels['localCopyDirection'].kernel_selection_parameters,
             ["gpuStream_t stream"] if is_gpu else []]
          | type_identifier_list -}}
   ) const
   {
      {% if not block_wise -%}
      WALBERLA_ASSERT_EQUAL(srcInterval.xSize(), dstInterval.xSize())
      WALBERLA_ASSERT_EQUAL(srcInterval.ySize(), dstInterval.ySize())
      WALBERLA_ASSERT_EQUAL(srcInterval.zSize(), dstInterval.zSize())

      {{kernels['localCopyDirection']
               | generate_call(cell_interval={src_field.name : 'srcInterval', dst_field.name : 'dstInterval'}, stream='stream')
               | indent(6) }}

      {%else%}
      WALBERLA_ABORT("Local communication is only implemented block wise")
      {%- endif %}
   }


   {% if nonuniform -%}
   void {{class_name}}::PackKernels::localCopyRedistribute(
      {{- [  "PdfField_T * " + src_field.name, "CellInterval & srcInterval",
             "PdfField_T * " + dst_field.name, "CellInterval & dstInterval", kernels['localCopyRedistribute'].kernel_selection_parameters,
             ["gpuStream_t stream"] if is_gpu else []]
          | type_identifier_list -}}
   ) const
   {
      {{kernels['localCopyRedistribute'] | generate_call(cell_interval={src_field.name : 'srcInterval', dst_field.name : 'dstInterval'}, stream='stream') | indent(6) }}
   }

   void {{class_name}}::PackKernels::localPartialCoalescence(
      {{- [  "PdfField_T * " + src_field.name, "MaskField_T * " + mask_field.name, "CellInterval & srcInterval",
             "PdfField_T * " + dst_field.name, "CellInterval & dstInterval", kernels['localPartialCoalescence'].kernel_selection_parameters,
             ["gpuStream_t stream"] if is_gpu else []]
          | type_identifier_list -}}
   ) const
   {
      {{kernels['localPartialCoalescence'] | generate_call(cell_interval={src_field.name : 'srcInterval', mask_field.name : 'srcInterval', dst_field.name : 'dstInterval'}, stream='stream') | indent(6) }}
   }

   void {{class_name}}::PackKernels::unpackRedistribute(
      {{- [ "PdfField_T * " + dst_field.name, "CellInterval & ci",
             "unsigned char * inBuffer", kernels['unpackDirection'].kernel_selection_parameters,
             ["gpuStream_t stream"] if is_gpu else []]
          | type_identifier_list -}}
   ) const
   {
      {{dtype}} * buffer = reinterpret_cast<{{dtype}}*>(inBuffer);
      {{kernels['unpackRedistribute'] | generate_call(cell_interval='ci', stream='stream') | indent(6) }}
   }

   void {{class_name}}::PackKernels::packPartialCoalescence(
      {{- [ "PdfField_T * " + src_field.name, "MaskField_T * " + mask_field.name, "CellInterval & ci",
             "unsigned char * outBuffer", kernels['packPartialCoalescence'].kernel_selection_parameters,
             ["gpuStream_t stream"] if is_gpu else []]
          | type_identifier_list -}}
   ) const
   {
      {{dtype}} * buffer = reinterpret_cast<{{dtype}}*>(outBuffer);
      {{kernels['packPartialCoalescence'] | generate_call(cell_interval='ci', stream='stream') | indent(6) }}
   }

   void {{class_name}}::PackKernels::zeroCoalescenceRegion(
      {{- [ "PdfField_T * " + dst_field.name, "CellInterval & ci",
             kernels['zeroCoalescenceRegion'].kernel_selection_parameters,
             ["gpuStream_t stream"] if is_gpu else []]
          | type_identifier_list -}}
   ) const
   {
      {{kernels['zeroCoalescenceRegion'] | generate_call(cell_interval='ci', stream='stream') | indent(6) }}
   }

   void {{class_name}}::PackKernels::unpackCoalescence(
      {{- [ "PdfField_T * " + dst_field.name, "CellInterval & ci",
             "unsigned char * inBuffer", kernels['unpackCoalescence'].kernel_selection_parameters,
             ["gpuStream_t stream"] if is_gpu else []]
          | type_identifier_list -}}
   ) const
   {
      {{dtype}} * buffer = reinterpret_cast<{{dtype}}*>(inBuffer);
      {{kernels['unpackCoalescence'] | generate_call(cell_interval='ci', stream='stream') | indent(6) }}
   }
   {%- endif %}
}  // namespace {{namespace}}
}  // namespace walberla