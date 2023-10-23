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

#if ( defined WALBERLA_CXX_COMPILER_IS_GNU ) || ( defined WALBERLA_CXX_COMPILER_IS_CLANG )
#   pragma GCC diagnostic push
#   pragma GCC diagnostic ignored "-Wfloat-equal"
#   pragma GCC diagnostic ignored "-Wshadow"
#   pragma GCC diagnostic ignored "-Wconversion"
#   pragma GCC diagnostic ignored "-Wunused-variable"
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
               | generate_call(cell_interval={src_field : 'srcInterval', dst_field : 'dstInterval'}, stream='stream')
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
      {{- [ "PdfField_T * " + src_field.name, "CellInterval & srcInterval",
             "PdfField_T * " + dst_field.name, "CellInterval & dstInterval",
             kernels['localCopyDirection'].kernel_selection_parameters,
             ["gpuStream_t stream"] if is_gpu else []]
          | type_identifier_list -}}
   ) const
   {
      WALBERLA_ASSERT_EQUAL(srcInterval.xSize(), dstInterval.xSize())
      WALBERLA_ASSERT_EQUAL(srcInterval.ySize(), dstInterval.ySize())
      WALBERLA_ASSERT_EQUAL(srcInterval.zSize(), dstInterval.zSize())

      {{kernels['localCopyDirection']
               | generate_call(cell_interval={src_field : 'srcInterval', dst_field : 'dstInterval'}, stream='stream')
               | indent(6) }}
   }

   {% if nonuniform -%}
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