//======================================================================================================================
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
//! \\file {{class_name}}.h
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

{% for kernel in kernels.values() %}
{{kernel |generate_definitions(target)}}
{% endfor %}

void {{class_name}}::packEqual({{-[
   "Field_T * " + src_field.name, "CellInterval & ci", "unsigned char * outBuffer",
   kernels['packEqual'].kernel_selection_parameters, ["gpuStream_t stream"] if is_gpu else[]
] | type_identifier_list -}}) const{
   {{dtype}} * buffer = reinterpret_cast<{{dtype}}*>(outBuffer);
   {{kernels['packEqual'] | generate_call(cell_interval='ci', stream='stream') | indent(3) }}
}

void {{class_name}}::unpackEqual({{-[
   "Field_T * " + dst_field.name, "CellInterval & ci", "unsigned char * inBuffer",
   kernels['unpackEqual'].kernel_selection_parameters, ["gpuStream_t stream"] if is_gpu else[]
] | type_identifier_list -}}) const{
   {{dtype}} * buffer = reinterpret_cast<{{dtype}}*>(inBuffer);
   {{kernels['unpackEqual'] | generate_call(cell_interval='ci', stream='stream') | indent(3) }}
}

void {{class_name}}::localCopyEqual(
   {{- [ "Field_T * " + src_field.name, "CellInterval & srcInterval",
          "Field_T * " + dst_field.name, "CellInterval & dstInterval",
          kernels['localCopyEqual'].kernel_selection_parameters,
          ["gpuStream_t stream"] if is_gpu else []]
       | type_identifier_list -}}) const{
   WALBERLA_ASSERT_EQUAL(srcInterval.xSize(), dstInterval.xSize())
   WALBERLA_ASSERT_EQUAL(srcInterval.ySize(), dstInterval.ySize())
   WALBERLA_ASSERT_EQUAL(srcInterval.zSize(), dstInterval.zSize())
   {{kernels['localCopyEqual']
            | generate_call(cell_interval={src_field.name : 'srcInterval', dst_field.name : 'dstInterval'}, stream='stream')
            | indent(3) }}
}

void {{class_name}}::packFineToCoarse({{-[
   "Field_T * " + src_field.name, "CellInterval & ci", "unsigned char * outBuffer",
   kernels['packFineToCoarse'].kernel_selection_parameters,
   ["gpuStream_t stream"] if is_gpu else[]
] | type_identifier_list -}}) const{
   {{dtype}} * buffer = reinterpret_cast<{{dtype}}*>(outBuffer);
   {{kernels['packFineToCoarse'] | generate_call(cell_interval='ci', stream='stream') | indent(3) }}
}

void {{class_name}}::unpackCoarseToFine({{-[
   "Field_T * " + dst_field.name, "CellInterval & ci", "unsigned char * inBuffer",
   kernels['unpackCoarseToFine'].kernel_selection_parameters, ["gpuStream_t stream"] if is_gpu else[]
] | type_identifier_list -}}) const{
   {{dtype}} * buffer = reinterpret_cast<{{dtype}}*>(inBuffer);
   {{kernels['unpackCoarseToFine'] | generate_call(cell_interval='ci', stream='stream') | indent(3) }}
}

}  // namespace {{namespace}}
}  // namespace walberla