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
{%- endif %}

#if ( defined WALBERLA_CXX_COMPILER_IS_GNU ) || ( defined WALBERLA_CXX_COMPILER_IS_CLANG )
#   pragma GCC diagnostic push
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

using namespace std;

namespace walberla {
namespace {{namespace}} {

{%if block_stream_collide -%}
{{block_stream_collide['kernel']|generate_definitions(target, block_stream_collide['max_threads'])}}
{{block_stream['kernel']|generate_definitions(target, block_stream['max_threads'])}}
{%endif%}

{% for kernel in kernels %}
{{kernel['kernel']|generate_definitions(target, kernel['max_threads'])}}
{% endfor %}

void {{class_name}}::blockStreamCollide({{- ["[[maybe_unused]] uint_t level", "[[maybe_unused]] uint8_t timestep", ["[[maybe_unused]] gpuStream_t stream"] if target == 'gpu' else []] | type_identifier_list -}})
{
   {%if block_stream_collide -%}

   {%if target is equalto 'gpu' -%}
   dim3 _grid = grid_[level];
   dim3 _block = block_[level];

   {%- for field in block_stream_collide['all_fields'] %}
   {{field.dtype.c_name}} ** {{block_stream_collide['indexed_to_field_name'][field.name]}} = {{field.name}}PointersGPU[level];
   {%- endfor %}

   {% else %}

   {%- for field in block_stream_collide['all_fields'] %}
   {{field.dtype.c_name}} ** {{block_stream_collide['indexed_to_field_name'][field.name]}} = {{field.name}}Pointers[level].data();
   {%- endfor %}

   {%- endif %}
   const int64_t _size_0 = size_0[level];
   int64_t _size_{{block_stream_collide['all_fields'][0].name}}_0 = size_1;
   int64_t _size_{{block_stream_collide['all_fields'][0].name}}_1 = size_2;
   int64_t _size_{{block_stream_collide['all_fields'][0].name}}_2 = size_3;

   {{block_stream_collide['kernel']|generate_field_strides()|indent(3)}}
   {{block_stream_collide['kernel']|generate_refs_for_kernel_parameters(prefix="this->", parameters_to_ignore=["_size_0"], ignore_fields=True, parameter_registration=parameter_scaling, level_known=True)|indent(3)}}
   {{block_stream_collide['kernel']|generate_call(stream='stream', plain_kernel_call=True)|indent(3)}}

   {%endif%}
}

void {{class_name}}::ghostLayerPropagation({{- ["[[maybe_unused]] uint_t level", "[[maybe_unused]] uint8_t timestep", ["[[maybe_unused]] gpuStream_t stream"] if target == 'gpu' else []] | type_identifier_list -}})
{
   {%if block_stream_collide -%}

   {{block_stream['kernel']|generate_field_strides()|indent(3)}}

   {%if target is equalto 'gpu' -%}
   auto parallelSection_ = parallelStreams_.parallelSection( stream );
   for (auto it = glPropagationPDFs[level].begin(); it != glPropagationPDFs[level].end(); it++){
      if(it->second.empty()){ continue;}

      int64_t _size_0 = int64_c(it->second.size());
      int64_t _size_{{pdf_field.name}}_0 = std::get<0>(it->first);
      int64_t _size_{{pdf_field.name}}_1 = std::get<1>(it->first);
      int64_t _size_{{pdf_field.name}}_2 = std::get<2>(it->first);

      {{pdf_field.dtype.c_name}} ** _data_{{pdf_field.name}}_dp = glPropagationPDFsGPU[level][it->first];
      dim3 _grid = glPropagationGrid_[level][it->first];
      dim3 _block = glPropagationBlock_[level][it->first];
      parallelSection_.run([&]( auto s ) {
      {{block_stream['kernel']|generate_call(stream='s', plain_kernel_call=True)|indent(9)}}
      });
   }

   {% else %}

   for (auto it = glPropagationPDFs[level].begin(); it != glPropagationPDFs[level].end(); it++){
      if(it->second.empty()){ continue;}

      int64_t _size_0 = int64_c(it->second.size());
      int64_t _size_{{pdf_field.name}}_0 = std::get<0>(it->first);
      int64_t _size_{{pdf_field.name}}_1 = std::get<1>(it->first);
      int64_t _size_{{pdf_field.name}}_2 = std::get<2>(it->first);

      {{pdf_field.dtype.c_name}} ** _data_{{pdf_field.name}}_dp = it->second.data();
      {{block_stream['kernel']|generate_call(stream='s', plain_kernel_call=True)|indent(6)}}
   }
   {%- endif %}

   {%endif%}
}

{% for kernel in kernels %}
void {{class_name}}::{{kernel['function_name']}}( {{kernel['kernel']|generate_plain_parameter_list(ghost_layers=True)}} )
{
   {{kernel['kernel']|generate_call(ghost_layers_to_include=kernel['ghost_layers_to_include'], stream='stream')|indent(3)}}
}
void {{class_name}}::{{kernel['function_name']}}CellInterval( {{kernel['kernel']|generate_plain_parameter_list(cell_interval='ci')}})
{
   {{kernel['kernel']|generate_call(stream='stream', cell_interval='ci')|indent(3)}}
}
{% endfor %}


} // namespace {{namespace}}
} // namespace walberla


#if ( defined WALBERLA_CXX_COMPILER_IS_GNU ) || ( defined WALBERLA_CXX_COMPILER_IS_CLANG )
#   pragma GCC diagnostic pop
#endif

#if ( defined WALBERLA_CXX_COMPILER_IS_INTEL )
#pragma warning pop
#endif
