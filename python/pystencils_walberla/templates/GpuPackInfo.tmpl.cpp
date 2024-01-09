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
//! \\author pystencils
//======================================================================================================================

#include "{{class_name}}.h"

{% if target is equalto 'cpu' -%}
#define FUNC_PREFIX
{%- elif target is equalto 'gpu' -%}
#define FUNC_PREFIX __global__
{%- endif %}


namespace walberla {
namespace {{namespace}} {

using walberla::cell::CellInterval;
using walberla::stencil::Direction;


{% for kernel in pack_kernels.values() %}
{{kernel|generate_definition(target)}}
{% endfor %}

{% for kernel in unpack_kernels.values() %}
{{kernel|generate_definition(target)}}
{% endfor %}



void {{class_name}}::pack(Direction dir, unsigned char * byte_buffer, IBlock * block, gpuStream_t stream)
{
    {{dtype}} * buffer = reinterpret_cast<{{dtype}}*>(byte_buffer);

    {{fused_kernel|generate_block_data_to_field_extraction(parameters_to_ignore=['buffer'])|indent(4)}}
    CellInterval ci;
    {% if gl_to_inner -%}
    {{field_name}}->getGhostRegion(dir, ci, 1, false);
    {%- else -%}
    {{field_name}}->getSliceBeforeGhostLayer(dir, ci, 1, false);
    {%- endif %}

    switch( dir )
    {
        {%- for direction_set, kernel in pack_kernels.items()  %}
        {%- for dir in direction_set %}
        case stencil::{{dir}}:
        {%- endfor %}
        {
            {{kernel|generate_call(cell_interval="ci", stream="stream")|indent(12)}}
            break;
        }
        {% endfor %}

        default:
            return;
    }
}


void {{class_name}}::unpack(Direction dir, unsigned char * byte_buffer, IBlock * block, gpuStream_t stream)
{
    {{dtype}} * buffer = reinterpret_cast<{{dtype}}*>(byte_buffer);

    {{fused_kernel|generate_block_data_to_field_extraction(parameters_to_ignore=['buffer'])|indent(4)}}
    CellInterval ci;
    {% if gl_to_inner -%}
    {{field_name}}->getSliceBeforeGhostLayer(dir, ci, 1, false);
    {%- else -%}
    {{field_name}}->getGhostRegion(dir, ci, 1, false);
    {%- endif %}
    auto communciationDirection = stencil::inverseDir[dir];

    switch( communciationDirection )
    {
        {%- for direction_set, kernel in unpack_kernels.items()  %}
        {%- for dir in direction_set %}
        case stencil::{{dir}}:
        {%- endfor %}
        {
            {{kernel|generate_call(cell_interval="ci", stream="stream")|indent(12)}}
            break;
        }
        {% endfor %}

        default:
            return;
    }
}


uint_t {{class_name}}::size(stencil::Direction dir, IBlock * block)
{
    {{fused_kernel|generate_block_data_to_field_extraction(parameters_to_ignore=['buffer'])|indent(4)}}
    CellInterval ci;
    {{field_name}}->getGhostRegion(dir, ci, 1, false);

    uint_t elementsPerCell = 0;

    switch( dir )
    {
        {%- for direction_set, elements in elements_per_cell.items()  %}
        {%- for dir in direction_set %}
        case stencil::{{dir}}:
        {%- endfor %}
            elementsPerCell = {{elements}};
            break;
        {% endfor %}
        default:
            elementsPerCell = 0;
    }
    return ci.numCells() * elementsPerCell * sizeof( {{dtype}} );
}



} // namespace {{namespace}}
} // namespace walberla