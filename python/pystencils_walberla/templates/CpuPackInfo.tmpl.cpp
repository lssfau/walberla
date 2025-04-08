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

#include "stencil/Directions.h"
#include "core/cell/CellInterval.h"
#include "core/DataTypes.h"
#include "{{class_name}}.h"

#if ( defined WALBERLA_CXX_COMPILER_IS_GNU ) || ( defined WALBERLA_CXX_COMPILER_IS_CLANG )
#   pragma GCC diagnostic push
#   pragma GCC diagnostic ignored "-Wfloat-equal"
#   pragma GCC diagnostic ignored "-Wshadow"
#   pragma GCC diagnostic ignored "-Wconversion"
#   pragma GCC diagnostic ignored "-Wunused-variable"
#endif

{% for header in headers %}
#include {{header}}
{% endfor %}

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


void {{class_name}}::pack(Direction dir, unsigned char * byte_buffer, IBlock * block) const
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
            {{kernel|generate_call(cell_interval="ci")|indent(12)}}
            break;
        }
        {% endfor %}

        default:
            WALBERLA_ASSERT(false);
    }
}


void {{class_name}}::unpack(Direction dir, unsigned char * byte_buffer, IBlock * block) const
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
            {{kernel|generate_call(cell_interval="ci")|indent(12)}}
            break;
        }
        {% endfor %}

        default:
            WALBERLA_ASSERT(false);
    }
}


uint_t {{class_name}}::size(stencil::Direction dir, const IBlock * block) const
{
    {{fused_kernel|generate_block_data_to_field_extraction(parameters_to_ignore=['buffer'])|indent(4)}}
    CellInterval ci;
    {{field_name}}->getGhostRegion(dir, ci, 1, false);

    uint_t elementsPerCell = uint_t{ 0u };

    switch( dir )
    {
        {%- for direction_set, elements in elements_per_cell.items()  %}
        {%- for dir in direction_set %}
        case stencil::{{dir}}:
        {%- endfor %}
            elementsPerCell = uint_t{ {{elements}}u };
            break;
        {% endfor %}
        default:
            elementsPerCell = uint_t{ 0u };
    }
    return ci.numCells() * elementsPerCell * sizeof( {{dtype}} );
}



} // namespace {{namespace}}
} // namespace walberla