#include "stencil/Directions.h"
#include "core/cell/CellInterval.h"
#include "cuda/GPUField.h"
#include "core/DataTypes.h"
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



void {{class_name}}::pack(Direction dir, unsigned char * byte_buffer, IBlock * block, cudaStream_t stream)
{
    {{dtype}} * buffer = reinterpret_cast<{{dtype}}*>(byte_buffer);

    {{fused_kernel|generate_block_data_to_field_extraction(parameters_to_ignore=['buffer'])|indent(4)}}
    CellInterval ci;
    {{field_name}}->getSliceBeforeGhostLayer(dir, ci, 1, false);

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
            WALBERLA_ASSERT(false);
    }
}


void {{class_name}}::unpack(Direction dir, unsigned char * byte_buffer, IBlock * block, cudaStream_t stream)
{
    {{dtype}} * buffer = reinterpret_cast<{{dtype}}*>(byte_buffer);

    {{fused_kernel|generate_block_data_to_field_extraction(parameters_to_ignore=['buffer'])|indent(4)}}
    CellInterval ci;
    {{field_name}}->getGhostRegion(dir, ci, 1, false);
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
            WALBERLA_ASSERT(false);
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