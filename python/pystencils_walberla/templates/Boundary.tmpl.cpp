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

#include <cmath>

#include "core/DataTypes.h"
#include "core/Macros.h"
#include "{{class_name}}.h"
{% if target == 'gpu' -%}
#include "cuda/ErrorChecking.h"
{%- endif %}


{% if target == 'cpu' -%}
#define FUNC_PREFIX
{%- elif target == 'gpu' -%}
#define FUNC_PREFIX __global__
{%- endif %}

using namespace std;

namespace walberla {
namespace {{namespace}} {

#ifdef __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wstrict-aliasing"
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wconversion"
#endif

#ifdef __CUDACC__
#pragma push
#pragma diag_suppress = declared_but_not_referenced
#endif


{{kernel|generate_definition}}

#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif

#ifdef __CUDACC__
#pragma pop
#endif


void {{class_name}}::run( IBlock * block, IndexVectors::Type type {% if target == 'gpu'%}, cudaStream_t stream {%endif%})
{
    auto * indexVectors = block->getData<IndexVectors>(indexVectorID);
    int64_t indexVectorSize = int64_c( indexVectors->indexVector(type).size() );
    if( indexVectorSize == 0)
        return;

    {% if target == 'gpu' -%}
    auto pointer = indexVectors->pointerGpu(type);
    {% else %}
    auto pointer = indexVectors->pointerCpu(type);
    {% endif %}

    uint8_t * _data_indexVector = reinterpret_cast<uint8_t*>(pointer);

    {{kernel|generate_block_data_to_field_extraction(['indexVector', 'indexVectorSize'])|indent(4)}}
    {{kernel|generate_call(spatial_shape_symbols=['indexVectorSize'], stream='stream')|indent(4)}}
}

void {{class_name}}::operator() ( IBlock * block{% if target == 'gpu'%}, cudaStream_t stream {%endif%} )
{
    run( block, IndexVectors::ALL{% if target == 'gpu'%}, stream {%endif%});
}

void {{class_name}}::inner( IBlock * block{% if target == 'gpu'%}, cudaStream_t stream {%endif%} )
{
    run( block, IndexVectors::INNER{% if target == 'gpu'%}, stream {%endif%} );
}

void {{class_name}}::outer( IBlock * block{% if target == 'gpu'%}, cudaStream_t stream {%endif%} )
{
    run( block, IndexVectors::OUTER{% if target == 'gpu'%}, stream {%endif%} );
}


} // namespace {{namespace}}
} // namespace walberla


