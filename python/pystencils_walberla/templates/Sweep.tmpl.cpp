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
//! \\ingroup lbm
//! \\author lbmpy
//======================================================================================================================

#include <cmath>

#include "core/DataTypes.h"
#include "core/Macros.h"
#include "{{class_name}}.h"
{% for header in headers %}
#include {{header}}
{% endfor %}


{% if target is equalto 'cpu' -%}
#define FUNC_PREFIX
{%- elif target is equalto 'gpu' -%}
#define FUNC_PREFIX __global__
{%- endif %}

#if ( defined WALBERLA_CXX_COMPILER_IS_GNU ) || ( defined WALBERLA_CXX_COMPILER_IS_CLANG )
#   pragma GCC diagnostic push
#   pragma GCC diagnostic ignored "-Wfloat-equal"
#   pragma GCC diagnostic ignored "-Wshadow"
#   pragma GCC diagnostic ignored "-Wconversion"
#   pragma GCC diagnostic ignored "-Wunused-variable"
#endif

#if ( defined WALBERLA_CXX_COMPILER_IS_INTEL )
#pragma warning push
#pragma warning( disable :  1599 )
#endif

using namespace std;

namespace walberla {
namespace {{namespace}} {


{{kernel|generate_definition(target)}}

void {{class_name}}::operator()( IBlock * block{%if target is equalto 'gpu'%} , cudaStream_t stream{% endif %} )
{
    {{kernel|generate_block_data_to_field_extraction|indent(4)}}
    {{kernel|generate_refs_for_kernel_parameters(prefix='this->', ignore_fields=True)|indent(4) }}
    {{kernel|generate_call(ghost_layers_to_include=ghost_layers_to_include, stream='stream')|indent(4)}}
    {{kernel|generate_swaps|indent(4)}}
}


void {{class_name}}::runOnCellInterval( const shared_ptr<StructuredBlockStorage> & blocks,
                                        const CellInterval & globalCellInterval,
                                        cell_idx_t ghostLayers,
                                        IBlock * block{%if target is equalto 'gpu'%} , cudaStream_t stream{% endif %} )
{
    CellInterval ci = globalCellInterval;
    CellInterval blockBB = blocks->getBlockCellBB( *block);
    blockBB.expand( ghostLayers );
    ci.intersect( blockBB );
    blocks->transformGlobalToBlockLocalCellInterval( ci, *block );
    if( ci.empty() )
        return;

    {{kernel|generate_block_data_to_field_extraction|indent(4)}}
    {{kernel|generate_refs_for_kernel_parameters(prefix='this->', ignore_fields=True)|indent(4) }}
    {{kernel|generate_call(stream='stream', cell_interval='ci')|indent(4)}}
    {{kernel|generate_swaps|indent(4)}}
}


} // namespace {{namespace}}
} // namespace walberla


#if ( defined WALBERLA_CXX_COMPILER_IS_GNU ) || ( defined WALBERLA_CXX_COMPILER_IS_CLANG )
#   pragma GCC diagnostic pop
#endif

#if ( defined WALBERLA_CXX_COMPILER_IS_INTEL )
#pragma warning pop
#endif
