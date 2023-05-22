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

#if ( defined WALBERLA_CXX_COMPILER_IS_GNU ) || ( defined WALBERLA_CXX_COMPILER_IS_CLANG )
#   pragma GCC diagnostic push
#   pragma GCC diagnostic ignored "-Wunused-variable"
#endif

#if ( defined WALBERLA_CXX_COMPILER_IS_INTEL )
#pragma warning push
#pragma warning( disable :  1599 )
#endif

using namespace std;

namespace walberla {
namespace {{namespace}} {

{% for kernel in kernels %}
{{kernel['kernel']|generate_definitions(target, kernel['max_threads'])}}
{% endfor %}


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
