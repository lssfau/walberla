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

#include "core/DataTypes.h"
#include "core/Macros.h"
#include "{{class_name}}.h"
{% if target == 'gpu' -%}
#include "gpu/ErrorChecking.h"
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
#ifdef __NVCC_DIAG_PRAGMA_SUPPORT__
#pragma nv_diag_suppress 177
#else
#pragma diag_suppress 177
#endif
#endif
//NOLINTBEGIN(readability-non-const-parameter*)
{{kernel|generate_definitions(target)}}
//NOLINTEND(readability-non-const-parameter*)
#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif

#ifdef __CUDACC__
#pragma pop
#endif


void {{class_name}}::run_impl(
   {{- ["IBlock * block", "IndexVectors::Type type",
        kernel.kernel_selection_parameters, ["gpuStream_t stream"] if target == 'gpu' else []]
       | type_identifier_list -}}
)
{
   auto * indexVectors = block->getData<IndexVectors>(indexVectorID);
   int32_t indexVectorSize = int32_c( indexVectors->indexVector(type).size() );
   if( indexVectorSize == 0)
      return;

   {% if target == 'gpu' -%}
   auto pointer = indexVectors->pointerGpu(type);
   {% else %}
   auto pointer = indexVectors->pointerCpu(type);
   {% endif %}

   uint8_t * _data_indexVector = reinterpret_cast<uint8_t*>(pointer);

   {{kernel|generate_block_data_to_field_extraction(['indexVector', 'indexVectorSize'])|indent(4)}}
   {{kernel|generate_timestep_advancements|indent(4)}}
   {{kernel|generate_refs_for_kernel_parameters(prefix='', parameters_to_ignore=['indexVectorSize'], ignore_fields=True)|indent(4) }}
   {{kernel|generate_call(spatial_shape_symbols=['indexVectorSize'], stream='stream')|indent(4)}}
}

void {{class_name}}::run(
   {{- ["IBlock * block", kernel.kernel_selection_parameters, ["gpuStream_t stream"] if target == 'gpu' else []] | type_identifier_list -}}
)
{
   run_impl( {{- ["block", "IndexVectors::ALL", kernel.kernel_selection_parameters, ["stream"] if target == 'gpu' else []] | identifier_list -}} );
}

void {{class_name}}::inner(
   {{- ["IBlock * block", kernel.kernel_selection_parameters, ["gpuStream_t stream"] if target == 'gpu' else []] | type_identifier_list -}}
)
{
   run_impl( {{- ["block", "IndexVectors::INNER", kernel.kernel_selection_parameters, ["stream"] if target == 'gpu' else []] | identifier_list -}} );
}

void {{class_name}}::outer(
   {{- ["IBlock * block", kernel.kernel_selection_parameters, ["gpuStream_t stream"] if target == 'gpu' else []] | type_identifier_list -}}
)
{
   run_impl( {{- ["block", "IndexVectors::OUTER", kernel.kernel_selection_parameters, ["stream"] if target == 'gpu' else []] | identifier_list -}} );
}

} // namespace {{namespace}}
} // namespace walberla


