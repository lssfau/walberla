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

#include <cmath>

#include "core/DataTypes.h"
#include "core/Macros.h"
#include "lbm/field/PdfField.h"
#include "lbm/sweeps/Streaming.h"
#include "{{class_name}}.h"

#ifdef _MSC_VER
#  pragma warning( disable : 4458 )
#endif

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
#   pragma GCC diagnostic ignored "-Wunused-parameter"
#endif

{% for header in headers %}
#include {{header}}
{% endfor %}


using namespace std;

namespace walberla {
namespace {{namespace}} {

{{stream_collide_kernel|generate_definition(target)}}
{{collide_kernel|generate_definition(target)}}
{{stream_kernel|generate_definition(target)}}


const {{dtype}} {{class_name}}::w[{{Q}}] = { {{weights}} };
const {{dtype}} {{class_name}}::wInv[{{Q}}] = { {{inverse_weights}} };

void {{class_name}}::Sweep::streamCollide( IBlock * block, const uint_t numberOfGhostLayersToInclude )
{
    {{stream_collide_kernel|generate_block_data_to_field_extraction(parameters=['pdfs', 'pdfs_tmp'])|indent(4)}}

    auto & lm = dynamic_cast< lbm::PdfField<{{class_name}}> * > (pdfs)->latticeModel();
    WALBERLA_ASSERT_EQUAL( *(lm.blockId_), block->getId() );

    {{stream_collide_kernel|generate_refs_for_kernel_parameters(prefix='lm.', parameters_to_ignore=['pdfs', 'pdfs_tmp'])|indent(4) }}
    {{stream_collide_kernel|generate_call('cell_idx_c(numberOfGhostLayersToInclude)')|indent(4)}}
    {{stream_collide_kernel|generate_swaps|indent(4)}}
}

void {{class_name}}::Sweep::collide( IBlock * block, const uint_t numberOfGhostLayersToInclude )
{
   {{collide_kernel|generate_block_data_to_field_extraction(parameters=['pdfs'])|indent(4)}}

    auto & lm = dynamic_cast< lbm::PdfField<{{class_name}}> * > (pdfs)->latticeModel();
    WALBERLA_ASSERT_EQUAL( *(lm.blockId_), block->getId() );

    {{collide_kernel|generate_refs_for_kernel_parameters(prefix='lm.', parameters_to_ignore=['pdfs', 'pdfs_tmp'])|indent(4) }}
    {{collide_kernel|generate_call('cell_idx_c(numberOfGhostLayersToInclude)')|indent(4)}}
}


void {{class_name}}::Sweep::stream( IBlock * block, const uint_t numberOfGhostLayersToInclude )
{
    {{stream_kernel|generate_block_data_to_field_extraction(parameters=['pdfs', 'pdfs_tmp'])|indent(4)}}

    {{stream_kernel|generate_call('cell_idx_c(numberOfGhostLayersToInclude)')|indent(4)}}

    {{stream_kernel|generate_swaps|indent(4)}}
}

// IMPORTANT REMARK:
// This is specifically implemented for using generated kernels in the waLBerla's free surface LBM and is
// implemented in rather unflexible fashion. Therefore, it should not be extended and in the long-term, the free
// surface implementation should be refactored such that the general generated stream() is applicable.
void {{class_name}}::Sweep::streamInCellInterval( {{stream_kernel|generate_field_type()}} * const pdfs,
                                                  {{stream_kernel|generate_field_type()}} * pdfs_tmp,
                                                  const CellInterval & ci )
{
    {{stream_kernel|generate_call(ghost_layers_to_include=0, cell_interval="ci", stream="stream")|indent(4)}}
}



} // namespace {{namespace}}
} // namespace walberla




// Buffer Packing

namespace walberla {
namespace mpi {

mpi::SendBuffer & operator<< (mpi::SendBuffer & buf, const ::walberla::{{namespace}}::{{class_name}} & lm)
{
    buf << lm.currentLevel;
    return buf;
}

mpi::RecvBuffer & operator>> (mpi::RecvBuffer & buf, ::walberla::{{namespace}}::{{class_name}} & lm)
{
    buf >> lm.currentLevel;
    return buf;
}


} // namespace mpi
} // namespace walberla

#if ( defined WALBERLA_CXX_COMPILER_IS_GNU ) || ( defined WALBERLA_CXX_COMPILER_IS_CLANG )
#   pragma GCC diagnostic pop
#endif