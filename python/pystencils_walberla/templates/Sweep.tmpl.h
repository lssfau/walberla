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
//! \\author pystencils
//======================================================================================================================

#pragma once
#include "core/DataTypes.h"

{% if target is equalto 'cpu' -%}
#include "field/GhostLayerField.h"
{%- elif target is equalto 'gpu' -%}
#include "cuda/GPUField.h"
{% if inner_outer_split -%}
#include "cuda/ParallelStreams.h"
{%- endif %}
{%- endif %}
#include "field/SwapableCompare.h"
#include "domain_decomposition/BlockDataID.h"
#include "domain_decomposition/IBlock.h"
#include "domain_decomposition/StructuredBlockStorage.h"
#include <set>

#ifdef __GNUC__
#define RESTRICT __restrict__
#elif _MSC_VER
#define RESTRICT __restrict
#else
#define RESTRICT
#endif

#if ( defined WALBERLA_CXX_COMPILER_IS_GNU ) || ( defined WALBERLA_CXX_COMPILER_IS_CLANG )
#   pragma GCC diagnostic push
#   pragma GCC diagnostic ignored "-Wunused-parameter"
#   pragma GCC diagnostic ignored "-Wreorder"
#endif

namespace walberla {
namespace {{namespace}} {


class {{class_name}}
{
public:
    {{class_name}}( {{kernel|generate_constructor_parameters}} {%if inner_outer_split%}, const Cell & outerWidth=Cell(1, 1, 1){% endif %})
        : {{ kernel|generate_constructor_initializer_list }}{%if inner_outer_split%}{% if kernel|generate_constructor_initializer_list|length %},{% endif %} outerWidth_(outerWidth){% endif %}
    {};

    {{ kernel| generate_destructor(class_name) |indent(4) }}

    void operator() ( IBlock * block{%if target is equalto 'gpu'%} , cudaStream_t stream = nullptr{% endif %} );
    void runOnCellInterval(const shared_ptr<StructuredBlockStorage> & blocks,
                           const CellInterval & globalCellInterval, cell_idx_t ghostLayers, IBlock * block
                           {%if target is equalto 'gpu'%} , cudaStream_t stream = nullptr{% endif %});



    static std::function<void (IBlock*)> getSweep(const shared_ptr<{{class_name}}> & kernel) {
        return [kernel](IBlock * b) { (*kernel)(b); };
    }

    static std::function<void (IBlock*{%if target is equalto 'gpu'%} , cudaStream_t {% endif %})>
            getSweepOnCellInterval(const shared_ptr<{{class_name}}> & kernel,
                                   const shared_ptr<StructuredBlockStorage> & blocks,
                                   const CellInterval & globalCellInterval,
                                   cell_idx_t ghostLayers=1 )
    {
        return [kernel, blocks, globalCellInterval, ghostLayers] (IBlock * b{%if target is equalto 'gpu'%} , cudaStream_t stream = nullptr{% endif %}) {
            kernel->runOnCellInterval(blocks, globalCellInterval, ghostLayers, b{%if target is equalto 'gpu'%}, stream {% endif %});
        };
    }

{% if inner_outer_split %}

    void inner( IBlock * block{%if target is equalto 'gpu'%} , cudaStream_t stream = nullptr{% endif %} );
    void outer( IBlock * block{%if target is equalto 'gpu'%} , cudaStream_t stream = nullptr{% endif %} );

    void setOuterPriority(int priority ) {
       {%if target is equalto 'gpu'%}
       parallelStreams_.setStreamPriority(priority);
       {%endif%}
    }
    {{kernel|generate_members|indent(4)}}

private:
    {%if target is equalto 'gpu'%}
    cuda::ParallelStreams parallelStreams_;
    {% endif %}

    Cell outerWidth_;
    std::vector<CellInterval> layers_;

{%- else -%}
    {{ kernel|generate_members|indent(4) }}
{% endif -%}

};


} // namespace {{namespace}}
} // namespace walberla


#if ( defined WALBERLA_CXX_COMPILER_IS_GNU ) || ( defined WALBERLA_CXX_COMPILER_IS_CLANG )
#   pragma GCC diagnostic pop
#endif
