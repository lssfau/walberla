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
#include "core/logging/Logging.h"
#include "core/Macros.h"

{% if target is equalto 'gpu' -%}
#include "gpu/GPUField.h"
#include "gpu/ParallelStreams.h"
{%- endif %}

#include "domain_decomposition/BlockDataID.h"
#include "domain_decomposition/IBlock.h"
#include "domain_decomposition/StructuredBlockStorage.h"

#include "field/SwapableCompare.h"
#include "field/GhostLayerField.h"

#include <set>
#include <cmath>

{% for header in headers %}
#include {{header}}
{% endfor %}

using namespace std::placeholders;

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
   enum Type { ALL = 0, INNER = 1, OUTER = 2 };

   {{class_name}}(const shared_ptr< StructuredBlockStorage > & blocks, {{kernel_list|generate_constructor_parameters}}, const Cell & outerWidth=Cell(1, 1, 1))
      : blocks_(blocks), {{ kernel_list|generate_constructor_initializer_list(parameter_registration=parameter_scaling) }}, outerWidth_(outerWidth)
   {
      {{kernel_list|generate_constructor(parameter_registration=parameter_scaling) |indent(6)}}

      validInnerOuterSplit_= true;

      for (auto& iBlock : *blocks)
      {
         if (int_c(blocks->getNumberOfXCells(iBlock)) <= outerWidth_[0] * 2 || int_c(blocks->getNumberOfYCells(iBlock)) <= outerWidth_[1] * 2 || int_c(blocks->getNumberOfZCells(iBlock)) <= outerWidth_[2] * 2)
            validInnerOuterSplit_ = false;
      }
   };

   {{ kernel_list| generate_destructor(class_name) |indent(4) }}

   /*************************************************************************************
   *                Internal Function Definitions with raw Pointer
   *************************************************************************************/

   {%- for kernel in kernels %}
   static void {{kernel['function_name']}} ({{kernel['kernel']|generate_plain_parameter_list(ghost_layers=0, stream="nullptr")}});
   static void {{kernel['function_name']}}CellInterval ({{kernel['kernel']|generate_plain_parameter_list(cell_interval='ci', stream="nullptr")}});
   {% endfor %}

   /*************************************************************************************
   *                Function Definitions for external Usage
   *************************************************************************************/

   {%- for kernel in kernels %}

   std::function<void (IBlock *)> {{kernel['function_name']}}()
   {
      return [{{- ["this", ] | type_identifier_list -}}](IBlock* block) { {{kernel['function_name']}}({{- ["block", ] | type_identifier_list -}}); };
   }

   std::function<void (IBlock *)> {{kernel['function_name']}}({{- ["Type type", ] | type_identifier_list -}})
   {
      if (!validInnerOuterSplit_ && type != Type::ALL)
         WALBERLA_ABORT_NO_DEBUG_INFO("innerOuterSplit too large - make it smaller, increase cellsPerBlock or avoid communication hiding")

      switch (type)
      {
      case Type::INNER:
         return [{{- ["this", ] | type_identifier_list -}}](IBlock* block) { {{kernel['function_name']}}Inner({{- ["block", ] | type_identifier_list -}}); };
      case Type::OUTER:
         return [{{- ["this", ] | type_identifier_list -}}](IBlock* block) { {{kernel['function_name']}}Outer({{- ["block", ] | type_identifier_list -}}); };
      default:
         return [{{- ["this", ] | type_identifier_list -}}](IBlock* block) { {{kernel['function_name']}}({{- ["block", ] | type_identifier_list -}}); };
      }
   }

   std::function<void (IBlock *)> {{kernel['function_name']}}({{- ["Type type", "const cell_idx_t ghost_layers"] | type_identifier_list -}})
   {
      if (!validInnerOuterSplit_ && type != Type::ALL)
         WALBERLA_ABORT_NO_DEBUG_INFO("innerOuterSplit too large - make it smaller, increase cellsPerBlock or avoid communication hiding")

      switch (type)
      {
      case Type::INNER:
         return [{{- ["this", ] | type_identifier_list -}}](IBlock* block) { {{kernel['function_name']}}Inner({{- ["block", ] | type_identifier_list -}}); };
      case Type::OUTER:
         return [{{- ["this", ] | type_identifier_list -}}](IBlock* block) { {{kernel['function_name']}}Outer({{- ["block", ] | type_identifier_list -}}); };
      default:
         return [{{- ["this", "ghost_layers"] | type_identifier_list -}}](IBlock* block) { {{kernel['function_name']}}({{- ["block", "ghost_layers"] | type_identifier_list -}}); };
      }
   }

   {% if target is equalto 'gpu' -%}
   std::function<void (IBlock *)> {{kernel['function_name']}}({{- ["Type type", "const cell_idx_t ghost_layers", "gpuStream_t gpuStream"] | type_identifier_list -}})
   {
      if (!validInnerOuterSplit_ && type != Type::ALL)
         WALBERLA_ABORT_NO_DEBUG_INFO("innerOuterSplit too large - make it smaller, increase cellsPerBlock or avoid communication hiding")

      switch (type)
      {
      case Type::INNER:
         return [{{- ["this", "gpuStream"] | type_identifier_list -}}](IBlock* block) { {{kernel['function_name']}}Inner({{- ["block", "gpuStream"] | type_identifier_list -}}); };
      case Type::OUTER:
         return [{{- ["this", "gpuStream"] | type_identifier_list -}}](IBlock* block) { {{kernel['function_name']}}Outer({{- ["block", "gpuStream"] | type_identifier_list -}}); };
      default:
         return [{{- ["this", "ghost_layers", "gpuStream"] | type_identifier_list -}}](IBlock* block) { {{kernel['function_name']}}({{- ["block", "ghost_layers", "gpuStream"] | type_identifier_list -}}); };
      }
   }

   std::function<void (IBlock *)> {{kernel['function_name']}}({{- ["Type type", "gpuStream_t gpuStream"] | type_identifier_list -}})
   {
      if (!validInnerOuterSplit_ && type != Type::ALL)
         WALBERLA_ABORT_NO_DEBUG_INFO("innerOuterSplit too large - make it smaller, increase cellsPerBlock or avoid communication hiding")

      switch (type)
      {
      case Type::INNER:
         return [{{- ["this", "gpuStream"] | type_identifier_list -}}](IBlock* block) { {{kernel['function_name']}}Inner({{- ["block", "gpuStream"] | type_identifier_list -}}); };
      case Type::OUTER:
         return [{{- ["this", "gpuStream"] | type_identifier_list -}}](IBlock* block) { {{kernel['function_name']}}Outer({{- ["block", "gpuStream"] | type_identifier_list -}}); };
      default:
         return [{{- ["this", "gpuStream"] | type_identifier_list -}}](IBlock* block) { {{kernel['function_name']}}({{- ["block", "cell_idx_c(0)", "gpuStream"] | type_identifier_list -}}); };
      }
   }
   {%- endif %}

   void {{kernel['function_name']}}({{- ["IBlock * block",] | type_identifier_list -}})
   {
      const cell_idx_t ghost_layers = 0;
      {% if target is equalto 'gpu' -%}
      gpuStream_t gpuStream = nullptr;
      {%- endif %}

      {{kernel['kernel']|generate_block_data_to_field_extraction|indent(6)}}
      {{kernel['kernel']|generate_refs_for_kernel_parameters(prefix='this->', ignore_fields=True, parameter_registration=parameter_scaling)|indent(6)}}
      {{kernel['kernel']|generate_timestep_advancements|indent(6)}}
      {{kernel['function_name']}}({{kernel['kernel']|generate_function_collection_call(ghost_layers='ghost_layers')}});
      {{kernel['kernel']|generate_swaps|indent(6)}}
   }

   void {{kernel['function_name']}}({{- ["IBlock * block", "const cell_idx_t ghost_layers"] | type_identifier_list -}})
   {
      {% if target is equalto 'gpu' -%}
      gpuStream_t gpuStream = nullptr;
      {%- endif %}

      {{kernel['kernel']|generate_block_data_to_field_extraction|indent(6)}}
      {{kernel['kernel']|generate_refs_for_kernel_parameters(prefix='this->', ignore_fields=True, parameter_registration=parameter_scaling)|indent(6)}}
      {{kernel['kernel']|generate_timestep_advancements|indent(6)}}
      {{kernel['function_name']}}({{kernel['kernel']|generate_function_collection_call(ghost_layers='ghost_layers')}});
      {{kernel['kernel']|generate_swaps|indent(6)}}
   }

   {% if target is equalto 'gpu' -%}
   void {{kernel['function_name']}}({{- ["IBlock * block", "const cell_idx_t ghost_layers", "gpuStream_t gpuStream"] | type_identifier_list -}})
   {
      {{kernel['kernel']|generate_block_data_to_field_extraction|indent(6)}}
      {{kernel['kernel']|generate_refs_for_kernel_parameters(prefix='this->', ignore_fields=True, parameter_registration=parameter_scaling)|indent(6)}}
      {{kernel['kernel']|generate_timestep_advancements|indent(6)}}
      {{kernel['function_name']}}({{kernel['kernel']|generate_function_collection_call(ghost_layers='ghost_layers')}});
      {{kernel['kernel']|generate_swaps|indent(6)}}
   }
   {%- endif %}

   void {{kernel['function_name']}}CellInterval({{- ["IBlock * block", "const CellInterval & ci", ["gpuStream_t gpuStream = nullptr"] if target == 'gpu' else []] | type_identifier_list -}})
   {
      {{kernel['kernel']|generate_block_data_to_field_extraction|indent(6)}}
      {{kernel['kernel']|generate_refs_for_kernel_parameters(prefix='this->', ignore_fields=True, parameter_registration=parameter_scaling)|indent(6)}}
      {{kernel['kernel']|generate_timestep_advancements|indent(6)}}
      {{kernel['function_name']}}CellInterval({{kernel['kernel']|generate_function_collection_call(cell_interval='ci')}});
      {{kernel['kernel']|generate_swaps|indent(6)}}
   }

   void {{kernel['function_name']}}Inner({{- ["IBlock * block", ["gpuStream_t gpuStream = nullptr"] if target == 'gpu' else []] | type_identifier_list -}})
   {
      {{kernel['kernel']|generate_block_data_to_field_extraction|indent(6)}}
      {{kernel['kernel']|generate_refs_for_kernel_parameters(prefix='this->', ignore_fields=True, parameter_registration=parameter_scaling)|indent(6)}}
      {{kernel['kernel']|generate_timestep_advancements(advance=False)|indent(6)}}

      CellInterval inner = {{kernel['field']}}->xyzSize();
      inner.expand(Cell(-outerWidth_[0], -outerWidth_[1], -outerWidth_[2]));

      {{kernel['function_name']}}CellInterval({{kernel['kernel']|generate_function_collection_call(cell_interval='inner')}});
   }

   void {{kernel['function_name']}}Outer({{- ["IBlock * block", ["gpuStream_t gpuStream = nullptr"] if target == 'gpu' else []] | type_identifier_list -}})
   {

      {{kernel['kernel']|generate_block_data_to_field_extraction|indent(6)}}
      {{kernel['kernel']|generate_refs_for_kernel_parameters(prefix='this->', ignore_fields=True, parameter_registration=parameter_scaling)|indent(6)}}
      {{kernel['kernel']|generate_timestep_advancements|indent(6)}}

      if( layers_.empty() )
      {
         CellInterval ci;

         {{kernel['field']}}->getSliceBeforeGhostLayer(stencil::T, ci, outerWidth_[2], false);
         layers_.push_back(ci);
         {{kernel['field']}}->getSliceBeforeGhostLayer(stencil::B, ci, outerWidth_[2], false);
         layers_.push_back(ci);

         {{kernel['field']}}->getSliceBeforeGhostLayer(stencil::N, ci, outerWidth_[1], false);
         ci.expand(Cell(0, 0, -outerWidth_[2]));
         layers_.push_back(ci);
         {{kernel['field']}}->getSliceBeforeGhostLayer(stencil::S, ci, outerWidth_[1], false);
         ci.expand(Cell(0, 0, -outerWidth_[2]));
         layers_.push_back(ci);

         {{kernel['field']}}->getSliceBeforeGhostLayer(stencil::E, ci, outerWidth_[0], false);
         ci.expand(Cell(0, -outerWidth_[1], -outerWidth_[2]));
         layers_.push_back(ci);
         {{kernel['field']}}->getSliceBeforeGhostLayer(stencil::W, ci, outerWidth_[0], false);
         ci.expand(Cell(0, -outerWidth_[1], -outerWidth_[2]));
         layers_.push_back(ci);
      }

      {%if target is equalto 'gpu'%}
      {
         auto parallelSection_ = parallelStreams_.parallelSection( gpuStream );
         for( auto & ci: layers_ )
         {
            parallelSection_.run([&]( auto s ) {
               {{kernel['function_name']}}CellInterval({{kernel['kernel']|generate_function_collection_call(cell_interval='ci')}});
            });
         }
      }
      {% else %}
      for( auto & ci: layers_ )
      {
         {{kernel['function_name']}}CellInterval({{kernel['kernel']|generate_function_collection_call(cell_interval='ci')}});
      }
      {% endif %}

      {{kernel['kernel']|generate_swaps|indent(9)}}
   }
   {% endfor %}

   {%if target is equalto 'gpu'%}
   void setOuterPriority(int priority)
   {
      parallelStreams_.setStreamPriority(priority);
   }
   {%endif%}

 private:
   shared_ptr< StructuredBlockStorage > blocks_;
   {{kernel_list|generate_members(parameter_registration=parameter_scaling)|indent(4)}}

   Cell outerWidth_;
   std::vector<CellInterval> layers_;
   bool validInnerOuterSplit_;

   {%if target is equalto 'gpu' -%}
   gpu::ParallelStreams parallelStreams_;
   // std::map<BlockID, gpuStream_t > streams_;
   {%- endif %}
};


} // namespace {{namespace}}
} // namespace walberla


#if ( defined WALBERLA_CXX_COMPILER_IS_GNU ) || ( defined WALBERLA_CXX_COMPILER_IS_CLANG )
#   pragma GCC diagnostic pop
#endif
