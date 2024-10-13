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

#include "blockforest/StructuredBlockForest.h"
#include "blockforest/BlockID.h"
#include "blockforest/Block.h"

#include "core/DataTypes.h"
#include "core/logging/Logging.h"
#include "core/Macros.h"

#include "field/AddToStorage.h"
#include "field/FlagField.h"
#include "field/iterators/FieldIterator.h"

{% if target is equalto 'gpu' -%}
#include "gpu/AddGPUFieldToStorage.h"
#include "gpu/GPUField.h"
#include "gpu/ParallelStreams.h"
{%- endif %}

#include "domain_decomposition/BlockDataID.h"
#include "domain_decomposition/IBlock.h"
#include "domain_decomposition/StructuredBlockStorage.h"

#include "field/SwapableCompare.h"
#include "field/GhostLayerField.h"

#include "stencil/Directions.h"
#include "stencil/{{stencil_name}}.h"

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
   using sizeTuple = std::tuple<int64_t, int64_t, int64_t>;

   {{class_name}}(const shared_ptr< StructuredBlockForest > & blocks, {{kernel_list|generate_constructor_parameters}}, const Cell & outerWidth=Cell(1, 1, 1))
      : blocks_(blocks), {{ kernel_list|generate_constructor_initializer_list(parameter_registration=parameter_scaling) }}, outerWidth_(outerWidth)
   {

      {{kernel_list|generate_constructor(parameter_registration=parameter_scaling) |indent(6)}}
      validInnerOuterSplit_ = true;

      for (auto& iBlock : *blocks)
      {
         if (int_c(blocks->getNumberOfXCells(iBlock)) <= outerWidth_[0] * 2 ||
             int_c(blocks->getNumberOfYCells(iBlock)) <= outerWidth_[1] * 2 ||
             int_c(blocks->getNumberOfZCells(iBlock)) <= outerWidth_[2] * 2)
            validInnerOuterSplit_ = false;
      }
   }

   void initialiseBlockPointer()
   {
      {%if block_stream_collide -%}
      blockWise_ = true;

      size_0.resize(blocks_->getNumberOfLevels());
      {%- for field in block_stream_collide['all_fields'] %}
      {{field.name}}Pointers.resize(blocks_->getNumberOfLevels());
      {%if target is equalto 'gpu' -%} {{field.name}}PointersGPU.resize(blocks_->getNumberOfLevels()); {% endif %}
      {%- endfor %}

      {%if target is equalto 'gpu' -%} block_.resize(blocks_->getNumberOfLevels()); {% endif %}
      {%if target is equalto 'gpu' -%} grid_.resize(blocks_->getNumberOfLevels()); {% endif %}

      glPropagationPDFs.resize(blocks_->getNumberOfLevels());
      {%if target is equalto 'gpu' -%} glPropagationPDFsGPU.resize(blocks_->getNumberOfLevels()); {% endif %}

      {%if target is equalto 'gpu' -%} glPropagationBlock_.resize(blocks_->getNumberOfLevels()); {% endif %}
      {%if target is equalto 'gpu' -%} glPropagationGrid_.resize(blocks_->getNumberOfLevels()); {% endif %}

      for( auto it = blocks_->begin(); it != blocks_->end(); ++it )
      {
         auto* local = dynamic_cast< Block* >(it.get());
         {%- for field in block_stream_collide['all_fields'] %}
         auto {{field.name}} = local->getData< {{field | field_type(is_gpu=is_gpu)}} >({{field.name}}ID);
         {%- endfor %}

         size_1 = int64_c({{block_stream_collide['all_fields'][0].name}}->xSize());
         size_2 = int64_c({{block_stream_collide['all_fields'][0].name}}->ySize());
         size_3 = int64_c({{block_stream_collide['all_fields'][0].name}}->zSize());

         {%- for field in block_stream_collide['all_fields'] %}

         stride_{{field.name}}_0 = int64_c({{field.name}}->xStride());
         stride_{{field.name}}_1 = int64_c({{field.name}}->yStride());
         stride_{{field.name}}_2 = int64_c({{field.name}}->zStride());
         stride_{{field.name}}_3 = int64_c(1 * int64_c({{field.name}}->fStride()));

         {%- endfor %}
         break;
      }

      for( auto it = blocks_->begin(); it != blocks_->end(); ++it )
      {
         auto* local = dynamic_cast< Block* >(it.get());
         const uint_t level = local->getLevel();
         {%- for field in block_stream_collide['all_fields'] %}
         auto {{field.name}} = local->getData< {{field | field_type(is_gpu=is_gpu)}} >({{field.name}}ID);
         {{field.name}}Pointers[level].emplace_back({{field.name}}->dataAt(0, 0, 0, 0));
         {%- endfor %}


         for(auto dir = stencil::{{stencil_name}}::beginNoCenter(); dir != stencil::{{stencil_name}}::end(); ++dir){
            uint_t nSecIdx = blockforest::getBlockNeighborhoodSectionIndex(*dir);
            // Propagate on ghost layers shadowing coarse or no blocks
            if(local->neighborhoodSectionHasLargerBlock(nSecIdx)){
               CellInterval ci;
               {{pdf_field.name}}->getGhostRegion(*dir, ci, 1);
               sizeTuple dirTuple = std::make_tuple(int64_c(ci.xSize()), int64_c(ci.ySize()), int64_c(ci.zSize()));
               glPropagationPDFs[level][dirTuple].emplace_back({{pdf_field.name}}->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0));
            }
         }
      }

      for (uint_t level = 0; level < blocks_->getNumberOfLevels(); level++) {
         size_0[level] = int64_c({{pdf_field.name}}Pointers[level].size());

         {%if target is equalto 'gpu' -%}

         int64_t indexingX = size_1 * size_0[level];
         int64_t indexingY = size_2;
         int64_t indexingZ = size_3;

         int64_t cudaBlockSize0 = 128;
         int64_t cudaBlockSize1 = 1;
         int64_t cudaBlockSize2 = 1;

         block_[level] = dim3((unsigned int)((1024 < ((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)) ? 1024 : ((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)), (unsigned int)((1024 < ((indexingY < cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)))) ? indexingY : cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0))))) ? 1024 : ((indexingY < cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)))) ? indexingY : cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0))))), (unsigned int)((64 < ((indexingZ < cudaBlockSize2*((int64_t)(cudaBlockSize0*cudaBlockSize1) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)*((indexingY < cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)))) ? indexingY : cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0))))))) ? indexingZ : cudaBlockSize2*((int64_t)(cudaBlockSize0*cudaBlockSize1) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)*((indexingY < cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)))) ? indexingY : cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)))))))) ? 64 : ((indexingZ < cudaBlockSize2*((int64_t)(cudaBlockSize0*cudaBlockSize1) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)*((indexingY < cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)))) ? indexingY : cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0))))))) ? indexingZ : cudaBlockSize2*((int64_t)(cudaBlockSize0*cudaBlockSize1) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)*((indexingY < cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)))) ? indexingY : cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)))))))));
         grid_[level]  = dim3((unsigned int)(( (indexingX) % (((1024 < ((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)) ? 1024 : ((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0))) == 0 ? (int64_t)(indexingX) / (int64_t)(((1024 < ((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)) ? 1024 : ((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0))) : ( (int64_t)(indexingX) / (int64_t)(((1024 < ((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)) ? 1024 : ((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0))) ) +1 )), (unsigned int)(( (indexingY) % (((1024 < ((indexingY < cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)))) ? indexingY : cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0))))) ? 1024 : ((indexingY < cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)))) ? indexingY : cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)))))) == 0 ? (int64_t)(indexingY) / (int64_t)(((1024 < ((indexingY < cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)))) ? indexingY : cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0))))) ? 1024 : ((indexingY < cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)))) ? indexingY : cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)))))) : ( (int64_t)(indexingY) / (int64_t)(((1024 < ((indexingY < cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)))) ? indexingY : cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0))))) ? 1024 : ((indexingY < cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)))) ? indexingY : cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)))))) ) +1 )), (unsigned int)(( (indexingZ) % (((64 < ((indexingZ < cudaBlockSize2*((int64_t)(cudaBlockSize0*cudaBlockSize1) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)*((indexingY < cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)))) ? indexingY : cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0))))))) ? indexingZ : cudaBlockSize2*((int64_t)(cudaBlockSize0*cudaBlockSize1) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)*((indexingY < cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)))) ? indexingY : cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)))))))) ? 64 : ((indexingZ < cudaBlockSize2*((int64_t)(cudaBlockSize0*cudaBlockSize1) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)*((indexingY < cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)))) ? indexingY : cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0))))))) ? indexingZ : cudaBlockSize2*((int64_t)(cudaBlockSize0*cudaBlockSize1) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)*((indexingY < cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)))) ? indexingY : cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0))))))))) == 0 ? (int64_t)(indexingZ) / (int64_t)(((64 < ((indexingZ < cudaBlockSize2*((int64_t)(cudaBlockSize0*cudaBlockSize1) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)*((indexingY < cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)))) ? indexingY : cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0))))))) ? indexingZ : cudaBlockSize2*((int64_t)(cudaBlockSize0*cudaBlockSize1) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)*((indexingY < cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)))) ? indexingY : cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)))))))) ? 64 : ((indexingZ < cudaBlockSize2*((int64_t)(cudaBlockSize0*cudaBlockSize1) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)*((indexingY < cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)))) ? indexingY : cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0))))))) ? indexingZ : cudaBlockSize2*((int64_t)(cudaBlockSize0*cudaBlockSize1) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)*((indexingY < cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)))) ? indexingY : cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0))))))))) : ( (int64_t)(indexingZ) / (int64_t)(((64 < ((indexingZ < cudaBlockSize2*((int64_t)(cudaBlockSize0*cudaBlockSize1) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)*((indexingY < cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)))) ? indexingY : cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0))))))) ? indexingZ : cudaBlockSize2*((int64_t)(cudaBlockSize0*cudaBlockSize1) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)*((indexingY < cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)))) ? indexingY : cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)))))))) ? 64 : ((indexingZ < cudaBlockSize2*((int64_t)(cudaBlockSize0*cudaBlockSize1) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)*((indexingY < cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)))) ? indexingY : cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0))))))) ? indexingZ : cudaBlockSize2*((int64_t)(cudaBlockSize0*cudaBlockSize1) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)*((indexingY < cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)))) ? indexingY : cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0))))))))) ) +1 )));

         {%- for field in block_stream_collide['all_fields'] %}

         WALBERLA_GPU_CHECK(gpuMalloc( (void**)&{{field.name}}PointersGPU[level], sizeof({{pdf_field.dtype.c_name}}* ) * {{field.name}}Pointers[level].size() ));
         WALBERLA_GPU_CHECK(gpuMemcpy( {{field.name}}PointersGPU[level], &{{field.name}}Pointers[level][0], sizeof({{pdf_field.dtype.c_name}} *) * {{field.name}}Pointers[level].size(), gpuMemcpyHostToDevice ));

         {%- endfor %}

         for (auto it = glPropagationPDFs[level].begin(); it != glPropagationPDFs[level].end(); it++){
            if(it->second.empty()){ continue;}

            indexingX = std::get<0>(it->first) * int64_c(it->second.size());
            indexingY = std::get<1>(it->first);
            indexingZ = std::get<2>(it->first);

            cudaBlockSize0 = 32;
            cudaBlockSize1 = 1;
            cudaBlockSize2 = 1;

            glPropagationBlock_[level][it->first] = dim3((unsigned int)((1024 < ((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)) ? 1024 : ((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)), (unsigned int)((1024 < ((indexingY < cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)))) ? indexingY : cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0))))) ? 1024 : ((indexingY < cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)))) ? indexingY : cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0))))), (unsigned int)((64 < ((indexingZ < cudaBlockSize2*((int64_t)(cudaBlockSize0*cudaBlockSize1) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)*((indexingY < cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)))) ? indexingY : cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0))))))) ? indexingZ : cudaBlockSize2*((int64_t)(cudaBlockSize0*cudaBlockSize1) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)*((indexingY < cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)))) ? indexingY : cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)))))))) ? 64 : ((indexingZ < cudaBlockSize2*((int64_t)(cudaBlockSize0*cudaBlockSize1) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)*((indexingY < cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)))) ? indexingY : cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0))))))) ? indexingZ : cudaBlockSize2*((int64_t)(cudaBlockSize0*cudaBlockSize1) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)*((indexingY < cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)))) ? indexingY : cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)))))))));
            glPropagationGrid_[level][it->first]  = dim3((unsigned int)(( (indexingX) % (((1024 < ((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)) ? 1024 : ((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0))) == 0 ? (int64_t)(indexingX) / (int64_t)(((1024 < ((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)) ? 1024 : ((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0))) : ( (int64_t)(indexingX) / (int64_t)(((1024 < ((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)) ? 1024 : ((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0))) ) +1 )), (unsigned int)(( (indexingY) % (((1024 < ((indexingY < cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)))) ? indexingY : cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0))))) ? 1024 : ((indexingY < cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)))) ? indexingY : cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)))))) == 0 ? (int64_t)(indexingY) / (int64_t)(((1024 < ((indexingY < cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)))) ? indexingY : cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0))))) ? 1024 : ((indexingY < cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)))) ? indexingY : cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)))))) : ( (int64_t)(indexingY) / (int64_t)(((1024 < ((indexingY < cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)))) ? indexingY : cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0))))) ? 1024 : ((indexingY < cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)))) ? indexingY : cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)))))) ) +1 )), (unsigned int)(( (indexingZ) % (((64 < ((indexingZ < cudaBlockSize2*((int64_t)(cudaBlockSize0*cudaBlockSize1) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)*((indexingY < cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)))) ? indexingY : cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0))))))) ? indexingZ : cudaBlockSize2*((int64_t)(cudaBlockSize0*cudaBlockSize1) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)*((indexingY < cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)))) ? indexingY : cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)))))))) ? 64 : ((indexingZ < cudaBlockSize2*((int64_t)(cudaBlockSize0*cudaBlockSize1) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)*((indexingY < cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)))) ? indexingY : cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0))))))) ? indexingZ : cudaBlockSize2*((int64_t)(cudaBlockSize0*cudaBlockSize1) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)*((indexingY < cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)))) ? indexingY : cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0))))))))) == 0 ? (int64_t)(indexingZ) / (int64_t)(((64 < ((indexingZ < cudaBlockSize2*((int64_t)(cudaBlockSize0*cudaBlockSize1) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)*((indexingY < cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)))) ? indexingY : cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0))))))) ? indexingZ : cudaBlockSize2*((int64_t)(cudaBlockSize0*cudaBlockSize1) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)*((indexingY < cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)))) ? indexingY : cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)))))))) ? 64 : ((indexingZ < cudaBlockSize2*((int64_t)(cudaBlockSize0*cudaBlockSize1) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)*((indexingY < cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)))) ? indexingY : cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0))))))) ? indexingZ : cudaBlockSize2*((int64_t)(cudaBlockSize0*cudaBlockSize1) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)*((indexingY < cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)))) ? indexingY : cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0))))))))) : ( (int64_t)(indexingZ) / (int64_t)(((64 < ((indexingZ < cudaBlockSize2*((int64_t)(cudaBlockSize0*cudaBlockSize1) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)*((indexingY < cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)))) ? indexingY : cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0))))))) ? indexingZ : cudaBlockSize2*((int64_t)(cudaBlockSize0*cudaBlockSize1) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)*((indexingY < cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)))) ? indexingY : cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)))))))) ? 64 : ((indexingZ < cudaBlockSize2*((int64_t)(cudaBlockSize0*cudaBlockSize1) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)*((indexingY < cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)))) ? indexingY : cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0))))))) ? indexingZ : cudaBlockSize2*((int64_t)(cudaBlockSize0*cudaBlockSize1) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)*((indexingY < cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0)))) ? indexingY : cudaBlockSize1*((int64_t)(cudaBlockSize0) / (int64_t)(((indexingX < cudaBlockSize0) ? indexingX : cudaBlockSize0))))))))) ) +1 )));

            WALBERLA_GPU_CHECK(gpuMalloc( (void**)&glPropagationPDFsGPU[level][it->first], sizeof({{pdf_field.dtype.c_name}}* ) * it->second.size() ))
            WALBERLA_GPU_CHECK(gpuMemcpy( glPropagationPDFsGPU[level][it->first], &it->second[0], sizeof({{pdf_field.dtype.c_name}}* ) * it->second.size(), gpuMemcpyHostToDevice ))
         }

         {% endif %}
      }
      {% endif %}

   };

   {%if block_stream_collide -%}
   ~{{class_name}}() {
      {%if target is equalto 'gpu' -%}
      for (uint_t level = 0; level < blocks_->getNumberOfLevels(); level++)
      {
         {%- for field in block_stream_collide['all_fields'] %}
         if(!{{field.name}}Pointers[level].empty()){
            WALBERLA_GPU_CHECK(gpuFree({{field.name}}PointersGPU[level]))
         }
         {%- endfor %}

         for (auto it = glPropagationPDFs[level].begin(); it != glPropagationPDFs[level].end(); it++){
            if(it->second.empty()){ continue;}
            WALBERLA_GPU_CHECK(gpuFree(glPropagationPDFsGPU[level][it->first]))
         }
      }
      {%- endif %}
   }
   {% else %}
   {{ kernel_list| generate_destructor(class_name) |indent(4) }}
   {% endif %}


   /*************************************************************************************
   *                Internal Function Definitions with raw Pointer
   *************************************************************************************/

   void blockStreamCollide({{- ["uint_t level", "uint8_t timestep", ["gpuStream_t stream = nullptr"] if target == 'gpu' else []] | type_identifier_list -}});
   void ghostLayerPropagation({{- ["uint_t level", "uint8_t timestep", ["gpuStream_t stream = nullptr"] if target == 'gpu' else []] | type_identifier_list -}});
   bool blockWise() {return blockWise_;};

   {%- for kernel in kernels %}
   static void {{kernel['function_name']}} ({{kernel['kernel']|generate_plain_parameter_list(ghost_layers=0, stream="nullptr")}});
   static void {{kernel['function_name']}}CellInterval ({{kernel['kernel']|generate_plain_parameter_list(cell_interval='ci', stream="nullptr")}});
   {% endfor %}

   /*************************************************************************************
   *                Function Definitions for external Usage
   *************************************************************************************/

   std::function<void ()> blockStreamCollideFck({{- ["uint_t level", "uint8_t timestep", ["gpuStream_t stream = nullptr"] if target == 'gpu' else []] | type_identifier_list -}} )
   {
      return [{{- ["this", "level", "timestep", ["stream"] if target == 'gpu' else []] | type_identifier_list -}}](){
         blockStreamCollide({{- ["level", "timestep", ["stream"] if target == 'gpu' else []] | type_identifier_list -}});
      };
   }

   void streamCollideOverBlocks({{- ["uint_t level", "uint8_t timestep", ["gpuStream_t stream = nullptr"] if target == 'gpu' else []] | type_identifier_list -}} )
   {
      blockStreamCollide({{- ["level", "timestep", ["stream"] if target == 'gpu' else []] | type_identifier_list -}});
   }


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

   {{kernel_list|generate_getter(parameter_registration=parameter_scaling)|indent(3)}}
   {{kernel_list|generate_setter(parameter_registration=parameter_scaling)|indent(3)}}

 private:
   shared_ptr< StructuredBlockForest > blocks_;
   {{kernel_list|generate_members(parameter_registration=parameter_scaling)|indent(4)}}

   Cell outerWidth_;
   std::vector<CellInterval> layers_;
   bool validInnerOuterSplit_;
   bool blockWise_{false};

   {%if target is equalto 'gpu' -%}
   gpu::ParallelStreams parallelStreams_;
   {%- endif %}

   {%if block_stream_collide -%}

   std::vector<int64_t> size_0;

   int64_t size_1;
   int64_t size_2;
   int64_t size_3;

   {%- for field in block_stream_collide['all_fields'] %}
   int64_t stride_{{field.name}}_0;
   int64_t stride_{{field.name}}_1;
   int64_t stride_{{field.name}}_2;
   int64_t stride_{{field.name}}_3;
   {% endfor %}

   {%- for field in block_stream_collide['all_fields'] %}
   std::vector<std::vector<{{field.dtype.c_name}} *>> {{field.name}}Pointers;
   {% endfor -%}

   std::vector<std::map<sizeTuple, std::vector<{{pdf_field.dtype.c_name}} *>>> glPropagationPDFs;
   {%if target is equalto 'gpu' -%}


   {%- for field in block_stream_collide['all_fields'] %}
   std::vector<{{field.dtype.c_name}} **> {{field.name}}PointersGPU;
   {% endfor -%}

   std::vector<dim3> block_;
   std::vector<dim3> grid_;

   std::vector<std::map<sizeTuple, {{pdf_field.dtype.c_name}} **>> glPropagationPDFsGPU;

   std::vector<std::map<sizeTuple, dim3>> glPropagationBlock_;
   std::vector<std::map<sizeTuple, dim3>> glPropagationGrid_;
   {%- endif %}
   {%- endif %}
};


} // namespace {{namespace}}
} // namespace walberla


#if ( defined WALBERLA_CXX_COMPILER_IS_GNU ) || ( defined WALBERLA_CXX_COMPILER_IS_CLANG )
#   pragma GCC diagnostic pop
#endif
