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

{% if target is equalto 'cpu' -%}
#include "field/GhostLayerField.h"
{%- elif target is equalto 'gpu' -%}
#include "gpu/FieldCopy.h"
#include "gpu/GPUField.h"
#include "gpu/GPUWrapper.h"
{%- endif %}
#include "domain_decomposition/BlockDataID.h"
#include "domain_decomposition/IBlock.h"
#include "blockforest/StructuredBlockForest.h"
#include "field/FlagField.h"
#include "core/debug/Debug.h"

#include <set>
#include <vector>

{% for header in interface_spec.headers %}
#include {{header}}
{% endfor %}

#ifdef __GNUC__
#define RESTRICT __restrict__
#elif _MSC_VER
#define RESTRICT __restrict
#else
#define RESTRICT
#endif

#ifdef WALBERLA_BUILD_WITH_HALF_PRECISION_SUPPORT
using walberla::half;
#endif

namespace walberla {
namespace {{namespace}} {


class {{class_name}}
{
public:
    {{StructDeclaration|indent(4)}}


    class IndexVectors
    {
    public:
        using CpuIndexVector = std::vector<{{StructName}}>;

        enum Type {
            ALL = 0,
            INNER = 1,
            OUTER = 2,
            NUM_TYPES = 3
        };

        IndexVectors() = default;
        bool operator==(IndexVectors const &other) const { return other.cpuVectors_ == cpuVectors_; }

        {% if target == 'gpu' -%}
        ~IndexVectors() {
            for( auto & gpuVec: gpuVectors_)
               WALBERLA_GPU_CHECK(gpuFree( gpuVec ));
        }
        {% endif -%}

        CpuIndexVector & indexVector(Type t) { return cpuVectors_[t]; }
        {{StructName}} * pointerCpu(Type t)  { return cpuVectors_[t].data(); }

        {% if target == 'gpu' -%}
        {{StructName}} * pointerGpu(Type t)  { return gpuVectors_[t]; }
        {% endif -%}

        void syncGPU()
        {
            {% if target == 'gpu' -%}
            for( auto & gpuVec: gpuVectors_)
               WALBERLA_GPU_CHECK(gpuFree( gpuVec ));
            gpuVectors_.resize( cpuVectors_.size() );

            WALBERLA_ASSERT_EQUAL(cpuVectors_.size(), NUM_TYPES);
            for(size_t i=0; i < cpuVectors_.size(); ++i )
            {
                auto & gpuVec = gpuVectors_[i];
                auto & cpuVec = cpuVectors_[i];
                WALBERLA_GPU_CHECK(gpuMalloc( &gpuVec, sizeof({{StructName}}) * cpuVec.size() ));
                WALBERLA_GPU_CHECK(gpuMemcpy( gpuVec, &cpuVec[0], sizeof({{StructName}}) * cpuVec.size(), gpuMemcpyHostToDevice ));
            }
            {%- endif %}
        }




    private:
        std::vector<CpuIndexVector> cpuVectors_{NUM_TYPES};

        {% if target == 'gpu' -%}
        using GpuIndexVector = {{StructName}} *;
        std::vector<GpuIndexVector> gpuVectors_;
        {%- endif %}
    };

    {% if calculate_force -%}

    struct ForceStruct {
       double F_0;
       double F_1;
       double F_2;
       ForceStruct() : F_0(double_c(0.0)), F_1(double_c(0.0)), F_2(double_c(0.0)) {}
       bool operator==(const ForceStruct & o) const {
          return floatIsEqual(F_0, o.F_0) && floatIsEqual(F_1, o.F_1) && floatIsEqual(F_2, o.F_2);
       }
    };

    class ForceVector
    {
     public:
       ForceVector() = default;
       bool operator==(ForceVector const &other) const { return other.cpuVector_ == cpuVector_; }

       {% if target == 'gpu' -%}
       ~ForceVector() {if(!gpuVector_.empty()){WALBERLA_GPU_CHECK(gpuFree( gpuVector_[0] ))}}
       {% endif -%}

       std::vector<ForceStruct> & forceVector() { return cpuVector_; }
       ForceStruct * pointerCpu()  { return cpuVector_.data(); }
       bool empty() {return cpuVector_.empty();}

       {% if target == 'gpu' -%}
       ForceStruct * pointerGpu()  { return gpuVector_[0]; }
       {% endif -%}

       Vector3<double> getForce()
       {
          syncCPU();
          Vector3<double> result(double_c(0.0));
          for(std::vector<ForceStruct>::iterator it = cpuVector_.begin(); it != cpuVector_.end(); ++it)
          {
             result[0] += it->F_0;
             result[1] += it->F_1;
             result[2] += it->F_2;
          }
          return result;
       }

       void syncGPU()
       {
          {% if target == 'gpu' -%}
          if(!gpuVector_.empty()){WALBERLA_GPU_CHECK(gpuFree( gpuVector_[0] ))}
          if(!cpuVector_.empty())
          {
             gpuVector_.resize(cpuVector_.size());
             WALBERLA_GPU_CHECK(gpuMalloc(&gpuVector_[0], sizeof(ForceStruct) * cpuVector_.size()))
             WALBERLA_GPU_CHECK(gpuMemcpy(gpuVector_[0], &cpuVector_[0], sizeof(ForceStruct) * cpuVector_.size(), gpuMemcpyHostToDevice))
          }
          {%- endif %}
       }

       void syncCPU()
       {
          {% if target == 'gpu' -%}
          WALBERLA_GPU_CHECK(gpuMemcpy( &cpuVector_[0], gpuVector_[0] , sizeof(ForceStruct) * cpuVector_.size(), gpuMemcpyDeviceToHost ))
          {%- endif %}
       }

     private:
       std::vector<ForceStruct> cpuVector_;
       {% if target == 'gpu' -%}
       std::vector<ForceStruct *> gpuVector_;
       {%- endif %}
    };

    {%- endif %}

    {{class_name}}( const shared_ptr<StructuredBlockForest> & blocks,
                   {{kernel|generate_constructor_parameters(['indexVector', 'indexVectorSize', 'forceVector', 'forceVectorSize'])}}{{additional_data_handler.constructor_arguments}})
        :{{additional_data_handler.initialiser_list}} {{ kernel|generate_constructor_initializer_list(['indexVector', 'indexVectorSize', 'forceVector', 'forceVectorSize']) }}
    {
        auto createIdxVector = []( IBlock * const , StructuredBlockStorage * const ) { return new IndexVectors(); };
        indexVectorID = blocks->addStructuredBlockData< IndexVectors >( createIdxVector, "IndexField_{{class_name}}");
        {% if calculate_force -%}
        auto createForceVector = []( IBlock * const , StructuredBlockStorage * const ) { return new ForceVector(); };
        forceVectorID = blocks->addStructuredBlockData< ForceVector >( createForceVector, "forceVector_{{class_name}}");
        {%- endif %}
    }

    void run (
        {{- ["IBlock * block", kernel.kernel_selection_parameters, ["gpuStream_t stream = nullptr"] if target == 'gpu' else []] | type_identifier_list -}}
    );

    {% if generate_functor -%}
    void operator() (
        {{- ["IBlock * block", kernel.kernel_selection_parameters, ["gpuStream_t stream = nullptr"] if target == 'gpu' else []] | type_identifier_list -}}
    )
    {
        run( {{- ["block", kernel.kernel_selection_parameters, ["stream"] if target == 'gpu' else []] | identifier_list -}} );
    }
    {%- endif %}

    void inner (
        {{- ["IBlock * block", kernel.kernel_selection_parameters, ["gpuStream_t stream = nullptr"] if target == 'gpu' else []] | type_identifier_list -}}
    );

    void outer (
        {{- ["IBlock * block", kernel.kernel_selection_parameters, ["gpuStream_t stream = nullptr"] if target == 'gpu' else []] | type_identifier_list -}}
    );

    Vector3<double> getForce(IBlock * {% if calculate_force -%}block{%else%}/*block*/{%- endif %})
    {
       {% if calculate_force -%}
       auto * forceVector = block->getData<ForceVector>(forceVectorID);
       if(forceVector->empty())
          return Vector3<double>(double_c(0.0));
       return forceVector->getForce();
       {% else %}
       WALBERLA_ABORT("Boundary condition was not generated including force calculation.")
       return Vector3<double>(double_c(0.0));
       {%- endif %}
    }

    std::function<void (IBlock *)> getSweep( {{- [interface_spec.high_level_args, ["gpuStream_t stream = nullptr"] if target == 'gpu' else []] | type_identifier_list -}} )
    {
        return [ {{- ["this", interface_spec.high_level_args, ["stream"] if target == 'gpu' else []] | identifier_list -}} ]
               (IBlock * b)
               { this->run( {{- [ ['b'], interface_spec.mapping_codes, ["stream"] if target == 'gpu' else [] ] | identifier_list -}} ); };
    }

    std::function<void (IBlock *)> getInnerSweep( {{- [interface_spec.high_level_args, ["gpuStream_t stream = nullptr"] if target == 'gpu' else []] | type_identifier_list -}} )
    {
        return [ {{- [ ['this'], interface_spec.high_level_args, ["stream"] if target == 'gpu' else [] ] | identifier_list -}} ]
               (IBlock * b)
               { this->inner( {{- [ ['b'], interface_spec.mapping_codes, ["stream"] if target == 'gpu' else [] ] | identifier_list -}} ); };
    }

    std::function<void (IBlock *)> getOuterSweep( {{- [interface_spec.high_level_args, ["gpuStream_t stream = nullptr"] if target == 'gpu' else []] | type_identifier_list -}} )
    {
        return [ {{- [ ['this'], interface_spec.high_level_args, ["stream"] if target == 'gpu' else [] ] | identifier_list -}} ]
               (IBlock * b)
               { this->outer( {{- [ ['b'], interface_spec.mapping_codes, ["stream"] if target == 'gpu' else [] ] | identifier_list -}} ); };
    }

    template<typename FlagField_T>
    void fillFromFlagField( const shared_ptr<StructuredBlockForest> & blocks, ConstBlockDataID flagFieldID,
                            FlagUID boundaryFlagUID, FlagUID domainFlagUID)
    {
        for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt )
            fillFromFlagField<FlagField_T>({{additional_data_handler.additional_arguments_for_fill_function}}&*blockIt, flagFieldID, boundaryFlagUID, domainFlagUID );
    }


    template<typename FlagField_T>
    void fillFromFlagField({{additional_data_handler.additional_parameters_for_fill_function}}IBlock * block, ConstBlockDataID flagFieldID,
                            FlagUID boundaryFlagUID, FlagUID domainFlagUID )
    {
        auto * indexVectors = block->getData< IndexVectors > ( indexVectorID );
        auto & indexVectorAll = indexVectors->indexVector(IndexVectors::ALL);
        auto & indexVectorInner = indexVectors->indexVector(IndexVectors::INNER);
        auto & indexVectorOuter = indexVectors->indexVector(IndexVectors::OUTER);
        {% if calculate_force -%}
        auto * forceVector = block->getData< ForceVector > ( forceVectorID );
        {%- endif %}

        auto * flagField = block->getData< FlagField_T > ( flagFieldID );
        {{additional_data_handler.additional_field_data|indent(4)}}

        if( !(flagField->flagExists(boundaryFlagUID) && flagField->flagExists(domainFlagUID) ))
            return;

        auto boundaryFlag = flagField->getFlag(boundaryFlagUID);
        auto domainFlag = flagField->getFlag(domainFlagUID);

        auto inner = flagField->xyzSize();
        inner.expand( cell_idx_t(-1) );

        indexVectorAll.clear();
        indexVectorInner.clear();
        indexVectorOuter.clear();

        {% if inner_or_boundary -%}
        {% if layout == "fzyx" -%}
        {%- for dirIdx, dirVec, offset in additional_data_handler.stencil_info %}
        for( auto it = flagField->beginWithGhostLayerXYZ( cell_idx_c( flagField->nrOfGhostLayers() - 1 ) ); it != flagField->end(); ++it )
        {
           if( ! isFlagSet(it, domainFlag) )
              continue;

           if ( isFlagSet( it.neighbor({{offset}} {%if dim == 3%}, 0 {%endif %}), boundaryFlag ) )
           {
              auto element = {{StructName}}(it.x(), it.y(), {%if dim == 3%} it.z(), {%endif %} {{dirIdx}} );
              {{additional_data_handler.data_initialisation(dirIdx)|indent(16)}}
              indexVectorAll.push_back( element );
              if( inner.contains( it.x(), it.y(), it.z() ) )
                 indexVectorInner.push_back( element );
              else
                 indexVectorOuter.push_back( element );
           }
        }
        {% endfor %}
        {%else%}
        for( auto it = flagField->beginWithGhostLayerXYZ( cell_idx_c( flagField->nrOfGhostLayers() - 1 ) ); it != flagField->end(); ++it )
        {
            if( ! isFlagSet(it, domainFlag) )
                continue;
            {%- for dirIdx, dirVec, offset in additional_data_handler.stencil_info %}
            if ( isFlagSet( it.neighbor({{offset}} {%if dim == 3%}, 0 {%endif %}), boundaryFlag ) )
            {
                auto element = {{StructName}}(it.x(), it.y(), {%if dim == 3%} it.z(), {%endif %} {{dirIdx}} );
                {{additional_data_handler.data_initialisation(dirIdx)|indent(16)}}
                indexVectorAll.push_back( element );
                if( inner.contains( it.x(), it.y(), it.z() ) )
                    indexVectorInner.push_back( element );
                else
                    indexVectorOuter.push_back( element );
            }
            {% endfor %}
        }
        {%endif%}
        {%else%}
        auto flagWithGLayers = flagField->xyzSizeWithGhostLayer();
        {% if single_link %}
        {{dtype}} dot = 0.0; {{dtype}} maxn = 0.0;
        cell_idx_t calculated_idx = 0;
        cell_idx_t dx = 0; cell_idx_t dy = 0; {%if dim == 3%}  cell_idx_t dz = 0; {% endif %}
        cell_idx_t sum_x = 0; cell_idx_t sum_y = 0; {%if dim == 3%} cell_idx_t sum_z = 0; {%endif %}
        {% endif -%}
        for( auto it = flagField->beginWithGhostLayerXYZ(); it != flagField->end(); ++it )
        {
            {% if single_link -%}
            sum_x = 0; sum_y = 0; {%if dim == 3%} sum_z = 0; {%endif %}
            {% endif %}
            if( ! isFlagSet(it, boundaryFlag) )
                continue;
            {%- for dirIdx, dirVec, offset in additional_data_handler.stencil_info %}
            if ( flagWithGLayers.contains(it.x() + cell_idx_c({{dirVec[0]}}), it.y() + cell_idx_c({{dirVec[1]}}), it.z() + cell_idx_c({{dirVec[2]}})) && isFlagSet( it.neighbor({{offset}} {%if dim == 3%}, 0 {%endif %}), domainFlag ) )
            {
                {% if single_link -%}
                sum_x += cell_idx_c({{dirVec[0]}}); sum_y += cell_idx_c({{dirVec[1]}}); {%if dim == 3%} sum_z += cell_idx_c({{dirVec[2]}}); {%endif %}
                {% else %}
                auto element = {{StructName}}(it.x(), it.y(), {%if dim == 3%} it.z(), {%endif %} {{dirIdx}} );
                {{additional_data_handler.data_initialisation(dirIdx)|indent(16)}}
                indexVectorAll.push_back( element );
                if( inner.contains( it.x(), it.y(), it.z() ) )
                    indexVectorInner.push_back( element );
                else
                    indexVectorOuter.push_back( element );
                {% endif %}
            }
            {% endfor %}

        {% if single_link %}
            dot = 0.0; maxn = 0.0; calculated_idx = 0;
            if(sum_x != 0 or sum_y !=0 {%if dim == 3%} or sum_z !=0 {%endif %})
            {
            {%- for dirIdx, dirVec, offset in additional_data_handler.stencil_info %}
                dx = {{dirVec[0]}}; dy = {{dirVec[1]}}; {%if dim == 3%} dz = {{dirVec[2]}}; {% endif %}
                dot = {{dtype}}( dx*sum_x + dy*sum_y {%if dim == 3%} + dz*sum_z {% endif %});
                if (dot > maxn)
                {
                    maxn = dot;
                    calculated_idx = {{dirIdx}};
                }
            {% endfor %}
                auto element = {{StructName}}(it.x(), it.y(), {%if dim == 3%} it.z(), {%endif %} calculated_idx );
                indexVectorAll.push_back( element );
                if( inner.contains( it.x(), it.y(), it.z() ) )
                indexVectorInner.push_back( element );
                else
                indexVectorOuter.push_back( element );
            }
        {% endif -%}

        }
        {% endif %}

        indexVectors->syncGPU();
        {% if calculate_force -%}
        forceVector->forceVector().resize(indexVectorAll.size());
        forceVector->syncGPU();
        {%- endif %}
    }

private:
    void run_impl(
        {{- ["IBlock * block", "IndexVectors::Type type",
             kernel.kernel_selection_parameters, ["gpuStream_t stream = nullptr"] if target == 'gpu' else []]
            | type_identifier_list -}}
   );

    BlockDataID indexVectorID;
    {% if calculate_force -%}
    BlockDataID forceVectorID;
    {%- endif %}
    {{additional_data_handler.additional_member_variable|indent(4)}}
public:
    {{kernel|generate_members(('indexVector', 'indexVectorSize', 'forceVector', 'forceVectorSize'))|indent(4)}}
};



} // namespace {{namespace}}
} // namespace walberla
