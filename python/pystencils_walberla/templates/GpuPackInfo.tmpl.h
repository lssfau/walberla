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
#include "core/cell/CellInterval.h"

#include "domain_decomposition/IBlock.h"

#include "stencil/Directions.h"

#include "gpu/GPUField.h"
#include "gpu/GPUWrapper.h"
#include "gpu/communication/GeneratedGPUPackInfo.h"

{% if target is equalto 'cpu' -%}
#define FUNC_PREFIX
{%- elif target is equalto 'gpu' -%}
#define FUNC_PREFIX __global__
{%- endif %}

#ifdef __GNUC__
#define RESTRICT __restrict__
#elif _MSC_VER
#define RESTRICT __restrict
#else
#define RESTRICT
#endif

namespace walberla {
namespace {{namespace}} {


class {{class_name}} : public ::walberla::gpu::GeneratedGPUPackInfo
{
public:
    {{class_name}}( {{fused_kernel|generate_constructor_parameters(parameters_to_ignore=['buffer'])}} )
        : {{ fused_kernel|generate_constructor_initializer_list(parameters_to_ignore=['buffer']) }}
    {};
    virtual ~{{class_name}}() {}

    void pack  (stencil::Direction dir, unsigned char * buffer, IBlock * block, gpuStream_t stream) override;
    void communicateLocal  ( stencil::Direction /*dir*/, const IBlock* /* sender */, IBlock* /* receiver */, gpuStream_t /* stream */ ) override
    {
       WALBERLA_ABORT("Local Communication not implemented yet for standard PackInfos. To run your application turn of local communication in the Communication class")
    }
    void unpack(stencil::Direction dir, unsigned char * buffer, IBlock * block, gpuStream_t stream) override;
    uint_t size  (stencil::Direction dir, IBlock * block) override;

private:
    {{fused_kernel|generate_members(parameters_to_ignore=['buffer'])|indent(4)}}
};


} // namespace {{namespace}}
} // namespace walberla
