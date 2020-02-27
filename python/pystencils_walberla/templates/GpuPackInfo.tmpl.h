#pragma once
#include "stencil/Directions.h"
#include "core/cell/CellInterval.h"
#include "cuda/GPUField.h"
#include "core/DataTypes.h"
#include "domain_decomposition/IBlock.h"
#include "cuda/communication/GeneratedGPUPackInfo.h"


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


class {{class_name}} : public ::walberla::cuda::GeneratedGPUPackInfo
{
public:
    {{class_name}}( {{fused_kernel|generate_constructor_parameters(parameters_to_ignore=['buffer'])}} )
        : {{ fused_kernel|generate_constructor_initializer_list(parameters_to_ignore=['buffer']) }}
    {};
    virtual ~{{class_name}}() {}

    virtual void pack  (stencil::Direction dir, unsigned char * buffer, IBlock * block, cudaStream_t stream);
    virtual void unpack(stencil::Direction dir, unsigned char * buffer, IBlock * block, cudaStream_t stream);
    virtual uint_t size  (stencil::Direction dir, IBlock * block);

private:
    {{fused_kernel|generate_members(parameters_to_ignore=['buffer'])|indent(4)}}
};


} // namespace {{namespace}}
} // namespace walberla
