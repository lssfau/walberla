#include "stencil/Directions.h"
#include "core/cell/CellInterval.h"
#include "cuda/GPUField.h"
#include "core/DataTypes.h"
#include "domain_decomposition/IBlock.h"
#include "cuda/communication/GeneratedGPUPackInfo.h"


#define FUNC_PREFIX __global__


namespace walberla {
namespace pystencils {


class UniformGridGPU_PackInfo : public ::walberla::cuda::GeneratedGPUPackInfo
{
public:
    UniformGridGPU_PackInfo( BlockDataID pdfsID_ )
        : pdfsID(pdfsID_)
    {};
    virtual ~UniformGridGPU_PackInfo() {}

    virtual void pack  (stencil::Direction dir, unsigned char * buffer, IBlock * block, cudaStream_t stream);
    virtual void unpack(stencil::Direction dir, unsigned char * buffer, IBlock * block, cudaStream_t stream);
    virtual uint_t size  (stencil::Direction dir, IBlock * block);

private:
    BlockDataID pdfsID;
};


} // namespace pystencils
} // namespace walberla