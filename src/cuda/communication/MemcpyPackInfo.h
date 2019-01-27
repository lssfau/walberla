#pragma once

#include "stencil/Directions.h"
#include "core/cell/CellInterval.h"
#include "cuda/GPUField.h"
#include "core/DataTypes.h"
#include "domain_decomposition/IBlock.h"
#include "cuda/communication/GeneratedGPUPackInfo.h"


namespace walberla {
namespace cuda {
namespace communication {

template<typename GPUFieldType>
class MemcpyPackInfo : public ::walberla::cuda::GeneratedGPUPackInfo
{
public:
    MemcpyPackInfo( BlockDataID pdfsID_ )
        : pdfsID(pdfsID_), numberOfGhostLayers_(0), communicateAllGhostLayers_(true)
    {};


    virtual void pack  (stencil::Direction dir, unsigned char * buffer, IBlock * block, cudaStream_t stream);
    virtual void unpack(stencil::Direction dir, unsigned char * buffer, IBlock * block, cudaStream_t stream);
    virtual uint_t size(stencil::Direction dir, IBlock * block);

private:
    BlockDataID pdfsID;
    uint_t numberOfGhostLayers_;
    bool communicateAllGhostLayers_;

    uint_t numberOfGhostLayersToCommunicate( const GPUFieldType * const field ) const;
};

} // namespace communication
} // namespace cuda
} // namespace walberla

#include "MemcpyPackInfo.impl.h"
