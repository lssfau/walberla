#pragma once

#include "core/DataTypes.h"
#include "core/cell/CellInterval.h"

#include "domain_decomposition/IBlock.h"

#include "stencil/Directions.h"

#include "gpu/GPUField.h"
#include "gpu/GPUWrapper.h"
#include "gpu/communication/GeneratedGPUPackInfo.h"

namespace walberla::gpu::communication {

template<typename GPUFieldType>
class MemcpyPackInfo : public ::walberla::gpu::GeneratedGPUPackInfo
{
public:
    MemcpyPackInfo( BlockDataID pdfsID_ ) : pdfsID(pdfsID_) {};
    ~MemcpyPackInfo() override = default;

    void pack  (stencil::Direction dir, unsigned char * buffer, IBlock * block, gpuStream_t stream) override;
    void communicateLocal  ( stencil::Direction dir, const IBlock *sender, IBlock *receiver, gpuStream_t stream ) override;
    void unpack(stencil::Direction dir, unsigned char * buffer, IBlock * block, gpuStream_t stream) override;
    uint_t size(stencil::Direction dir, IBlock * block) override;

private:
    BlockDataID pdfsID;
    uint_t numberOfGhostLayers_{0};
    bool communicateAllGhostLayers_{true};

    uint_t numberOfGhostLayersToCommunicate( const GPUFieldType * const field ) const;
};

} // namespace walberla::gpu::communication

#include "MemcpyPackInfo.impl.h"
