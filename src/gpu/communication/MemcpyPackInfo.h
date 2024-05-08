#pragma once

#include "core/DataTypes.h"
#include "core/cell/CellInterval.h"

#include "domain_decomposition/IBlock.h"

#include "stencil/Directions.h"

#include "gpu/GPUField.h"
#include "gpu/GPUWrapper.h"
#include "gpu/communication/GeneratedGPUPackInfo.h"

namespace walberla::gpu::communication {

/**
 * \brief Data packing/unpacking for ghost layer based communication of a \ref GPUField.
 *
 * Encapsulate information on how to extract data from blocks that should be
 * communicated to neighboring blocks (see \ref pack())
 * and how to inject this data in a receiving block (see \ref unpack()).
 * This involves a device memory buffer and two device-to-device memory copy operations.
 *
 * A special method exists for communication between two blocks which are
 * allocated on the same process (see \ref communicateLocal()).
 * In this case the data does not have be communicated via a device buffer,
 * but can be sent directly. This involves a single device-to-device memory
 * copy operation.
 *
 * Data that is packed in direction "dir" at one block is unpacked in
 * direction "stencil::inverseDir[dir]" at the neighboring block.
 * This behavior must be implemented in \ref communicateLocal()!
 *
 * \ingroup gpu
 * \tparam GPUFieldType   A fully qualified \ref GPUField.
 */
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
