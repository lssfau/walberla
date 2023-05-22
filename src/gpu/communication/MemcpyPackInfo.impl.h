#include "blockforest/Block.h"
#include "field/GhostRegions.h"
#include "field/Layout.h"
#include "stencil/Directions.h"
#include "core/cell/CellInterval.h"
#include "gpu/GPUField.h"
#include "gpu/GPUCopy.h"
#include "core/DataTypes.h"
#include "MemcpyPackInfo.h"


namespace walberla {
namespace gpu
{
namespace communication {

template<typename GPUFieldType>
void MemcpyPackInfo< GPUFieldType >::pack(stencil::Direction dir, unsigned char * byte_buffer,
                                          IBlock * block, gpuStream_t stream)
{
   // Extract field data pointer from the block
   const GPUFieldType * fieldPtr = block->getData< GPUFieldType >( pdfsID );
   WALBERLA_ASSERT_NOT_NULLPTR( fieldPtr )
   // 
   cell_idx_t nrOfGhostLayers = cell_idx_c( numberOfGhostLayersToCommunicate( fieldPtr ) );
   CellInterval fieldCi = field::getSliceBeforeGhostLayer( *fieldPtr, dir, nrOfGhostLayers, false );

   // Base offsets into the buffer and GPUField, respectively
   auto dstOffset = std::make_tuple( uint_c(0), uint_c(0), uint_c(0), uint_c(0) );
   auto srcOffset = std::make_tuple( uint_c(fieldCi.xMin() + nrOfGhostLayers),
                                     uint_c(fieldCi.yMin() + nrOfGhostLayers),
                                     uint_c(fieldCi.zMin() + nrOfGhostLayers),
                                     uint_c(0) );

   // Size of data to pack, in terms of elements of the field
   auto intervalSize = std::make_tuple( fieldCi.xSize(), fieldCi.ySize(),
                                        fieldCi.zSize(), fieldPtr->fSize() );

   if ( fieldPtr->layout() == field::fzyx )
   {
      const uint_t dstAllocSizeZ = fieldCi.zSize();
      const uint_t srcAllocSizeZ = fieldPtr->zAllocSize();

      gpuPitchedPtr byteBufferPitchedPtr = make_gpuPitchedPtr( byte_buffer,
                                                                 fieldCi.xSize() * sizeof(typename GPUFieldType::value_type),
                                                                 fieldCi.xSize() * sizeof(typename GPUFieldType::value_type),
                                                                 fieldCi.ySize() );

      copyDevToDevFZYX( byteBufferPitchedPtr, fieldPtr->pitchedPtr(), dstOffset, srcOffset,
                        dstAllocSizeZ, srcAllocSizeZ, sizeof(typename GPUFieldType::value_type),
                        intervalSize, stream );
   }
   else
   {
      const uint_t dstAllocSizeZ = fieldCi.ySize();
      const uint_t srcAllocSizeZ = fieldPtr->yAllocSize();

      gpuPitchedPtr byteBufferPitchedPtr = make_gpuPitchedPtr( byte_buffer,
                                                                 fieldPtr->fSize() * sizeof(typename GPUFieldType::value_type),
                                                                 fieldPtr->fSize() * sizeof(typename GPUFieldType::value_type),
                                                                 fieldCi.xSize() );
      copyDevToDevZYXF( byteBufferPitchedPtr, fieldPtr->pitchedPtr(), dstOffset, srcOffset,
                        dstAllocSizeZ, srcAllocSizeZ, sizeof(typename GPUFieldType::value_type),
                        intervalSize, stream );
   }
}

template<typename GPUFieldType>
void MemcpyPackInfo< GPUFieldType >::communicateLocal( stencil::Direction dir, const IBlock* sender, IBlock* receiver, gpuStream_t stream )
{
   // WALBERLA_ABORT("The MemcpyPackInfo does not provide a thread safe local communication. Thus is can not be used in local mode. To use it set local useLocalCommunication to false in the communication scheme")


   // Extract field data pointer from the block
   const GPUFieldType * senderFieldPtr = sender->getData< GPUFieldType >( pdfsID );
   const GPUFieldType * receiverFieldPtr = receiver->getData< GPUFieldType >( pdfsID );
   WALBERLA_ASSERT_NOT_NULLPTR( senderFieldPtr )
   WALBERLA_ASSERT_NOT_NULLPTR( receiverFieldPtr )

   //
   cell_idx_t nrOfGhostLayers = cell_idx_c( numberOfGhostLayersToCommunicate( senderFieldPtr ) );
   WALBERLA_ASSERT_EQUAL(nrOfGhostLayers, cell_idx_c( numberOfGhostLayersToCommunicate( receiverFieldPtr )))
   WALBERLA_ASSERT_EQUAL(senderFieldPtr->layout(), receiverFieldPtr->layout() )
   WALBERLA_ASSERT_EQUAL(senderFieldPtr->fSize(), receiverFieldPtr->fSize() )

   CellInterval senderCi = field::getSliceBeforeGhostLayer( *senderFieldPtr, dir, nrOfGhostLayers, false );
   CellInterval receiverCi = field::getGhostRegion( *receiverFieldPtr, stencil::inverseDir[dir], nrOfGhostLayers, false );

   // Base offsets into the buffer and GPUField, respectively
   auto srcOffset = std::make_tuple( uint_c(senderCi.xMin() + nrOfGhostLayers),
                                     uint_c(senderCi.yMin() + nrOfGhostLayers),
                                     uint_c(senderCi.zMin() + nrOfGhostLayers),
                                     uint_c(0) );

   auto dstOffset = std::make_tuple( uint_c(receiverCi.xMin() + nrOfGhostLayers),
                                     uint_c(receiverCi.yMin() + nrOfGhostLayers),
                                     uint_c(receiverCi.zMin() + nrOfGhostLayers),
                                     uint_c(0) );


   // Size of data to pack, in terms of elements of the field
   auto intervalSize = std::make_tuple( senderCi.xSize(), senderCi.ySize(),
                                        senderCi.zSize(), senderFieldPtr->fSize() );

   WALBERLA_ASSERT_EQUAL(intervalSize, std::make_tuple( receiverCi.xSize(), receiverCi.ySize(), receiverCi.zSize(), receiverFieldPtr->fSize() ))

   if ( senderFieldPtr->layout() == field::fzyx )
   {
      const uint_t dstAllocSizeZ = receiverFieldPtr->zAllocSize();
      const uint_t srcAllocSizeZ = senderFieldPtr->zAllocSize();

      copyDevToDevFZYX( receiverFieldPtr->pitchedPtr(), senderFieldPtr->pitchedPtr(), dstOffset, srcOffset,
                       dstAllocSizeZ, srcAllocSizeZ, sizeof(typename GPUFieldType::value_type),
                       intervalSize, stream );
   }
   else
   {
      const uint_t dstAllocSizeZ = receiverFieldPtr->yAllocSize();
      const uint_t srcAllocSizeZ = senderFieldPtr->yAllocSize();

      copyDevToDevZYXF( receiverFieldPtr->pitchedPtr(), senderFieldPtr->pitchedPtr(), dstOffset, srcOffset,
                       dstAllocSizeZ, srcAllocSizeZ, sizeof(typename GPUFieldType::value_type),
                       intervalSize, stream );
   }
}

template<typename GPUFieldType>
void MemcpyPackInfo< GPUFieldType >::unpack(stencil::Direction dir, unsigned char * byte_buffer,
                                            IBlock * block, gpuStream_t stream)
{
   GPUFieldType * fieldPtr = block->getData< GPUFieldType >( pdfsID );
   WALBERLA_ASSERT_NOT_NULLPTR(fieldPtr)

   cell_idx_t nrOfGhostLayers = cell_idx_c( numberOfGhostLayersToCommunicate( fieldPtr ) );

   CellInterval fieldCi = field::getGhostRegion( *fieldPtr, dir, nrOfGhostLayers, false );
   auto dstOffset = std::make_tuple( uint_c(fieldCi.xMin() + nrOfGhostLayers),
                                     uint_c(fieldCi.yMin() + nrOfGhostLayers),
                                     uint_c(fieldCi.zMin() + nrOfGhostLayers),
                                     uint_c(0) );
   auto srcOffset = std::make_tuple( uint_c(0), uint_c(0), uint_c(0), uint_c(0) );

   auto intervalSize = std::make_tuple( fieldCi.xSize(), fieldCi.ySize(), fieldCi.zSize(), fieldPtr->fSize() );

   if ( fieldPtr->layout() == field::fzyx )
   {
      const uint_t dstAllocSizeZ = fieldPtr->zAllocSize();
      const uint_t srcAllocSizeZ = fieldCi.zSize();

      gpuPitchedPtr byteBufferPitchedPtr = make_gpuPitchedPtr( byte_buffer,
                                                                 fieldCi.xSize() * sizeof(typename GPUFieldType::value_type),
                                                                 fieldCi.xSize() * sizeof(typename GPUFieldType::value_type),
                                                                 fieldCi.ySize() );

      copyDevToDevFZYX( fieldPtr->pitchedPtr(), byteBufferPitchedPtr, dstOffset, srcOffset,
                        dstAllocSizeZ, srcAllocSizeZ, sizeof(typename GPUFieldType::value_type),
                        intervalSize, stream );
   }
   else
   {
      const uint_t dstAllocSizeY = fieldPtr->yAllocSize();
      const uint_t srcAllocSizeY = fieldCi.ySize();
      gpuPitchedPtr byteBufferPitchedPtr = make_gpuPitchedPtr( byte_buffer,
                                                                 fieldPtr->fSize() * sizeof(typename GPUFieldType::value_type),
                                                                 fieldPtr->fSize() * sizeof(typename GPUFieldType::value_type),
                                                                 fieldCi.xSize() );
      copyDevToDevZYXF( fieldPtr->pitchedPtr(), byteBufferPitchedPtr, dstOffset, srcOffset,
                        dstAllocSizeY, srcAllocSizeY, sizeof(typename GPUFieldType::value_type),
                        intervalSize, stream );
   }
}

template<typename GPUFieldType>
uint_t MemcpyPackInfo< GPUFieldType >::size(stencil::Direction dir, IBlock * block)
{
    auto pdfs = block->getData< GPUFieldType >(pdfsID);

    CellInterval ci;
    cell_idx_t nrOfGhostLayers = cell_idx_c( numberOfGhostLayersToCommunicate( pdfs ) );
    pdfs->getGhostRegion(dir, ci, nrOfGhostLayers, false);

    /*
    uint_t elementsPerCell = 0;

    switch( dir )
    {
        case stencil::SW:
            elementsPerCell = 1;
            break;
        
        case stencil::S:
            elementsPerCell = 5;
            break;
        
        case stencil::W:
            elementsPerCell = 5;
            break;
        
        case stencil::B:
            elementsPerCell = 5;
            break;
        
        case stencil::T:
            elementsPerCell = 5;
            break;
        
        case stencil::BN:
            elementsPerCell = 1;
            break;
        
        case stencil::N:
            elementsPerCell = 5;
            break;
        
        case stencil::TE:
            elementsPerCell = 1;
            break;
        
        case stencil::E:
            elementsPerCell = 5;
            break;
        
        case stencil::BE:
            elementsPerCell = 1;
            break;
        
        case stencil::SE:
            elementsPerCell = 1;
            break;
        
        case stencil::C:
            elementsPerCell = 1;
            break;
        
        case stencil::TN:
            elementsPerCell = 1;
            break;
        
        case stencil::TS:
            elementsPerCell = 1;
            break;
        
        case stencil::NE:
            elementsPerCell = 1;
            break;
        
        case stencil::BW:
            elementsPerCell = 1;
            break;
        
        case stencil::NW:
            elementsPerCell = 1;
            break;
        
        case stencil::BS:
            elementsPerCell = 1;
            break;
        
        case stencil::TW:
            elementsPerCell = 1;
            break;
        
        default:
            elementsPerCell = 0;
    }

    return ci.numCells() * elementsPerCell * sizeof(typename GPUFieldType::value_type);
    */
    uint_t totalCells = ci.numCells() * pdfs->fSize() * sizeof(typename GPUFieldType::value_type);
    return totalCells;
}

template<typename GPUFieldType>
uint_t MemcpyPackInfo< GPUFieldType >::numberOfGhostLayersToCommunicate( const GPUFieldType * const field ) const
{
   if( communicateAllGhostLayers_ )
   {
      return field->nrOfGhostLayers();
   }
   else
   {
      WALBERLA_ASSERT_LESS_EQUAL( numberOfGhostLayers_, field->nrOfGhostLayers() )
      return numberOfGhostLayers_;
   }
}

} // namespace communication
} // namespace gpu
} // namespace walberla