#include "../gpu/01_GameOfLife_kernels.h"


namespace walberla {


__global__ void gameOfLifeKernel( gpu::FieldAccessor<real_t> src, gpu::FieldAccessor<real_t> dst  )
{
   src.set( blockIdx, threadIdx );
   dst.set( blockIdx, threadIdx );

   // Count number of living neighbors
   int liveNeighbors = 0;
   if ( src.getNeighbor(  1, 0,0 ) > 0.5 ) ++liveNeighbors;
   if ( src.getNeighbor( -1, 0,0 ) > 0.5 ) ++liveNeighbors;
   if ( src.getNeighbor(  0,+1,0 ) > 0.5 ) ++liveNeighbors;
   if ( src.getNeighbor(  0,-1,0 ) > 0.5 ) ++liveNeighbors;

   if ( src.getNeighbor(  -1, -1, 0 ) > 0.5 ) ++liveNeighbors;
   if ( src.getNeighbor(  -1, +1, 0 ) > 0.5 ) ++liveNeighbors;
   if ( src.getNeighbor(  +1, -1,0 ) > 0.5 ) ++liveNeighbors;
   if ( src.getNeighbor(  +1, +1,0 ) > 0.5 ) ++liveNeighbors;


   // cell dies because of under- or over-population
   if ( liveNeighbors < 2 || liveNeighbors > 3 )
      dst.get() = 0.0;
   else if ( liveNeighbors == 3 ) // cell comes alive
      dst.get() = 1.0;
   else
      dst.get() = src.get();
}

void GameOfLifeSweepCUDA::operator()(IBlock * block)
{
   auto srcCudaField = block->getData< gpu::GPUField<real_t> > ( gpuFieldSrcID_ );
   auto dstCudaField = block->getData< gpu::GPUField<real_t> > ( gpuFieldDstID_ );

   auto srcIndexing = gpu::FieldIndexing<real_t>::xyz( *srcCudaField );
   auto dstIndexing = gpu::FieldIndexing<real_t>::xyz( *dstCudaField );

   auto srcAccess = srcIndexing.gpuAccess();
   auto dstAccess = dstIndexing.gpuAccess();

   const dim3 gridDim = srcIndexing.gridDim();
   const dim3 blockDim = srcIndexing.blockDim();

   gameOfLifeKernel<<<gridDim, blockDim, 0, nullptr >>>(srcAccess, dstAccess );

   srcCudaField->swapDataPointers( dstCudaField );
}




} // namespace walberla
