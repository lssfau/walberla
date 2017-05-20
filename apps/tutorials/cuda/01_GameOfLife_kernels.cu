#include "../cuda/01_GameOfLife_kernels.h"

#include <iostream>



namespace walberla {


__global__ void gameOfLifeKernel( cuda::FieldAccessor<double> src, cuda::FieldAccessor<double> dst  )
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




} // namespace walberla
