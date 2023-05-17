#pragma once
#include "core/DataTypes.h"
#include "core/logging/Logging.h"

#include "gpu/FieldIndexing.h"

#include "field/SwapableCompare.h"
#include "domain_decomposition/BlockDataID.h"
#include "domain_decomposition/IBlock.h"
#include "domain_decomposition/StructuredBlockStorage.h"

namespace walberla {

class GameOfLifeSweepCUDA
{
 public:
   GameOfLifeSweepCUDA( BlockDataID gpuFieldSrcID, BlockDataID gpuFieldDstID )
      : gpuFieldSrcID_( gpuFieldSrcID ), gpuFieldDstID_( gpuFieldDstID ){}

   void operator() ( IBlock * block );

 private:
   BlockDataID gpuFieldSrcID_;
   BlockDataID gpuFieldDstID_;
};


__global__ void gameOfLifeKernel(gpu::FieldAccessor<real_t> src, gpu::FieldAccessor<real_t> dst  );


} // namespace walberla
