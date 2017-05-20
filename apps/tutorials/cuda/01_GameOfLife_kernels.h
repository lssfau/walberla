#include <iostream>

#include "cuda/FieldIndexing.h"


namespace walberla {


__global__ void gameOfLifeKernel( cuda::FieldAccessor<double> src, cuda::FieldAccessor<double> dst  );


} // namespace walberla
