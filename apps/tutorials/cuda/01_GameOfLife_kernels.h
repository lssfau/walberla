#pragma once

#include <iostream>

#include "cuda/FieldAccessor.h"


namespace walberla {


__global__ void gameOfLifeKernel( cuda::FieldAccessor<double> src, cuda::FieldAccessor<double> dst  );


} // namespace walberla
