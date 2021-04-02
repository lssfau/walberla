//======================================================================================================================
//
//  This file is part of waLBerla. waLBerla is free software: you can
//  redistribute it and/or modify it under the terms of the GNU General Public
//  License as published by the Free Software Foundation, either version 3 of
//  the License, or (at your option) any later version.
//
//  waLBerla is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
//  for more details.
//
//  You should have received a copy of the GNU General Public License along
//  with waLBerla (see COPYING.txt). If not, see <http://www.gnu.org/licenses/>.
//
//! \file
//! \author Igor Ostanin <i.ostanin@skoltech.ru>
//! \author Grigorii Drozdov <drozd013@umn.edu>
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/mpi/Reduce.h"

#include <string>
#include <vector>

namespace walberla {
namespace mesa_pd {

struct Statistics
{
   real_t kineticEnergy   = 0_r;
   real_t vdwEnergy       = 0_r;
   real_t stretchEnergy   = 0_r;
   real_t tensileEnergy   = 0_r;
   real_t shearEnergy     = 0_r;
   real_t bendingEnergy   = 0_r;
   real_t twistingEnergy  = 0_r;
   int64_t numAlignedPairs = 0;
};

void saveStatistics(const std::string& filename,
                    const std::vector<Statistics>& stats);

} //namespace mesa_pd

namespace mpi {
void reduceInplace( mesa_pd::Statistics& stats, mpi::Operation operation, int recvRank = 0, MPI_Comm comm = MPI_COMM_WORLD );
} //namespace mpi

} //namespace walberla
