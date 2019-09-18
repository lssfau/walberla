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
//! \file   LinkedCells.cpp
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include <mesa_pd/data/LinkedCells.h>

#include <core/Environment.h>
#include <core/logging/Logging.h>

#include <algorithm>
#include <iostream>

namespace walberla {
namespace mesa_pd {

int main( int argc, char ** argv )
{
   Environment env(argc, argv);
   WALBERLA_UNUSED(env);
   mpi::MPIManager::instance()->useWorldComm();

   data::LinkedCells lc(math::AABB(real_t(0),real_t(0),real_t(0),
                                   real_t(4.5),real_t(5.5),real_t(6.5)),
                        real_t(1));

   WALBERLA_CHECK_EQUAL(lc.numCellsPerDim_[0], 5);
   WALBERLA_CHECK_EQUAL(lc.numCellsPerDim_[1], 6);
   WALBERLA_CHECK_EQUAL(lc.numCellsPerDim_[2], 7);

   auto cellIdx = data::getCellIdx(lc, 2, 1, 3);
   WALBERLA_CHECK_EQUAL(cellIdx, 90 + 5 + 2);
   int64_t x = 0;
   int64_t y = 0;
   int64_t z = 0;
   data::getCellCoordinates(lc, cellIdx, x, y, z);
   WALBERLA_CHECK_EQUAL(x, 2);
   WALBERLA_CHECK_EQUAL(y, 1);
   WALBERLA_CHECK_EQUAL(z, 3);

   return EXIT_SUCCESS;
}

} //namespace mesa_pd
} //namespace walberla

int main( int argc, char ** argv )
{
   return walberla::mesa_pd::main(argc, argv);
}
