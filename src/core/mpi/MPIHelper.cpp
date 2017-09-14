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
//! \file MPIHelper.cpp
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include "MPIHelper.h"

namespace walberla {
namespace mpi {

int translateRank(const MPI_Comm srcComm, const MPI_Comm destComm, const int srcRank)
{
   int destRank = -1;
   MPI_Group srcGroup, destGroup;
   MPI_Comm_group(srcComm, &srcGroup);
   MPI_Comm_group(destComm, &destGroup);
   MPI_Group_translate_ranks(srcGroup, 1, const_cast<int*>(&srcRank), destGroup, &destRank);
   MPI_Group_free(&srcGroup);
   MPI_Group_free(&destGroup);
   return destRank;
}

std::vector<int> translateRank(const MPI_Comm srcComm, const MPI_Comm destComm, const std::vector<int>& srcRank)
{
   std::vector<int> destRank(srcRank.size(), -1);
   MPI_Group srcGroup, destGroup;
   MPI_Comm_group(srcComm, &srcGroup);
   MPI_Comm_group(destComm, &destGroup);
   MPI_Group_translate_ranks(srcGroup, int_c(srcRank.size()), const_cast<int*>(&srcRank[0]), destGroup, &destRank[0]);
   MPI_Group_free(&srcGroup);
   MPI_Group_free(&destGroup);
   return destRank;
}

} // namespace mpi
} // namespace walberla
