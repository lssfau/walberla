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

#include <core/debug/CheckFunctions.h>

namespace walberla {
namespace mpi {

//!
//! \brief This functions maps the rank in one communicator to the rank in another communicator.
//! \param srcComm source communicator
//! \param destComm destination communicator
//! \param srcRank rank in the source communicator
//! \return rank in the destination communicator or -1 if not available
//!
int translateRank(const MPI_Comm srcComm, const MPI_Comm destComm, const int srcRank)
{
   int destRank = -1;
   MPI_Group srcGroup, destGroup;
   MPI_Comm_group(srcComm, &srcGroup);
   MPI_Comm_group(destComm, &destGroup);
   MPI_Group_translate_ranks(srcGroup, 1, const_cast<int*>(&srcRank), destGroup, &destRank);
   int size;
   MPI_Comm_size(destComm, &size);
   if (destRank < 0 || destRank >= size)
   {
      WALBERLA_CHECK_EQUAL( destRank,  MPI_UNDEFINED );
      destRank = -1;
   }
   MPI_Group_free(&srcGroup);
   MPI_Group_free(&destGroup);
   return destRank;
}

//!
//! \brief This functions converts a array of ranks in one communicator to an array of ranks in another communicator.
//! \param srcComm source communicator
//! \param destComm destination communicator
//! \param srcRank source ranks
//! \return converted ranks, -1 if not available
//!
std::vector<int> translateRank(const MPI_Comm srcComm, const MPI_Comm destComm, const std::vector<int>& srcRank)
{
   std::vector<int> destRank(srcRank.size(), -1);
   MPI_Group srcGroup, destGroup;
   MPI_Comm_group(srcComm, &srcGroup);
   MPI_Comm_group(destComm, &destGroup);
   MPI_Group_translate_ranks(srcGroup, int_c(srcRank.size()), const_cast<int*>(&srcRank[0]), destGroup, &destRank[0]);
   int size;
   MPI_Comm_size(destComm, &size);
   for (auto& dstRnk : destRank)
   {
      if (dstRnk < 0 || dstRnk >= size)
      {
         WALBERLA_CHECK_EQUAL( dstRnk,  MPI_UNDEFINED );
         dstRnk = -1;
      }
   }
   MPI_Group_free(&srcGroup);
   MPI_Group_free(&destGroup);
   return destRank;
}

} // namespace mpi
} // namespace walberla
