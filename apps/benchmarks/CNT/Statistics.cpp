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

#include "Statistics.h"
#include "TerminalColors.h"

#include "core/logging/Logging.h"

#include <fstream>

namespace walberla {
namespace mesa_pd {

void saveStatistics(const std::string& filename,
                    const std::vector<Statistics>& stats)
{
   WALBERLA_LOG_INFO_ON_ROOT(CYAN << "saving statistics: " << filename << RESET);
   std::ofstream fout(filename);
   fout.precision(10);
   fout << "# kineticEnergy vdwEnergy stretchEnergy tensileEnergy shearEnergy bendingEnergy twistingEnergy numAlignedPairs" << std::endl;
   for (const auto& stat : stats)
   {
      fout << stat.kineticEnergy << " "
           << stat.vdwEnergy << " "
           << stat.stretchEnergy << " "
           << stat.tensileEnergy << " "
           << stat.shearEnergy << " "
           << stat.bendingEnergy << " "
           << stat.twistingEnergy << " "
           << stat.numAlignedPairs << std::endl;
   }
   fout.close();
}

} //namespace mesa_pd

namespace mpi {
void reduceInplace( mesa_pd::Statistics& stats, mpi::Operation operation, int recvRank, MPI_Comm comm )
{
   reduceInplace(stats.kineticEnergy, operation, recvRank, comm);
   reduceInplace(stats.vdwEnergy, operation, recvRank, comm);
   reduceInplace(stats.stretchEnergy, operation, recvRank, comm);
   reduceInplace(stats.tensileEnergy, operation, recvRank, comm);
   reduceInplace(stats.shearEnergy, operation, recvRank, comm);
   reduceInplace(stats.bendingEnergy, operation, recvRank, comm);
   reduceInplace(stats.twistingEnergy, operation, recvRank, comm);
   reduceInplace(stats.numAlignedPairs, operation, recvRank, comm);
}
} //namespace mpi

} //namespace walberla