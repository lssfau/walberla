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

#include "Config.h"
#include "TerminalColors.h"

#include "core/mpi/Gatherv.h"
#include "core/mpi/MPIIO.h"
#include "core/mpi/RecvBuffer.h"
#include "core/mpi/SendBuffer.h"

namespace walberla {
namespace mesa_pd {

void saveConfig(const shared_ptr<data::ParticleStorage>& ps,
                const std::string& filename,
                const IO_MODE io_mode)
{
   mpi::SendBuffer sb;

   for (const auto& p : *ps)
   {
      using namespace walberla::mesa_pd::data::particle_flags;
      if (isSet(p.getFlags(), GHOST)) continue; //skip ghosts

      sb << static_cast<int64_t> (p.getSegmentID());
      sb << static_cast<int64_t> (p.getClusterID());
//      sb << static_cast<int64_t> (p.getGroupID());

      sb << static_cast<double> (p.getPosition()[0]);
      sb << static_cast<double> (p.getPosition()[1]);
      sb << static_cast<double> (p.getPosition()[2]);

      Vec3 orient = p.getRotation().getMatrix() * Vec3(1,0,0);
      // Now orient has the shape (sin(th)cos(ph), sin(th)sin(ph), cos(th))
      sb << static_cast<double> (std::acos(orient[2]));
      sb << static_cast<double> (std::atan2(orient[1], orient[0]));

      sb << static_cast<double> (p.getLinearVelocity()[0]);
      sb << static_cast<double> (p.getLinearVelocity()[1]);
      sb << static_cast<double> (p.getLinearVelocity()[2]);

      sb << static_cast<double> (p.getAngularVelocity()[0]);
      sb << static_cast<double> (p.getAngularVelocity()[1]);
      sb << static_cast<double> (p.getAngularVelocity()[2]);
   }

   switch (io_mode)
   {
      case IO_MODE::REDUCE_TO_ROOT:
      {
         WALBERLA_LOG_INFO_ON_ROOT(CYAN << "Saving configuration (reduce to root): " << filename << RESET);
         walberla::mpi::RecvBuffer rb;
         walberla::mpi::gathervBuffer(sb, rb);
         WALBERLA_ROOT_SECTION()
         {
            uint_t dataSize = sizeof(mpi::RecvBuffer::ElementType) * rb.size();
            std::ofstream binfile;
            binfile.open(filename + ".sav",
                         std::ios::out | std::ios::binary);
            binfile.write(reinterpret_cast< const char * >(rb.ptr()), numeric_cast<std::streamsize>(dataSize));
            binfile.close();
         }
         break;
      }
      case IO_MODE::USE_MPI_IO:
      {
         WALBERLA_LOG_INFO_ON_ROOT(CYAN << "Saving configuration (MPIIO): " << filename << RESET);
         walberla::mpi::writeMPIIO(filename + ".mpiio", sb);
         break;
      }
      case IO_MODE::ONE_FILE_PER_RANK:
      {
         WALBERLA_LOG_INFO_ON_ROOT(CYAN << "Saving configuration (one file per process): " << filename << RESET);
         uint_t dataSize = sizeof(mpi::SendBuffer::ElementType) * sb.size();
         std::ofstream binfile;
         binfile.open(filename + "_" + std::to_string(mpi::MPIManager::instance()->rank()) + ".sav",
                      std::ios::out | std::ios::binary);
         binfile.write(reinterpret_cast< const char * >(sb.ptr()), numeric_cast<std::streamsize>(dataSize));
         binfile.close();
      }
         break;
      default:
      WALBERLA_ABORT("Unknown IO_MODE");
   }
}

} //namespace mesa_pd
} //namespace walberla
