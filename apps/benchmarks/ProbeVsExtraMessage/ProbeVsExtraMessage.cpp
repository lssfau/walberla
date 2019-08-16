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
//! \file ProbeVsExtraMessage.h
//! \author Martin Bauer <martin.bauer@fau.de>
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//! \brief Micro Benchmark, measuring time for different variable sized communications
//
//======================================================================================================================

#include "core/DataTypes.h"
#include "core/Environment.h"
#include "core/math/Vector3.h"
#include "core/mpi/BufferSystem.h"
#include "core/mpi/MPIManager.h"
#include "core/timing/TimingPool.h"
#include "sqlite/SQLite.h"
#include "stencil/D3Q27.h"
#include "stencil/D3Q19.h"
#include "stencil/D3Q7.h"

#include <array>
#include <iostream>
#include <sstream>

namespace walberla {

class CustomBufferSystem : public mpi::BufferSystem
{
public:
   explicit CustomBufferSystem( const MPI_Comm & communicator, int tag = 0 )
      : mpi::BufferSystem(communicator, tag)
   {}
   auto& recvBuffer ( walberla::mpi::MPIRank rank )
   {
      auto it = recvInfos_.find(rank);
      WALBERLA_CHECK_UNEQUAL(it, recvInfos_.end(), recvInfos_.size());
      return it->second.buffer;
   }

   auto& getRecvInfos() {return recvInfos_;}
};

class MPIInfo
{
public:
   MPIInfo( const Vector3<uint_t>& procs, const Vector3<bool>& periodicity );

   int getNeighborRank(const stencil::Direction& dir);
private:
   shared_ptr<mpi::MPIManager> manager_;
   Vector3<uint_t> procs_;
   Vector3<bool>   periodicity_;
   Vector3<int>    pos_;
   std::array<int, 27> ranks_;
};

MPIInfo::MPIInfo( const Vector3<uint_t>& procs, const Vector3<bool>& periodicity )
   : manager_(mpi::MPIManager::instance())
   , procs_(procs)
   , periodicity_(periodicity)
{
   manager_->createCartesianComm(procs[0], procs[1], procs[2], periodicity[0], periodicity[1], periodicity[2]);
   manager_->cartesianCoord(pos_.data());

   for (auto dirIt = stencil::D3Q27::beginNoCenter(); dirIt != stencil::D3Q27::end(); ++dirIt)
   {
      auto neighborCoord = pos_;
      neighborCoord[0] += stencil::cx[*dirIt];
      neighborCoord[1] += stencil::cy[*dirIt];
      neighborCoord[2] += stencil::cz[*dirIt];
      if (!periodicity_[0] && (neighborCoord[0] < 0)) ranks_[*dirIt] = -1;
      if (!periodicity_[1] && (neighborCoord[1] < 0)) ranks_[*dirIt] = -1;
      if (!periodicity_[2] && (neighborCoord[2] < 0)) ranks_[*dirIt] = -1;
      if (!periodicity_[0] && (neighborCoord[0] >= int_c(procs_[0]))) ranks_[*dirIt] = -1;
      if (!periodicity_[1] && (neighborCoord[1] >= int_c(procs_[1]))) ranks_[*dirIt] = -1;
      if (!periodicity_[2] && (neighborCoord[2] >= int_c(procs_[2]))) ranks_[*dirIt] = -1;
      ranks_[*dirIt] = manager_->cartesianRank(uint_c(neighborCoord[0]), uint_c(neighborCoord[1]), uint_c(neighborCoord[2]));
   }
}

inline
int MPIInfo::getNeighborRank( const stencil::Direction& dir )
{
   return ranks_[dir];
}

template <typename Stencil>
void communicate( MPIInfo& mpiInfo,
                  const uint_t iterations,
                  const uint_t messageSize,
                  const bool iProbe,
                  WcTimingPool& tp)
{
   std::vector<char> sendBuf(messageSize);
   std::vector<char> recvBuf(messageSize);

   CustomBufferSystem bs( mpi::MPIManager::instance()->comm() );
   bs.useIProbe(iProbe);

   for( uint_t i =0; i < iterations; ++i )
   {
      WALBERLA_MPI_BARRIER();
      tp["pack"].start();
      for (auto dirIt = Stencil::beginNoCenter(); dirIt != Stencil::end(); ++dirIt)
      {
         auto recvRank = mpiInfo.getNeighborRank( *dirIt );
         if (recvRank == -1) continue;
         bs.sendBuffer(recvRank) << sendBuf;
         WALBERLA_ASSERT_EQUAL(bs.sendBuffer(recvRank).size(), messageSize + sizeof(size_t));
      }
      tp["pack"].end();

      WALBERLA_MPI_BARRIER();
      tp["communicate"].start();
      bs.setReceiverInfoFromSendBufferState(false, true);
      bs.sendAll();
      for( auto it = bs.begin(); it != bs.end(); ++it )
      {
         WALBERLA_ASSERT_EQUAL(it.buffer().size(), messageSize + sizeof(size_t));
         WALBERLA_ASSERT_EQUAL(recvBuf.size(), messageSize);
      }
      tp["communicate"].end();

      WALBERLA_MPI_BARRIER();
      tp["unpack"].start();
      auto& recvInfos = bs.getRecvInfos();
      for (auto recvIt = recvInfos.begin(); recvIt != recvInfos.end(); ++recvIt)
      {
         auto& rb = recvIt->second.buffer;
         rb >> recvBuf;
         WALBERLA_ASSERT(rb.isEmpty());
      }
      tp["unpack"].end();
   }
}

std::string envToString(const char* env)
{
   return env != nullptr ? std::string(env) : "";
}

int main( int argc, char ** argv )
{
   mpi::Environment mpiEnv( argc, argv );

   if ( argc != 10 )
   {
      WALBERLA_ROOT_SECTION()
      {
         std::cout << "Usage ./probeVsExtraMessage x y z px py pz iterations messageSize stencil " << std::endl;
         std::cout << std::endl;
         std::cout << "x: number of processes in x direction" << std::endl;
         std::cout << "y: number of processes in y direction" << std::endl;
         std::cout << "z: number of processes in z direction" << std::endl;
         std::cout << "px: periodic in x direction?" << std::endl;
         std::cout << "py: periodic in y direction?" << std::endl;
         std::cout << "pz: periodic in z direction?" << std::endl;
         std::cout << "iterations: number of communications per case" << std::endl;
         std::cout << "messageSize: size of the SendBuffer in bytes" << std::endl;
         std::cout << "stencil: communication stencil (D3Q27, D3Q19, D3Q7)" << std::endl;

      }
      return EXIT_FAILURE;
   }

   Vector3<uint_t> procs;
   procs[0] = std::stoul(argv[1]);
   procs[1] = std::stoul(argv[2]);
   procs[2] = std::stoul(argv[3]);
   WALBERLA_CHECK_EQUAL(procs[0]*procs[1]*procs[2], mpi::MPIManager::instance()->numProcesses());
   Vector3<bool> periodicity;
   periodicity[0] = std::stoul(argv[4]);
   periodicity[1] = std::stoul(argv[5]);
   periodicity[2] = std::stoul(argv[6]);

   uint_t iterations  = std::stoul(argv[7]);
   uint_t messageSize = std::stoul(argv[8]);
   std::string stencil = argv[9];
   WALBERLA_CHECK_EQUAL(stencil, "D3Q27", "only D3Q27 is supported!");

   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(procs);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(periodicity);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(iterations);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(messageSize);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(stencil);

   MPIInfo mpiInfo(procs, periodicity);

   WcTimingPool tp_twoMessages;
   WcTimingPool tp_probe;
   WALBERLA_MPI_BARRIER();
   if (stencil == "D3Q27")
   {
      communicate<stencil::D3Q27>(mpiInfo, iterations, messageSize, false, tp_twoMessages);
      communicate<stencil::D3Q27>(mpiInfo, iterations, messageSize, true, tp_probe);
   } else if (stencil == "D3Q19")
   {
      communicate<stencil::D3Q19>(mpiInfo, iterations, messageSize, false, tp_twoMessages);
      communicate<stencil::D3Q19>(mpiInfo, iterations, messageSize, true, tp_probe);
   } else if (stencil == "D3Q7")
   {
      communicate<stencil::D3Q7>(mpiInfo, iterations, messageSize, false, tp_twoMessages);
      communicate<stencil::D3Q7>(mpiInfo, iterations, messageSize, true, tp_probe);
   } else
   {
      WALBERLA_ABORT("stencil not supported: " << stencil);
   }
   WALBERLA_LOG_INFO_ON_ROOT(tp_twoMessages);
   WALBERLA_LOG_INFO_ON_ROOT(tp_probe);

   WALBERLA_ROOT_SECTION()
   {
      std::map< std::string, walberla::int64_t > integerProperties;
      std::map< std::string, double >            realProperties;
      std::map< std::string, std::string >       stringProperties;

      integerProperties["procs_x"]      = int64_c(procs[0]);
      integerProperties["procs_y"]      = int64_c(procs[1]);
      integerProperties["procs_z"]      = int64_c(procs[2]);
      integerProperties["priodicity_x"] = int64_c(periodicity[0]);
      integerProperties["priodicity_y"] = int64_c(periodicity[1]);
      integerProperties["priodicity_z"] = int64_c(periodicity[2]);
      integerProperties["iterations"]   = int64_c(iterations);
      integerProperties["messageSize"]  = int64_c(messageSize);
      stringProperties["stencil"]       = stencil;
      stringProperties["SLURM_CLUSTER_NAME"]       = envToString(std::getenv( "SLURM_CLUSTER_NAME" ));
      stringProperties["SLURM_CPUS_ON_NODE"]       = envToString(std::getenv( "SLURM_CPUS_ON_NODE" ));
      stringProperties["SLURM_CPUS_PER_TASK"]      = envToString(std::getenv( "SLURM_CPUS_PER_TASK" ));
      stringProperties["SLURM_JOB_ACCOUNT"]        = envToString(std::getenv( "SLURM_JOB_ACCOUNT" ));
      stringProperties["SLURM_JOB_ID"]             = envToString(std::getenv( "SLURM_JOB_ID" ));
      stringProperties["SLURM_JOB_CPUS_PER_NODE"]  = envToString(std::getenv( "SLURM_JOB_CPUS_PER_NODE" ));
      stringProperties["SLURM_JOB_NAME"]           = envToString(std::getenv( "SLURM_JOB_NAME" ));
      stringProperties["SLURM_JOB_NUM_NODES"]      = envToString(std::getenv( "SLURM_JOB_NUM_NODES" ));
      stringProperties["SLURM_NTASKS"]             = envToString(std::getenv( "SLURM_NTASKS" ));
      stringProperties["SLURM_NTASKS_PER_CORE"]    = envToString(std::getenv( "SLURM_NTASKS_PER_CORE" ));
      stringProperties["SLURM_NTASKS_PER_NODE"]    = envToString(std::getenv( "SLURM_NTASKS_PER_NODE" ));
      stringProperties["SLURM_NTASKS_PER_SOCKET"]  = envToString(std::getenv( "SLURM_NTASKS_PER_SOCKET" ));
      stringProperties["SLURM_TASKS_PER_NODE"]     = envToString(std::getenv( "SLURM_TASKS_PER_NODE" ));

      auto runId = sqlite::storeRunInSqliteDB( "ProbeVsTwoMessages.sqlite", integerProperties, stringProperties, realProperties );
      sqlite::storeTimingPoolInSqliteDB( "ProbeVsTwoMessages.sqlite", runId, tp_twoMessages, "twoMessages" );
      sqlite::storeTimingPoolInSqliteDB( "ProbeVsTwoMessages.sqlite", runId, tp_probe, "probe" );
   }

   return 0;
}
} // namespace walberla

int main( int argc, char* argv[] )
{
   return walberla::main( argc, argv );
}
