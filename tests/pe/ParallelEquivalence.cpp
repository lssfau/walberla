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
//! \file ParallelEquivalence.cpp
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================


#include "blockforest/all.h"
#include "core/all.h"
#include "domain_decomposition/all.h"

#include "pe/basic.h"
#include "pe/ccd/SimpleCCDDataHandling.h"
#include "pe/synchronization/SyncNextNeighbors.h"

#include "CheckVitalParameters.h"

#include "core/debug/TestSubsystem.h"

#include <tuple>

#include <algorithm>
#include <vector>

namespace walberla {
using namespace walberla::pe;

int runs = 1000;

typedef std::tuple<Sphere, Plane> BodyTuple ;

struct BodyData
{
   BodyData(const walberla::id_t uid, const Vec3& pos, const Vec3& linVel)
      : uid_(uid)
      , pos_(pos)
      , linVel_(linVel)
   {}
   BodyData(mpi::RecvBuffer& recv)
   {
      recv >> uid_;
      recv >> pos_;
      recv >> linVel_;
   }

   void pack(mpi::SendBuffer& buf)
   {
      buf << uid_;
      buf << pos_;
      buf << linVel_;
   }

   walberla::id_t uid_;
   Vec3 pos_;
   Vec3 linVel_;
};

void checkBodyData(BodyData& a, BodyData& b)
{
   WALBERLA_CHECK_EQUAL(a.uid_, b.uid_);
   WALBERLA_CHECK_FLOAT_EQUAL(a.pos_, b.pos_, "uid: " << a.uid_ << "\nuid: " << b.uid_);
   WALBERLA_CHECK_FLOAT_EQUAL(a.linVel_, b.linVel_, "uid: " << a.uid_ << "\nuid: " << b.uid_);
}

bool comp(const BodyData& a, const BodyData& b)
{
   return a.uid_ < b.uid_;
}

void sim(shared_ptr< StructuredBlockForest > forest, std::vector<BodyData>& res, bool useHashGrids = true)
{
   res.clear();

   MaterialID granular = Material::find( "granular" );

   shared_ptr<BodyStorage> globalStorage = make_shared<BodyStorage> ();

//   WALBERLA_LOG_DEVEL("process: " << mpi::MPIManager::instance()->rank() << "\nnumber of blocks: " << forest->size());

   SetBodyTypeIDs<BodyTuple>::execute();

   auto storageID_           = forest->addBlockData(createStorageDataHandling<BodyTuple>(), "Storage");
   BlockDataID ccdID_;
   if (useHashGrids)
   {
      ccdID_                 = forest->addBlockData(ccd::createHashGridsDataHandling( globalStorage, storageID_ ), "CCD");
   } else
   {
      ccdID_                 = forest->addBlockData(ccd::createSimpleCCDDataHandling( globalStorage, storageID_ ), "CCD");
   }
   auto fcdID_               = forest->addBlockData(fcd::createGenericFCDDataHandling<BodyTuple, fcd::AnalyticCollideFunctor>(), "FCD");

   cr::DEM cr( globalStorage, forest->getBlockStoragePointer(), storageID_, ccdID_, fcdID_, nullptr);

//   auto vtkOutput_   = make_shared<DefaultBodyVTKOutput>(storageID_, *forest) ;
//   auto vtkWriter_   = vtk::createVTKOutput_PointData(vtkOutput_, "Bodies", 1);

   const int dim = 10;

   pe::createPlane( *globalStorage, 0, Vec3(0, 0, +1), Vec3(5,5, 0), granular);
   pe::createPlane( *globalStorage, 0, Vec3(0, 0, -1), Vec3(5,5,10), granular);

   pe::createPlane( *globalStorage, 0, Vec3(0, +1, 0), Vec3(5, 0,5), granular);
   pe::createPlane( *globalStorage, 0, Vec3(0, -1, 0), Vec3(5,10,5), granular);

   pe::createPlane( *globalStorage, 0, Vec3(+1, 0, 0), Vec3( 0,5,5), granular);
   pe::createPlane( *globalStorage, 0, Vec3(-1, 0, 0), Vec3(10,5,5), granular);

   walberla::id_t counter = 0;

   const real_t dv = real_c(0.1);

   math::seedRandomGenerator(1337);
   for (int z = 0; z < dim; ++z)
      for (int y = 0; y < dim; ++y)
         for (int x = 0; x < dim; ++x)
         {
            SphereID sp = pe::createSphere( *globalStorage, forest->getBlockStorage(), storageID_,
                              ++counter, Vec3(real_c(x) + real_c(0.5), real_c(y) + real_c(0.5), real_c(z) + real_c(0.5)), real_c(0.3));
            Vec3 randVel = Vec3(math::realRandom<real_t>(-dv, dv), math::realRandom<real_t>(-dv, dv), math::realRandom<real_t>(-dv, dv));
            if (sp != nullptr) sp->setLinearVel(randVel);
         }

   WALBERLA_CHECK_EQUAL(globalStorage->size(), 6);

   syncNextNeighbors<BodyTuple>( forest->getBlockForest(), storageID_);
   syncNextNeighbors<BodyTuple>( forest->getBlockForest(), storageID_);
   syncNextNeighbors<BodyTuple>( forest->getBlockForest(), storageID_);

   for (int i = 0; i < runs; ++i)
   {
//      WALBERLA_LOG_DEVEL_ON_ROOT("step: " << i);
      cr.timestep( real_c(0.01) );
      syncNextNeighbors<BodyTuple>( forest->getBlockForest(), storageID_);
//      if ((i%10) == 0)
//         vtkWriter_->write(true);
   }

   for (auto it = forest->begin(); it != forest->end(); ++it){
      IBlock & currentBlock = *it;
      Storage * storage = currentBlock.getData< Storage >( storageID_ );
      BodyStorage& localStorage = (*storage)[0];
      for (auto bodyIt = localStorage.begin(); bodyIt != localStorage.end(); ++bodyIt)
      {
         BodyID b = bodyIt.getBodyID();
         res.emplace_back(b->getID(), b->getPosition(), b->getLinearVel());
      }
   }
   mpi::SendBuffer sendBuf;
   mpi::RecvBuffer recvBuf;
   for (unsigned long i = 0; i < res.size(); ++i)
   {
      res.at(i).pack(sendBuf);
   }
//   WALBERLA_LOG_DEVEL("on process: " << res.size() << "\t" << sendBuf.size());
   mpi::gathervBuffer(sendBuf, recvBuf);

   res.clear();
   while (!recvBuf.isEmpty())
   {
      res.emplace_back(recvBuf );
   }

   forest.reset();
}



uint_t balance( SetupBlockForest & forest, const uint_t, const memory_t )
{
   for( auto block = forest.begin(); block != forest.end(); ++block )
      block->assignTargetProcess( uint_t(0) );

   return uint_t(1);
}



static shared_ptr< StructuredBlockForest > createBlockStructure( const uint_t numberOfXBlocks,        const uint_t numberOfYBlocks,        const uint_t numberOfZBlocks,
                                                                 const uint_t numberOfXCellsPerBlock, const uint_t numberOfYCellsPerBlock, const uint_t numberOfZCellsPerBlock,
                                                                 const real_t dx,
                                                                 const bool   xPeriodic = false,      const bool   yPeriodic = false,      const bool   zPeriodic = false,
                                                                 const bool keepGlobalBlockInformation = false )
{
   // initialize SetupBlockForest = determine domain decomposition
   SetupBlockForest sforest;

   sforest.init( math::AABB( real_c(0), real_c(0), real_c(0), dx * real_c( numberOfXBlocks * numberOfXCellsPerBlock ),
                                                              dx * real_c( numberOfYBlocks * numberOfYCellsPerBlock ),
                                                              dx * real_c( numberOfZBlocks * numberOfZCellsPerBlock ) ),
                 numberOfXBlocks, numberOfYBlocks, numberOfZBlocks, xPeriodic, yPeriodic, zPeriodic );

   // calculate process distribution
   const memory_t memoryLimit = math::Limits< memory_t >::inf();

   sforest.balanceLoad( balance, uint_c( MPIManager::instance()->numProcesses() ), real_t(0), memoryLimit, true );

   WALBERLA_LOG_INFO_ON_ROOT( sforest );

   MPIManager::instance()->useWorldComm();

   // create StructuredBlockForest (encapsulates a newly created BlockForest)
   shared_ptr< StructuredBlockForest > sbf =
         make_shared< StructuredBlockForest >( make_shared< BlockForest >( uint_c( MPIManager::instance()->rank() ), sforest, keepGlobalBlockInformation ),
                                               numberOfXCellsPerBlock, numberOfYCellsPerBlock, numberOfZCellsPerBlock );
   sbf->createCellBoundingBoxes();
   return sbf;
}

int main( int argc, char ** argv )
{
   debug::enterTestMode();

   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );

   for( int i = 1; i < argc; ++i )
      if( std::strcmp( argv[i], "-runs" ) == 0 )
      {
         runs = std::atoi(argv[i+1]);
         WALBERLA_LOG_DEVEL_ON_ROOT("runs: " << runs);
      }

   createMaterial( "granular", real_c(7830.0), real_c(0.25), real_c(0.4), real_c(0.4), real_c(0.35), real_c(1.39e11), real_c(5.18e7), real_c(1.07e2), real_c(1.07e2) );

//   logging::Logging::instance()->setFileLogLevel( logging::Logging::DETAIL );
//   logging::Logging::instance()->includeLoggingToFile("SyncLog");

   std::vector<BodyData> OneBlockOneProcess;
   std::vector<BodyData> ManyBlocksOneProcess;
   std::vector<BodyData> ManyBlocksManyProcesses;
   std::vector<BodyData> ManyBlocksManyProcessesSimple;

   sim(
      createBlockStructure(
               uint_c( 1), uint_c( 1), uint_c( 1), // number of blocks in x,y,z direction
               uint_c( 1), uint_c( 1), uint_c( 1), // how many cells per block (x,y,z)
               real_c(10),                         // dx: length of one cell in physical coordinates
               false, false, false ),              // no periodicity
      OneBlockOneProcess);

   MPIManager::instance()->resetMPI();

   sim(
      createBlockStructure(
               uint_c( 2), uint_c( 2), uint_c( 2), // number of blocks in x,y,z direction
               uint_c( 1), uint_c( 1), uint_c( 1), // how many cells per block (x,y,z)
               real_c( 5),                         // dx: length of one cell in physical coordinates
               false, false, false ),              // no periodicity
      ManyBlocksOneProcess);

   sim(
      blockforest::createUniformBlockGrid(
               uint_c( 2), uint_c( 2), uint_c( 2), // number of blocks in x,y,z direction
               uint_c( 1), uint_c( 1), uint_c( 1), // how many cells per block (x,y,z)
               real_c( 5),                         // dx: length of one cell in physical coordinates
               0,                                  // max blocks per process
               false, false,                       // include metis / force metis
               false, false, false ),              // no periodicity
      ManyBlocksManyProcesses);

   sim(
      blockforest::createUniformBlockGrid(
               uint_c( 2), uint_c( 2), uint_c( 2), // number of blocks in x,y,z direction
               uint_c( 1), uint_c( 1), uint_c( 1), // how many cells per block (x,y,z)
               real_c( 5),                         // dx: length of one cell in physical coordinates
               0,                                  // max blocks per process
               false, false,                       // include metis / force metis
               false, false, false ),              // no periodicity
      ManyBlocksManyProcessesSimple,
            false);

   WALBERLA_ROOT_SECTION()
   {
      WALBERLA_LOG_DEVEL_ON_ROOT("OneBlockOneProcess: " << OneBlockOneProcess.size());
      WALBERLA_LOG_DEVEL_ON_ROOT("ManyBlocksOneProcess: " << ManyBlocksOneProcess.size());
      WALBERLA_LOG_DEVEL_ON_ROOT("ManyBlocksManyProcesses: " << ManyBlocksManyProcesses.size());
      WALBERLA_LOG_DEVEL_ON_ROOT("ManyBlocksManyProcessesSimple: " << ManyBlocksManyProcessesSimple.size());

      std::sort(OneBlockOneProcess.begin(), OneBlockOneProcess.end(), comp);
      std::sort(ManyBlocksOneProcess.begin(), ManyBlocksOneProcess.end(), comp);
      std::sort(ManyBlocksManyProcesses.begin(), ManyBlocksManyProcesses.end(), comp);
      std::sort(ManyBlocksManyProcessesSimple.begin(), ManyBlocksManyProcessesSimple.end(), comp);

      WALBERLA_CHECK_EQUAL(OneBlockOneProcess.size(), ManyBlocksOneProcess.size());
      WALBERLA_CHECK_EQUAL(OneBlockOneProcess.size(), ManyBlocksManyProcesses.size());
      WALBERLA_CHECK_EQUAL(OneBlockOneProcess.size(), ManyBlocksManyProcessesSimple.size());

      WALBERLA_LOG_DEVEL("checking OneBlockOneProcess against ManyBlocksOneProcess");
      for (unsigned long i = 0; i < OneBlockOneProcess.size(); ++i)
      {
         checkBodyData(OneBlockOneProcess.at(i), ManyBlocksOneProcess.at(i));
      }

      WALBERLA_LOG_DEVEL("checking OneBlockOneProcess against ManyBlocksManyProcesses");
      for (unsigned long i = 0; i < OneBlockOneProcess.size(); ++i)
      {
         checkBodyData(OneBlockOneProcess.at(i), ManyBlocksManyProcesses.at(i));
      }

      WALBERLA_LOG_DEVEL("checking OneBlockOneProcess against ManyBlocksManyProcessesSimple");
      for (unsigned long i = 0; i < OneBlockOneProcess.size(); ++i)
      {
         checkBodyData(OneBlockOneProcess.at(i), ManyBlocksManyProcessesSimple.at(i));
      }
   }

   return EXIT_SUCCESS;
}
} // namespace walberla

int main( int argc, char* argv[] )
{
  return walberla::main( argc, argv );
}