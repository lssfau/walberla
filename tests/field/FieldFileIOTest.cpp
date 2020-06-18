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
//! \file FieldFileIOTest.cpp
//! \ingroup field
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//
//======================================================================================================================

#include "blockforest/SetupBlockForest.h"
#include "blockforest/StructuredBlockForest.h"
#include "blockforest/loadbalancing/StaticCurve.h"

#include "core/debug/TestSubsystem.h"
#include "core/math/IntegerFactorization.h"
#include "core/math/Random.h"
#include "core/mpi/Environment.h"
#include "core/stringToNum.h"
#include "core/timing/Timer.h"

#include "field/AddToStorage.h"
#include "field/Field.h"
#include "field/FileIO.h"

namespace mpi_file_io_test
{
using namespace walberla;
using walberla::uint8_t;

static void refinementSelectionFunction(SetupBlockForest& forest)
{
   const uint_t numRootBlocks = forest.getNumberOfRootBlocks();
   for (uint_t i = 0; i < numRootBlocks; i += 8)
   {
      SetupBlock* block = forest.getRootBlock(i);

      if (!block->hasChildren()) block->setMarker(true);
   }
}

static void workloadMemorySUIDAssignmentFunction(SetupBlockForest& forest)
{
   std::vector< SetupBlock* > blocks;
   forest.getBlocks(blocks);

   for (uint_t i = 0; i != blocks.size(); ++i)
   {
      blocks[i]->setMemory(1.0);
      blocks[i]->setWorkload(1.0);
   }
}

void test(int argc, char* argv[])
{
   typedef field::GhostLayerField< double, 3 > FieldType;

   std::vector< std::string > args(argv, argv + argc);

   uint_t numBlocks  = 8;
   uint_t xBlockSize = 3;
   uint_t yBlockSize = 5;
   uint_t zBlockSize = 7;

   if (args.size() == 5)
   {
      numBlocks  = stringToNum< uint_t >(args[1]);
      xBlockSize = stringToNum< uint_t >(args[2]);
      yBlockSize = stringToNum< uint_t >(args[3]);
      zBlockSize = stringToNum< uint_t >(args[4]);
   }
   else if (args.size() > 5)
   {
      WALBERLA_ABORT("USAGE:\n\n"
                     << args[0] << " <NUMBER_OF_COARSE_BLOCKS> <X_BLOCK_SIZE> <Y_BLOCK_SIZE> <Z_BLOCK_SIZE>");
   }

   SetupBlockForest sforest;

   sforest.addRefinementSelectionFunction(refinementSelectionFunction);
   sforest.addWorkloadMemorySUIDAssignmentFunction(workloadMemorySUIDAssignmentFunction);

   AABB domain(0, 0, 0, 100, 100, 100);

   auto factors = math::getFactors3D(numBlocks);

   sforest.init(domain, factors[0], factors[1], factors[2], true, false, false);

   sforest.balanceLoad(blockforest::StaticLevelwiseCurveBalance(true), uint_c(MPIManager::instance()->numProcesses()));

   auto sbf = make_shared< StructuredBlockForest >(
      make_shared< BlockForest >(uint_c(MPIManager::instance()->rank()), sforest, true), xBlockSize, yBlockSize,
      zBlockSize);

   auto originalFieldId = field::addToStorage< FieldType >(sbf, "OriginalField");
   auto readFieldId     = field::addToStorage< FieldType >(sbf, "ReadField");

   math::seedRandomGenerator(numeric_cast< std::mt19937::result_type >(MPIManager::instance()->rank()));

   for (auto it = sbf->begin(); it != sbf->end(); ++it)
   {
      auto field = it->getData< FieldType >(originalFieldId);

      for (auto dataIt = field->begin(); dataIt != field->end(); ++dataIt)
         *dataIt = math::realRandom< FieldType::value_type >();
   }

   WcTimer timer;

   WALBERLA_MPI_BARRIER();
   timer.start();
   field::writeToFile< FieldType >("mpiFile.wlb", sbf->getBlockStorage(), originalFieldId);
   WALBERLA_MPI_BARRIER();
   timer.end();
   WALBERLA_LOG_INFO_ON_ROOT("Writing took " << timer.last() << "s");

   WALBERLA_MPI_BARRIER();
   timer.start();
   field::readFromFile< FieldType >("mpiFile.wlb", sbf->getBlockStorage(), readFieldId);
   WALBERLA_MPI_BARRIER();
   timer.end();
   WALBERLA_LOG_INFO_ON_ROOT("Reading took " << timer.last() << "s");

   for (auto it = sbf->begin(); it != sbf->end(); ++it)
   {
      auto originalField = it->getData< FieldType >(originalFieldId);
      auto readField     = it->getData< FieldType >(readFieldId);

      auto readIt = readField->begin();
      for (auto origIt = originalField->begin(); origIt != originalField->end(); ++origIt, ++readIt)
         WALBERLA_CHECK_IDENTICAL(*origIt, *readIt);
   }
}

} // namespace mpi_file_io_test

int main(int argc, char* argv[])
{
   walberla::debug::enterTestMode();

   walberla::mpi::Environment mpiEnv(argc, argv);

   // test with MPI_WORLD_COMM
   walberla::MPIManager::instance()->useWorldComm();
   mpi_file_io_test::test(argc, argv);

   // test with Cartesian MPI communicator
   // this is tested additionally because some versions of OpenMPI are known to produce segmentation faults when using
   // MPI-IO with a 3D Cartesian MPI communicator; for those OpenMPI versions, serial I/O is used instead
   if (walberla::MPIManager::instance()->numProcesses() == 16)
   {
      walberla::MPIManager::instance()->resetMPI();

      walberla::MPIManager::instance()->createCartesianComm(walberla::uint_c(4), walberla::uint_c(2),
                                                            walberla::uint_c(2), false, false, false);
      mpi_file_io_test::test(argc, argv);
   }
}
