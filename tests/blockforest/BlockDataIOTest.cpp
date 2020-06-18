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
//! \file BlockDataIOTest.cpp
//! \ingroup field
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//
//======================================================================================================================

#include "blockforest/BlockForestEvaluation.h"
#include "blockforest/SetupBlockForest.h"
#include "blockforest/StructuredBlockForest.h"
#include "blockforest/loadbalancing/StaticCurve.h"

#include "core/debug/TestSubsystem.h"
#include "core/math/Random.h"
#include "core/mpi/Environment.h"

#include "field/AddToStorage.h"
#include "field/Field.h"

namespace block_data_io_test
{
using namespace walberla;
using walberla::uint8_t;

const SUID Empty("empty");
const Set< SUID > None(Set< SUID >::emptySet());

static void refinementSelectionFunction(SetupBlockForest& forest)
{
   for (auto block = forest.begin(); block != forest.end(); ++block)
      if (block->getAABB().contains(Vector3< real_t >(real_t(75))))
         if (!block->hasFather()) block->setMarker(true);
}

static void workloadMemorySUIDAssignmentFunction(SetupBlockForest& forest)
{
   for (auto block = forest.begin(); block != forest.end(); ++block)
   {
      block->setMemory(memory_t(1));
      block->setWorkload(workload_t(1));
      if (block->getAABB().contains(Vector3< real_t >(real_t(25)))) block->addState(Empty);
   }
}

void test()
{
   typedef field::GhostLayerField< double, 2 > FieldType;

   SetupBlockForest sforest;

   sforest.addRefinementSelectionFunction(refinementSelectionFunction);
   sforest.addWorkloadMemorySUIDAssignmentFunction(workloadMemorySUIDAssignmentFunction);

   AABB domain(0, 0, 0, 100, 100, 100);

   sforest.init(domain, uint_t(2), uint_t(2), uint_t(2), true, false, false);

   sforest.balanceLoad(blockforest::StaticLevelwiseCurveBalance(true), uint_c(MPIManager::instance()->numProcesses()));

   auto sbf = make_shared< StructuredBlockForest >(
      make_shared< BlockForest >(uint_c(MPIManager::instance()->rank()), sforest, true), uint_t(10), uint_t(8),
      uint_t(14));

   blockforest::BlockForestEvaluation evaluation(sbf->getBlockForest());
   WALBERLA_LOG_INFO_ON_ROOT("BlockForest:\n" << evaluation.toString());

   // auto originalFieldId = field::addToStorage< FieldType >( sbf, "OriginalField", 0.0, field::zyxf, uint_t(3), false,
   // None, Empty );
   auto dataHandling    = make_shared< field::DefaultBlockDataHandling< FieldType > >(sbf, uint_t(3), 0.0, field::zyxf);
   auto originalFieldId = sbf->addBlockData(dataHandling, "OriginalField", None, Empty);

   math::seedRandomGenerator(numeric_cast< std::mt19937::result_type >(MPIManager::instance()->rank()));

   for (auto it = sbf->begin(None, Empty); it != sbf->end(); ++it)
   {
      auto field = it->getData< FieldType >(originalFieldId);

      for (auto dataIt = field->begin(); dataIt != field->end(); ++dataIt)
         *dataIt = math::realRandom< FieldType::value_type >();
   }

   sbf->saveBlockData("block.data", originalFieldId);

   WALBERLA_MPI_BARRIER()

   auto readFieldId = sbf->loadBlockData("block.data", dataHandling, "ReadField", None, Empty);

   for (auto it = sbf->begin(None, Empty); it != sbf->end(); ++it)
   {
      auto originalField = it->getData< FieldType >(originalFieldId);
      auto readField     = it->getData< FieldType >(readFieldId);

      auto readIt = readField->begin();
      for (auto origIt = originalField->begin(); origIt != originalField->end(); ++origIt, ++readIt)
         WALBERLA_CHECK_IDENTICAL(*origIt, *readIt);
   }
}

} // namespace block_data_io_test

int main(int argc, char* argv[])
{
   walberla::debug::enterTestMode();

   walberla::mpi::Environment mpiEnv(argc, argv);

   // test with MPI_WORLD_COMM
   walberla::MPIManager::instance()->useWorldComm();
   block_data_io_test::test();

   // test with Cartesian MPI communicator
   // this is tested additionally because some versions of OpenMPI are known to produce segmentation faults when using
   // MPI-IO with a 3D Cartesian MPI communicator; for those OpenMPI versions, serial I/O is used instead
   if (walberla::MPIManager::instance()->numProcesses() == 8)
   {
      walberla::MPIManager::instance()->resetMPI();

      walberla::MPIManager::instance()->createCartesianComm(walberla::uint_c(2), walberla::uint_c(2),
                                                            walberla::uint_c(2), false, false, false);
      block_data_io_test::test();
   }
}
