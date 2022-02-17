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
//! \file TimeloopAndSweepRegisterTest.cpp
//! \ingroup timeloop
//! \author Markus Holzer <markus.holzer@fau.de>
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//! \brief Test if registering and execution of sweeps in a timeloop with different selectors works correctly.
//
//======================================================================================================================

#include "blockforest/Initialization.h"

#include "core/DataTypes.h"
#include "core/Environment.h"
#include "core/debug/TestSubsystem.h"

#include "field/Field.h"

#include "timeloop/SweepTimeloop.h"

#include <vector>

namespace walberla
{
namespace timeloop
{
namespace TimeloopAndSweepRegisterTest
{

using Field_T = Field< uint_t, 1 >;

auto FieldAdder = [](IBlock* const block, StructuredBlockStorage* const storage) {
   return new Field_T(storage->getNumberOfXCells(*block), storage->getNumberOfYCells(*block),
                      storage->getNumberOfZCells(*block), uint_t(0), field::fzyx,
                      make_shared< field::AllocateAligned< uint_t, uint_t(64) > >());
};

class Sweep1
{
 public:
   Sweep1(BlockDataID fieldID) : fieldID_(fieldID) {}

   void operator()(IBlock* const block)
   {
      Field_T* field = block->getData< Field_T >(fieldID_);

      WALBERLA_FOR_ALL_CELLS(fieldIt, field, { *fieldIt += uint_c(1); }) // WALBERLA_FOR_ALL_CELLS
   }

 private:
   BlockDataID fieldID_;
}; // class Sweep1

class Sweep2
{
 public:
   Sweep2(BlockDataID fieldID) : fieldID_(fieldID) {}

   void operator()(IBlock* const block)
   {
      Field_T* field = block->getData< Field_T >(fieldID_);
      WALBERLA_FOR_ALL_CELLS(fieldIt, field, { *fieldIt += uint_c(2); }) // WALBERLA_FOR_ALL_CELLS
   }

 private:
   BlockDataID fieldID_;
}; // class Sweep2

int main(int argc, char** argv)
{
   debug::enterTestMode();
   mpi::Environment env(argc, argv);

   std::vector< std::string > expectedSequence;
   std::vector< std::string > sequence;

   const SUID sweepSelect1("Sweep1");
   const SUID sweepSelect2("Sweep2");

   const shared_ptr< StructuredBlockForest > blockForest = blockforest::createUniformBlockGrid(
      uint_c(4), uint_c(2), uint_c(2), uint_c(10), uint_c(10), uint_c(10), real_c(1), false, false, false, false);

   const BlockDataID fieldID = blockForest->addStructuredBlockData< Field_T >(FieldAdder, "Test Field");

   for (auto& block : *blockForest)
   {
      if (block.getAABB().min()[0] < real_c(20)) { block.setState(sweepSelect1); }
      else
      {
         block.setState(sweepSelect2);
      }
   }

   const uint_t timesteps = uint_c(10);
   SweepTimeloop timeloop(blockForest, timesteps);

   timeloop.add() << Sweep(Sweep1(fieldID), "Sweep 1", sweepSelect1, sweepSelect2);
   timeloop.add() << Sweep(Sweep2(fieldID), "Sweep 2", sweepSelect2, sweepSelect1);

   timeloop.run();

   for (const auto& block : *blockForest)
   {
      const Field_T* field = block.getData< Field_T >(fieldID);

      if (block.getAABB().min()[0] < real_c(20))
      {
         WALBERLA_FOR_ALL_CELLS(fieldIt, field,
                                { WALBERLA_CHECK_EQUAL(*fieldIt, timesteps); }) // WALBERLA_FOR_ALL_CELLS
      }
      else
      {
         WALBERLA_FOR_ALL_CELLS(fieldIt, field,
                                { WALBERLA_CHECK_EQUAL(*fieldIt, timesteps * uint_c(2)); }) // WALBERLA_FOR_ALL_CELLS
      }
   }

   return EXIT_SUCCESS;
}

} // namespace TimeloopAndSweepRegisterTest
} // namespace timeloop
} // namespace walberla

int main(int argc, char** argv) { return walberla::timeloop::TimeloopAndSweepRegisterTest::main(argc, argv); }