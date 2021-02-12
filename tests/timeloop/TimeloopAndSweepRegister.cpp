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
//! \file TimeloopAndSweepRegister.cpp
//! \ingroup timeloop
//! \author Markus Holzer <markus.holzer@fau.de>
//! \brief test cases that test the registering of Sweeps at timeloop
//
//======================================================================================================================

#include "blockforest/Initialization.h"

#include "core/DataTypes.h"
#include "core/Environment.h"
#include "core/debug/TestSubsystem.h"

#include "field/Field.h"

#include "timeloop/SweepTimeloop.h"

#include <string>
#include <vector>

using namespace walberla;

using Field_T = Field< uint_t, 1 >;

auto FieldAdder = [](IBlock* const block, StructuredBlockStorage* const storage) {
   return new Field_T(storage->getNumberOfXCells(*block), storage->getNumberOfYCells(*block),
                      storage->getNumberOfZCells(*block), uint_t(0.0), field::fzyx,
                      make_shared< field::AllocateAligned< uint_t, 64 > >());
};

class Sweep1
{
 public:
   Sweep1(BlockDataID fieldID) : fieldID_(fieldID) {}

   void operator()(IBlock* block)
   {
      auto field = block->getData< Field_T >(fieldID_);

      for (auto iter = field->begin(); iter != field->end(); ++iter)
         *iter += 1;
   }

 private:
   BlockDataID fieldID_;
};

class Sweep2
{
 public:
   Sweep2(BlockDataID fieldID) : fieldID_(fieldID) {}

   void operator()(IBlock* block)
   {
      auto field = block->getData< Field_T >(fieldID_);

      for (auto iter = field->begin(); iter != field->end(); ++iter)
         *iter += 2;
   }

 private:
   BlockDataID fieldID_;
};

int main(int argc, char** argv)
{
   debug::enterTestMode();
   mpi::Environment env(argc, argv);

   std::vector<std::string> expectedSequence;
   std::vector<std::string> sequence;

   SUID sweepSelect1("Sweep1");
   SUID sweepSelect2("Sweep2");

   shared_ptr< StructuredBlockForest > blocks = blockforest::createUniformBlockGrid(
      uint_c(4), uint_c(2), uint_c(2), uint_c(10), uint_c(10), uint_c(10), real_c(1), false, false, false, false);

   BlockDataID fieldID = blocks->addStructuredBlockData< Field_T >(FieldAdder, "Test Field");

   for (auto& block : *blocks)
   {
      if (block.getAABB().min()[0] < 20)
         block.setState(sweepSelect1);
      else
         block.setState(sweepSelect2);
   }

   uint_t timesteps = 10;
   SweepTimeloop timeloop(blocks->getBlockStorage(), timesteps);

   timeloop.add() << Sweep(Sweep1(fieldID), "Sweep 1", sweepSelect1, sweepSelect2);
   timeloop.add() << Sweep(Sweep2(fieldID), "Sweep 2", sweepSelect2, sweepSelect1);

   WcTimingPool timingPool;

   timeloop.run(timingPool);
   for (auto& block : *blocks)
   {
      auto field = block.getData< Field_T >(fieldID);
      if (block.getAABB().min()[0] < 20)
      {
         for (auto iter = field->begin(); iter != field->end(); ++iter)
            WALBERLA_CHECK_EQUAL(*iter, timesteps)
      }
      else
      {
         for (auto iter = field->begin(); iter != field->end(); ++iter)
            WALBERLA_CHECK_EQUAL(*iter, timesteps * 2)
      }
   }
   
   return EXIT_SUCCESS;
}
