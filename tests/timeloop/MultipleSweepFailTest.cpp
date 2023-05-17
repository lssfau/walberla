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
//! \file MultipleSweepFailTest.cpp
//! \ingroup timeloop
//! \author Markus Holzer <markus.holzer@fau.de>
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//! \brief Test if multiple sweeps (with the same selector) can not be added to the timeloop at once.
//!
//! THIS TEST MUST FAIL!
//
//======================================================================================================================

#include "blockforest/Initialization.h"

#include "core/DataTypes.h"
#include "core/debug/TestSubsystem.h"
#include "core/mpi/Environment.h"

#include "timeloop/SweepTimeloop.h"

namespace walberla
{
namespace timeloop
{
namespace MultipleSweepFailTest
{

int main(int argc, char** argv)
{
   debug::enterTestMode();
   mpi::Environment const env(argc, argv);

   const std::shared_ptr< StructuredBlockForest > blockForest = blockforest::createUniformBlockGrid(
      uint_c(1), uint_c(1), uint_c(1), uint_c(1), uint_c(1), uint_c(1), real_c(1), false, false, false, false);

   SweepTimeloop timeloop(blockForest, uint_c(1));

   // empty sweep that does nothing
   const auto emptySweep = [](IBlock*) {};

   // this must fail, as two sweeps are added at once with the same selectors (default selectors: Set<SUID>::emptySet())
   timeloop.add() << Sweep(emptySweep, "Sweep 1") << Sweep(emptySweep, "Sweep 2");
   timeloop.singleStep();

   return EXIT_SUCCESS;
}
} // namespace MultipleSweepFailTest
} // namespace timeloop
} // namespace walberla

int main(int argc, char** argv) { return walberla::timeloop::MultipleSweepFailTest::main(argc, argv); }
