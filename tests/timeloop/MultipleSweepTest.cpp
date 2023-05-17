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
//! \file MultipleSweepTest.cpp
//! \ingroup timeloop
//! \author Markus Holzer <markus.holzer@fau.de>
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//! \brief Test if multiple sweeps (with different selectors) can be added to the timeloop at once.
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
namespace MultipleSweepTest
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

   const SUID selector1("selector1");
   const SUID selector2("selector2");
   const SUID selector3("selector3");
   const SUID selector4("selector4");

   // this must not fail, as two sweeps with the different (required) selectors are added
   timeloop.add() << Sweep(emptySweep, "Sweep 1", selector1, Set< SUID >::emptySet())
                  << Sweep(emptySweep, "Sweep 2", selector2, Set< SUID >::emptySet());
   timeloop.add() << Sweep(emptySweep, "Sweep 3", selector3, selector4)
                  << Sweep(emptySweep, "Sweep 4", selector4, selector3);

   timeloop.singleStep();

   return EXIT_SUCCESS;
}
} // namespace MultipleSweepTest
} // namespace timeloop
} // namespace walberla

int main(int argc, char** argv) { return walberla::timeloop::MultipleSweepTest::main(argc, argv); }
