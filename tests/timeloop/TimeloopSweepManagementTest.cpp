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
//! \file TimeloopSweepManagementTest.cpp
//! \ingroup timeloop
//! \author Markus Holzer <markus.holzer@fau.de>
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//! \brief Test if a sweep in the timeloop works on its members rather can re-initializing them in every time step.
//
//======================================================================================================================

#include "blockforest/Initialization.h"

#include "core/DataTypes.h"
#include "core/Environment.h"
#include "core/debug/TestSubsystem.h"

#include "timeloop/SweepTimeloop.h"

namespace walberla
{
namespace timeloop
{
namespace TimeloopSweepManagementTest
{

class CounterSweep
{
 public:
   CounterSweep(const std::shared_ptr< uint_t >& externalCounter)
      : externalCounter_(externalCounter), internalCounter_(uint_c(0))
   {}

   void operator()(IBlock*)
   {
      ++(*externalCounter_);
      ++internalCounter_;

      WALBERLA_CHECK_EQUAL(*externalCounter_, internalCounter_);
   }

 private:
   std::shared_ptr< uint_t > externalCounter_;
   uint_t internalCounter_;

}; // class CounterSweep

int main(int argc, char** argv)
{
   debug::enterTestMode();
   mpi::Environment env(argc, argv);

   shared_ptr< StructuredBlockForest > blockForest =
      blockforest::createUniformBlockGrid(uint_c(1), uint_c(1), uint_c(1), uint_c(1), uint_c(1), uint_c(1), real_c(1));

   SweepTimeloop timeloop(blockForest, uint_c(10));

   const std::shared_ptr< uint_t > counter = std::make_shared< uint_t >(uint_c(0));

   auto TestSweep = Sweep(CounterSweep(counter), "Counter sweep");
   timeloop.add() << TestSweep;
   timeloop.run();

   return EXIT_SUCCESS;
}
} // namespace TimeloopSweepManagementTest
} // namespace timeloop
} // namespace walberla

int main(int argc, char** argv) { return walberla::timeloop::TimeloopSweepManagementTest::main(argc, argv); }
