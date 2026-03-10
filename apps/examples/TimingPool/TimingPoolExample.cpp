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
//! \file TimingPoolExample.cpp
//! \author Michael Zikeli <michael.zikeli@fau.de>
//
//======================================================================================================================

#include "blockforest/Initialization.h"
#include "blockforest/communication/UniformBufferedScheme.h"

#include "core/DataTypes.h"
#include "core/Environment.h"
#include "core/math/Random.h"
#include "core/timing/RemainingTimeLogger.h"
#include "core/timing/TimingPool.h"

#include "field/AddToStorage.h"
#include "field/communication/PackInfo.h"

#include "stencil/D3Q27.h"

#include "timeloop/SweepTimeloop.h"

#include <chrono>
#include <thread>

namespace walberla::ExampleTimingPool
{

using ScalarField = GhostLayerField< real_t, 1 >;

// function sweep
void myTrivialKernel(IBlock* block, BlockDataID fieldBlockDataID, const real_t delay_factor)
{
   const auto numProcesses = mpi::MPIManager::instance()->numProcesses();
   // This weight controls the delay time and factors in:
   const real_t delay_weight =
      500 * delay_factor                                // emulation for sweeps of different costly executions
      * real_c(mpi::MPIManager::instance()->rank() + 1) // emulation heterogeneous compute power per rank
      * (1.0 / real_c(numProcesses)); // ensure similar total execution independet of number of processes

   // retrieve the field-data from the blockforest
   auto field = block->getData< Field< real_t, 1 > >(fieldBlockDataID);

   const auto threshold = 50.0 * real_c(numProcesses) / real_c(field->allocSize());
   // some bogus "algorithm"
   for (auto cell = field->begin(); cell != field->end(); ++cell)
   {
      // some arbitrary write access to the field
      (*cell) = walberla::math::realRandom(real_t(0), real_t(1));

      // introduce some random delay to simulate non-uniform computation time
      if ((*cell) <= threshold)
      {
         // the delay is proportional to the MPI rank to showcase and test the thread based timing pool reduction.
         const uint_t delay = uint_c(std::round(delay_weight * walberla::math::realRandom(real_t(0), real_t(10))));
         std::this_thread::sleep_for(std::chrono::microseconds(delay));
      }
   }
} // myTrivialKernel()

void run(int argc, char** argv)
{
   Environment env(argc, argv);
   auto config = env.config();

   WALBERLA_LOG_INFO_ON_ROOT(*config);

   ///////////////////////////
   /// BLOCK STORAGE SETUP ///
   ///////////////////////////

   auto blocks = blockforest::createUniformBlockGridFromConfig(config);

   //////////////
   /// FIELDS ///
   //////////////

   BlockDataID fieldId = field::addToStorage< ScalarField >(blocks,      // BlockStorage pointer
                                                            "field",     // field identifier
                                                            real_c(1.0), // initial value
                                                            field::fzyx, // layout type ( fzyx == SoA | xyzf == AoS )
                                                            uint_c(1),   // number of ghost layers
                                                            true);       // force field initialization

   /////////////////////
   /// COMMUNICATION ///
   /////////////////////

   blockforest::communication::UniformBufferedScheme< stencil::D3Q27 > commScheme(blocks);
   commScheme.addPackInfo(make_shared< field::communication::PackInfo< ScalarField > >(fieldId));

   /////////////////////////////
   /// SIMULATION PARAMETERS ///
   /////////////////////////////

   Config::BlockHandle simParams = config->getBlock("Parameters");
   const uint_t timesteps{ simParams.getParameter< uint_t >("timesteps", 10) };
   const real_t delayFactorKernel1{ simParams.getParameter< real_t >("delayFactorKernel1", 1) };
   const real_t delayFactorKernel2{ simParams.getParameter< real_t >("delayFactorKernel2", 1) };
   const real_t delayFactorKernel3{ simParams.getParameter< real_t >("delayFactorKernel3", 1) };

   // Begin-Timing-Pool-Definition
   ////////////////
   /// TIMELOOP ///
   ////////////////

   SweepTimeloop timeloop(blocks, timesteps);

   timeloop.add() << BeforeFunction(commScheme, "Communication")
                  << Sweep([fieldId, delay_factor = delayFactorKernel1](
                              IBlock* block) { myTrivialKernel(block, fieldId, delay_factor); },
                           "Kernel 1");
   timeloop.add() << Sweep([fieldId, delay_factor = delayFactorKernel2](
                              IBlock* block) { myTrivialKernel(block, fieldId, delay_factor); },
                           "Kernel 2")
                  << AfterFunction(commScheme, "Communication");
   timeloop.add() << Sweep(
      [fieldId, delay_factor = delayFactorKernel3](IBlock* block) { myTrivialKernel(block, fieldId, delay_factor); },
      "Kernel 3");
   // End-Timing-Pool-Definition

   // Adds a logger to print the simulation progress.
   auto RemainingTimeLogger = timing::RemainingTimeLogger(timeloop.getNrOfTimeSteps());
   timeloop.addFuncAfterTimeStep(RemainingTimeLogger, "remaining time logger");

   ///////////////////
   /// Measure Run ///
   ///////////////////

   //! [Marker-Timing-Pool-Execution]
   WcTimingPool timeloopTiming; // collects timing data for each registered sweep
   WcTimer simTimer;            // measures the overall simulation time

   WALBERLA_MPI_WORLD_BARRIER()
   WALBERLA_LOG_INFO_ON_ROOT("Starting simulation with " << timesteps << " time steps")

   simTimer.start();
   timeloop.run(timeloopTiming);
   simTimer.end();

   WALBERLA_LOG_INFO_ON_ROOT("Simulation finished")

   const auto reducedTimeloopTiming = timeloopTiming.getReduced();
   WALBERLA_LOG_RESULT_ON_ROOT("Time loop timing:\n" << *reducedTimeloopTiming);
   //! [Marker-Timing-Pool-Execution]

   double time = simTimer.max();
   WALBERLA_MPI_SECTION() { walberla::mpi::reduceInplace(time, walberla::mpi::MAX); }
   WALBERLA_LOG_RESULT_ON_ROOT("Execution time: " << time);

   WALBERLA_ROOT_SECTION() { WALBERLA_CHECK_GREATER(time, 0.0); }

   WALBERLA_LOG_INFO_ON_ROOT("Execution Successful");

} // run()
} // namespace walberla::ExampleTimingPool

int main(int argc, char** argv)
{
   walberla::ExampleTimingPool::run(argc, argv);
   return EXIT_SUCCESS;
}
