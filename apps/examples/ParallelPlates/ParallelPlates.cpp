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
//! \file ParallelPlates.cpp
//! \author Frederik Hennig <frederik.hennig@fau.de>
//
//======================================================================================================================

#include "blockforest/all.h"

#include "core/DataTypes.h"
#include "core/all.h"

#include "domain_decomposition/SharedSweep.h"

#include "stencil/all.h"

#include "timeloop/all.h"

#include "vtk/all.h"

#include <limits>

#include "gen/ParallelPlatesSweeps.hpp"
#include "walberla/V8.hpp"
#include "walberla/v8/Testing.hpp"
#include "walberla/v8/memory/MemoryTraits.hpp"
#include "walberla/v8/testing/Testutils.hpp"

namespace ParallelPlates
{
using namespace walberla;
using namespace walberla::v8;

using MemoryTag = v8::memtag::automatic;

using ScalarField_T = memory::Field< real_t, 1, MemoryTag >;
using VectorField_T = memory::Field< real_t, 3, MemoryTag >;

using LbStencil  = stencil::D3Q19;
using PdfField_T = memory::Field< real_t, LbStencil::Q, MemoryTag >;

struct ChannelType
{
   enum class Type { COUETTE, POISEUILLE };
   using enum Type;

   static ChannelType::Type fromStr(const std::string& channelType)
   {
      if (channelType == "couette") return COUETTE;
      if (channelType == "poiseuille") return POISEUILLE;
      throw std::invalid_argument{ channelType };
   }
};

void run(int argc, char** argv)
{
   Environment env{ argc, argv };
   auto config = env.config();

   auto blocks = blockforest::createUniformBlockGridFromConfig(config);

   ScalarField_T rho{ *blocks };
   VectorField_T u{ *blocks };
   PdfField_T pdfs{ *blocks };

   Config::BlockHandle simParams = config->getBlock("Parameters");
   std::string channelTypeStr    = simParams.getParameter< std::string >("channelType");
   ChannelType::Type channelType = ChannelType::fromStr(channelTypeStr);
   const real_t latticeViscosity = simParams.getParameter< real_t >("nu");
   const real_t channelVelocity  = simParams.getParameter< real_t >("u_max");
   const real_t errorThreshold   = simParams.getParameter< real_t >("errorThreshold");

   //  Prepare sweep functors

   std::function< void(IBlock*) > setAnalytical;
   std::function< void(IBlock*) > initializePdfs;
   std::function< void(IBlock*) > streamCollide;

   switch (channelType)
   {
   case ChannelType::COUETTE: {
      setAnalytical  = gen::Couette::SetAnalyticalSolution{ blocks, rho, u, channelVelocity };
      initializePdfs = gen::Couette::LBM::InitPdfs{ pdfs, rho, u };
      streamCollide =
         makeSharedSweep(std::make_shared< gen::Couette::LBM::StreamCollide >(pdfs, rho, u, latticeViscosity));
   }
   break;
   case ChannelType::POISEUILLE: {
      setAnalytical  = gen::Poiseuille::SetAnalyticalSolution{ blocks, rho, u, channelVelocity };
      initializePdfs = gen::Poiseuille::LBM::InitPdfs{ blocks, pdfs, rho, u, latticeViscosity, channelVelocity };
      streamCollide  = makeSharedSweep(std::make_shared< gen::Poiseuille::LBM::StreamCollide >(
         blocks, pdfs, rho, u, latticeViscosity, channelVelocity));
   }
   break;
   }

   //  Set up initial state

   for (auto& b : *blocks)
   {
      setAnalytical(&b);
      initializePdfs(&b);
   }

   //  Set up ghost layer communication
   auto haloExchange = HaloExchange::create< LbStencil, MemoryTag >(blocks)
                          .sync(halo_exchange::streamPullSync< LbStencil >(pdfs))
                          .makeShared();

   //  Set up boundary conditions
   auto intersectsUpperWall = [&](auto link) -> bool {
      blocks->transformBlockLocalToGlobalCell(link.wallCell, link.block);
      return link.wallCell.z() > blocks->getDomainCellBB().zMax();
   };

   auto intersectsLowerWall = [&](auto link) -> bool {
      blocks->transformBlockLocalToGlobalCell(link.wallCell, link.block);
      return link.wallCell.z() < blocks->getDomainCellBB().zMin();
   };

   auto noSlip = gen::NoSlipFactory{ blocks, pdfs }.fromLinks([&](auto link) {
      return intersectsLowerWall(link) || (channelType == ChannelType::POISEUILLE && intersectsUpperWall(link));
   });

   auto ubb = gen::UBBFactory{ blocks, pdfs, channelVelocity }.fromLinks(
      [&](auto link) { return channelType == ChannelType::COUETTE && intersectsUpperWall(link); });

   //  Timeloop
   const uint_t numTimesteps{ simParams.getParameter< uint_t >("timesteps") };
   SweepTimeloop loop{ blocks->getBlockStorage(), numTimesteps };

   loop.add() << Sweep(streamCollide) << AfterFunction(SharedFunctor(haloExchange));
   loop.add() << Sweep(noSlip);
   loop.add() << Sweep(ubb);

   RemainingTimeLogger logger{ numTimesteps };
   loop.addFuncAfterTimeStep(logger);

   //  VTK Output

   Config::BlockHandle outputParams = config->getBlock("Output");

   const uint_t vtkWriteFrequency = outputParams.getParameter< uint_t >("vtkWriteFrequency", 0);
   if (vtkWriteFrequency > 0)
   {
      auto vtkOutput =
         vtk::createVTKOutput_BlockData(*blocks, "vtk", vtkWriteFrequency, 0, false, "vtk_out_" + channelTypeStr,
                                        "simulation_step", false, true, true, false, 0);

      auto densityWriter = make_shared< memory::FieldVtkWriter< ScalarField_T, float32 > >(rho, "density");
      vtkOutput->addCellDataWriter(densityWriter);

      auto velWriter = make_shared< memory::FieldVtkWriter< VectorField_T, float32 > >(u, "velocity");
      vtkOutput->addCellDataWriter(velWriter);

      loop.addFuncAfterTimeStep(vtk::writeFiles(vtkOutput), "VTK Output");
   }

   //  Run the Simulation

   WALBERLA_LOG_INFO_ON_ROOT("Commencing simulation with " << numTimesteps << " timesteps")

   loop.run();

   //  Check solution
   WALBERLA_LOG_INFO_ON_ROOT("Checking for convergence...")

   auto velocityErrorLmax =
      std::allocate_shared< real_t, v8::memory::MemoryTraits< MemoryTag, real_t >::AllocatorType >(
         {}, -std::numeric_limits< real_t >::infinity());
   std::function< void(IBlock*) > computeVelocityError;

   switch (channelType)
   {
   case ChannelType::COUETTE: {
      computeVelocityError = gen::Couette::VelocityErrorLmax{ blocks, u, velocityErrorLmax.get(), channelVelocity };
   }
   break;
   case ChannelType::POISEUILLE: {
      computeVelocityError = gen::Poiseuille::VelocityErrorLmax{ blocks, u, velocityErrorLmax.get(), channelVelocity };
   }
   break;
   }

   for (auto& b : *blocks)
   {
      computeVelocityError(&b);
   }

   mpi::reduceInplace(*velocityErrorLmax, mpi::MAX);

   WALBERLA_LOG_INFO_ON_ROOT("Lmax error of x-velocity: " << *velocityErrorLmax);

   WALBERLA_ROOT_SECTION() { testing::assert_less(*velocityErrorLmax, errorThreshold); }
}
} // namespace ParallelPlates

int main(int argc, char** argv)
{
   ParallelPlates::run(argc, argv);
   return EXIT_SUCCESS;
}
