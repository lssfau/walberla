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
//! \file DoubleShearLayer.cpp
//! \author Frederik Hennig <frederik.hennig@fau.de>
//
//======================================================================================================================

#include "blockforest/all.h"

#include "core/all.h"

#include "domain_decomposition/SharedSweep.h"

#include "stencil/all.h"

#include "timeloop/all.h"

#include "vtk/all.h"

#include "gen/DoubleShearLayerSweeps.hpp"
#include "walberla/v8/HaloExchange.hpp"
#include "walberla/v8/Memory.hpp"
#include "walberla/v8/Sweep.hpp"

namespace DoubleShearLayer
{

//! [begin aliases]
using namespace walberla;
using namespace walberla::v8;

using MemoryTag = memtag::automatic;

using ScalarField_T = memory::Field< real_t, 1, MemoryTag >;
using VectorField_T = memory::Field< real_t, 3, MemoryTag >;

using LbStencil  = stencil::D3Q19;
using PdfField_T = memory::Field< real_t, LbStencil::Q, MemoryTag >;
//! [end aliases]

void run(int argc, char** argv)
{
   Environment env{ argc, argv };
   auto config = env.config();

   Config::BlockHandle domainParams = config->getBlock("Domain");

   Vector3< uint_t > numBlocks     = domainParams.getParameter< Vector3< uint_t > >("blocks");
   Vector3< uint_t > cellsPerBlock = domainParams.getParameter< Vector3< uint_t > >("cellsPerBlock");

   WALBERLA_CHECK_EQUAL(cellsPerBlock[0], cellsPerBlock[1], "Number of cells in x- and y- direction must be the same");

   AABB domainAabb{ 0., 0., 0., 1., 1., real_c(cellsPerBlock[2]) / real_c(cellsPerBlock[0]) };
   std::array< bool, 3 > periodic{ true, true, true };

   auto blocks = blockforest::createUniformBlockGrid(
      domainAabb, numBlocks[0], numBlocks[1], numBlocks[2], cellsPerBlock[0], cellsPerBlock[1], cellsPerBlock[2],
      /*oneBlockPerProcess*/ true, periodic[0], periodic[1], periodic[2]);

   //! [end domain setup]

   //! [begin fields]
   ScalarField_T rho{ *blocks };
   VectorField_T u{ *blocks };
   PdfField_T pdfs{ *blocks };
   ScalarField_T vorticity{ *blocks };
   //! [end fields]

   //  Simulation Parameters
   Config::BlockHandle simParams = config->getBlock("Parameters");

   const real_t reynolds{ simParams.getParameter< real_t >("Reynolds") };
   const real_t delta{ simParams.getParameter< real_t >("delta") };
   const real_t kappa{ simParams.getParameter< real_t >("kappa") };
   const real_t u_0{ simParams.getParameter< real_t >("u_0") };

   //  Initial State

   gen::SetInitialState setInitialState{ blocks, rho, u, delta, kappa, u_0 };
   gen::LBM::InitPdfs lbInit{ pdfs, rho, u };

   for (auto& b : *blocks)
   {
      setInitialState(&b);
      lbInit(&b);
   }
   //!  [end initial state]

   //  Set up LB stream/collide sweep

   //  Compute relaxation rate
   const real_t N = real_c(blocks->getDomainCellBB().xSize());
   const real_t nu{ (u_0 * N) / reynolds };
   const real_t theta{
      1. / 3,
   };
   const real_t tau{ nu / theta };
   const real_t omega{ 2. / (2. * tau + 1.) };

   auto streamCollide = std::make_shared< gen::LBM::StreamCollide >(pdfs, rho, u, omega);

   //! [begin halo exchange]
   //  Set up ghost layer communication
   auto haloExchange = HaloExchange::create< LbStencil, MemoryTag >(blocks)
                          .sync(halo_exchange::streamPullSync< LbStencil >(pdfs))
                          .sync(u)
                          .makeShared();

   //! [end halo exchange]

   //  Timeloop
   const uint_t numTimesteps{ simParams.getParameter< uint_t >("timesteps") };
   SweepTimeloop loop{ blocks->getBlockStorage(), numTimesteps };

   loop.add() << Sweep(makeSharedSweep(streamCollide)) << AfterFunction(SharedFunctor(haloExchange));
   loop.add() << Sweep(gen::ComputeVorticity{ blocks, u, vorticity });

   RemainingTimeLogger logger{ numTimesteps };
   loop.addFuncAfterTimeStep(logger);

   //  VTK Output

   Config::BlockHandle outputParams = config->getBlock("Output");

   const uint_t vtkWriteFrequency = outputParams.getParameter< uint_t >("vtkWriteFrequency", 0);
   if (vtkWriteFrequency > 0)
   {
      auto vtkOutput = vtk::createVTKOutput_BlockData(*blocks, "vtk", vtkWriteFrequency, 0, false, "vtk_out",
                                                      "simulation_step", false, true, true, false, 0);

      auto densityWriter = make_shared< memory::FieldVtkWriter< ScalarField_T, float32 > >(rho, "density");
      vtkOutput->addCellDataWriter(densityWriter);

      auto velWriter = make_shared< memory::FieldVtkWriter< VectorField_T, float32 > >(u, "velocity");
      vtkOutput->addCellDataWriter(velWriter);

      auto vorticityWriter = make_shared< memory::FieldVtkWriter< ScalarField_T, float32 > >(vorticity, "vorticity");
      vtkOutput->addCellDataWriter(vorticityWriter);

      loop.addFuncAfterTimeStep(vtk::writeFiles(vtkOutput), "VTK Output");
   }
   //! [end vtk output]

   //  Run the Simulation

   WALBERLA_LOG_INFO_ON_ROOT("Commencing simulation with " << numTimesteps << " timesteps")

   loop.run();
}
} // namespace DoubleShearLayer

int main(int argc, char** argv)
{
   DoubleShearLayer::run(argc, argv);
   return EXIT_SUCCESS;
}
//! [end file]
