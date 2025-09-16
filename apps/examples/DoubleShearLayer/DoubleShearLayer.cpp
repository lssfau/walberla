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

#include "core/all.h"
#include "blockforest/all.h"
#include "stencil/all.h"
#include "field/all.h"
#include "timeloop/all.h"

#include "blockforest/communication/UniformBufferedScheme.h"
#include "field/communication/StencilRestrictedPackInfo.h"
#include "field/communication/PackInfo.h"
#include "domain_decomposition/SharedSweep.h"

#include "walberla/experimental/Sweep.hpp"

#include "gen/DoubleShearLayerSweeps.hpp"

namespace DoubleShearLayer
{
    using namespace walberla;

    using ScalarField_T = field::GhostLayerField<real_t, 1>;
    using VectorField_T = field::GhostLayerField<real_t, 3>;

    using LbStencil = stencil::D3Q19;
    using PdfField_T = field::GhostLayerField<real_t, LbStencil::Q>;

    using CommScheme = blockforest::communication::UniformBufferedScheme<LbStencil>;
    using PdfsPackInfo = field::communication::StencilRestrictedPackInfo<PdfField_T, LbStencil>;

    void run(int argc, char **argv)
    {
        Environment env{argc, argv};
        auto config = env.config();

        Config::BlockHandle domainParams = config->getBlock("Domain");

        Vector3< uint_t > numBlocks = domainParams.getParameter< Vector3< uint_t > >("blocks");
        Vector3< uint_t > cellsPerBlock = domainParams.getParameter< Vector3< uint_t > >("cellsPerBlock");

        WALBERLA_CHECK_EQUAL(cellsPerBlock[0], cellsPerBlock[1], "Number of cells in x- and y- direction must be the same");

        AABB domainAabb { 0., 0., 0., 1., 1., real_c(cellsPerBlock[2]) / real_c(cellsPerBlock[0]) };
        std::array< bool, 3 > periodic { true, true, true };

        auto blocks = blockforest::createUniformBlockGrid(
            domainAabb,
            numBlocks[0], numBlocks[1], numBlocks[2],
            cellsPerBlock[0], cellsPerBlock[1], cellsPerBlock[2],
            /*oneBlockPerProcess*/ true,
            periodic[0], periodic[1], periodic[2]
        );

        //! [end domain setup]

        BlockDataID pdfsId = field::addToStorage<PdfField_T>(blocks, "pdfs", real_c(0.0), field::fzyx, 1);
        BlockDataID rhoId = field::addToStorage<ScalarField_T>(blocks, "rho", real_c(1.0), field::fzyx, 0);
        BlockDataID uId = field::addToStorage<VectorField_T>(blocks, "u", real_c(0.0), field::fzyx, 1);
        BlockDataID vorticityId = field::addToStorage<ScalarField_T>(blocks, "vorticity", real_c(1.0), field::fzyx, 0);
        
        //  Simulation Parameters
        Config::BlockHandle simParams = config->getBlock("Parameters");
        
        const real_t reynolds{ simParams.getParameter< real_t >("Reynolds") };
        const real_t delta{ simParams.getParameter< real_t >("delta") };
        const real_t kappa{ simParams.getParameter< real_t >("kappa") };
        const real_t u_0{ simParams.getParameter< real_t >("u_0") };
        
        //  Initial State

        gen::SetInitialState setInitialState{ blocks, rhoId, uId, delta, kappa, u_0};
        gen::LBM::InitPdfs lbInit{pdfsId, rhoId, uId};

        for (auto &b : *blocks)
        {
            setInitialState(&b);
            lbInit(&b);
        }
        //!  [end initial state]

        //  Set up LB stream/collide sweep

        //  Compute relaxation rate
        const real_t N = real_c(blocks->getDomainCellBB().xSize());
        const real_t nu { (u_0 * N) / reynolds };
        const real_t theta { 1. / 3, };
        const real_t tau { nu / theta };
        const real_t omega{ 2. / (2. * tau + 1.) };

        auto streamCollide = std::make_shared< gen::LBM::StreamCollide >(pdfsId, rhoId, uId, omega);

        //  Set up ghost layer communication
        CommScheme comm{blocks};
        comm.addPackInfo(std::make_shared<PdfsPackInfo>(pdfsId));
        comm.addPackInfo(std::make_shared<field::communication::PackInfo< VectorField_T > >(uId)); // Finite differences in vorticity need neighbor information

        //  Timeloop
        const uint_t numTimesteps{simParams.getParameter<uint_t>("timesteps")};
        SweepTimeloop loop{blocks->getBlockStorage(), numTimesteps};

        loop.add() << Sweep(makeSharedSweep(streamCollide)) << AfterFunction(comm);
        loop.add() << Sweep(gen::ComputeVorticity{blocks, uId, vorticityId});

        RemainingTimeLogger logger{numTimesteps};
        loop.addFuncAfterTimeStep(logger);

        //  VTK Output

        Config::BlockHandle outputParams = config->getBlock("Output");

        const uint_t vtkWriteFrequency = outputParams.getParameter<uint_t>("vtkWriteFrequency", 0);
        if (vtkWriteFrequency > 0)
        {
            auto vtkOutput = vtk::createVTKOutput_BlockData(*blocks, "vtk", vtkWriteFrequency, 0, false, "vtk_out",
                                                            "simulation_step", false, true, true, false, 0);

            auto densityWriter = make_shared<field::VTKWriter<ScalarField_T, float32>>(rhoId, "density");
            vtkOutput->addCellDataWriter(densityWriter);

            auto velWriter = make_shared<field::VTKWriter<VectorField_T, float32>>(uId, "velocity");
            vtkOutput->addCellDataWriter(velWriter);

            auto vorticityWriter = make_shared<field::VTKWriter<ScalarField_T, float32>>(vorticityId, "vorticity");
            vtkOutput->addCellDataWriter(vorticityWriter);

            loop.addFuncAfterTimeStep(vtk::writeFiles(vtkOutput), "VTK Output");
        }

        //  Run the Simulation

        WALBERLA_LOG_INFO_ON_ROOT("Commencing simulation with " << numTimesteps << " timesteps")

        loop.run();
    }
}  // namespace DoubleShearLayer

int main(int argc, char **argv)
{
    DoubleShearLayer::run(argc, argv);
    return EXIT_SUCCESS;
}
//! [end file]
