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
//! \file FlowAroundSphere.cpp
//! \author Philipp Suffa <philipp.suffa@fau.de>
//
//======================================================================================================================

#include "core/all.h"
#include "blockforest/all.h"
#include "geometry/all.h"
#include "lbm/lattice_model/CollisionModel.h"
#include "stencil/all.h"
#include "field/all.h"
#include "timeloop/all.h"

#include "blockforest/communication/UniformBufferedScheme.h"
#include "field/communication/StencilRestrictedPackInfo.h"
#include "domain_decomposition/SharedSweep.h"

#include "walberla/experimental/Sweep.hpp"

#include "gen/FlowAroundSphereExample.hpp"

namespace FlowAroundSphere
{
    using namespace walberla;

    using ScalarField_T = field::GhostLayerField<real_t, 1>;
    using VectorField_T = field::GhostLayerField<real_t, 3>;

    using LbStencil = stencil::D3Q19;
    using PdfField_T = field::GhostLayerField<real_t, LbStencil::Q>;


    FlowAroundSphereExample::gen::QBBData sqSignedDistanceToSphere(geometry::Sphere sphere, Vector3<real_t> point) {
       real_t distance = sqrt(pow(sphere.midpoint()[0] - point[0], 2) + pow(sphere.midpoint()[1] - point[1], 2) + pow(sphere.midpoint()[2] - point[2], 2)) - sphere.radius();
       real_t sign = 1.0 ? distance >= 0 : -1.0;
       return FlowAroundSphereExample::gen::QBBData(sign * distance * distance);
    }

    void run(int argc, char **argv)
    {
        Environment env{argc, argv};
        auto config = env.config();

        auto blocks = blockforest::createUniformBlockGridFromConfig(config);
        auto dx = blocks->dx();
        auto domainAABB = blocks->getDomain();
        WALBERLA_LOG_INFO_ON_ROOT("dx " << dx << " domain " << domainAABB)
        BlockDataID pdfsId = field::addToStorage<PdfField_T>(blocks, "pdfs", real_c(0.0), field::fzyx, 1);
        BlockDataID rhoId = field::addToStorage<ScalarField_T>(blocks, "rho", real_c(1.0), field::fzyx, 0);
        BlockDataID uId = field::addToStorage<VectorField_T>(blocks, "u", real_c(0.0), field::fzyx, 0);
        
        Config::BlockHandle simParams = config->getBlock("Parameters");
        const real_t reynoldsNumber{ simParams.getParameter< real_t >("re") };
        const real_t vel_inflow = simParams.getParameter< real_t >("u_max");
        const real_t refLen = domainAABB.xSize()*0.1;
        const real_t viscosity = vel_inflow * refLen / reynoldsNumber;
        real_t omega = lbm::collision_model::omegaFromViscosity(viscosity);

        // Initialize macroscopic fields and PDF field
        FlowAroundSphereExample::gen::IntializeMacroFields initFields{rhoId, uId, vel_inflow};
        FlowAroundSphereExample::gen::LBM::InitPdfs lbInit{pdfsId, rhoId, uId};

        for (auto &b : *blocks)
        {
           initFields(&b);
           lbInit(&b);
        }

        //  Set up LB stream/collide sweep
        auto streamCollide = std::make_shared< FlowAroundSphereExample::gen::LBM::StreamCollide >(pdfsId, rhoId, uId, omega);

        //  Set up ghost layer communication
        blockforest::communication::UniformBufferedScheme<LbStencil> comm{blocks};
        comm.addPackInfo(std::make_shared<field::communication::StencilRestrictedPackInfo<PdfField_T, LbStencil>>(pdfsId));

        //  Set up boundary conditions
        geometry::Sphere sphere(Vector3<real_t> (domainAABB.xSize()*0.35,domainAABB.ySize()*0.5,domainAABB.zSize()*0.5), refLen*0.5);

        auto inflowLinks = [&](auto link) -> bool {
           blocks->transformBlockLocalToGlobalCell(link.wallCell, link.block);
           return link.wallCell.x() < blocks->getDomainCellBB().xMin();
        };

        auto outflowLinks = [&](auto link) -> bool {
           blocks->transformBlockLocalToGlobalCell(link.wallCell, link.block);
           return link.wallCell.x() > blocks->getDomainCellBB().xMax();
        };

        auto sideWallLinks = [&](auto link) -> bool {
           blocks->transformBlockLocalToGlobalCell(link.wallCell, link.block);
           return link.wallCell.y() < blocks->getDomainCellBB().yMin() || link.wallCell.y() > blocks->getDomainCellBB().yMax() ||
                  link.wallCell.z() < blocks->getDomainCellBB().zMin() || link.wallCell.z() > blocks->getDomainCellBB().zMax();
        };

        auto sphereLinks = [&](auto link) -> std::optional< FlowAroundSphereExample::gen::QBBData >  {
           blocks->transformBlockLocalToGlobalCell(link.wallCell, link.block);
           blocks->transformBlockLocalToGlobalCell(link.fluidCell, link.block);
           auto cellCenterSolid = blocks->getCellCenter(link.wallCell);
           auto cellCenterFluid = blocks->getCellCenter(link.fluidCell);
           if (contains(sphere, cellCenterSolid) && !contains(sphere, cellCenterFluid))
              return sqSignedDistanceToSphere(sphere, cellCenterFluid);
           else
              return std::nullopt;
        };


        auto qbb = FlowAroundSphereExample::gen::QBBFactory{ blocks, pdfsId, omega }.fromLinks(sphereLinks);
        auto noSlip = FlowAroundSphereExample::gen::NoSlipFactory{ blocks, pdfsId }.fromLinks(sideWallLinks);
        auto inflow = FlowAroundSphereExample::gen::UBBFactory{ blocks, pdfsId, vel_inflow }.fromLinks(inflowLinks);
        auto outflow = FlowAroundSphereExample::gen::OutflowFactory{ blocks, pdfsId }.fromLinks(outflowLinks);
        
        //  Timeloop
        const uint_t numTimesteps{simParams.getParameter<uint_t>("timesteps")};
        SweepTimeloop loop{blocks->getBlockStorage(), numTimesteps};

        loop.add() << Sweep(makeSharedSweep(streamCollide)) << AfterFunction(comm);
        loop.add() << Sweep(noSlip);
        loop.add() << Sweep(inflow);
        loop.add() << Sweep(outflow);
        loop.add() << Sweep(qbb);


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

            loop.addFuncAfterTimeStep(vtk::writeFiles(vtkOutput), "VTK Output");
        }

        //  Run the Simulation
        WALBERLA_LOG_INFO_ON_ROOT("Commencing simulation with " << numTimesteps << " timesteps")
        loop.run();
    }
}

int main(int argc, char **argv)
{
    FlowAroundSphere::run(argc, argv);
    return EXIT_SUCCESS;
}
