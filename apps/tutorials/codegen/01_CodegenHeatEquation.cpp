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
//! \file 01_CodegenHeatEquation.cpp
//! \author Frederik Hennig <frederik.hennig@fau.de>
//
//======================================================================================================================

#include "blockforest/Initialization.h"
#include "blockforest/communication/UniformBufferedScheme.h"

#include "core/Environment.h"
#include "core/math/Constants.h"

#include "field/AddToStorage.h"
#include "field/communication/PackInfo.h"
#include "field/vtk/VTKWriter.h"

#include "pde/boundary/Neumann.h"

#include "stencil/D2Q9.h"

#include "timeloop/SweepTimeloop.h"

#include "HeatEquationKernel.h"

namespace walberla
{
typedef GhostLayerField< real_t, 1 > ScalarField;

void swapFields(StructuredBlockForest& blocks, BlockDataID uID, BlockDataID uTmpID)
{
   for (auto block = blocks.begin(); block != blocks.end(); ++block)
   {
      ScalarField* u     = block->getData< ScalarField >(uID);
      ScalarField* u_tmp = block->getData< ScalarField >(uTmpID);

      u->swapDataPointers(u_tmp);
   }
}

void initDirichletBoundaryNorth(shared_ptr< StructuredBlockForest > blocks, BlockDataID uId, BlockDataID uTmpId)
{
   for (auto block = blocks->begin(); block != blocks->end(); ++block)
   {
      if (blocks->atDomainYMaxBorder(*block))
      {
         ScalarField* u     = block->getData< ScalarField >(uId);
         ScalarField* u_tmp = block->getData< ScalarField >(uTmpId);

         CellInterval xyz = u->xyzSizeWithGhostLayer();

         xyz.yMin() = xyz.yMax();
         for (auto cell = xyz.begin(); cell != xyz.end(); ++cell)
         {
            const Vector3< real_t > p = blocks->getBlockLocalCellCenter(*block, *cell);
            //  Set the dirichlet boundary to f(x) = 1 + sin(x) * x^2
            real_t v          = real_c(1.0 + std::sin(2 * math::pi * p[0]) * p[0] * p[0]);
            u->get(*cell)     = v;
            u_tmp->get(*cell) = v;
         }
      }
   }
}

int main(int argc, char** argv)
{
   mpi::Environment env(argc, argv);

   /////////////////////////////
   /// SIMULATION PARAMETERS ///
   /////////////////////////////

   //  Ensure matching aspect ratios of cells and domain.
   const uint_t xCells = uint_c(25);
   const uint_t yCells = uint_c(25);

   const real_t xSize = real_c(1.0);
   const real_t ySize = real_c(1.0);

   const uint_t xBlocks = uint_c(1);
   const uint_t yBlocks = uint_c(1);

   const uint_t processes = uint_c(MPIManager::instance()->numProcesses());

   if (processes != xBlocks * yBlocks)
   { WALBERLA_ABORT("The number of processes must be equal to the number of blocks!"); }

   const real_t dx = xSize / real_c(xBlocks * xCells + uint_t(1));
   const real_t dy = ySize / real_c(yBlocks * yCells + uint_t(1));

   WALBERLA_CHECK_FLOAT_EQUAL(dx, dy);

   const real_t dt    = real_c(1e-4);
   const real_t kappa = real_c(1.0);

   ///////////////////////////
   /// BLOCK STORAGE SETUP ///
   ///////////////////////////

   auto aabb = math::AABB(real_c(0.5) * dx, real_c(0.5) * dy, real_c(0.0), xSize - real_c(0.5) * dx,
                          ySize - real_c(0.5) * dy, dx);

   shared_ptr< StructuredBlockForest > blocks = blockforest::createUniformBlockGrid(
      aabb, xBlocks, yBlocks, uint_c(1), xCells, yCells, 1, true, false, false, false);

   //////////////
   /// FIELDS ///
   //////////////

   BlockDataID uFieldId    = field::addToStorage< ScalarField >(blocks, "u", real_c(0.0), field::fzyx, uint_c(1));
   BlockDataID uTmpFieldId = field::addToStorage< ScalarField >(blocks, "u_tmp", real_c(0.0), field::fzyx, uint_c(1));

   /////////////////////
   /// COMMUNICATION ///
   /////////////////////

   blockforest::communication::UniformBufferedScheme< stencil::D2Q9 > commScheme(blocks);
   commScheme.addPackInfo(make_shared< field::communication::PackInfo< ScalarField > >(uFieldId));

   //////////////////////////
   /// DIRICHLET BOUNDARY ///
   //////////////////////////

   initDirichletBoundaryNorth(blocks, uFieldId, uTmpFieldId);

   ////////////////////////
   /// NEUMANN BOUNDARY ///
   ////////////////////////

   pde::NeumannDomainBoundary< ScalarField > neumann(*blocks, uFieldId);

   neumann.excludeBoundary(stencil::N);
   neumann.excludeBoundary(stencil::B);
   neumann.excludeBoundary(stencil::T);

   ////////////////
   /// TIMELOOP ///
   ////////////////

   SweepTimeloop timeloop(blocks, uint_c(2e4));

   timeloop.add() << BeforeFunction(commScheme, "Communication") << BeforeFunction(neumann, "Neumann Boundaries")
                  << Sweep(pystencils::HeatEquationKernel(uFieldId, uTmpFieldId, dt, dx, kappa), "HeatEquationKernel")
                  << AfterFunction([blocks, uFieldId, uTmpFieldId]() { swapFields(*blocks, uFieldId, uTmpFieldId); },
                                   "Swap");

   auto vtkWriter = field::createVTKOutput< ScalarField, float >( uFieldId, *blocks, "temperature", uint_c(200), uint_c(0) );
   vtkWriter();
   timeloop.addFuncAfterTimeStep(vtkWriter, "VTK");

   timeloop.run();

   return EXIT_SUCCESS;
}
}

int main(int argc, char** argv) { walberla::main(argc, argv); }
