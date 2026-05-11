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
//! \author Frederik Hennig
//
//======================================================================================================================

#pragma once

#include "blockforest/all.h"

#include "core/all.h"

#include "field/all.h"

#include "geometry/all.h"

#include "stencil/all.h"

#include "vtk/all.h"

#include "gen/LbmAlgorithms.hpp"
#include "walberla/V8.hpp"
#include "walberla/v8/Testing.hpp"
#include "walberla/v8/domain/GridGeometry.hpp"
#include "walberla/v8/sweep/DomainSlices.hpp"
#include "walberla/v8/testing/Testutils.hpp"

namespace lbm_scenarios
{

using namespace walberla;
using namespace walberla::v8;

template< memory::MemTag MemoryTag >
struct FreeSlipPipe
{
   using ScalarField_T = memory::Field< real_t, 1, MemoryTag >;
   using VectorField_T = memory::Field< real_t, 3, MemoryTag >;
   using LbStencil     = stencil::D3Q19;
   using PdfField_T    = memory::Field< real_t, LbStencil::Q, MemoryTag >;

   using FlagField_T = FlagField< uint8_t >;

   static void run()
   {
      uint_t numProcesses = uint_c(MPIManager::instance()->numProcesses());

      Vector3< uint_t > numBlocks{ numProcesses, 1, 1 };
      Vector3< uint_t > cellsPerBlock{ 4, 32, 32 };
      Vector3< bool > periodic{ true, false, false };

      auto blocks = blockforest::createUniformBlockGrid(numBlocks[0], numBlocks[1], numBlocks[2],             //
                                                        cellsPerBlock[0], cellsPerBlock[1], cellsPerBlock[2], //
                                                        1.0, true,                                            //
                                                        periodic[0], periodic[1], periodic[2]);

      ScalarField_T rho{ *blocks, 0UL, real_c(1.0) };
      VectorField_T u{ *blocks, 1UL, real_c(0.0) };
      PdfField_T pdfs{ *blocks };

      const BlockDataID flagFieldId = field::addFlagFieldToStorage< FlagField_T >(blocks, "flagField");

      const FlagUID fluidFlagUid{ "Fluid" };
      const FlagUID freeSlipFlagUID{ "FreeSlip" };

      const real_t pipeRadius{ 14.0 };
      const Vector3< real_t > pipeAnchor{ 0.0, 16.0, 16.0 };
      const real_t maxVelocity{ 0.02 };
      const Vector3< real_t > force{ 0., 0., 0. };
      const real_t omega{ 1.0 };

      gen::LBM::InitPdfs initPdfs{ pdfs, rho, u, force };
      for (auto& b : *blocks)
      {
         FieldView uView{ u, b };
         domain::BlockGeometry bGeom{ *blocks, b };

         FlagField_T& flagField = *b.getData< FlagField_T >(flagFieldId);
         const uint8_t freeSlipFlag{ flagField.getOrRegisterFlag(freeSlipFlagUID) };

         sweep::forAllCells(exectag::Serial{}, uView.slices().withGhostLayers(1), [&](Cell c) {
            auto cellCenter = bGeom.localGrid().cellCenter(c);
            cellCenter[0] = 0.0;

            Vector3< real_t > initVelocity;
            real_t radialDistance = (cellCenter - pipeAnchor).length();
            if (radialDistance > pipeRadius)
            {
               flagField.addFlag(c, freeSlipFlag);
               initVelocity = Vector3< real_t >{ NAN };
            }
            else
            {
               initVelocity = Vector3< real_t >{ maxVelocity, 0., 0. };
            }

            uView(c, 0) = initVelocity[0];
            uView(c, 2) = initVelocity[1];
            uView(c, 1) = initVelocity[2];
         });

         initPdfs(&b);
      }

      geometry::setNonBoundaryCellsToDomain< FlagField_T >(*blocks, flagFieldId, fluidFlagUid);

      gen::LBM::StreamCollide streamCollide(pdfs, rho, u, force, omega);
      auto freeSlip = gen::bc_sparse::FreeSlipFactory(blocks, pdfs)
                         .fromFlagField< FlagField_T >(flagFieldId, freeSlipFlagUID, fluidFlagUid);

      auto haloExchange = HaloExchange::create< LbStencil, MemoryTag >(blocks)
                             .sync(halo_exchange::streamPullSync< LbStencil >(pdfs))
                             .build();

      for (uint_t t = 0; t < 10; ++t)
      {
         for (auto& b : *blocks)
         {
            streamCollide(&b);
         }
         haloExchange();

         sweep::sync();

         for (auto& block : *blocks)
         {
            FieldView uView{ u, block };

            FlagField_T& flagField = *block.getData< FlagField_T >(flagFieldId);
            uint8_t fluidFlag      = flagField.getFlag(fluidFlagUid);

            sweep::forAllCells(exectag::Serial{}, *blocks, [&](Cell c) {
               if (flagField.isFlagSet(c, fluidFlag))
               {
                  v8::testing::assert_close(uView(c, 0), maxVelocity);
                  v8::testing::with_tolerance{ 1e-12, 1e-5 }.assert_close(uView(c, 1), 0.);
                  v8::testing::with_tolerance{ 1e-12, 1e-5 }.assert_close(uView(c, 2), 0.);
               }
            });
         }

         for (auto& b : *blocks)
         {
            freeSlip(&b);
         }
      }
   }
};
} // namespace lbm_scenarios