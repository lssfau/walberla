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

#include "geometry/all.h"

#include "stencil/all.h"

#include "vtk/all.h"

#include "gen/LbmAlgorithms.hpp"
#include "walberla/V8.hpp"
#include "walberla/v8/Testing.hpp"

namespace lbm_scenarios
{

using namespace walberla;
using namespace walberla::v8;

template< memory::MemTag MemoryTag >
struct FullyPeriodic
{
   using ScalarField_T = memory::Field< real_t, 1, MemoryTag >;
   using VectorField_T = memory::Field< real_t, 3, MemoryTag >;
   using LbStencil     = stencil::D3Q19;
   using PdfField_T    = memory::Field< real_t, LbStencil::Q, MemoryTag >;

   static void run()
   {
      uint_t numProcesses = uint_c(MPIManager::instance()->numProcesses());
      Vector3< uint_t > numBlocks{ math::getFactors3D(numProcesses) };
      Vector3< uint_t > cellsPerBlock{ 4, 4, 4 };
      Vector3< bool > periodic{ true, true, true };

      auto blocks = blockforest::createUniformBlockGrid(numBlocks[0], numBlocks[1], numBlocks[2], cellsPerBlock[0],
                                                        cellsPerBlock[1], cellsPerBlock[2], 1.0, true, periodic[0],
                                                        periodic[1], periodic[2]);

      ScalarField_T rho{ *blocks, 0UL, real_c(1.0) };
      VectorField_T u{ *blocks, 0UL, real_c(0.0) };
      PdfField_T pdfs{ *blocks };

      const Vector3< real_t > force{ 0.005, 0., 0. };
      const real_t omega{ 1.0 };

      gen::LBM::InitPdfs init{ pdfs, rho, u, force };
      for (auto& b : *blocks)
      {
         init(&b);
      }

      gen::LBM::StreamCollide streamCollide(pdfs, rho, u, force, omega);

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

            sweep::forAllCells(exectag::Serial{}, *blocks, [&](Cell c) {
               real_t expected{ real_c(t) * force[0] };
               real_t actual{ uView(c, 0) };
               testing::with_tolerance(1e-10, 1e-6).assert_close(actual, expected);
            });
         }
      }
   }
};
} // namespace lbm_scenario_tests