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
#include "walberla/v8/sweep/DomainSlices.hpp"

namespace lbm_scenarios
{

using namespace walberla;
using namespace walberla::v8;

template< memory::MemTag MemoryTag >
struct MirroredHalfChannel
{
   using ScalarField_T = memory::Field< real_t, 1, MemoryTag >;
   using VectorField_T = memory::Field< real_t, 3, MemoryTag >;
   using LbStencil     = stencil::D3Q19;
   using PdfField_T    = memory::Field< real_t, LbStencil::Q, MemoryTag >;

   static void run()
   {
      uint_t numProcesses = uint_c(MPIManager::instance()->numProcesses());
      std::vector< uint_t > numBlocksXY{ math::getFactors(numProcesses, 2u) };

      Vector3< uint_t > cellsPerBlock{ 4, 4, 64 };
      Vector3< bool > periodic{ true, true, false };

      auto blocks = blockforest::createUniformBlockGrid(numBlocksXY[0], numBlocksXY[1], 1,                    //
                                                        cellsPerBlock[0], cellsPerBlock[1], cellsPerBlock[2], //
                                                        1.0, true,                                            //
                                                        periodic[0], periodic[1], periodic[2]);

      ScalarField_T rho{ *blocks, 0UL, real_c(1.0) };
      VectorField_T u{ *blocks, 0UL, real_c(0.0) };
      PdfField_T pdfs{ *blocks };

      /* Hagen-Poiseuille-law in lattice units */
      const real_t u_max{ 0.025 };
      const real_t reynolds{ 10.0 };
      const real_t L_z{ real_c(2 * cellsPerBlock[2]) };
      const real_t radius{ L_z / 2.0 };
      const real_t r_squared{ radius * radius };
      const real_t lattice_viscosity{ L_z * u_max / reynolds };
      const real_t omega{ 2. / (6. * lattice_viscosity + 1.) };
      const real_t acceleration{ (u_max * 2.0 * lattice_viscosity) / r_squared };

      const Vector3< real_t > force{ acceleration, 0., 0. };

      auto velocityProfile = [&](real_t z) -> real_t {
         return acceleration / (2.0 * lattice_viscosity) * (r_squared - z * z);
      };

      gen::LBM::InitPdfs initPdfs{ pdfs, rho, u, force };
      for (auto& b : *blocks)
      {
         FieldView uView{ u, b };
         domain::BlockGeometry bGeom{ *blocks, b };

         sweep::forAllCells(exectag::Serial{}, *blocks, [&](Cell c) {
            auto cellCenter = bGeom.localGrid().cellCenter(c);
            uView(c, 0)     = velocityProfile(cellCenter[2]);
         });

         initPdfs(&b);
      }

      gen::LBM::StreamCollide streamCollide(pdfs, rho, u, force, omega);
      sweep::SweepFactory sweepFactory{ blocks };

      auto noSlipTop = sweepFactory.atDomainBorder< stencil::Direction::T >(gen::bc_grid_aligned::NoSlipTop(pdfs));
      auto freeSlipBottom =
         sweepFactory.atDomainBorder< stencil::Direction::B >(gen::bc_grid_aligned::FreeSlipBottom(pdfs));

      auto haloExchange = HaloExchange::create< LbStencil, MemoryTag >(blocks)
                             .sync(halo_exchange::streamPullSync< LbStencil >(pdfs))
                             .build();

      for (uint_t t = 0; t < 50; ++t)
      {
         for (auto& b : *blocks)
         {
            streamCollide(&b);
         }
         haloExchange();

         sweep::sync();

         for (auto& block : *blocks)
         {
            domain::BlockGeometry bGeom{ *blocks, block };
            FieldView uView{ u, block };

            sweep::forAllCells(exectag::Serial{}, *blocks, [&](Cell c) {
               auto cellCenter = bGeom.localGrid().cellCenter(c);
               real_t expected{ velocityProfile(cellCenter[2]) };
               real_t actual{ uView(c, 0) };
               testing::with_tolerance(1e-5, 1e-1).assert_close(actual, expected);
            });
         }

         for(auto& b: *blocks){
            noSlipTop(&b);
            freeSlipBottom(&b);
         }
      }
   }
};
} // namespace lbm_scenarios