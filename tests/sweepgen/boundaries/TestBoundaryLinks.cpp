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
//! \file BoundaryLinks.cpp
//! \author Philipp Suffa <philipp.suffa@fau.de>
//
//======================================================================================================================

#include "blockforest/all.h"
#include "blockforest/communication/UniformBufferedScheme.h"

#include "core/all.h"

#include "domain_decomposition/SharedSweep.h"

#include "field/all.h"
#include "field/communication/StencilRestrictedPackInfo.h"

#include "geometry/all.h"

#include "lbm/lattice_model/CollisionModel.h"

#include "stencil/all.h"

#include "timeloop/all.h"

#include "walberla/experimental/Sweep.hpp"

#if defined(WALBERLA_BUILD_WITH_OPENMESH)
#   include "mesh_common/DistanceFunction.h"
#   include "mesh_common/MeshIO.h"
#   include "mesh_common/TriangleMeshes.h"
#   include "mesh_common/distance_octree/DistanceOctree.h"
#endif

#include "gen/TestBoundaryLinks.hpp"

namespace TestBoundaryLinks
{
using namespace walberla;

using LbStencil  = stencil::D3Q7;
using PdfField_T = field::GhostLayerField< real_t, LbStencil::Q >;

template< typename Link >
std::vector< Link > createSolutionIndexVector()
{
   if constexpr (std::same_as< Link, gen::QBBLink >)
   {
      return std::vector< Link >(
         { Link(1, 1, 0, 5, 0.390625), Link(2, 1, 0, 5, 0.390625), Link(1, 2, 0, 5, 0.390625),
           Link(2, 2, 0, 5, 0.390625), Link(1, 0, 1, 1, 0.390625), Link(2, 0, 1, 1, 0.390625),
           Link(0, 1, 1, 4, 0.390625), Link(3, 1, 1, 3, 0.390625), Link(0, 2, 1, 4, 0.390625),
           Link(3, 2, 1, 3, 0.390625), Link(1, 3, 1, 2, 0.390625), Link(2, 3, 1, 2, 0.390625),
           Link(1, 0, 2, 1, 0.390625), Link(2, 0, 2, 1, 0.390625), Link(0, 1, 2, 4, 0.390625),
           Link(3, 1, 2, 3, 0.390625), Link(0, 2, 2, 4, 0.390625), Link(3, 2, 2, 3, 0.390625),
           Link(1, 3, 2, 2, 0.390625), Link(2, 3, 2, 2, 0.390625), Link(1, 1, 3, 6, 0.390625),
           Link(2, 1, 3, 6, 0.390625), Link(1, 2, 3, 6, 0.390625), Link(2, 2, 3, 6, 0.390625) });
   }
   else
   {
      return std::vector< Link >({ Link(1, 1, 0, 5), Link(2, 1, 0, 5), Link(1, 2, 0, 5), Link(2, 2, 0, 5),
                                   Link(1, 0, 1, 1), Link(2, 0, 1, 1), Link(0, 1, 1, 4), Link(3, 1, 1, 3),
                                   Link(0, 2, 1, 4), Link(3, 2, 1, 3), Link(1, 3, 1, 2), Link(2, 3, 1, 2),
                                   Link(1, 0, 2, 1), Link(2, 0, 2, 1), Link(0, 1, 2, 4), Link(3, 1, 2, 3),
                                   Link(0, 2, 2, 4), Link(3, 2, 2, 3), Link(1, 3, 2, 2), Link(2, 3, 2, 2),
                                   Link(1, 1, 3, 6), Link(2, 1, 3, 6), Link(1, 2, 3, 6), Link(2, 2, 3, 6) });
   }
}

template< typename Link >
void testIndexVector(
   shared_ptr< StructuredBlockForest > blocks,
   walberla::experimental::sweep::SparseIndexList< Link, walberla::experimental::memory::memtag::host > indexVector)
{
   auto reference_vec = createSolutionIndexVector< Link >();
   for (auto& block : *blocks)
   {
      auto idx_vec = indexVector.view(block);
      WALBERLA_CHECK_EQUAL(reference_vec.size(), idx_vec.size())
      for (const auto idx : idx_vec)
      {
         WALBERLA_CHECK(std::ranges::any_of(reference_vec, [&](const Link& link) {
            if constexpr (std::same_as< Link, gen::QBBLink >)
               return link.x == idx.x && link.y == idx.y && link.z == idx.z && link.dir == idx.dir &&
                      abs(link.q - idx.q) < 0.01;
            else
               return link.x == idx.x && link.y == idx.y && link.z == idx.z && link.dir == idx.dir;
         }))
      }
   }
}

void run(int argc, char** argv)
{
   Environment env{ argc, argv };
   AABB domainAABB(0.0, 0.0, 0.0, 3.0, 3.0, 3.0);
   auto blocks = blockforest::createUniformBlockGrid(domainAABB, 1, 1, 1, 4, 4, 4, true, true, true, true, true);

   BlockDataID pdfsId = field::addToStorage< PdfField_T >(blocks, "pdfs", real_c(0.0), field::fzyx, 1);

   real_t omega = 1.0;

   //  Set up boundary conditions
   AABB aabbBody(1.0, 1.0, 1.0, 2.0, 2.0, 2.0);

   auto cubeLinks = [&](auto link) -> bool {
      blocks->transformBlockLocalToGlobalCell(link.wallCell, link.block);
      blocks->transformBlockLocalToGlobalCell(link.fluidCell, link.block);
      auto cellCenterSolid = blocks->getCellCenter(link.wallCell);
      auto cellCenterFluid = blocks->getCellCenter(link.fluidCell);
      return walberla::geometry::contains(aabbBody, cellCenterSolid) &&
             !walberla::geometry::contains(aabbBody, cellCenterFluid);
   };

   auto cubeLinksData = [&](auto link) -> std::optional< gen::QBBData > {
      blocks->transformBlockLocalToGlobalCell(link.wallCell, link.block);
      blocks->transformBlockLocalToGlobalCell(link.fluidCell, link.block);
      auto cellCenterSolid = blocks->getCellCenter(link.wallCell);
      auto cellCenterFluid = blocks->getCellCenter(link.fluidCell);
      if (walberla::geometry::contains(aabbBody, cellCenterSolid) &&
          !walberla::geometry::contains(aabbBody, cellCenterFluid))
         return aabbBody.sqSignedDistance(cellCenterFluid);
      else
         return std::nullopt;
   };

   auto cubeDistanceFunc = [&](auto link) -> gen::QBBData {
      auto cellCenterFluid = blocks->getGlobalCellCenterFromBlockLocalCell(link.fluidCell, link.block);
      return aabbBody.sqSignedDistance(cellCenterFluid);
   };

   using flag_t            = walberla::uint8_t;
   using FlagField_T       = FlagField< flag_t >;
   BlockDataID flagFieldId = field::addFlagFieldToStorage< FlagField_T >(blocks, "flag field");

   const FlagUID fluidFlagUID("Fluid");
   const FlagUID bodyFlagUID("Body");

   for (auto& b : *blocks)
   {
      FlagField_T* flagField = b.template getData< FlagField_T >(flagFieldId);
      auto fluidFlag         = flagField->registerFlag(fluidFlagUID);
      auto bodyFlag          = flagField->registerFlag(bodyFlagUID);

      for (auto cIt = flagField->beginXYZ(); cIt != flagField->end(); ++cIt)
      {
         Cell c{ cIt.cell() };
         if (aabbBody.contains(blocks->getBlockLocalCellCenter(b, c))) { field::addFlag(cIt, bodyFlag); }
         else
         {
            field::addFlag(cIt, fluidFlag);
         }
      }
   }

   auto noSlip_link      = gen::NoSlipFactory{ blocks, pdfsId }.fromLinks(cubeLinks);
   auto noSlip_body      = gen::NoSlipFactory{ blocks, pdfsId }.fromBody(aabbBody);
   auto noSlip_flagField = gen::NoSlipFactory{ blocks, pdfsId }.fromFlagField< FlagField_T >(
      flagFieldId, bodyFlagUID, fluidFlagUID);

   auto qbb_link = gen::QBBFactory{ blocks, pdfsId, omega }.fromLinks(cubeLinksData);
   auto qbb_body = gen::QBBFactory{ blocks, pdfsId, omega }.fromBody(aabbBody, cubeDistanceFunc);

   testIndexVector(blocks, noSlip_link.indexVector());
   testIndexVector(blocks, noSlip_body.indexVector());
   testIndexVector(blocks, noSlip_flagField.indexVector());
   testIndexVector(blocks, qbb_link.indexVector());
   testIndexVector(blocks, qbb_body.indexVector());

   //  Timeloop
   SweepTimeloop loop{ blocks->getBlockStorage(), 1 };

   loop.add() << Sweep(noSlip_link);
   loop.add() << Sweep(noSlip_body);

   loop.add() << Sweep(qbb_link);
   loop.add() << Sweep(qbb_body);

#if defined(WALBERLA_BUILD_WITH_OPENMESH)
   auto meshFile = "cube.obj";
   auto mesh     = make_shared< mesh::TriangleMesh >();
   mesh::readAndBroadcast(meshFile, *mesh);
   mesh::translate(*mesh, Vector3< real_t >(1.0));

   auto triDist              = make_shared< mesh::TriangleDistance< mesh::TriangleMesh > >(mesh);
   auto distanceOctree       = make_shared< mesh::DistanceOctree< mesh::TriangleMesh > >(triDist);
   auto meshDistanceFunction = makeMeshDistanceFunction(distanceOctree);

   auto meshDistance = [&](auto link) -> gen::QBBData {
      auto cellCenterFluid = blocks->getGlobalCellCenterFromBlockLocalCell(link.fluidCell, link.block);
      return gen::QBBData(meshDistanceFunction(cellCenterFluid));
   };

   auto noSlip_mesh = gen::NoSlipFactory{ blocks, pdfsId }.fromDistanceFunction(meshDistanceFunction);
   auto qbb_mesh =
      gen::QBBFactory{ blocks, pdfsId, omega }.fromDistanceFunction(meshDistanceFunction, meshDistance);

   testIndexVector(blocks, noSlip_mesh.indexVector());
   testIndexVector(blocks, qbb_mesh.indexVector());

   loop.add() << Sweep(noSlip_mesh);
   loop.add() << Sweep(qbb_mesh);
#endif
   loop.run();
}
} // namespace TestBoundaryLinks

int main(int argc, char** argv)
{
   TestBoundaryLinks::run(argc, argv);
   return EXIT_SUCCESS;
}
