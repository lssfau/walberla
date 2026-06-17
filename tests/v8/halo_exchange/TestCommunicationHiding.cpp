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
//! \file TestCommunicationHiding.cpp
//! \author Philipp Suffa <philipp.suffa@fau.de>
//
//======================================================================================================================

#include "blockforest/all.h"

#include "core/all.h"

#include "stencil/all.h"

#include "timeloop/all.h"

#include "gen/TestCommunicationHiding.hpp"
#include "walberla/v8/HaloExchange.hpp"
#include "walberla/v8/Sweep.hpp"
#include "walberla/v8/Testing.hpp"

namespace TestCommunicationHiding
{

using namespace walberla;
using namespace walberla::v8;

using MemoryTag = memtag::automatic;

using ScalarField_T = memory::Field< real_t, 1, MemoryTag >;
using VectorField_T = memory::Field< real_t, 3, MemoryTag >;

using LbStencil  = stencil::D3Q19;
using PdfField_T = memory::Field< real_t, LbStencil::Q, MemoryTag >;

void TestCommunicationHiding()
{
   Vector3< uint_t > cellsPerBlock = Vector3< uint_t >(10, 10, 10);

   const uint_t numProcs = uint_t(mpi::MPIManager::instance()->numProcesses());
   const Vector3< uint_t > numBlocks = math::getFactors3D(numProcs);

   AABB domainAabb{ 0., 0., 0., 1., 1., real_c(cellsPerBlock[2]) / real_c(cellsPerBlock[0]) };
   std::array< bool, 3 > periodic{ true, true, true };

   auto blocks = blockforest::createUniformBlockGrid(
      domainAabb, numBlocks[0], numBlocks[1], numBlocks[2], cellsPerBlock[0], cellsPerBlock[1], cellsPerBlock[2],
      true, periodic[0], periodic[1], periodic[2]);

   ScalarField_T rho{ *blocks };
   VectorField_T u{ *blocks };
   PdfField_T pdfs{ *blocks };

   const real_t maxVelocity{0.01};
   for (auto& block : *blocks){
      FieldView uv{ u, block };
      FieldView rhov{ rho, block };
      sweep::forAllCells(exectag::Serial{}, uv.slices().interior(), [&](Cell c){
         uv(c, 0) = maxVelocity;
         uv(c, 1) = 0.0;
         uv(c, 2) = 0.0;
         rhov(c) = 1.0;
      });
   }

   gen::LBM::InitPdfs lbInit{ pdfs, rho, u };
   for (auto& b : *blocks)
   {
      lbInit(&b);
   }
   const real_t omega{0.8};
   auto haloExchange = HaloExchange::create< LbStencil, MemoryTag >(blocks)
                          .sync(halo_exchange::streamPullSync< LbStencil >(pdfs))
                          .makeShared();


   auto sweepFactory = sweep::SweepFactory(blocks);
   auto innerOuterSplit = Cell(1,1,1);
   auto streamCollide = gen::LBM::StreamCollide(pdfs, rho, u, omega);

   auto streamCollideInner = sweepFactory.commHidingSweep<sweep::INNER>( streamCollide, innerOuterSplit);
   auto streamCollideOuter = sweepFactory.commHidingSweep<sweep::OUTER>( streamCollide, innerOuterSplit);

   // One stream Collide needed up front to bring the pdfs into the correct state
   for (auto& b : *blocks)
   {
      (streamCollide)(&b);
   }
   SweepTimeloop loop{ blocks->getBlockStorage(), 10 };
   loop.add() << BeforeFunction([haloExchange](){haloExchange->startCommunication();}) << Sweep(streamCollideInner) ;
   loop.add() << BeforeFunction([haloExchange](){haloExchange->wait();}) << Sweep(streamCollideOuter) ;

   loop.run();
   sweep::sync();
   for (auto& block : *blocks)
   {
      FieldView uv{ u, block };
      sweep::forAllCells(exectag::Serial{}, *blocks, [&](Cell c) {
         testing::assert_close(uv(c, 0), maxVelocity);
         testing::with_tolerance{ 1e-12, 1e-5 }.assert_close(uv(c, 1), real_c(0.));
         testing::with_tolerance{ 1e-12, 1e-5 }.assert_close(uv(c, 2), real_c(0.));
      });
   }
}

void TestInnerOuterSweep()
{
   const Vector3< uint_t > cellsPerBlock(5, 6, 7);

   const uint_t numProcs = uint_t(mpi::MPIManager::instance()->numProcesses());
   const Vector3< uint_t > numBlocks = math::getFactors3D(numProcs);

   auto blocks = blockforest::createUniformBlockGrid(numBlocks[0], numBlocks[1], numBlocks[2], cellsPerBlock[0],
                                                     cellsPerBlock[1], cellsPerBlock[2], 1., true, true, true, true);

   const real_t unsetValue = real_c(-1);
   VectorField_T field{ *blocks, 1, unsetValue };

   auto sweepFactory = sweep::SweepFactory(blocks);
   const Cell innerOuterSplit(1, 1, 1);
   gen::SetCellCenter setCenterSweep(blocks, field);
   auto innerSweep = sweepFactory.commHidingSweep< sweep::INNER >(setCenterSweep, innerOuterSplit);
   auto outerSweep = sweepFactory.commHidingSweep< sweep::OUTER >(setCenterSweep, innerOuterSplit);

   const auto isInnerCell = [innerOuterSplit, cellsPerBlock](Cell c) {
      return c.x() >= innerOuterSplit.x() && c.x() < cell_idx_c(cellsPerBlock[0]) - innerOuterSplit.x() &&
             c.y() >= innerOuterSplit.y() && c.y() < cell_idx_c(cellsPerBlock[1]) - innerOuterSplit.y() &&
             c.z() >= innerOuterSplit.z() && c.z() < cell_idx_c(cellsPerBlock[2]) - innerOuterSplit.z();
   };

   for (auto& block : *blocks)
   {
      innerSweep(&block);
   }

   for (auto& block : *blocks)
   {
      FieldView fView{ field, block };
      sweep::forAllCells(exectag::Serial{}, *blocks, [&](Cell c) {
         if (isInnerCell(c))
         {
            const auto center = blocks->getBlockLocalCellCenter(block, c);
            testing::assert_close(fView(c, 0), center[0]);
            testing::assert_close(fView(c, 1), center[1]);
            testing::assert_close(fView(c, 2), center[2]);
         }
         else
         {
            testing::assert_close(fView(c, 0), unsetValue);
            testing::assert_close(fView(c, 1), unsetValue);
            testing::assert_close(fView(c, 2), unsetValue);
         }
      });
   }

   for (auto& block : *blocks)
   {
      outerSweep(&block);
   }

   for (auto& block : *blocks)
   {
      FieldView fView{ field, block };
      sweep::forAllCells(exectag::Serial{}, *blocks, [&](Cell c) {
         const auto center = blocks->getBlockLocalCellCenter(block, c);
         testing::assert_close(fView(c, 0), center[0]);
         testing::assert_close(fView(c, 1), center[1]);
         testing::assert_close(fView(c, 2), center[2]);
      });
   }
}

} // namespace

int main(int argc, char** argv)
{
   walberla::mpi::Environment env{ argc, argv };

   return walberla::v8::testing::TestsRunner({
      { "testCommunicationHiding", &TestCommunicationHiding::TestCommunicationHiding },
      { "testInnerOuterSplit", &TestCommunicationHiding::TestInnerOuterSweep }}).run(argc, argv);
}
