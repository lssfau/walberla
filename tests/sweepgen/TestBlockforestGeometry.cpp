#include "blockforest/all.h"

#include "core/all.h"

#include "field/all.h"

#include "gen/GeometryKernels.hpp"
#include "walberla/experimental/Sweep.hpp"

using namespace walberla;
namespace wex = walberla::experimental;

using Int64Field = field::GhostLayerField< int64_t, 3 >;
using RealField  = field::GhostLayerField< real_t, 3 >;

int main(int argc, char** argv)
{
   Environment env{ argc, argv };

   auto blocks = blockforest::createUniformBlockGrid(2, 2, 2, 4, 4, 4, 1.2);

   auto intFieldId  = field::addToStorage< Int64Field >(blocks, "intField");
   auto realFieldId = field::addToStorage< RealField >(blocks, "realField");

   wex::sweep::SerialSweeper sweeper{ blocks };

   {
      CellCentersGlobal cellCenters{ blocks, realFieldId };
      sweeper.sweep(cellCenters);

      sweeper.sweep([&](IBlock& block) {
         RealField& outp = *block.getData< RealField >(realFieldId);
         sweeper.forAllCells([&](Cell c) {
            Vector3< real_t > desired = blocks->getBlockLocalCellCenter(block, c);

            WALBERLA_CHECK_FLOAT_EQUAL(outp.get(c, 0), desired[0]);
            WALBERLA_CHECK_FLOAT_EQUAL(outp.get(c, 1), desired[1]);
            WALBERLA_CHECK_FLOAT_EQUAL(outp.get(c, 2), desired[2]);
         });
      });
   }

   {
      CellCentersLocal cellCenters{ blocks, realFieldId };
      sweeper.sweep(cellCenters);

      sweeper.sweep([&](IBlock& block) {
         RealField& outp = *block.getData< RealField >(realFieldId);
         AABB blockAbb = block.getAABB();
         sweeper.forAllCells([&](Cell c) {
            Vector3< real_t > desired = blocks->getBlockLocalCellCenter(block, c) - blockAbb.min();

            WALBERLA_CHECK_FLOAT_EQUAL(outp.get(c, 0), desired[0]);
            WALBERLA_CHECK_FLOAT_EQUAL(outp.get(c, 1), desired[1]);
            WALBERLA_CHECK_FLOAT_EQUAL(outp.get(c, 2), desired[2]);
         });
      });
   }

   {
    CellIndicesLocal localIndices{ intFieldId };
      sweeper.sweep(localIndices);

      sweeper.sweep([&](IBlock& block) {
         Int64Field& outp = *block.getData< Int64Field >(intFieldId);
         sweeper.forAllCells([&](Cell c) {
            WALBERLA_CHECK_EQUAL(outp.get(c, 0), c[0]);
            WALBERLA_CHECK_EQUAL(outp.get(c, 1), c[1]);
            WALBERLA_CHECK_EQUAL(outp.get(c, 2), c[2]);
         });
      });
   }

   {
    CellIndicesGlobal glocalIndices{ blocks, intFieldId };
      sweeper.sweep(glocalIndices);

      sweeper.sweep([&](IBlock& block) {
         Int64Field& outp = *block.getData< Int64Field >(intFieldId);
         sweeper.forAllCells([&](Cell c) {
            Cell globalCell{ c };
            blocks->transformBlockLocalToGlobalCell(globalCell, block);

            WALBERLA_CHECK_EQUAL(outp.get(c, 0), globalCell[0]);
            WALBERLA_CHECK_EQUAL(outp.get(c, 1), globalCell[1]);
            WALBERLA_CHECK_EQUAL(outp.get(c, 2), globalCell[2]);
         });
      });
   }

   /* The same with cell intervals */

   CellInterval ci { Cell {1, 2, 1}, Cell (3, 3, 2) };

   {
      CellCentersGlobal cellCenters{ blocks, realFieldId };
      sweeper.sweep([&](IBlock * b){ cellCenters.runOnCellInterval(b, ci); });

      sweeper.sweep([&](IBlock& block) {
         RealField& outp = *block.getData< RealField >(realFieldId);
         sweeper.forAllCells(ci, [&](Cell c) {
            Vector3< real_t > desired = blocks->getBlockLocalCellCenter(block, c);

            WALBERLA_CHECK_FLOAT_EQUAL(outp.get(c, 0), desired[0]);
            WALBERLA_CHECK_FLOAT_EQUAL(outp.get(c, 1), desired[1]);
            WALBERLA_CHECK_FLOAT_EQUAL(outp.get(c, 2), desired[2]);
         });
      });
   }

   {
      CellCentersLocal cellCenters{ blocks, realFieldId };
      sweeper.sweep([&](IBlock * b){ cellCenters.runOnCellInterval(b, ci); });

      sweeper.sweep([&](IBlock& block) {
         RealField& outp = *block.getData< RealField >(realFieldId);
         AABB blockAbb = block.getAABB();
         sweeper.forAllCells(ci, [&](Cell c) {
            Vector3< real_t > desired = blocks->getBlockLocalCellCenter(block, c) - blockAbb.min();

            WALBERLA_CHECK_FLOAT_EQUAL(outp.get(c, 0), desired[0]);
            WALBERLA_CHECK_FLOAT_EQUAL(outp.get(c, 1), desired[1]);
            WALBERLA_CHECK_FLOAT_EQUAL(outp.get(c, 2), desired[2]);
         });
      });
   }

   {
    CellIndicesLocal localIndices{ intFieldId };
      sweeper.sweep([&](IBlock * b){ localIndices.runOnCellInterval(b, ci); });

      sweeper.sweep([&](IBlock& block) {
         Int64Field& outp = *block.getData< Int64Field >(intFieldId);
         sweeper.forAllCells(ci, [&](Cell c) {
            WALBERLA_CHECK_EQUAL(outp.get(c, 0), c[0]);
            WALBERLA_CHECK_EQUAL(outp.get(c, 1), c[1]);
            WALBERLA_CHECK_EQUAL(outp.get(c, 2), c[2]);
         });
      });
   }

   {
    CellIndicesGlobal glocalIndices{ blocks, intFieldId };
    sweeper.sweep([&](IBlock * b){ glocalIndices.runOnCellInterval(b, ci); });

      sweeper.sweep([&](IBlock& block) {
         Int64Field& outp = *block.getData< Int64Field >(intFieldId);
         sweeper.forAllCells(ci, [&](Cell c) {
            Cell globalCell{ c };
            blocks->transformBlockLocalToGlobalCell(globalCell, block);

            WALBERLA_CHECK_EQUAL(outp.get(c, 0), globalCell[0]);
            WALBERLA_CHECK_EQUAL(outp.get(c, 1), globalCell[1]);
            WALBERLA_CHECK_EQUAL(outp.get(c, 2), globalCell[2]);
         });
      });
   }

   return EXIT_SUCCESS;
}