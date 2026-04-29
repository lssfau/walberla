#include "gen/TestSweepGeneration.hpp"

#include "blockforest/all.h"

#include "core/all.h"

#include "field/all.h"

#include "walberla/v8/Domain.hpp"
#include "walberla/v8/Memory.hpp"
#include "walberla/v8/Sweep.hpp"
#include "walberla/v8/Testing.hpp"

namespace
{

using namespace walberla;
using namespace walberla::v8;

void test1DField()
{
   using Field_T = field::GhostLayerField< real_t, 3 >;

   auto blocks = blockforest::createUniformBlockGrid(1, 1, 1, 32, 1, 1, 1.0, true);

   auto fId_fzyx = field::addToStorage< Field_T >(blocks, "f", -1., field::fzyx);
   auto fId_zyxf = field::addToStorage< Field_T >(blocks, "f", -1., field::zyxf);

   gen::Kernels1D::Set_fzyx set_fzyx{ blocks, fId_fzyx };
   gen::Kernels1D::Set_zyxf set_zyxf{ blocks, fId_zyxf };

   for (auto& b : *blocks)
   {
      auto check = [&](BlockDataID fId) {
         Field_T& f = *b.getData< Field_T >(fId);

         CellInterval interior = f.xyzSize();

         for (auto fIt = f.beginWithGhostLayer(); fIt != f.end(); ++fIt)
         {
            Vector3< real_t > cc{ blocks->getBlockLocalCellCenter(b, fIt.cell()) };
            if (interior.contains(fIt.cell())) { WALBERLA_CHECK_FLOAT_EQUAL(*fIt, cc[uint_c(fIt.f())]); }
            else
            {
               WALBERLA_CHECK_FLOAT_EQUAL(*fIt, -1.);
            }
         }
      };

      set_fzyx(&b);
      check(fId_fzyx);

      set_zyxf(&b);
      check(fId_zyxf);
   }
}

void test2DField()
{
   using Field_T = field::GhostLayerField< real_t, 3 >;

   auto blocks = blockforest::createUniformBlockGrid(1, 1, 1, 32, 32, 1, 1.0, true);

   auto fId_fzyx = field::addToStorage< Field_T >(blocks, "f", -1., field::fzyx);
   auto fId_zyxf = field::addToStorage< Field_T >(blocks, "f", -1., field::zyxf);

   gen::Kernels2D::Set_fzyx set_fzyx{ blocks, fId_fzyx };
   gen::Kernels2D::Set_zyxf set_zyxf{ blocks, fId_zyxf };

   for (auto& b : *blocks)
   {
      auto check = [&](BlockDataID fId) {
         Field_T& f = *b.getData< Field_T >(fId);

         CellInterval interior = f.xyzSize();

         for (auto fIt = f.beginWithGhostLayer(); fIt != f.end(); ++fIt)
         {
            Vector3< real_t > cc{ blocks->getBlockLocalCellCenter(b, fIt.cell()) };
            if (interior.contains(fIt.cell())) { WALBERLA_CHECK_FLOAT_EQUAL(*fIt, cc[uint_c(fIt.f())]); }
            else
            {
               WALBERLA_CHECK_FLOAT_EQUAL(*fIt, -1.);
            }
         }
      };

      set_fzyx(&b);
      check(fId_fzyx);

      set_zyxf(&b);
      check(fId_zyxf);
   }
}

void testV8Fields()
{
   auto blocks = blockforest::createUniformBlockGrid(1, 1, 1, 16, 16, 16, 1.0, true);

   {
      const memory::Field< double, 3, memtag::host > f{ *blocks };
      const memory::Field< double, 1, memtag::host > g{ *blocks };
      const memory::Field< double, 1, memtag::host > gDesired{ *blocks };

      gen::TestExperimentalFields::InitializeF initF{ blocks, f };
      gen::TestExperimentalFields::PointwiseSum pointwiseSum{ f, g };

      const sweep::SerialSweeper sweeper{ blocks };

      for (auto& b : *blocks)
      {
         initF(&b);
         pointwiseSum(&b);

         FieldView gDesiredView{ gDesired, b };

         sweeper.forAllCells([&](auto cell) {
            const Vector3< double > cellCenter = blocks->getBlockLocalCellCenter(b, cell);
            gDesiredView(cell)                 = cellCenter[0] + cellCenter[1] + cellCenter[2];
         });

         FieldView gView{ g, b };

         testing::assert_allclose(gView, gDesiredView);
      }
   }

   {
      const memory::Field< double, 3, memtag::host > f{ *blocks, 0 };
      const memory::Field< double, 1, memtag::host > g{ *blocks, 0, -1. };
      const memory::Field< double, 1, memtag::host > gDesired{ *blocks, 0 };

      gen::TestExperimentalFields::InitializeF initF{ blocks, f };
      gen::TestExperimentalFields::PointwiseSum pointwiseSum{ f, g };

      const sweep::SerialSweeper sweeper{ blocks };
      CellInterval ci{ { 3, 4, 5 }, { 12, 13, 14 } };

      for (auto& b : *blocks)
      {
         initF.runOnCellInterval(&b, ci);
         pointwiseSum.runOnCellInterval(&b, ci);

         FieldView gDesiredView{ gDesired, b };

         sweeper.forAllCells([&](auto cell) {
            if (ci.contains(cell))
            {
               const Vector3< double > cellCenter = blocks->getBlockLocalCellCenter(b, cell);
               gDesiredView(cell)                 = cellCenter[0] + cellCenter[1] + cellCenter[2];
            }
            else
            {
               gDesiredView(cell) = -1.;
            }
         });

         FieldView gView{ g, b };

         testing::assert_allclose(gView, gDesiredView);
      }
   }
}

void testV8FieldSwaps()
{
   auto blocks = blockforest::createUniformBlockGrid(1, 1, 1, 16, 16, 16, 1.0, true);
   constexpr double initValue{ 1.3 };
   const memory::Field< double, 1, memtag::host > f{ *blocks, 1, initValue };
   const memory::Field< double, 1, memtag::host > fExpected{ *blocks, 1 };

   gen::TestFieldSwaps::TestStencil stencilSweep{ f };
   sweep::SerialSweeper sweeper{ blocks };

   for (auto& block : *blocks)
   {
      auto& origBuffer     = f.bufferSystem().buffer(block);
      double* origPtr      = origBuffer.data();
      double* origAllocPtr = origBuffer.allocData();

      stencilSweep(&block);

      double* shadowPtr      = origBuffer.data();
      double* shadowAllocPtr = origBuffer.allocData();

      testing::assert_inequal(origPtr, shadowPtr);
      testing::assert_inequal(origAllocPtr, shadowAllocPtr);

      FieldView fExpectedView{ fExpected, block };

      sweeper.forAllCells([&](auto cell) { fExpectedView(cell) = 6.0 * initValue; });

      testing::assert_allclose(FieldView(f, block), fExpectedView);

      stencilSweep(&block);

      // After two swaps we're back to the original

      testing::assert_equal(origPtr, origBuffer.data());
      testing::assert_equal(origAllocPtr, origBuffer.allocData());
   }
}

void test1DFieldV8()
{
   using Field_T = memory::Field< real_t, 3, memtag::host >;

   auto blocks = blockforest::createUniformBlockGrid(1, 1, 1, 32, 1, 1, 1.0, true);

   Field_T f_fzyx{ *blocks, 1, -1., field::fzyx };
   Field_T f_zyxf{ *blocks, 1, -1., field::zyxf };

   gen::v8::Kernels1D::Set_fzyx set_fzyx{ blocks, f_fzyx };
   gen::v8::Kernels1D::Set_zyxf set_zyxf{ blocks, f_zyxf };

   for (auto& b : *blocks)
   {
      const domain::BlockGeometry bGeom{ *blocks, b };

      auto check = [&](Field_T& f) {
         FieldView fView{ f, b };
         auto interior = f.slices().interior();

         sweep::forAllCells(sweep::exectag::Serial{}, f.slices().withGhostLayers(), [&](Cell c) {
            auto cc = bGeom.localGrid().cellCenter(c);
            for (size_t i = 0; i < 3; ++i)
            {
               if (interior.contains(c)) { testing::assert_close(fView(c, cell_idx_c(i)), cc[i]); }
               else
               {
                  testing::assert_close(fView(c, cell_idx_c(i)), -1.);
               }
            }
         });
      };

      set_fzyx(&b);
      check(f_fzyx);

      set_zyxf(&b);
      check(f_zyxf);
   }
}

void test2DFieldV8()
{
   using Field_T = memory::Field< real_t, 3, memtag::host >;

   auto blocks = blockforest::createUniformBlockGrid(1, 1, 1, 32, 32, 1, 1.0, true);

   Field_T f_fzyx{ *blocks, 1, -1., field::fzyx };
   Field_T f_zyxf{ *blocks, 1, -1., field::zyxf };

   gen::v8::Kernels2D::Set_fzyx set_fzyx{ blocks, f_fzyx };
   gen::v8::Kernels2D::Set_zyxf set_zyxf{ blocks, f_zyxf };

   for (auto& b : *blocks)
   {
      const domain::BlockGeometry bGeom{ *blocks, b };

      auto check = [&](Field_T& f) {
         FieldView fView{ f, b };
         auto interior = f.slices().interior();

         sweep::forAllCells(sweep::exectag::Serial{}, f.slices().withGhostLayers(), [&](Cell c) {
            auto cc = bGeom.localGrid().cellCenter(c);
            for (size_t i = 0; i < 3; ++i)
            {
               if (interior.contains(c)) { testing::assert_close(fView(c, cell_idx_c(i)), cc[i]); }
               else
               {
                  testing::assert_close(fView(c, cell_idx_c(i)), -1.);
               }
            }
         });
      };

      set_fzyx(&b);
      check(f_fzyx);

      set_zyxf(&b);
      check(f_zyxf);
   }
}

} // namespace

int main(int argc, char** argv)
{
   walberla::mpi::Environment env{ argc, argv };

   return walberla::v8::testing::TestsRunner( //
             {
                { "test1DField", &test1DField },
                { "test2DField", &test2DField },
                { "testV8Fields", &testV8Fields },
                { "testV8FieldSwaps", &testV8FieldSwaps },
                { "test1DFieldV8", &test1DFieldV8 },
                { "test2DFieldV8", &test2DFieldV8 },
             })
      .run(argc, argv);
}
