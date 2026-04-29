#include "blockforest/all.h"

#include "core/all.h"

#include <span>

#include "walberla/v8/Memory.hpp"
#include "walberla/v8/Sweep.hpp"
#include "walberla/v8/Testing.hpp"

namespace
{
using namespace walberla;
using namespace walberla::v8;

void testFieldCreation()
{
   auto blocks = blockforest::createUniformBlockGrid(1, 1, 1, 2, 3, 4, 1.,
                                                     /* oneBlockPerProcess */ false);

   {
      const Field< double, 3, memtag::stdmem > f{ *blocks, 1, 0., field::fzyx };

      testing::assert_equal(f.numGhostLayers(), 1UL);
      testing::assert_equal(f.layout(), field::fzyx);

      auto& bsys     = f.bufferSystem();
      auto& indexing = bsys.indexing();

      testing::assert_range_equal(indexing.shape(), std::array< size_t, 4 >{ 2, 3, 4, 3 });
      testing::assert_range_equal(indexing.allocShape(), std::array< size_t, 4 >{ 4, 5, 6, 3 });
      testing::assert_range_equal(indexing.allocOffsets(), std::array< size_t, 4 >{ 1, 1, 1, 0 });
      testing::assert_range_equal(indexing.strides(), std::array< size_t, 4 >{ 1, 4, 20, 120 });

      for (auto& block : *blocks)
      {
         auto& buf = bsys.buffer(block);
         testing::assert_allclose(buf.linearAllocView(), 0.);
      }
   }

   {
      const Field< double, 3, memtag::stdmem > f{ *blocks, 1, 2.3, field::zyxf };

      testing::assert_equal(f.numGhostLayers(), 1UL);
      testing::assert_equal(f.layout(), field::zyxf);

      auto& bsys     = f.bufferSystem();
      auto& indexing = bsys.indexing();

      testing::assert_range_equal(indexing.shape(), std::array< size_t, 4 >{ 2, 3, 4, 3 });
      testing::assert_range_equal(indexing.allocShape(), std::array< size_t, 4 >{ 4, 5, 6, 3 });
      testing::assert_range_equal(indexing.allocOffsets(), std::array< size_t, 4 >{ 1, 1, 1, 0 });
      testing::assert_range_equal(indexing.strides(), std::array< size_t, 4 >{ 3, 12, 60, 1 });

      for (auto& block : *blocks)
      {
         auto& buf = bsys.buffer(block);
         testing::assert_allclose(buf.linearAllocView(), 2.3);
      }
   }
}

template< typename FieldType, typename BlockRef, typename ViewType >
concept CanCreateView = requires(const FieldType& field, BlockRef blockref) {
   { ViewType{ field, blockref } } -> std::same_as< ViewType >;
};

void testConstCorrectness()
{
   using FieldType      = Field< double, 3, memory::memtag::stdmem >;
   using ConstFieldType = Field< const double, 3, memory::memtag::stdmem >;

   using ViewType      = FieldView< double, 3, memory::memtag::stdmem >;
   using ConstViewType = FieldView< const double, 3, memory::memtag::stdmem >;

   using BlockRef      = IBlock&;
   using ConstBlockRef = const IBlock&;

   static_assert(CanCreateView< FieldType, BlockRef, ViewType >);
   static_assert(CanCreateView< FieldType, BlockRef, ConstViewType >);
   static_assert(!CanCreateView< FieldType, ConstBlockRef, ViewType >);
   static_assert(CanCreateView< FieldType, ConstBlockRef, ConstViewType >);

   static_assert(!CanCreateView< ConstFieldType, BlockRef, ViewType >);
   static_assert(CanCreateView< ConstFieldType, BlockRef, ConstViewType >);
   static_assert(!CanCreateView< ConstFieldType, ConstBlockRef, ViewType >);
   static_assert(CanCreateView< ConstFieldType, ConstBlockRef, ConstViewType >);

   static_assert(std::convertible_to< FieldType, ConstFieldType >);
   static_assert(!std::convertible_to< ConstFieldType, ConstViewType >);

   static_assert(std::convertible_to< ViewType, ConstViewType >);
   static_assert(!std::convertible_to< ConstViewType, ViewType >);
}

void testFieldView()
{
   static_assert(memory::IFieldView< FieldView< double, 3, memtag::stdmem > >);

   auto blocks = blockforest::createUniformBlockGrid(1, 1, 1, 2, 3, 4, 1.,
                                                     /* oneBlockPerProcess */ false);

   {
      const Field< double, 3, memtag::stdmem > f{ *blocks, 1, 0., field::fzyx };

      for (auto& block : *blocks)
      {
         FieldView fView{ f, block };

         testing::assert_equal(fView.numGhostLayers(), 1UL);
         testing::assert_equal(fView.layout(), field::fzyx);

         sweep::forAllCells(exectag::Serial{}, fView.slices().interior(), [&](Cell cell) {
            const Vector3< double > cellCenter{ blocks->getBlockLocalCellCenter(block, cell) };
            fView(cell, 0) = cellCenter[0];
            fView(cell, 1) = cellCenter[1];
            fView(cell, 2) = cellCenter[2];
         });
      }

      for (auto& block : *blocks)
      {
         ConstFieldView fView{ f, block };

         testing::assert_equal(fView.numGhostLayers(), 1UL);
         testing::assert_equal(fView.layout(), field::fzyx);

         sweep::forAllCells(exectag::Serial{}, fView.slices().interior(), [&](Cell cell) {
            const Vector3< double > cellCenter{ blocks->getBlockLocalCellCenter(block, cell) };
            const Vector3< double > actual{ fView(cell, 0), fView(cell, 1), fView(cell, 2) };
            testing::assert_allclose(std::span{ actual.data(), 3 }, std::span{ cellCenter.data(), 3 });
         });
      }
   }
}

void testScalarFieldView()
{
   auto blocks = blockforest::createUniformBlockGrid(1, 1, 1, 2, 3, 4, 1.,
                                                     /* oneBlockPerProcess */ false);

   {
      const Field< double, 1, memtag::stdmem > f{ *blocks, 1, 0., field::fzyx };

      for (auto& block : *blocks)
      {
         FieldView fView{ f, block };

         testing::assert_equal(fView.numGhostLayers(), 1UL);
         testing::assert_equal(fView.layout(), field::fzyx);

         sweep::forAllCells(exectag::Serial{}, fView.slices().interior(), [&](Cell cell) {
            const Vector3< double > cellCenter{ blocks->getBlockLocalCellCenter(block, cell) };
            fView(cell) = cellCenter[0] * cellCenter[1] * cellCenter[2];
         });
      }

      for (auto& block : *blocks)
      {
         ConstFieldView fView{ f, block };

         testing::assert_equal(fView.numGhostLayers(), 1UL);
         testing::assert_equal(fView.layout(), field::fzyx);

         sweep::forAllCells(exectag::Serial{}, fView.slices().interior(), [&](Cell cell) {
            const Vector3< double > cellCenter{ blocks->getBlockLocalCellCenter(block, cell) };
            const double expected{ cellCenter[0] * cellCenter[1] * cellCenter[2] };
            const double actual{ fView(cell) };
            testing::assert_close(actual, expected);
         });
      }
   }
}

void testHostLineAlignment()
{
   auto blocks = blockforest::createUniformBlockGrid(1, 1, 1, 2, 3, 4, 1.,
                                                     /* oneBlockPerProcess */ false);

   constexpr size_t expectedAlignment{ memory::MemoryTraits< memtag::host, double >::alignment() };

   for (auto numGhostLayers : std::views::iota(size_t(0), size_t(4)))
   {
      // FZYX: For all values of y, z, and q, address at `x == 0` must be aligned
      {
         const Field< double, 3, memtag::host > f{ *blocks, numGhostLayers, 0., field::fzyx };

         for (auto& block : *blocks)
         {
            FieldView fView{ f, block };

            for (cell_idx_t z = 0; z < cell_idx_c(blocks->getNumberOfZCellsPerBlock()); ++z)
               for (cell_idx_t y = 0; y < cell_idx_c(blocks->getNumberOfYCellsPerBlock()); ++y)
                  for (cell_idx_t q = 0; q < 3; ++q)
                  {
                     const size_t lineStart = (size_t) std::addressof(fView(0, y, z, q));
                     testing::assert_equal(lineStart % expectedAlignment, size_t(0));
                  }
         }
      }
   }

   // ZYXF: No line alignment requirements, but entries of each cell must be packed
   {
      const Field< double, 3, memtag::host > f{ *blocks, 1, 0., field::zyxf };

      for (auto& block : *blocks)
      {
         FieldView fView{ f, block };

         for (cell_idx_t z = 0; z < 4; ++z)
            for (cell_idx_t y = 0; y < 4; ++y)
            {
               const std::ptrdiff_t cellSize{ std::addressof(fView(1, y, z, 0)) - std::addressof(fView(0, y, z, 0)) };

               testing::assert_equal(cellSize, 3L);
            }
      }
   }
}

} // namespace

int main(int argc, char** argv)
{
   walberla::mpi::Environment env{ argc, argv };

   return walberla::v8::testing::TestsRunner({
                                                          { "testFieldCreation", &testFieldCreation },
                                                          { "testConstCorrectness", &testConstCorrectness },
                                                          { "testFieldView", &testFieldView },
                                                          { "testScalarFieldView", &testScalarFieldView },
                                                          { "testHostLineAlignment", &testHostLineAlignment },
                                                       })
      .run(argc, argv);
}
