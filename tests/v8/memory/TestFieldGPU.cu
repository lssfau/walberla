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

template< typename FieldType, typename BlockRef, typename ViewType >
concept CanCreateView = requires(const FieldType& field, BlockRef blockref) {
   { ViewType{ field, blockref } } -> std::same_as< ViewType >;
};

void testConstCorrectness()
{
   using FieldType      = Field< double, 3, memory::memtag::unified >;
   using ConstFieldType = Field< const double, 3, memory::memtag::unified >;

   using ViewType      = FieldView< double, 3, memory::memtag::unified >;
   using ConstViewType = FieldView< const double, 3, memory::memtag::unified >;

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

void testFieldViewUnified()
{
   static_assert(memory::IFieldView< FieldView< double, 3, memtag::unified > >);

   auto blocks = blockforest::createUniformBlockGrid(1, 1, 1, 2, 3, 4, 1.,
                                                     /* oneBlockPerProcess */ false);

   {
      const Field< double, 3, memtag::unified > f{ *blocks, 1, 0., field::fzyx };

      for (auto& block : *blocks)
      {
         FieldView fView{ f, block };

         testing::assert_equal(fView.numGhostLayers(), 1UL);
         testing::assert_equal(fView.layout(), field::fzyx);

         sweep::forAllCells(exectag::GPU{}, fView.slices().interior(), 
         [fView] WALBERLA_HOST_DEVICE (Cell cell) {
            fView(cell, 0) = cell.x();
            fView(cell, 1) = cell.y();
            fView(cell, 2) = cell.z();
         });

         sweep::sync();
      }

      for (auto& block : *blocks)
      {
         ConstFieldView fView{ f, block };

         testing::assert_equal(fView.numGhostLayers(), 1UL);
         testing::assert_equal(fView.layout(), field::fzyx);

         sweep::forAllCells(exectag::Serial{}, fView.slices().interior(), 
         [&](Cell cell) {
            testing::assert_equal(cell_idx_c(fView(cell, 0)), cell.x());
            testing::assert_equal(cell_idx_c(fView(cell, 1)), cell.y());
            testing::assert_equal(cell_idx_c(fView(cell, 2)), cell.z());
         });
      }
   }
}

void testCopyDeviceToPinned()
{
   static_assert(memory::IFieldView< FieldView< double, 3, memtag::device > >);
   static_assert(memory::IFieldView< FieldView< double, 3, memtag::pinned > >);

   auto blocks = blockforest::createUniformBlockGrid(1, 1, 1, 5, 5, 5, 1.,
                                                     /* oneBlockPerProcess */ false);

   {
      const Field< double, 3, memtag::device > fDev{ *blocks, 1, 42.42, field::fzyx };
      const Field< double, 3, memtag::pinned > fPin{ *blocks, 1, 0.,    field::fzyx };

      for (auto& block : *blocks)
      {
         FieldView fDevView{ fDev, block };

         sweep::forAllCells(exectag::GPU{}, fDevView.slices().interior(), 
         [fDevView] WALBERLA_HOST_DEVICE (Cell cell) {
            fDevView(cell, 0) += cell.x();
            fDevView(cell, 1) += cell.y();
            fDevView(cell, 2) += cell.z();
         });

         sweep::forAllCells(exectag::GPU{}, fDevView.slices().ghostSlice(stencil::Direction::N, 1, true), 
         [fDevView] WALBERLA_HOST_DEVICE (Cell cell) {
            fDevView(cell, 0) += cell.x() + 7.0;
            fDevView(cell, 1) += cell.y() + 8.0;
            fDevView(cell, 2) += cell.z() + 9.0;
         });
         sweep::sync();
      }

      for (auto& block : *blocks)
      {
         ConstFieldView fDevView{ fDev, block };
         FieldView fPinView{ fPin, block };
         fieldCpy(fPinView, fDevView);
      }

      for (auto& block : *blocks)
      {
         ConstFieldView fPinView{ fPin, block };

         sweep::forAllCells(exectag::Serial{}, fPinView.slices().interior(), 
         [&](Cell cell) {
            testing::assert_close(fPinView(cell, 0), cell.x() + 42.42);
            testing::assert_close(fPinView(cell, 1), cell.y() + 42.42);
            testing::assert_close(fPinView(cell, 2), cell.z() + 42.42);
         });

         sweep::forAllCells(exectag::Serial{}, fPinView.slices().ghostSlice(stencil::Direction::N, 1, true), 
         [&](Cell cell) {
            testing::assert_close(fPinView(cell, 0), cell.x() + 42.42 + 7.0);
            testing::assert_close(fPinView(cell, 1), cell.y() + 42.42 + 8.0);
            testing::assert_close(fPinView(cell, 2), cell.z() + 42.42 + 9.0);
         });
      }
   }
}

void testScalarFieldView()
{
   auto blocks = blockforest::createUniformBlockGrid(1, 1, 1, 2, 3, 4, 1.,
                                                     /* oneBlockPerProcess */ false);

   {
      const Field< double, 1, memtag::unified > f{ *blocks, 1, 42.42, field::fzyx };

      for (auto& block : *blocks)
      {
         FieldView fView{ f, block };

         testing::assert_equal(fView.numGhostLayers(), 1UL);
         testing::assert_equal(fView.layout(), field::fzyx);

         sweep::forAllCells(exectag::GPU{}, fView.slices().interior(), 
         [fView] WALBERLA_HOST_DEVICE (Cell cell) {
            fView(cell) += cell[0] * cell[1] * cell[2];
         });

         sweep::sync();
      }

      for (auto& block : *blocks)
      {
         ConstFieldView fView{ f, block };

         testing::assert_equal(fView.numGhostLayers(), 1UL);
         testing::assert_equal(fView.layout(), field::fzyx);

         sweep::forAllCells(exectag::Serial{}, fView.slices().interior(), 
         [&](Cell cell) {
            const double expected{ cell[0] * cell[1] * cell[2] + 42.42 };
            const double actual{ fView(cell) };
            testing::assert_close(actual, expected);
         });
      }
   }
}

template< IFieldView FVT >
__global__ void checkAlignmentX(FVT fView, size_t expectedAlignment, bool* flag){
   if (blockIdx.x != 0 || threadIdx.x != 0) { return; } 
   *flag = true;
   
   for (cell_idx_t q = 0; q < fView.shape()[3]; ++q)
      for (cell_idx_t z = 0; z < fView.shape()[2]; ++z)
         for (cell_idx_t y = 0; y < fView.shape()[1]; ++y)
         {
            const size_t lineStart = (size_t) &(fView(0, y, z, q));
            if (lineStart % expectedAlignment) { *flag = false; }
         }
}

template< IFieldView FVT >
__global__ void checkPackQ(FVT fView, bool* flag){
   if (blockIdx.x != 0 || threadIdx.x != 0) { return; } 
   *flag = true;

   for (cell_idx_t z = 0; z < fView.shape()[2]; ++z)
      for (cell_idx_t y = 0; y < fView.shape()[1]; ++y)
      {
         const std::ptrdiff_t cellSize{ &(fView(1, y, z, 0)) - &(fView(0, y, z, 0)) };
         if (cellSize !=  fView.shape()[3]) { *flag = false; }
      }
}

void testDeviceLineAlignment()
{
   auto blocks = blockforest::createUniformBlockGrid(1, 1, 1, 2, 3, 4, 1.,
                                                     /* oneBlockPerProcess */ false);

   constexpr size_t expectedAlignment{ memory::MemoryTraits< memtag::unified, double >::alignment() };

   for (auto numGhostLayers : std::views::iota(size_t(0), size_t(4)))
   {
      // FZYX: For all values of y, z, and q, address at `x == 0` must be aligned
      {
         const Field< double, 3, memtag::unified > f{ *blocks, numGhostLayers, 0., field::fzyx };

         for (auto& block : *blocks)
         {
            FieldView fView{ f, block };

            bool* flag = nullptr;
            WALBERLA_GPU_CHECK(gpuMallocManaged(&flag, sizeof(bool)));
            checkAlignmentX<<<1,1>>>(fView, expectedAlignment, flag);
            WALBERLA_GPU_CHECK(gpuDeviceSynchronize());
            testing::assert_true(*flag);
         }
      }
   }

   // ZYXF: No line alignment requirements, but entries of each cell must be packed
   {
      const Field< double, 3, memtag::unified > f{ *blocks, 1, 0., field::zyxf };

      for (auto& block : *blocks)
      {
         FieldView fView{ f, block };

         bool* flag = nullptr;
         WALBERLA_GPU_CHECK(gpuMallocManaged(&flag, sizeof(bool)));
         checkPackQ<<<1,1>>>(fView, flag);
         WALBERLA_GPU_CHECK(gpuDeviceSynchronize());
         testing::assert_true(*flag);
      }
   }
}


} // namespace

int main(int argc, char** argv)
{
   walberla::mpi::Environment env{ argc, argv };

   return walberla::v8::testing::TestsRunner({
                                                          { "testConstCorrectness", &testConstCorrectness },
                                                          { "testFieldViewUnified", &testFieldViewUnified },
                                                          { "testCopyDeviceToPinned", &testCopyDeviceToPinned },
                                                          { "testScalarFieldView", &testScalarFieldView },
                                                          { "testDeviceLineAlignment", &testDeviceLineAlignment },
                                                       })
      .run(argc, argv);
}
