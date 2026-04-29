
#include "blockforest/all.h"

#include "core/all.h"

#include "walberla/v8/Memory.hpp"
#include "walberla/v8/Sweep.hpp"
#include "walberla/v8/Testing.hpp"

namespace
{
using namespace walberla;
using namespace walberla::v8;

void testBufferIndexing()
{
   auto blocks = blockforest::createUniformBlockGrid(1, 1, 1, 4, 4, 4, 1.,
                                                     /* oneBlockPerProcess */ false);

   auto buffers = memory::BufferSystem< double, 3, memory::memtag::stdmem >::create(*blocks, { 4, 4, 4 })
                     .linearizationOrder({ 2, 1, 0 })
                     .build();

   auto& indexing = buffers.indexing();

   testing::assert_range_equal(indexing.shape(), std::array< size_t, 3 >{ 4, 4, 4 });
   testing::assert_range_equal(indexing.allocShape(), std::array< size_t, 3 >{ 4, 4, 4 });
   testing::assert_range_equal(indexing.allocOffsets(), std::array< size_t, 3 >{ 0, 0, 0 });
   testing::assert_range_equal(indexing.strides(), std::array< size_t, 3 >{ 16, 4, 1 });

   for (auto& block : *blocks)
   {
      auto bufView = buffers.view(block);

      testing::assert_range_equal(bufView.shape(), std::array< size_t, 3 >{ 4, 4, 4 });
      testing::assert_range_equal(bufView.allocShape(), std::array< size_t, 3 >{ 4, 4, 4 });
      testing::assert_range_equal(bufView.allocOffsets(), std::array< size_t, 3 >{ 0, 0, 0 });
      testing::assert_range_equal(bufView.strides(), std::array< size_t, 3 >{ 16, 4, 1 });

      testing::assert_equal(std::addressof(bufView(0, 0, 1)) - std::addressof(bufView(0, 0, 0)), std::ptrdiff_t(1));

      testing::assert_equal(std::addressof(bufView(0, 1, 0)) - std::addressof(bufView(0, 0, 0)), std::ptrdiff_t(4));

      testing::assert_equal(std::addressof(bufView(1, 0, 0)) - std::addressof(bufView(0, 0, 0)), std::ptrdiff_t(16));

      const size_t bufAllocSize = buffers.buffer(block).linearAllocSize();
      testing::assert_equal(bufAllocSize, indexing.linearBufferAllocSize());

      double* bufData      = buffers.buffer(block).data();
      double* bufAllocData = buffers.buffer(block).allocData();

      testing::assert_equal(bufView.data(), bufData);
      testing::assert_equal(bufView.data(), bufAllocData);
   }

   sweep::SerialSweeper sweeper{ blocks };
   for (auto& block : *blocks)
   {
      auto bView    = buffers.view(block);
      double* bData = buffers.buffer(block).data();

      sweeper.forAllCells([&](auto cell) {
         bView(cell.x(), cell.y(), cell.z()) = double(cell.x()) + 10. * double(cell.y()) + 100. * double(cell.z());
      });

      sweeper.forAllCells([&](auto cell) {
         const size_t linearIndex{ size_t(cell.x()) * 16UL + size_t(cell.y()) * 4UL + size_t(cell.z()) };
         const double expectedValue{ double(cell.x()) + 10. * double(cell.y()) + 100. * double(cell.z()) };
         testing::assert_close(bData[linearIndex], expectedValue);
      });
   }
}

void testBufferIndexingWithPadding()
{
   auto blocks = blockforest::createUniformBlockGrid(1, 1, 1, 4, 4, 4, 1.,
                                                     /* oneBlockPerProcess */ false);

   auto buffers = memory::BufferSystem< double, 3, memory::memtag::stdmem >::create(*blocks, { 4, 4, 4 })
                     .allocShape({ 6, 6, 8 })
                     .allocOffsets({ 1, 1, 2 })
                     .linearizationOrder({ 2, 1, 0 })
                     .build();

   auto& indexing = buffers.indexing();

   testing::assert_range_equal(indexing.shape(), std::array< size_t, 3 >{ 4, 4, 4 });
   testing::assert_range_equal(indexing.allocShape(), std::array< size_t, 3 >{ 6, 6, 8 });
   testing::assert_range_equal(indexing.allocOffsets(), std::array< size_t, 3 >{ 1, 1, 2 });
   testing::assert_range_equal(indexing.strides(), std::array< size_t, 3 >{ 48, 8, 1 });

   for (auto& block : *blocks)
   {
      auto bufView = buffers.view(block);
      auto& buf    = buffers.buffer(block);

      testing::assert_equal(bufView.data(), buf.data());
      testing::assert_equal(std::addressof(bufView(0, 0, 0)), buf.data());
      testing::assert_equal(std::addressof(bufView(-1, -1, -2)), buf.allocData());
   }
}

template< typename BSys, typename BlockRef, typename BView >
concept CanGetView = requires(const BSys& bsys, BlockRef blockref) {
   { bsys.view(blockref) } -> std::convertible_to< BView >;
};

template< typename BSys, typename BlockRef, typename Buf >
concept CanGetBufferRef = requires(const BSys& bsys, BlockRef blockref) {
   { bsys.buffer(blockref) } -> std::convertible_to< Buf& >;
};

void testConstCorrectness()
{
   using BSysType      = memory::BufferSystem< double, 3, memory::memtag::stdmem >;
   using ConstBSysType = memory::BufferSystem< const double, 3, memory::memtag::stdmem >;

   using BufferType =
      memory::buffer_system::LinearBuffer< double,
                                           memory::MemoryTraits< memory::memtag::stdmem, double >::AllocatorType >;

   using BlockRef      = IBlock&;
   using ConstBlockRef = const IBlock&;

   using BViewType      = memory::buffer_system::BufferView< double, 3, memory::memtag::stdmem >;
   using ConstBViewType = memory::buffer_system::BufferView< const double, 3, memory::memtag::stdmem >;

   static_assert(std::convertible_to< BViewType, ConstBViewType >);
   static_assert(std::convertible_to< BSysType, ConstBSysType >);
   static_assert(!std::convertible_to< ConstBViewType, BViewType >);
   static_assert(!std::convertible_to< ConstBSysType, BSysType >);

   static_assert(CanGetView< BSysType, BlockRef, BViewType >);
   static_assert(CanGetView< BSysType, BlockRef, ConstBViewType >);
   static_assert(!CanGetView< BSysType, ConstBlockRef, BViewType >);
   static_assert(CanGetView< BSysType, ConstBlockRef, ConstBViewType >);

   static_assert(!CanGetView< ConstBSysType, BlockRef, BViewType >);
   static_assert(CanGetView< ConstBSysType, BlockRef, ConstBViewType >);
   static_assert(!CanGetView< ConstBSysType, ConstBlockRef, BViewType >);
   static_assert(CanGetView< ConstBSysType, ConstBlockRef, ConstBViewType >);

   static_assert(CanGetBufferRef< BSysType, BlockRef, BufferType >);
   static_assert(CanGetBufferRef< BSysType, BlockRef, const BufferType >);
   static_assert(!CanGetBufferRef< BSysType, ConstBlockRef, BufferType >);
   static_assert(CanGetBufferRef< BSysType, ConstBlockRef, const BufferType >);

   static_assert(!CanGetBufferRef< ConstBSysType, BlockRef, BufferType >);
   static_assert(CanGetBufferRef< ConstBSysType, BlockRef, const BufferType >);
   static_assert(!CanGetBufferRef< ConstBSysType, ConstBlockRef, BufferType >);
   static_assert(CanGetBufferRef< ConstBSysType, ConstBlockRef, const BufferType >);

   auto blocks = blockforest::createUniformBlockGrid(1, 1, 1, 4, 4, 4, 1.,
                                                     /* oneBlockPerProcess */ false);

   auto buffers  = BSysType::create(*blocks, { 4, 4, 4 }).build();
   auto buffers2 = BSysType::create(*blocks, { 4, 4, 4 }).build();

   // non-const to const copy init
   [[maybe_unused]] ConstBSysType cbuffers{ buffers };

   testing::assert_equal(uint_t(buffers.buffersId()), uint_t(cbuffers.buffersId()));

   // non-const to const assignment
   cbuffers = buffers2;
   testing::assert_equal(uint_t(buffers2.buffersId()), uint_t(cbuffers.buffersId()));

   for (auto& block : *blocks)
   {
      auto bView  = buffers.view(block);
      auto b2View = buffers2.view(block);

      // non-const to const copy init
      [[maybe_unused]] ConstBViewType cbView{ bView };

      testing::assert_equal((const double*) bView.data(), cbView.data());

      // non-const to const assignment
      cbView = b2View;

      testing::assert_equal((const double*) b2View.data(), cbView.data());
   }
}

} // namespace

int main(int argc, char** argv)
{
   walberla::mpi::Environment env{ argc, argv };

   return walberla::v8::testing::TestsRunner(
             {
                { "testBufferIndexing", &testBufferIndexing },
                { "testBufferIndexingWithPadding", &testBufferIndexingWithPadding },
                { "testConstCorrectness", &testConstCorrectness },
             })
      .run(argc, argv);
}
