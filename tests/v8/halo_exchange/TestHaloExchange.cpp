#include "core/Stdlib.hpp"
#include WALBERLA_STDLIB(span)

#include "blockforest/all.h"

#include "core/all.h"

#include "stencil/all.h"

#include <algorithm>
#include <atomic>
#include <ranges>
#include <span>

#include "walberla/v8/HaloExchange.hpp"
#include "walberla/v8/Memory.hpp"
#include "walberla/v8/StencilRanges.hpp"
#include "walberla/v8/Sweep.hpp"
#include "walberla/v8/Testing.hpp"

#if defined(WALBERLA_BUILD_WITH_GPU_SUPPORT)

namespace device
{
/**
 * Functions defined next door in `TestHaloExchangeDevice.cu`
 */
void testHaloExchangeInterfacesDevice();
void testDeviceFieldPackInfo();
void testDeviceStreamPullPackInfo();
} // namespace device
#endif

namespace
{

using namespace walberla;
using namespace walberla::v8;

class TestHostPackInfo
{
 public:
   static inline std::atomic< size_t > valuesSent     = 0;
   static inline std::atomic< size_t > valuesReceived = 0;
   static void resetCounts()
   {
      valuesSent     = 0;
      valuesReceived = 0;
   }

   using value_type = int;

   void pack(const IBlock& /*srcBlock*/, stencil::Direction commDir, stdlib::span< value_type > buffer,
             exectag::Serial) const
   {
      testing::assert_equal(buffer.size(), packetSize(commDir));

      buffer[0] = stencil::cx[commDir];
      buffer[1] = stencil::cy[commDir];
      buffer[2] = stencil::cz[commDir];

      TestHostPackInfo::valuesSent += buffer.size();
   };

   void unpack(IBlock& /*dstBlock*/, stencil::Direction commDir, stdlib::span< const value_type > buffer,
               exectag::Serial) const
   {
      testing::assert_equal(buffer.size(), packetSize(commDir));

      testing::assert_range_equal(
         buffer, std::array< int, 3 >{ stencil::cx[commDir], stencil::cy[commDir], stencil::cz[commDir] });

      TestHostPackInfo::valuesReceived += buffer.size();
   };

   size_t packetSize(stencil::Direction /*commDir*/) const { return 3; };
};

void testHaloExchangeInterfacesHost()
{
   static_assert(halo_exchange::packinfo_interface::IPackInfo< TestHostPackInfo, exectag::Serial >);

   mpi::MPIManager::instance()->useWorldComm();

   Vector3< size_t > numBlocks = math::getFactors3D(size_t(mpi::MPIManager::instance()->numProcesses()));
   auto blocks = blockforest::createUniformBlockGrid(numBlocks[0], numBlocks[1], numBlocks[2], 4, 4, 4, 1.,
                                                     /* oneBlockPerProcess */
                                                     true,
                                                     /* periodic */ true, true, true);
   auto he     = HaloExchange::create< stencil::D3Q6, memtag::stdmem >(blocks) //
                .sync(TestHostPackInfo())
                .build();

   he.startCommunication();
   he.wait();

   testing::assert_equal((size_t) TestHostPackInfo::valuesSent, size_t(6 * 3));
   testing::assert_equal((size_t) TestHostPackInfo::valuesReceived, size_t(6 * 3));
}

void testHostFieldPackInfo()
{
   using FieldType = memory::Field< double, 3, memtag::stdmem >;
   static_assert(halo_exchange::packinfo_interface::IPackInfo<                              //
                 halo_exchange::fields::GenericFieldPackInfo< FieldType, exectag::Serial >, //
                 exectag::Serial >);

   mpi::MPIManager::instance()->useWorldComm();

   Vector3< size_t > numBlocks = math::getFactors3D(size_t(mpi::MPIManager::instance()->numProcesses()));
   auto blocks = blockforest::createUniformBlockGrid(numBlocks[0], numBlocks[1], numBlocks[2], 4, 4, 4, 1.,
                                                     /* oneBlockPerProcess */
                                                     true,
                                                     /* periodic */ true, true, true);

   FieldType f{ *blocks };

   auto he = HaloExchange::create< stencil::D3Q6, memtag::stdmem >(blocks) //
                .sync(f)
                .build();

   for (auto& block : *blocks)
   {
      FieldView fView{ f, block };
      sweep::forAllCells(exectag::Serial{}, fView.slices().interior(), [&](Cell cell) {
         auto cc        = blocks->getBlockLocalCellCenter(block, cell);
         fView(cell, 0) = cc[0];
         fView(cell, 1) = cc[1];
         fView(cell, 2) = cc[2];
      });
   }

   he.communicate();

   for (auto& block : *blocks)
   {
      const FieldView fView{ f, block };

      for (auto dir : stencil_ranges::noCenter< stencil::D3Q6 >())
      {
         const CellInterval dstRegion{ fView.slices().ghostSlice(dir) };
         sweep::forAllCells(exectag::Serial{}, dstRegion, [&](Cell cell) {
            auto cc = blocks->getBlockLocalCellCenter(block, cell);
            blocks->mapToPeriodicDomain(cc);

            testing::assert_close(fView(cell, 0), cc[0]);
            testing::assert_close(fView(cell, 1), cc[1]);
            testing::assert_close(fView(cell, 2), cc[2]);
         });
      }
   }
}

void testStreamPullPackInfo()
{
   using Stencil   = stencil::D3Q27;
   using FieldType = memory::Field< double, Stencil::Q, memtag::stdmem >;
   static_assert(halo_exchange::packinfo_interface::IPackInfo<
                 halo_exchange::fields::GenericStreamPullPackInfo< FieldType, Stencil, exectag::Serial >, //
                 exectag::Serial >);

   mpi::MPIManager::instance()->useWorldComm();

   Vector3< size_t > numBlocks = math::getFactors3D(size_t(mpi::MPIManager::instance()->numProcesses()));
   auto blocks = blockforest::createUniformBlockGrid(numBlocks[0], numBlocks[1], numBlocks[2], 4, 4, 4, 1.,
                                                     /* oneBlockPerProcess */
                                                     true,
                                                     /* periodic */ true, true, true);

   FieldType f{ *blocks, 1, -1. };

   auto he =
      HaloExchange::create< Stencil, memtag::stdmem >(blocks).sync(halo_exchange::streamPullSync< Stencil >(f)).build();

   for (auto& block : *blocks)
   {
      FieldView fView{ f, block };
      sweep::forAllCells(exectag::Serial{}, fView.slices().interior(), [&](Cell cell) {
         auto cc = blocks->getBlockLocalCellCenter(block, cell);
         for (auto dir : stencil_ranges::all< Stencil >())
         {
            fView(cell, cell_idx_t(dir)) = cc.sqrLength() * double(dir);
         }
      });
   }

   he.communicate();

   for (auto& block : *blocks)
   {
      const FieldView fView{ f, block };

      for (auto dir : stencil_ranges::noCenter< stencil::D3Q6 >())
      {
         const CellInterval dstRegion{ fView.slices().ghostSlice(dir) };
         sweep::forAllCells(exectag::Serial{}, dstRegion, [&](Cell cell) {
            auto cc = blocks->getBlockLocalCellCenter(block, cell);
            blocks->mapToPeriodicDomain(cc);

            auto subdirs = stencil_ranges::subdirections< Stencil >(stencil::inverseDir[dir]);

            for (auto streamingDir : stencil_ranges::all< Stencil >())
            {
               if (std::ranges::find(subdirs, streamingDir) != std::ranges::end(subdirs))
               {
                  testing::assert_close(fView(cell, cell_idx_t(streamingDir)), cc.sqrLength() * double(streamingDir));
               }
               else
               {
                  testing::assert_close(fView(cell, cell_idx_t(streamingDir)), -1.);
               }
            }
         });
      }
   }
}

} // namespace

int main(int argc, char** argv)
{
   walberla::mpi::Environment env{ argc, argv };

   return walberla::v8::testing::TestsRunner( //
             {
                { "testHaloExchangeInterfacesHost", &testHaloExchangeInterfacesHost },
                { "testHostFieldPackInfo", &testHostFieldPackInfo },
                { "testStreamPullPackInfo", &testStreamPullPackInfo },
#if defined(WALBERLA_BUILD_WITH_GPU_SUPPORT)
                { "testHaloExchangeInterfacesDevice", &device::testHaloExchangeInterfacesDevice },
                { "testDeviceFieldPackInfo", &device::testDeviceFieldPackInfo },
                { "testDeviceStreamPullPackInfo", &device::testDeviceStreamPullPackInfo },
#endif
             })
      .run(argc, argv);
}
