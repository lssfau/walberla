#include "core/Stdlib.hpp"
#include WALBERLA_STDLIB(span)

#include "blockforest/all.h"

#include "core/all.h"

#include "gpu/GPUWrapper.h"

#include "stencil/all.h"

#include <algorithm>
#include <atomic>
#include <ranges>
#include <span>

#include "walberla/v8/Domain.hpp"
#include "walberla/v8/HaloExchange.hpp"
#include "walberla/v8/Memory.hpp"
#include "walberla/v8/StencilRanges.hpp"
#include "walberla/v8/Sweep.hpp"
#include "walberla/v8/Testing.hpp"

namespace device
{

using namespace walberla;
using namespace walberla::v8;

namespace testpackinfo
{
static __global__ void doPack(stdlib::span< int > buffer, stencil::Direction commDir)
{
   constexpr stencil_ranges::Velocities cs;

   buffer[0] = cs.x(commDir);
   buffer[1] = cs.y(commDir);
   buffer[2] = cs.z(commDir);
}

static __global__ void doUnpack(stdlib::span< const int > buffer, stdlib::span< int > dst)
{
   dst[0] = buffer[0];
   dst[1] = buffer[1];
   dst[2] = buffer[2];
}
} // namespace testpackinfo

class TestDevicePackInfo
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
             exectag::GPU xtag) const
   {
      testing::assert_equal(buffer.size(), packetSize(commDir));

      // clang-format off
      testpackinfo::doPack<<< dim3(1), dim3(1), 0, xtag.stream() >>>(buffer, commDir);
      // clang-format on

      TestDevicePackInfo::valuesSent += buffer.size();
   };

   void unpack(IBlock& /*dstBlock*/, stencil::Direction commDir, stdlib::span< const value_type > buffer,
               exectag::GPU xtag) const
   {
      testing::assert_equal(buffer.size(), packetSize(commDir));

      std::vector< value_type, memory::UnifiedAllocator< value_type > > results(3, value_type{});
      stdlib::span< value_type > resultsView{ results.data(), results.size() };

      // clang-format off
      testpackinfo::doUnpack<<< dim3(1), dim3(1), 0, xtag.stream() >>>(buffer, resultsView);
      // clang-format on

      xtag.sync();

      testing::assert_range_equal(
         results, std::array< int, 3 >{ stencil::cx[commDir], stencil::cy[commDir], stencil::cz[commDir] });

      TestDevicePackInfo::valuesReceived += buffer.size();
   };

   size_t packetSize(stencil::Direction /*commDir*/) const { return 3; };
};

void testHaloExchangeInterfacesDevice()
{
   static_assert(halo_exchange::packinfo_interface::IPackInfo< TestDevicePackInfo, exectag::GPU >);

   mpi::MPIManager::instance()->useWorldComm();

   Vector3< size_t > numBlocks = math::getFactors3D(size_t(mpi::MPIManager::instance()->numProcesses()));
   auto blocks = blockforest::createUniformBlockGrid(numBlocks[0], numBlocks[1], numBlocks[2], 4, 4, 4, 1.,
                                                     /* oneBlockPerProcess */
                                                     true,
                                                     /* periodic */ true, true, true);
   auto he     = HaloExchange::create< stencil::D3Q6, memtag::device >(blocks) //
                .sync(TestDevicePackInfo())
                .build();

   he.startCommunication();
   he.wait();

   testing::assert_equal((size_t) TestDevicePackInfo::valuesSent, size_t(6 * 3));
   testing::assert_equal((size_t) TestDevicePackInfo::valuesReceived, size_t(6 * 3));
}

void testDeviceFieldPackInfo()
{
   using Stencil   = stencil::D3Q27;
   using FieldType = memory::Field< double, Stencil::Q, memtag::unified >;
   static_assert(halo_exchange::packinfo_interface::IPackInfo<
                 halo_exchange::fields::GenericStreamPullPackInfo< FieldType, Stencil, exectag::GPU >, //
                 exectag::GPU >);

   mpi::MPIManager::instance()->useWorldComm();

   Vector3< size_t > numBlocks = math::getFactors3D(size_t(mpi::MPIManager::instance()->numProcesses()));
   auto blocks = blockforest::createUniformBlockGrid(numBlocks[0], numBlocks[1], numBlocks[2], 4, 4, 4, 1.,
                                                     /* oneBlockPerProcess */
                                                     true,
                                                     /* periodic */ true, true, true);

   FieldType f{ *blocks };

   auto he = HaloExchange::create< stencil::D3Q6, memtag::device >(blocks) //
                .sync(f)
                .build();

   for (auto& block : *blocks)
   {
      FieldView fView{ f, block };
      domain::BlockGeometry bGeom{ *blocks, block };

      sweep::forAllCells(exectag::GPU{}, fView.slices().interior(), [=] WALBERLA_HOST_DEVICE(Cell cell) {
         auto cc        = bGeom.localGrid().cellCenter(cell);
         fView(cell, 0) = cc[0];
         fView(cell, 1) = cc[1];
         fView(cell, 2) = cc[2];
      });
   }

   sweep::sync();

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

void testDeviceStreamPullPackInfo()
{
   using Stencil   = stencil::D3Q27;
   using FieldType = memory::Field< double, Stencil::Q, memtag::unified >;
   static_assert(halo_exchange::packinfo_interface::IPackInfo<
                 halo_exchange::fields::GenericStreamPullPackInfo< FieldType, Stencil, exectag::GPU >, //
                 exectag::GPU >);

   mpi::MPIManager::instance()->useWorldComm();

   Vector3< size_t > numBlocks = math::getFactors3D(size_t(mpi::MPIManager::instance()->numProcesses()));
   auto blocks = blockforest::createUniformBlockGrid(numBlocks[0], numBlocks[1], numBlocks[2], 4, 4, 4, 1.,
                                                     /* oneBlockPerProcess */
                                                     true,
                                                     /* periodic */ true, true, true);

   FieldType f{ *blocks, 1, -1. };

   auto he =
      HaloExchange::create< Stencil, memtag::device >(blocks).sync(halo_exchange::streamPullSync< Stencil >(f)).build();

   for (auto& block : *blocks)
   {
      FieldView fView{ f, block };
      domain::BlockGeometry bGeom{ *blocks, block };

      sweep::forAllCells(exectag::GPU{}, fView.slices().interior(), [=] WALBERLA_HOST_DEVICE(Cell cell) {
         auto cc = bGeom.localGrid().cellCenter(cell);
         for (auto dir : stencil_ranges::all< Stencil >())
         {
            fView(cell, cell_idx_t(dir)) = cc.sqrLength() * double(dir);
         }
      });
   }

   sweep::sync();

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

} // namespace device