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
//! \file TestGPUSweeper.cu
//! \author Behzad Safaei <behzad.safaei@fau.de>
//
//======================================================================================================================

#include "blockforest/all.h"
#include "core/all.h"
#include "gpu/GPURAII.h"

#include "walberla/v8/Memory.hpp"
#include "walberla/v8/Sweep.hpp"
#include "walberla/v8/Testing.hpp"



namespace test_gpu_sweeper
{
using namespace walberla;
using namespace walberla::v8;

template< IFieldView FView >
struct IncrFunctor {
   IncrFunctor(FView fv, FView::element_type val) : 
      fv_(fv), val_(val) {}

   WALBERLA_HOST_DEVICE inline void operator()(Cell c){
      fv_(c, 0) += val_;
   }

   FView fv_;
   FView::element_type val_;
};

void testForAllCells(){
   auto blocks = blockforest::createUniformBlockGrid(1, 1, 1, 2, 3, 4, 1.,
                                                     /* oneBlockPerProcess */ false);
   
   const double initVal = 2.2;
   const double incrVal = 42.42;
   const Field< double, 1, memtag::unified > f{ *blocks, 1, initVal, field::fzyx };
   
   for (auto& block : *blocks){
      FieldView fv{ f, block };
      IncrFunctor myFtor{ fv, incrVal };
      // Sweep on device
      sweep::forAllCells(exectag::GPU{}, fv.slices().interior(), myFtor);
      sweep::sync();
   }

   for (auto& block : *blocks){
      ConstFieldView fv{ f, block };
      // Sweep on host
      sweep::forAllCells(exectag::Serial{}, fv.slices().interior(), [&](Cell c){
         testing::assert_close(fv(c, 0), initVal + incrVal);
      });
   }
}

void testForAllCellsStreams() {
   auto blocks = blockforest::createUniformBlockGrid(1, 1, 1, 8, 8, 8, 1.,
                                                      /* oneBlockPerProcess */ false);

   const double initVal1 = 1.1;
   const double initVal2 = 2.2;
   const double incrVal1 = 10.1;
   const double incrVal2 = 20.2;

   Field<double, 1, memtag::unified> field1{ *blocks, 1, initVal1, field::fzyx };
   Field<double, 1, memtag::unified> field2{ *blocks, 1, initVal2, field::fzyx };
   CellInterval ci{ field1.slices().interior() };

   // Create two GPU execution spaces with a separate stream for each
   auto stream1 = gpu::StreamRAII::newStream();
   auto stream2 = gpu::StreamRAII::newStream();
   exectag::GPU gex1{ stream1 };
   exectag::GPU gex2{ stream2 };

   // Launch first sweep on field1
   for (auto& block : *blocks){
      FieldView fv1{ field1, block };
      IncrFunctor ftor1{ fv1, incrVal1 };
      sweep::forAllCells(gex1, ci, ftor1);
   }

   // Launch second sweep on field2 concurrently
   for (auto& block : *blocks){
      FieldView fv2{ field2, block };
      IncrFunctor ftor2{ fv2, incrVal2 };
      sweep::forAllCells(gex2, ci, ftor2);
   }

   // Synchronize both streams
   gex1.sync();
   gex2.sync();

   // Verify 
   for (auto& block : *blocks){
      ConstFieldView cfv1{ field1, block };
      ConstFieldView cfv2{ field2, block };
      // Sweep on host
      sweep::forAllCells(exectag::Serial{}, ci, [&](Cell c){
         testing::assert_close(cfv1(c, 0), initVal1 + incrVal1);
         testing::assert_close(cfv2(c, 0), initVal2 + incrVal2);
      });
   }

}

template< IFieldView FView >
struct SFunctor {
   SFunctor(FView fv, FView::element_type val) : 
      fv_(fv), val_(val) {}

   WALBERLA_HOST_DEVICE inline void operator()(Cell c){
      for (int i = 0; i < 100000; ++i) {
            fv_(c, 0) += val_ / 1.0001 * (11111 % int(fv_(c, 0))); 
            fv_(c, 0) -= val_ / 1.0002 * (22222 % int(fv_(c, 0))); 
            fv_(c, 0) *= val_ / 1.0003 * (33333 % int(fv_(c, 0))); 
      }
   }

   FView fv_;
   FView::element_type val_;
};

void testForAllCellsStreamsConcurrency() {
   auto blocks = blockforest::createUniformBlockGrid(
      1, 1, 1, 
      32, 16, 8,  // i.e. 4096 threads per stream
      1., false
   );

   const double initVal1 = 1.1;
   const double initVal2 = 2.2;
   const double incrVal1 = 10.1;
   const double incrVal2 = 20.2;

   Field<double, 1, memtag::unified> field1{ *blocks, 1, initVal1, field::fzyx };
   Field<double, 1, memtag::unified> field2{ *blocks, 1, initVal2, field::fzyx };
   CellInterval ci{ field1.slices().interior() };

   // Two GPU streams
   auto stream1 = gpu::StreamRAII::newStream();
   auto stream2 = gpu::StreamRAII::newStream();
   exectag::GPU gex1{ stream1 };
   exectag::GPU gex2{ stream2 };

   // Create GPU events for time measurements
   gpu::EventRAII start1, end1, start2, end2;
   gpu::EventRAII startTotalCunc, endTotalCunc;
   gpu::EventRAII startTotalSeq, endTotalSeq;

   // Concurrent mode --------------------------------------------
   // Measure total time on the default stream
   WALBERLA_GPU_CHECK(gpuEventRecord(startTotalCunc, 0));
   {  
      // Launch first kernel on field1 and measure the time spent by stream 1
      WALBERLA_GPU_CHECK(gpuEventRecord(start1, gex1.stream()));
      FieldView fv1{ field1, *(blocks->begin()) };
      SFunctor ftor1{ fv1, incrVal1 };
      sweep::forAllCells(gex1, ci, ftor1);
      WALBERLA_GPU_CHECK(gpuEventRecord(end1, gex1.stream()));

      // Launch second kernel on field2 concurrently and measure the time spent by stream 2 
      WALBERLA_GPU_CHECK(gpuEventRecord(start2, gex2.stream()));
      FieldView fv2{ field2, *(blocks->begin()) };
      SFunctor ftor2{ fv2, incrVal2 };
      sweep::forAllCells(gex2, ci, ftor2);
      WALBERLA_GPU_CHECK(gpuEventRecord(end2, gex2.stream()));

      // Synchronize both streams and end their timers
      WALBERLA_GPU_CHECK(gpuEventSynchronize(end1));
      WALBERLA_GPU_CHECK(gpuEventSynchronize(end2));
   }
   WALBERLA_GPU_CHECK(gpuEventRecord(endTotalCunc, 0));
   WALBERLA_GPU_CHECK(gpuEventSynchronize(endTotalCunc));

   // Sequential mode --------------------------------------------
   // Repeat the same experiment but without concurrency
   WALBERLA_GPU_CHECK(gpuEventRecord(startTotalSeq, 0));
   {
      // Launch first stream and synchronize
      FieldView fv1{ field1, *(blocks->begin()) };
      SFunctor ftor1{ fv1, incrVal1 };
      sweep::forAllCells(gex1, ci, ftor1);
      gex1.sync();

      // Launch second kernel and synchronzie
      FieldView fv2{ field2, *(blocks->begin()) };
      SFunctor ftor2{ fv2, incrVal2 };
      sweep::forAllCells(gex2, ci, ftor2);
      gex2.sync();
   }
   WALBERLA_GPU_CHECK(gpuEventRecord(endTotalSeq, 0));
   WALBERLA_GPU_CHECK(gpuEventSynchronize(endTotalSeq));

   // Measure times
   float ms1 = 0.0, ms2 = 0.0, totalConc = 0.0,totalSeq = 0.0;
   WALBERLA_GPU_CHECK(gpuEventElapsedTime(&ms1, start1, end1));
   WALBERLA_GPU_CHECK(gpuEventElapsedTime(&ms2, start2, end2));
   WALBERLA_GPU_CHECK(gpuEventElapsedTime(&totalConc, startTotalCunc, endTotalCunc));
   WALBERLA_GPU_CHECK(gpuEventElapsedTime(&totalSeq, startTotalSeq, endTotalSeq));

   std::cout << "Kernel1 time: " << ms1 << " ms" << std::endl;
   std::cout << "Kernel2 time: " << ms2 << " ms" << std::endl;
   std::cout << "Total time (concurrent mode): "   << totalConc  << " ms" << std::endl;  // Should be less than (ms1 + ms2)
   std::cout << "Total time (sequential mode): "   << totalSeq   << " ms" << std::endl;  // Should roughly be twice totalConc

   testing::assert_less(totalConc, ms1 + ms2);
}



} // namespace test_gpu_sweeper

int main(int argc, char** argv){
   walberla::mpi::Environment env{ argc, argv };

   using namespace test_gpu_sweeper;

   return walberla::v8::testing::TestsRunner(
             {
                { "testForAllCells",                        &testForAllCells},
                { "testForAllCellsStreams",                 &testForAllCellsStreams},
                { "testForAllCellsStreamsConcurrency",      &testForAllCellsStreamsConcurrency},
             })
      .run(argc, argv);
}
