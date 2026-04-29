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
//! \file Sweeper.hpp
//! \author Frederik Hennig <frederik.hennig@fau.de>
//! \author Behzad Safaei <behzad.safaei@fau.de>
//
//======================================================================================================================

#pragma once

#include "blockforest/StructuredBlockForest.h"

#include "core/OpenMP.h"
#include "core/cell/all.h"

#include <memory>

#include "./ExecutionTags.hpp"

namespace walberla::v8::sweep
{

/**
 * @defgroup v8core-sweep-dispatchers Cellwise Kernel Dispatchers
 * @ingroup v8core-sweep
 * @brief Run an operation for every cell in the given interval.
 *
 * The functions in this module iterate the given cell interval `ci` in x-fastest, z-slowest order
 * and call `func` for every cell. The iteration is dispatched
 * according to the given `ExecutionTag`:
 *
 *  - exectag::Serial : Sequential iteration on the host.
 *
 *  - exectag::OpenMP : Parallel iteration using OpenMP threads.
 *
 *  - exectag::GPU : CUDA/HIP kernel exectuion with one GPU thread per cell in the given 'ci'.
 *                   Launches the kernel asynchronously on the GPU stream specified by `exTag` (defaults to the default
 * stream). Synchronization must be performed by the caller.
 */

/**
 * @brief Cellwise kernel function
 * @ingroup
 *
 * Kernel function taking a `Cell` as its single argument, as required by the `forAllCells` kernel dispatchers.
 */
template< typename Func >
concept CellFunction = requires(Func func, Cell cell) { func(cell); };

template< typename Func >
concept CellPairFunction = requires(Func func, Cell cell) { func(cell, cell); };

// ------------ Serial ------------

/**
 * @brief Single-thread CPU cellwise kernel dispatcher
 * @ingroup v8core-sweep-dispatchers
 *
 * Execute a cellwise kernel function on all cells of a given interval in sequence on the CPU.
 */
template< CellFunction Func >
void forAllCells(const exectag::Serial& /* exTag */, const CellInterval& ci, Func&& func)
{
   for (cell_idx_t z = ci.zMin(); z <= ci.zMax(); ++z)
      for (cell_idx_t y = ci.yMin(); y <= ci.yMax(); ++y)
         for (cell_idx_t x = ci.xMin(); x <= ci.xMax(); ++x)
            func({ x, y, z });
}

// ------------ OpenMP ------------

/**
 * @brief Multi-threaded CPU cellwise kernel dispatcher using OpenMP
 * @ingroup v8core-sweep-dispatchers
 *
 * Execute a cellwise kernel function on all cells of a given interval in parallel on the CPU using OpenMP.
 *
 * @warning If OpenMP is not enabled in waLBerla's build configuration, falls back to sequential iteration.
 */
#if defined(_OPENMP)
template< CellFunction Func >
void forAllCells(const exectag::OpenMP& exTag, const CellInterval& ci, Func&& func)
{
#   pragma omp parallel for collapse(3) schedule(static) \
      num_threads(exTag.numThreads() > 0 ? exTag.numThreads() : omp_get_max_threads())
   for (cell_idx_t z = ci.zMin(); z <= ci.zMax(); ++z)
      for (cell_idx_t y = ci.yMin(); y <= ci.yMax(); ++y)
         for (cell_idx_t x = ci.xMin(); x <= ci.xMax(); ++x)
            func({ x, y, z });
}
#else
template< CellFunction Func >
void forAllCells(const exectag::OpenMP& /*exTag*/, const CellInterval& ci, Func&& func)
{
   forAllCells(exectag::Serial{}, ci, std::forward< Func >(func));
}
#endif

// ------------ GPU ------------
#if defined(__CUDACC__) || defined(__HIPCC__)
namespace detail
{
template< CellFunction Func >
__global__ void sweepKernel(CellInterval ci, Func func)
{
   const cell_idx_t x = blockIdx.x * blockDim.x + threadIdx.x;
   const cell_idx_t y = blockIdx.y * blockDim.y + threadIdx.y;
   const cell_idx_t z = blockIdx.z * blockDim.z + threadIdx.z;

   if (x < ci.xSize() && y < ci.ySize() && z < ci.zSize()) func({ ci.xMin() + x, ci.yMin() + y, ci.zMin() + z });
}
} // namespace detail
#endif

#if defined(__CUDACC__) || defined(__HIPCC__)
template< CellFunction Func >
void forAllCells(const exectag::GPU& exTag, const CellInterval& ci, Func&& func)
{
   dim3 block{ std::min< unsigned int >(exTag.block().x, ci.xSize()),
               std::min< unsigned int >(exTag.block().y, ci.ySize()),
               std::min< unsigned int >(exTag.block().z, ci.zSize()) };

   dim3 grid{ ((unsigned int) ci.xSize() + block.x - 1) / block.x, ((unsigned int) ci.ySize() + block.y - 1) / block.y,
              ((unsigned int) ci.zSize() + block.z - 1) / block.z };

   // clang-format off
   detail::sweepKernel<<< grid, block, 0, exTag.stream() >>>(ci, std::forward< Func >(func));
   // clang-format on
}
#else

/**
 * @brief GPU-accelerated kernel dispatcher
 * @ingroup v8core-sweep-dispatchers
 *
 * Execute a cellwise kernel function on all cells of a given interval in parallel on the GPU
 *
 * @warning Only available if CUDA or HIP are enabled
 */
template< CellFunction Func >
void forAllCells(const exectag::GPU&, const CellInterval&, Func&&)
{
   static_assert(never_true< Func >::value, "forAllCells(exectag::GPU, ...) requires GPU compiler.");
}

#endif

inline void sync() { WALBERLA_GPU_CHECK(gpuDeviceSynchronize()) }

/**
 * @brief Run an operation for every pair of cells in the given intervals
 *
 * Iterates the given cell intervals sequentially and in lockstep in z-y-x order
 * (`x` is the fastest, `z` the slowest coordinate)
 * and calls `func` for every pair of cells.
 *
 * @throws std::invalid_argument In debug mode, throws if the given cell intervals are of different size.
 */
template< ExecTag XTag, CellPairFunction Func >
inline void forAllCellPairs(XTag xtag, const CellInterval& firstCi, const CellInterval& secondCi, Func func)
{
   WALBERLA_DEBUG_SECTION()
   {
      if (firstCi.xSize() != secondCi.xSize() || firstCi.ySize() != secondCi.ySize() ||
          firstCi.zSize() != secondCi.zSize())
      {
         throw std::invalid_argument{ "Cell intervals must have the same size" };
      }
   }

   const CellInterval referenceInterval{ { cell_idx_t(0), cell_idx_t(0), cell_idx_t(0) },
                                         {
                                            firstCi.xSize() - cell_idx_t(1),
                                            firstCi.ySize() - cell_idx_t(1),
                                            firstCi.zSize() - cell_idx_t(1),
                                         } };

   forAllCells(xtag, referenceInterval, [firstCi, secondCi, func] WALBERLA_HOST_DEVICE(Cell offsets) {
      const Cell firstCell{ firstCi.min() + offsets };
      const Cell secondCell{ secondCi.min() + offsets };
      func(firstCell, secondCell);
   });
}

class SerialSweeper
{
 private:
   std::shared_ptr< StructuredBlockForest > blocks_;

 public:
   SerialSweeper(const std::shared_ptr< StructuredBlockForest >& blocks) : blocks_{ blocks } {}

   void sweep(std::function< void(IBlock&) > func) const
   {
      for (auto& block : *blocks_)
         func(block);
   }

   void sweep(std::function< void(IBlock*) > func) const
   {
      for (auto& block : *blocks_)
         func(&block);
   }

   void forAllBlocks(std::function< void(IBlock&) > func) const { sweep(func); }

   void forAllBlocks(std::function< void(IBlock*) > func) const { sweep(func); }

   void forAllCells(std::function< void(Cell) > func) const
   {
      for (cell_idx_t z = 0; z < cell_idx_c(blocks_->getNumberOfZCellsPerBlock()); ++z)
         for (cell_idx_t y = 0; y < cell_idx_c(blocks_->getNumberOfYCellsPerBlock()); ++y)
            for (cell_idx_t x = 0; x < cell_idx_c(blocks_->getNumberOfXCellsPerBlock()); ++x)
               func({ x, y, z });
   }
};

} // namespace walberla::v8::sweep
