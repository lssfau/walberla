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
//! \file ExecutionTags.hpp
//! \author Behzad Safaei <behzad.safaei@fau.de>
//
//======================================================================================================================

#pragma once

#include <concepts>

#include "waLBerlaDefinitions.h"
#include "gpu/ErrorChecking.h"


namespace walberla::v8
{
namespace sweep
{
namespace exectag
{

/**
 * @defgroup v8core-sweep-exectags Execution Tags
 * @ingroup v8core-sweep
 * @brief Control kernel execution between CPU and accelerators
 * 
 * Execution tags are used to control the execution of kernels (see @ref v8core-sweep-dispatchers).
 */


struct _exec_tag {};

/**
 * @brief Serial execution on a single CPU thread
 * @ingroup v8core-sweep-exectags
 */
struct Serial : public _exec_tag {
    void sync() const {}
};


/**
 * @brief Thread-parallel execution using OpenMP
 * @ingroup v8core-sweep-exectags
 * 
 * Use this execution tag to dispatch a kernel using OpenMP, with a specified number of threads.
 */
struct OpenMP : public _exec_tag {
    OpenMP() = default;
    OpenMP(int t) : numThreads_(t) {}

    int numThreads() const { return numThreads_; }
    void sync() const {}

private:
    int numThreads_ = 0;
};

/**
 * @brief Accelerated execution on GPU
 * @ingroup v8core-sweep-exectags
 * 
 * Use this execution tag to dispatch a kernel to a GPU accelerator.
 * Allows setting the desired thread block size and GPU stream.
 */
struct GPU : public _exec_tag {
#if defined(WALBERLA_BUILD_WITH_CUDA) || defined(WALBERLA_BUILD_WITH_HIP)
    GPU() = default;
    GPU(gpuStream_t s) : stream_(s) {}
    GPU(dim3 b) : block_(b) {}
    GPU(dim3 b, gpuStream_t s): block_(b), stream_(s) {}

    const gpuStream_t& stream() const { return stream_; }
    dim3 block() const { return block_; }
    void sync() const { WALBERLA_GPU_CHECK(gpuStreamSynchronize(stream_)) }

private:
    dim3 block_{ 32, 4, 4 };
    gpuStream_t stream_{ gpuStreamDefault };
#else
    void sync() const {}
#endif
};

} // namespace exectag

/**
 * @brief Concept identifying valid execution tags.
 * @ingroup v8core-sweep-exectags
 */
template< typename T >
concept ExecTag = 
    std::derived_from< T, exectag::_exec_tag > && 
    requires(const T& t) {{ t.sync() } -> std::same_as<void>;};

} // namepsace sweep

namespace exectag = sweep::exectag;

} // namespace walberla::v8
