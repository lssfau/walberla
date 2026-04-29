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
//! \file Stdlib.hpp
//! \ingroup core
//! \author Behzad Safaei <behzad.safaei@fau.de>
//! \author Frederik Hennig <frederik.hennig@fau.de>
//
//======================================================================================================================

#pragma once

#include "waLBerlaDefinitions.h"

#if defined(WALBERLA_ENABLE_V8CORE)

/**
 * Host/Device function qualifiers
 */

#   if defined(__CUDACC__) || defined(__HIPCC__)
// If this file is being compiled by nvcc or hipcc
#      define WALBERLA_HOST_DEVICE __host__ __device__
#      define WALBERLA_HOST __host__
#      define WALBERLA_DEVICE __device__
#   else
#      define WALBERLA_HOST_DEVICE
#      define WALBERLA_HOST
#      define WALBERLA_DEVICE
#   endif

/**
 * Portable standard library include proxies
 */

#   if defined(WALBERLA_BUILD_WITH_CUDA)
#      define WALBERLA_STDLIB(HEADER) <cuda/std/HEADER>
#   elif defined(WALBERLA_BUILD_WITH_HIP)
#      define WALBERLA_STDLIB(HEADER) <hip/std/HEADER>
#   else
#      define WALBERLA_STDLIB(HEADER) <HEADER>
#   endif

#   include WALBERLA_STDLIB(cstddef) // Minimal header to declare cuda:: or hip::

/**
 * Portable standard library aliases
 */

namespace walberla
{

#   if defined(WALBERLA_BUILD_WITH_CUDA)
namespace stdlib = ::cuda::std;
#   elif defined(WALBERLA_BUILD_WITH_HIP)
namespace stdlib = ::hip::std;
#   else
namespace stdlib = ::std;
#   endif

} // namespace walberla

#else

/**
 * If v8core is not enabled, don't do anything special
 */

#   define WALBERLA_STDLIB(HEADER) <HEADER>
#   define WALBERLA_HOST_DEVICE
#   define WALBERLA_HOST
#   define WALBERLA_DEVICE

#   include <cstddef>

namespace walberla
{
namespace stdlib = ::std;
}
#endif
