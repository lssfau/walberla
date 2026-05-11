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
//! \author Frederik Hennig
//
//======================================================================================================================

#include "core/all.h"

#include "walberla/V8.hpp"
#include "walberla/v8/Testing.hpp"

// Scenarios
#include "FreeSlipPipe.hpp"
#include "FullyPeriodic.hpp"
#include "MirroredHalfChannel.hpp"

using namespace walberla;
using namespace walberla::v8;

// clang-format off
#define TEST_LANGUAGE_@_testLanguage@
// clang-format on

#if defined(TEST_LANGUAGE_CUDA) || defined(TEST_LANGUAGE_HIP)
using MemoryTag = memtag::unified;
#else
using MemoryTag = memtag::host;
#endif

using namespace lbm_scenarios;

int main(int argc, char** argv)
{
   mpi::Environment env{ argc, argv };

   return walberla::v8::testing::TestsRunner( //
             {
                { "FullyPeriodic", &FullyPeriodic< MemoryTag >::run },
                { "MirroredHalfChannel", &MirroredHalfChannel< MemoryTag >::run },
                { "FreeSlipPipe", &FreeSlipPipe< MemoryTag >::run },
             })
      .run(argc, argv);
}
