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
//! \file AlignmentTest.h
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#include "cuda/AlignedAllocation.h"
#include "core/mpi/Environment.h"
#include "core/debug/TestSubsystem.h"
#include "core/logging/Logging.h"


using namespace walberla;
using namespace cuda;


int main( int argc, char ** argv )
{
   debug::enterTestMode();
   mpi::Environment env( argc, argv );

   size_t pitch = 0;
   size_t width = 7;
   size_t height = 20;
   size_t alignment = 512;
   size_t offset = 16;
   void *ptr = allocate_pitched_with_offset( pitch, width, height, alignment, offset );
   WALBERLA_LOG_INFO("Pitch " << pitch);

   char * cptr = reinterpret_cast<char*>( ptr );
   WALBERLA_CHECK_EQUAL( size_t(cptr + offset) % alignment, 0 );

   free_aligned_with_offset( ptr );

   return 0;
}
