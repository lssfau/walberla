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
//! \file RandomUUID.cpp
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include <core/debug/TestSubsystem.h>
#include <core/logging/Logging.h>
#include <core/mpi/BufferDataTypeExtensions.h>
#include <core/mpi/Environment.h>
#include <core/mpi/RecvBuffer.h>
#include <core/mpi/SendBuffer.h>
#include <core/RandomUUID.h>

namespace walberla{

int main( int /*argc*/, char** /*argv*/ )
{
   debug::enterTestMode();

   RandomUUID uuid1;
   RandomUUID uuid2;

   WALBERLA_CHECK_EQUAL(uuid1, uuid1);
   WALBERLA_CHECK_EQUAL(uuid2, uuid2);

   WALBERLA_CHECK_UNEQUAL(uuid1, uuid2);
   WALBERLA_CHECK_UNEQUAL(uuid2, uuid1);

   WALBERLA_LOG_DEVEL(uuid1);

   WALBERLA_MPI_SECTION()
   {
      //pack&unpack test
      mpi::SendBuffer sb;
      sb << uuid1 <<uuid2;
      mpi::RecvBuffer rb(sb);
      RandomUUID uuid3, uuid4;
      rb >> uuid3 >> uuid4;
      WALBERLA_CHECK_EQUAL(uuid1, uuid3);
      WALBERLA_CHECK_EQUAL(uuid2, uuid4);
   }

   return EXIT_SUCCESS;
}

} //namespace walberla

int main( int argc, char** argv )
{
   walberla::mpi::Environment env(argc, argv);
   return walberla::main(argc, argv);
}
