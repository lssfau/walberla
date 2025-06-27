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
//! \file BufferIOTest.cpp
//! \ingroup core
//! \author Markus Holzer <markus.holzer@fau.de>
//! \brief Tests for BufferIO: test writting and reading of Buffers
//
//======================================================================================================================

#include "core/Abort.h"
#include "core/debug/TestSubsystem.h"
#include "core/logging/Logging.h"
#include "core/math/Random.h"
#include "core/mpi/BufferIO.h"
#include "core/mpi/RecvBuffer.h"
#include "core/mpi/SendBuffer.h"
#include "core/mpi/Environment.h"

#include <vector>

using namespace walberla;

/**
 * Writes random values from a SendBuffer to a file and reads the values back in to a RecvBuffer.
 * The test checks if the values match
 */
void bufferIO(const bool forceSerialIO)
{
   constexpr uint_t s = 10;

   std::vector<double> values;
   values.reserve(s);
   for(uint_t i = 0; i < s; ++i){
      values.emplace_back( math::realRandom(real_c(0.0), real_c(1.0)) );
   }

   mpi::SendBuffer sb;
   mpi::RecvBuffer rb;

   for (const auto v : values)
   {
      sb << v;
   }

   const std::string extension = forceSerialIO ? "serialIO" : "mpiIO";
   const std::string filename = "bufferIO" + extension + ".data";
   walberla::mpi::writeBuffer(filename, sb, forceSerialIO);
   walberla::mpi::readBuffer(filename, rb, forceSerialIO);

   // Check that the same size of data was read
   WALBERLA_CHECK_EQUAL(sb.size(), rb.size());

   // Check that data is identical
   for(uint_t i = 0; i < s; ++i){
      double read;
      rb >> read;

      WALBERLA_CHECK_FLOAT_EQUAL(values[i], read);
   }

   // Check that all values have been checked
   WALBERLA_CHECK_EQUAL(rb.size(), 0);
   
}



int main(int argc, char**argv)
{
   debug::enterTestMode();
   mpi::Environment mpiEnv( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();
   
   // enforce serial I/O
   bufferIO(true);
   // use MPIIO if mpi version is compatible
   bufferIO(false);

   return EXIT_SUCCESS;
}
