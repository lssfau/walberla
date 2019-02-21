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
//! \file OpenMPWrapperTest.cpp
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include <core/debug/TestSubsystem.h>
#include <core/logging/Logging.h>
#include <core/OpenMP.h>
#include <core/mpi/Environment.h>

using namespace walberla;

int main( int argc, char** argv )
{
   mpi::Environment env(argc, argv);
   debug::enterTestMode();
   
   WALBERLA_OPENMP_SECTION()
   {
      WALBERLA_LOG_DEVEL_VAR(omp_get_num_threads  ());
      WALBERLA_LOG_DEVEL_VAR(omp_get_dynamic      ());
      WALBERLA_LOG_DEVEL_VAR(omp_get_nested       ());
      WALBERLA_LOG_DEVEL_VAR(omp_get_max_threads  ());
      WALBERLA_LOG_DEVEL_VAR(omp_get_thread_num   ());
      WALBERLA_LOG_DEVEL_VAR(omp_get_num_procs    ());
   }
}
