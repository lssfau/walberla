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
//! \file PrintStacktraceTest.cpp
//! \ingroup core
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//
//======================================================================================================================

#include "core/Abort.h"
#include "core/debug/PrintStacktrace.h"
#include "core/mpi/MPIManager.h"

#include <iostream>


void function2( int b )
{
   walberla::debug::printStacktrace();
   std::cout << b << "\n\n";

   std::cerr << "\n\n";
   WALBERLA_ABORT("My Error message");
}

void function1( int a )
{
   function2( a );
}


int main( int argc, char**argv )
{
   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );

   function1(5);
   return 0;
}
