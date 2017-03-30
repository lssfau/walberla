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
//! \file ConfigFromPythonTest.cpp
//! \ingroup core
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#include "core/Abort.h"
#include "core/debug/TestSubsystem.h"
#include "core/mpi/MPIManager.h"
#include "core/mpi/Environment.h"
#include "python_coupling/CreateConfig.h"

using namespace walberla;


int main( int argc, char** argv )
{
   debug::enterTestMode();

   mpi::Environment env( argc, argv );

   auto config =  python_coupling::createConfigFromPythonScript( argv[1] );

   WALBERLA_CHECK_EQUAL( argc, 2 );

   config->listParameters();
   WALBERLA_CHECK_EQUAL      ( int        ( config->getParameter<int>        ("testInt")  ),    4          );
   WALBERLA_CHECK_EQUAL      ( std::string( config->getParameter<std::string>("testString")), "someString" );
   WALBERLA_CHECK_FLOAT_EQUAL( double     ( config->getParameter<real_t>     ("testDouble")), real_t(42.42));


   auto subBlock = config->getBlock("subBlock");
   WALBERLA_CHECK_EQUAL      ( std::string( subBlock.getParameter<std::string>("subKey1") ), std::string("abc") );
   WALBERLA_CHECK_EQUAL      ( std::string( subBlock.getParameter<std::string>("subKey2") ), std::string("def") );

   return 0;
}
