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
//! \file ConfigTest.cpp
//! \ingroup core
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//
//======================================================================================================================

#include "core/Abort.h"
#include "core/config/Config.h"
#include "core/debug/TestSubsystem.h"
#include "core/mpi/MPIManager.h"


using namespace walberla;



int main( int argc, char** argv )
{
	debug::enterTestMode();

	MPIManager::instance()->initializeMPI( &argc, &argv );

	WALBERLA_CHECK_EQUAL( argc, 2 );

	Config config;

	config.readParameterFile( argv[1] );

	if(!config)
	   WALBERLA_ABORT(config.error());

	// Test case insensitivity

	Config::BlockHandle block  = config.getOneBlock( "Block" );
	Config::BlockHandle block1 = config.getOneBlock( "BLOCK" );
	Config::BlockHandle block2 = config.getOneBlock( "block" );

	WALBERLA_CHECK_EQUAL( block.getKey(),  "Block" );
	WALBERLA_CHECK_EQUAL( block1.getKey(), "Block" );
	WALBERLA_CHECK_EQUAL( block2.getKey(), "Block" );

	std::string value  = block.getParameter<std::string>("Key");
	std::string value0 = block.getParameter<std::string>("KEY");
	std::string value1 = block.getParameter<std::string>("key");

	WALBERLA_CHECK_EQUAL( value,  "Value" );
	WALBERLA_CHECK_EQUAL( value0, "Value" );
	WALBERLA_CHECK_EQUAL( value1, "Value" );

	return 0;
}
