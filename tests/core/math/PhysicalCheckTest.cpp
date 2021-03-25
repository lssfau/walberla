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
//! \file PhysicalCheckTest.cpp
//! \ingroup core
//! \author David Staubach <david.staubach@fau.de>
//
//======================================================================================================================

#include "core/Abort.h"
#include "core/DataTypes.h"
#include "core/Environment.h"
#include "core/debug/TestSubsystem.h"
#include "core/logging/Logging.h"
#include "core/math/PhysicalCheck.h"

#include <iostream>
#include <string>
#include <vector>


using namespace walberla;
using namespace walberla::math;


int testPhysicalCheck1( shared_ptr<Config> & config )
{
   Config::BlockHandle pcConfigBlock = config->getBlock( "Physical_Check" );

   if( !pcConfigBlock )
      WALBERLA_ABORT( "You have to specify a \"Physical_Check\" block in the configuration file!" );

   PhysicalCheck pc( pcConfigBlock );

   WALBERLA_CHECK( pc.isDefined(std::string("parameter1")) );
   WALBERLA_CHECK( pc.isDefined(std::string("parameter2")) );
   WALBERLA_CHECK( pc.isDefined(std::string("var1")) );
   WALBERLA_CHECK( pc.isDefined(std::string("var2")) );
   // TODO: further checks

   // Check for functionality is within the function
   pc.completeConfig( config );

   return 0;
}

int testPhysicalCheck2( shared_ptr<Config> & config )
{
   Config::BlockHandle pcConfigBlock = config->getBlock( "Physical_Check" );

   if( !pcConfigBlock )
      WALBERLA_ABORT( "You have to specify a \"Physical_Check\" block in the configuration file!" );

   PhysicalCheck pc;

   pc.addBlock(pcConfigBlock);

   WALBERLA_CHECK( pc.isDefined(std::string("parameter1")) );
   WALBERLA_CHECK( pc.isDefined(std::string("parameter2")) );
   WALBERLA_CHECK( pc.isDefined(std::string("var1")) );
   WALBERLA_CHECK( pc.isDefined(std::string("var2")) );
   // TODO: further checks

   // Check for functionality is within the function
   pc.completeConfig( config );

   return 0;
}

int main( int argc, char** argv )
{
   walberla::Environment env( argc, argv );

   debug::enterTestMode();

   try {
      shared_ptr<Config> config = env.config();

      int value;
      value = testPhysicalCheck1( config );

      value = testPhysicalCheck2( config );

      return value;
   }
   catch( std::exception & e )
   {
      WALBERLA_LOG_INFO( "Unhandled exception raised: " << e.what() );
   }
}
