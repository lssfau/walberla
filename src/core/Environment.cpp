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
//! \file Environment.cpp
//! \ingroup core
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#include "Environment.h"

#include "core/config/Create.h"
#include "core/logging/Initialization.h"
#include "core/logging/Logging.h"
#include "core/logging/Tracing.h"
#include "core/mpi/MPIManager.h"
#include "core/StringUtility.h"
#include "core/uid/GlobalState.h"
#include "core/uid/SUID.h"

#include <algorithm>
#include <functional>
#include <sstream>
#include <string>


namespace walberla
{



Environment::Environment( int & argc, char ** & argv, bool mpiAbortOnException )
     : mpiEnv_( argc, argv, mpiAbortOnException )
{
  // Read configuration file
  if ( argc > 1 )
     config_ = config::create( argc, argv );

  WALBERLA_ROOT_SECTION() { std::cout << logging::Logging::getHeaderFooter( true ) << std::endl; }
  WALBERLA_MPI_WORLD_BARRIER();

  // Use configuration file to setup logging
  logging::configureLogging( config_ );

  // Use configuration file to setup the global state
  configureGlobalState( config_ );
}



Environment::~Environment()
{
   WALBERLA_MPI_WORLD_BARRIER();
   WALBERLA_ROOT_SECTION() { std::cout << logging::Logging::getHeaderFooter( false ) << std::endl; }
}



//======================================================================================================================
//
//  Global State Setup
//
//======================================================================================================================

void configureGlobalState( const shared_ptr<Config> & config ) {

   if( !!config && config->isDefined( "GlobalState" ) ) {

      std::string suids = config->getParameter< std::string >( "GlobalState" );

      std::vector< std::string > states = string_split( suids, ", \t" );
      states.erase( std::remove_if( states.begin(), states.end(), std::bind( &std::string::empty, std::placeholders::_1 ) ), states.end() );

      Set<SUID> state;
      for( auto it = states.begin(); it != states.end(); ++it )
         state += SUID( *it );

      uid::GlobalState::instance()->configure( state );
   }
}



} // namespace walberla
