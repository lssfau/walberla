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
//! \file Environment.h
//! \ingroup core
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//! \author Martin Bauer <martin.bauer@fau.de>
//! \brief Initialization of waLBerla
//
//======================================================================================================================

#pragma once

#include "DataTypes.h"
#include "core/config/Config.h"
#include "core/mpi/Environment.h"


namespace walberla {


//**********************************************************************************************************************
/*! RAII Object to initialize waLBerla using command line parameters
*/
//**********************************************************************************************************************
class Environment
{
public:

   /****************************************************************************************************************//**
   * Constructor
   *
   * \param argc,argv  If command line parameters are present they have to contain at least
   *                   the path to the configuration file and optionally pairs of param-value:
   *                   "-blockName.parameterName=parameterValue"
   *                   These values are then substituted in configuration files.
   *                   It is also possible to pass no command line options (see description below )
   *
   * If command line arguments are present the constructor initializes static objects
   * and singleton classes using information of configuration file
   *  - MPIManager ( via mpi::Environment )
   *  - reads configuration file
   *  - configures Logging
   *  - configures GlobalState
   *
   * If no command line arguments are present only MPI subsystem is initialized and
   *  - logging is configured to default state
   *  - the global state is configured to a default state (the global state still can be configured
   *    by manually calling GlobalState::instance()->configure())
   *
   ********************************************************************************************************************/
   Environment( int & argc, char ** & argv, bool mpiAbortOnException = true );

   ~Environment();

   /****************************************************************************************************************//**
   * Returns configuration object, or null if no configuration object exists
   *
   * \return handle to configuration object, if a path to configuration file was
   *         in the command line arguments as passed to constructor, otherwise a null pointer
   *
   ********************************************************************************************************************/
   shared_ptr<Config> config() { return config_; }

protected:
   shared_ptr<Config> config_;
   mpi::Environment mpiEnv_;
};



//======================================================================================================================
//
//  Function for manual initialization, if possible use Environment
//
//======================================================================================================================

//**********************************************************************************************************************
/*! Configures the global state given a config object
*/
//**********************************************************************************************************************
void configureGlobalState( const shared_ptr<Config> & config );



} // namespace walberla
