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
//! \file CreateConfig.h
//! \ingroup python
//! \author Martin Bauer <martin.bauer@fau.de>
//! \author Markus Holzer <markus.holzer@fau.de>
//! \brief Creates a walberla::Config object from a python script
//
//======================================================================================================================

#pragma once

#include "core/config/Config.h"
#include "core/config/Iterator.h"
#include "core/DataTypes.h"
#include <string>
#include <vector>

namespace walberla {
namespace python_coupling {


//** General Functions - fallback to non-python *****************************************************************
/*! \name General Functions - fallback to non-python  */
//@{

   shared_ptr<Config> createConfig( int argc, char ** argv );



          config::Iterator configBegin( int argc, char ** argv);
   inline config::Iterator configEnd()                              { return config::Iterator(); }

//@}
//****************************************************************************************************************





//** Functions to load config directly from script    ************************************************************
/*! \name Functions to load config directly from script  */
//@{

   shared_ptr<Config> createConfigFromPythonScript( const std::string & scriptFile,
                                                    const std::string & pythonFunctionName = "config",
                                                    const std::vector<std::string> & argv = std::vector<std::string>() );


   config::Iterator createConfigIteratorFromPythonScript( const std::string & scriptFile,
                                                          const std::string & pythonFunctionName = "config",
                                                          const std::vector<std::string> & argv = std::vector<std::string>() );

//@}
//****************************************************************************************************************



} // namespace python_coupling
} // namespace walberla


