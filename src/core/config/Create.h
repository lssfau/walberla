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
//! \file Create.h
//! \ingroup core
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#pragma once

#include "Config.h"
#include "Iterator.h"


namespace walberla {
namespace config {


   //** Load configuration file(s) using command line arguments        *************************************************
   /*! \name Load configuration file using command line arguments  */
   //@{



   //*******************************************************************************************************************
   /*! Loads configuration file using command line arguments
   *
   * The last command line parameter which does not begin with '-' is taken as the configuration file name.
   * Parameters starting with '-' are used to overwrite individual configuration values.
   *
   * For the syntax of the configuration file see config::Config
   *
   * For details on Parameter substitution using additional command line parameters see
   * config::substituteCommandLineArgs
   */
   //*******************************************************************************************************************
   shared_ptr<Config> create(int argc, char ** argv);






   //*******************************************************************************************************************
   /*! Loads multiple numbered configuration files.
   *
   * In order to do a parameter study i.e. running multiple simulations in a row this function can be used.
   * The configuration file name has to have the following structure:  someName000.someExtension
   * If there are zeros before the extension dot, it is interpreted as the first file of a numbered series.
   *
   * If there are no zeros in the filename before the dot, the iterator goes over the single file i.e. its save
   * to always use the iterator variant.
   *
   * Command line substitution does not work when iterating over multiple files.
   */
   //*******************************************************************************************************************
          Iterator begin( int argc, char ** argv);
   inline Iterator end()                           { return config::Iterator(); }



   //@}
   //*******************************************************************************************************************











   //** Helper functions  ********************************************************************************************
   /*! \name Helper functions  */
   //@{
   void createFromTextFile( Config & out, const std::string & pathToTextFile );

   Iterator createConfigIteratorFromTextFiles ( const std::string & baseName,
                                                const std::string & extension=".prm",
                                                int numberOfDigits = 3 );


   std::string usageString(const std::string & executableName);
   void substituteCommandLineArgs( Config & config, int argc, char**argv );
   void substituteCommandLineArgs( Config & config, const std::vector<std::string> & params );
   //@}
   //****************************************************************************************************************

} // namespace config
} // namespace walberla




