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
//! \file Initialization.cpp
//! \ingroup core
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#include "Initialization.h"
#include "Logging.h"
#include "core/StringUtility.h"

#include <algorithm>
#include <sstream>
#include <string>


namespace walberla {
namespace logging {



static void parseIgnoreBlocks( const Config::Blocks & ignoreBlocks, std::vector< walberla::regex > & regexes )
{
   for( auto ignoreBlock = ignoreBlocks.begin(); ignoreBlock != ignoreBlocks.end(); ++ignoreBlock )
   {
      if( ignoreBlock->isDefined( "callerPathPattern" ) )
      {
         std::string regexString = ignoreBlock->getParameter<std::string>( "callerPathPattern" );
         // Replace slashes with a regex to allow for windows/linux compatibility
         string_replace_all( regexString, "/", "(\\\\|/)" );
         try
         {
            walberla::regex regex( regexString );
            regexes.push_back( regex );
         }
         catch( const walberla::regex_error & e )
         {
            std::ostringstream oss;
            oss << __FILE__ << ":" << __LINE__ << " - Error parsing regular Expression \"" << regexString
                << "\" from ignore block at position " << e.code() << ".";
            throw std::runtime_error( oss.str() );
         }
      }
      else
      {
         std::ostringstream oss;
         oss << __FILE__ << ":" << __LINE__ << " - You specified an ignore block without an \"callerPathPattern\" key!";
         throw std::runtime_error( oss.str() );
      }
   }
}



//**********************************************************************************************************************
/*!
*   \brief Function for initializing the logging singleton from file
*
*   This initialization function reads data stored in a configuration file and uses thi data to configure the logging
*   singleton.
*   The structure of the configuration file which is recognized by this function must look like as follows:
*   \code
*   Logging { // the whole 'Logging' block is optional
*
*      logLevel       warning|info|progress|detail|tracing; // optional, default=info
*      streamLogLevel warning|info|progress|detail|tracing; // optional, default=info
*      fileLogLevel   warning|info|progress|detail|tracing; // optional, default=info
*
*      logfile [filename]; // optional, default=no logfiles are created, logging only uses stdout & stderr
*      append true|false;  // optional, only valid if 'logfile' is provided, default=false
*
*      time true|false; // optional, default=true
*
*      logCallerPath true|false; // optional, default=false
*
*      ignore { callerPathPattern [regex]; } // optional
*      ignore { callerPathPattern [regex]; }
*      ignore ...
*
*      ignoreWarning  { callerPathPattern [regex]; } // optional
*      ignoreWarning  { callerPathPattern [regex]; }
*      ignoreWarning  ...
*
*   } // Logging
*   \endcode
*
*   Ignoring of log messages:
*
*   Log messages can be ignored by specifying a regular expression on the filename and line number of the source
*   file the log message originates from. The regex has to be specified in an ignore block as the parameter
*   callerPathPattern. If the regular expression is found somewhere in the string <filename>:<line number>,
*   the log message is ignored. To divide components of your path always uses slashes!
*   They will automatically be converted to the regex (/|\\), matching slashes and back-slashes. For fancy regexes
*   you can use perl regex syntax.
*   To ignore warnings use the spcecial ignoreWarning block.
*   Note that you cannot ignore Errors since they abort your program!
*
*   Examples:
*   ignore { callerPathPattern /src/core/; } // ignores log messages from core
*   ignore { callerPathPattern /src/vtk/; } // ignores log messages from vtk
*   ignore { callerPathPattern (/src/core/|/src/vtk/); } // ignores log messages from core & vtk
*   ignore { callerPathPattern /src/core/FILENAME.h:416; } // ignores a specific log message
*   ignoreWarning { callerPathPattern /src/core/FILENAME.h:212; } // ignores a specific warning
*
*   \param config The configuration
*/
//**********************************************************************************************************************
void configureLogging( const shared_ptr< Config > & config )
{
   if( !!config )
   {
      auto block = config->getGlobalBlock();
      if( block )
         configureLogging( block.getBlock( "Logging" ) );
   }
}



void configureLogging( const Config::BlockHandle & loggingBlock )
{
   if( !loggingBlock )
      return;

   if( loggingBlock.isDefined( "logLevel" ) )
   {
      std::string type = loggingBlock.getParameter< std::string >( "logLevel" );
      string_to_lower( type );

      if( type == "warning" ){
         logging::Logging::instance()->setLogLevel( logging::Logging::WARNING );
      } else if( type == "info" ){
         logging::Logging::instance()->setLogLevel( logging::Logging::INFO );
      } else if( type == "progress" ){
         logging::Logging::instance()->setLogLevel( logging::Logging::PROGRESS );
      } else if( type == "detail" ){
         logging::Logging::instance()->setLogLevel( logging::Logging::DETAIL );
      } else if( type == "tracing" ){
         logging::Logging::instance()->setLogLevel( logging::Logging::TRACING );
      } else
         throw std::runtime_error("Error: Unknown parameter for 'logLevel'. Possible parameters are: warning|info|progress|detail|tracing");
   }

   if( loggingBlock.isDefined( "streamLogLevel" ) )
   {
      std::string type = loggingBlock.getParameter< std::string >( "streamLogLevel" );
      string_to_lower( type );

      if( type == "warning" ){
         logging::Logging::instance()->setStreamLogLevel( logging::Logging::WARNING );
      } else if( type == "info" ){
         logging::Logging::instance()->setStreamLogLevel( logging::Logging::INFO );
      } else if( type == "progress" ){
         logging::Logging::instance()->setStreamLogLevel( logging::Logging::PROGRESS );
      } else if( type == "detail" ){
         logging::Logging::instance()->setStreamLogLevel( logging::Logging::DETAIL );
      } else if( type == "tracing" ){
         logging::Logging::instance()->setStreamLogLevel( logging::Logging::TRACING );
      } else
         throw std::runtime_error("Error: Unknown parameter for 'streamLogLevel'. Possible parameters are: warning|info|progress|detail|tracing");
   }

   if( loggingBlock.isDefined( "fileLogLevel" ) )
   {
      std::string type = loggingBlock.getParameter< std::string >( "fileLogLevel" );
      string_to_lower( type );

      if( type == "warning" ){
         logging::Logging::instance()->setFileLogLevel( logging::Logging::WARNING );
      } else if( type == "info" ){
         logging::Logging::instance()->setFileLogLevel( logging::Logging::INFO );
      } else if( type == "progress" ){
         logging::Logging::instance()->setFileLogLevel( logging::Logging::PROGRESS );
      } else if( type == "detail" ){
         logging::Logging::instance()->setFileLogLevel( logging::Logging::DETAIL );
      } else if( type == "tracing" ){
         logging::Logging::instance()->setFileLogLevel( logging::Logging::TRACING );
      } else
         throw std::runtime_error("Error: Unknown parameter for 'fileLogLevel'. Possible parameters are: warning|info|progress|detail|tracing");
   }

   if( loggingBlock.isDefined( "logfile" ) )
   {
      const std::string filename = loggingBlock.getParameter< std::string >( "logfile" );
      const bool          append = loggingBlock.getParameter< bool >( "append", false );

      logging::Logging::instance()->stopLoggingToFile();
      logging::Logging::instance()->includeLoggingToFile( filename, append );
   }
   else if( loggingBlock.isDefined( "append" ) )
   {
      throw std::runtime_error("Error: Found parameter 'append' without a corresponding parameter 'logfile'.");
   }

   logging::Logging::instance()->showTimeStamp( loggingBlock.getParameter< bool >( "time", true ) );
   logging::Logging::instance()->logCallerPath( loggingBlock.getParameter< bool >( "logCallerPath", false ) );

   Config::Blocks ignoreBlocks;
   Config::Blocks ignoreWarningBlocks;

   loggingBlock.getBlocks( "ignore", ignoreBlocks );
   std::vector< walberla::regex > regexes;
   parseIgnoreBlocks( ignoreBlocks, regexes );
   for( auto regex = regexes.begin(); regex != regexes.end(); ++regex )
      logging::Logging::instance()->addIgnoreRegex( *regex );

   regexes.clear();
   loggingBlock.getBlocks( "ignoreWarning", ignoreWarningBlocks );
   parseIgnoreBlocks( ignoreWarningBlocks, regexes );
   for( auto regex = regexes.begin(); regex != regexes.end(); ++regex )
      logging::Logging::instance()->addIgnoreWarningRegex( *regex );
}



} // namespace logging
} // namespace walberla
