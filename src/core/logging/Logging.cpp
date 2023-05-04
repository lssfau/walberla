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
//! \file Logging.cpp
//! \ingroup core
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#include "Logging.h"

#include "core/Filesystem.h"
#include "core/Regex.h"
#include <ctime>


namespace walberla {
namespace logging {



const std::string Logging::ERROR_TAG    = std::string( "[ERROR   ]" );
const std::string Logging::DEVEL_TAG    = std::string( "[DEVEL   ]" );
const std::string Logging::RESULT_TAG   = std::string( "[RESULT  ]" );
const std::string Logging::WARNING_TAG  = std::string( "[WARNING ]" );
const std::string Logging::INFO_TAG     = std::string( "[INFO    ]" );
const std::string Logging::PROGRESS_TAG = std::string( "[PROGRESS]" );
const std::string Logging::DETAIL_TAG   = std::string( "[DETAIL  ]" );
const std::string Logging::TRACING_TAG  = std::string( "[TRACING ]" );

const uint_t Logging::TAG_WIDTH       = uint_t(10);
const uint_t Logging::TIMESTAMP_WIDTH = uint_t(17);

void Logging::setStreamLogLevel( LogLevel logLevel )
{
#ifndef WALBERLA_LOGLEVEL_INFO
   if( logLevel == INFO )
      logWarning( "You are trying to set the stream log level to INFO, but INFO logs are deactivated by CMake!"
                  "The current WALBERLA_LOGLEVEL is: " WALBERLA_LOGLEVEL_STRING "!"
                  "Set WALBERLA_LOGLEVEL=INFO to activate INFO logs.",
                  "Logging::setStreamLogLevel", -1 );
#endif
#ifndef WALBERLA_LOGLEVEL_PROGRESS
   if( logLevel == PROGRESS )
      logWarning( "You are trying to set the stream log level to PROGRESS, but PROGRESS logs are deactivated by CMake!"
                  "The current WALBERLA_LOGLEVEL is: " WALBERLA_LOGLEVEL_STRING "!"
                  "Set WALBERLA_LOGLEVEL=PROGRESS to activate PROGRESS logs.",
                 "Logging::setStreamLogLevel", -1 );
#endif
#ifndef WALBERLA_LOGLEVEL_DETAIL
   if( logLevel == DETAIL )
      logWarning( "You are trying to set the stream log level to DETAIL, but DETAIL logs are deactivated by CMake!"
                  "The current WALBERLA_LOGLEVEL is: " WALBERLA_LOGLEVEL_STRING "!"
                  "Set WALBERLA_LOGLEVEL=DETAIL to activate DETAIL logs.",
                  "Logging::setStreamLogLevel", -1 );
#endif
#ifndef WALBERLA_LOGLEVEL_TRACING
   if( logLevel == TRACING )
      logWarning( "You are trying to set the stream log level to TRACING, but TRACING logs are deactivated by CMake!"
                  "The current WALBERLA_LOGLEVEL is: " WALBERLA_LOGLEVEL_STRING "!"
                  "Set WALBERLA_LOGLEVEL=TRACING to activate TRACING logs.",
                  "Logging::setStreamLogLevel", -1 );
#endif
   streamLogLevel_ = logLevel;
}



void Logging::setFileLogLevel( LogLevel logLevel )
{
#ifndef WALBERLA_LOGLEVEL_INFO
   if( logLevel == INFO )
      logWarning( "You are trying to set the file log level to INFO, but INFO logs are deactivated by CMake!"
                  "The current WALBERLA_LOGLEVEL is: " WALBERLA_LOGLEVEL_STRING "!"
                  "Set WALBERLA_LOGLEVEL=INFO to activate INFO logs.",
                  "Logging::setFileLogLevel", -1 );
#endif
#ifndef WALBERLA_LOGLEVEL_PROGRESS
   if( logLevel == PROGRESS )
      logWarning( "You are trying to set the file log level to PROGRESS, but PROGRESS logs are deactivated by CMake!"
                  "The current WALBERLA_LOGLEVEL is: " WALBERLA_LOGLEVEL_STRING "!"
                  "Set WALBERLA_LOGLEVEL=PROGRESS to activate PROGRESS logs.",
                  "Logging::setFileLogLevel", -1 );
#endif
#ifndef WALBERLA_LOGLEVEL_DETAIL
   if( logLevel == DETAIL )
      logWarning( "You are trying to set the file log level to DETAIL, but DETAIL logs are deactivated by CMake!"
                  "The current WALBERLA_LOGLEVEL is: " WALBERLA_LOGLEVEL_STRING "!"
                  "Set WALBERLA_LOGLEVEL=DETAIL to activate DETAIL logs.",
                  "Logging::setFileLogLevel", -1 );
#endif
#ifndef WALBERLA_LOGLEVEL_TRACING
   if( logLevel == TRACING )
      logWarning( "You are trying to set the file log level to TRACING, but TRACING logs are deactivated by CMake!"
                  "The current WALBERLA_LOGLEVEL is: " WALBERLA_LOGLEVEL_STRING "!"
                  "Set WALBERLA_LOGLEVEL=TRACING to activate TRACING logs.",
                  "Logging::setFileLogLevel", -1 );
#endif
   fileLogLevel_ = logLevel;
}



void Logging::includeLoggingToFile( const std::string & file, bool append )
{
   if( file_.is_open() )
   {
      std::cerr << "Error: Trying to open a file for logging while another file is already opened and used for logging!" << std::endl;
      throw std::runtime_error( "Error: Trying to open a file for logging while another file is already opened and used for logging!" );
   }

   std::ostringstream rank;
   rank << std::setfill('0') << std::setw( int_c( std::ceil( std::log10( real_c( numberOfProcesses_ ) ) ) ) ) << processId_;

   std::ostringstream filename;
   if( file.empty() )
   {
      filesystem::path logDir( "logging" );
      if( !filesystem::exists( logDir ) )
         filesystem::create_directory( logDir );

      if( numberOfProcesses_ == 1 )
         filename << "logging/logfile.txt";
      else
         filename << "logging/process_" << rank.str() << ".log";
   }
   else
   {
      if( numberOfProcesses_ == 1 )
         filename << file;
      else
      {
         filesystem::path filePath( file );
         const std::string stem      = filePath.stem().string();
         const std::string extension = filePath.extension().string();
         const std::string modifiedFilename = stem + std::string("_") + rank.str() + extension;
         filePath.remove_filename();
         filePath /= modifiedFilename;
         filename << filePath.string();
      }
   }

   if( append )
      file_.open( filename.str().c_str(), std::ofstream::out | std::ofstream::app );
   else
      file_.open( filename.str().c_str(), std::ofstream::out );

   if( file_.fail() || file_.bad() )
   {
      std::cerr << "Error: Opening file '" << filename.str() << "' failed!" << std::endl;
      throw std::runtime_error( std::string("Error: Opening file '") + filename.str() + std::string("' failed!") );
   }

   file_ << getHeaderFooter( true ) << std::endl;
}



std::string Logging::getHeaderFooter( bool header )
{
   std::time_t t;
   std::time( &t );

   char cTimeString[64];
   std::strftime( cTimeString, 64, "%A, %d.%B %Y, %H:%M:%S", std::localtime( &t ) );
   std::string timeString( cTimeString );

   std::string beginEnd( header ? "BEGIN LOGGING" : "  END LOGGING" );

   std::ostringstream oss;
   oss << "================================================================================\n"
       << "|       " << beginEnd << " - " << timeString
       << std::setw( 70 - int_c( timeString.length() ) - int_c( beginEnd.length() ) ) << std::setfill(' ') << std::right << "|\n"
       << "================================================================================";

   return oss.str();
}



std::string Logging::createLog( const std::string & type, const std::string & message,
                                const std::string & callerPath, const int line ) const
{
   std::ostringstream log;
   log << getRankStamp() << type;

   if( showTimeStamp_ )
   {
      std::ostringstream timeStamp;
      timeStamp << std::setw( int_c(TIMESTAMP_WIDTH) ) << std::setfill('-') << std::right << getTimeStamp();
      log << timeStamp.str();
   }

   if( additionalStamp_ )
   {
      std::ostringstream customStamp;
      customStamp << std::setw( int_c(additionalStamp_->maxStampWidth()) ) << std::setfill('-') << std::right << additionalStamp_->stamp();
      log << customStamp.str();
   }

   log << " " << message;
   if( logCallerPath_ )
      log << "\n\n(from: " << callerPath << ":" << line << ")";

   return resetLinebreaks( log.str() );
}


bool Logging::isInIgnoreCallerPaths( const std::vector< walberla::regex > & regexes,
                                            const std::string & callerPath, const int line ) const
{
   if( !regexes.empty() )
   {
      std::stringstream callerPathAndLine;
      callerPathAndLine << callerPath << ":" << line;

      for( auto regex = regexes.begin(); regex != regexes.end(); ++regex )
         if( walberla::regex_search( callerPathAndLine.str(), *regex ) )
            return true;
   }

   return false;
}



} // namespace logging
} // namespace walberla