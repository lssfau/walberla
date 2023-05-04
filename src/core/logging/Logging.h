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
//! \file Logging.h
//! \ingroup core
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "waLBerlaDefinitions.h"

#include "core/DataTypes.h"
#include "core/logging/CMakeDefs.h"
#include "core/mpi/MPIManager.h"
#include "core/singleton/Singleton.h"
#include "core/timing/WcPolicy.h"
#include "core/Regex.h"
#include "core/StringUtility.h"

#include <functional>

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>


namespace walberla {

class Abort;       // forward declaration
class Environment; // forward declaration

namespace logging {


namespace internal { class LoggingMacroDispatcher; } // forward declaration
class Tracer;                                        // forward declaration


class Logging : public singleton::Singleton<Logging>
{
   WALBERLA_BEFRIEND_SINGLETON;

public:

   friend class walberla::Abort;
   friend class walberla::Environment;
   friend class internal::LoggingMacroDispatcher;
   friend class Tracer;

   enum LogLevel { WARNING = 0, INFO = 1, PROGRESS = 2, DETAIL = 3, TRACING = 4 };

   static const std::string ERROR_TAG;
   static const std::string DEVEL_TAG;
   static const std::string RESULT_TAG;
   static const std::string WARNING_TAG;
   static const std::string INFO_TAG;
   static const std::string PROGRESS_TAG;
   static const std::string DETAIL_TAG;
   static const std::string TRACING_TAG;

   static const uint_t TAG_WIDTH;
   static const uint_t TIMESTAMP_WIDTH;

   class CustomStamp
   {
   public:
      virtual ~CustomStamp() = default;
      virtual std::string stamp() = 0;
      virtual uint_t      maxStampWidth() = 0;
   };



   inline ~Logging();

   void setLogLevel( LogLevel logLevel ) { setStreamLogLevel(logLevel); setFileLogLevel(logLevel); }

   void     setStreamLogLevel( LogLevel logLevel );
   LogLevel getStreamLogLevel() const { return streamLogLevel_; }

   void     setFileLogLevel( LogLevel logLevel );
   LogLevel getFileLogLevel() const { return fileLogLevel_; }

          void includeLoggingToFile( const std::string & file = std::string(""), bool append = false );
   inline void    stopLoggingToFile();
   inline bool        loggingToFile() const { return file_.is_open(); }

   void showTimeStamp( bool show ) { showTimeStamp_ = show; }
   bool showTimeStamp() const { return showTimeStamp_; }

   inline void   addCustomStamp( const shared_ptr< CustomStamp > & stamp  ) { additionalStamp_ = stamp; }
   inline void clearCustomStamp() { additionalStamp_.reset(); }
   inline bool  usesCustomStamp() const { return !!additionalStamp_; }

   void logCallerPath( bool log ) { logCallerPath_ = log; }
   bool logCallerPath() const { return logCallerPath_; }

   void                                addIgnoreRegex( const walberla::regex & regex ) { ignoreRegexes_.push_back( regex ); }
   const std::vector< walberla::regex > & getIgnoreRegexes() const { return ignoreRegexes_; }
   void                              clearIgnoreRegexes() { ignoreRegexes_.clear(); }

   void                                addIgnoreWarningRegex( const walberla::regex & regex ) { ignoreWarningRegexes_.push_back( regex ); }
   const std::vector< walberla::regex > & getIgnoreWarningRegexes() const { return ignoreWarningRegexes_; }
   void                              clearIgnoreWarningRegexes() { ignoreWarningRegexes_.clear(); }



   static inline void printHeaderOnStream();
   static inline void printFooterOnStream();

   bool onRoot() const { return processId_ == 0; }

   inline bool logDevel   ( const std::string & callerPath, const int line );
   inline bool logResult  ( const std::string & callerPath, const int line );
   inline bool logWarning ( const std::string & callerPath, const int line );
   inline bool logInfo    ( const std::string & callerPath, const int line );
   inline bool logProgress( const std::string & callerPath, const int line );
   inline bool logDetail  ( const std::string & callerPath, const int line );
   inline bool logTracing ( const std::string & callerPath, const int line );

private:

   inline Logging();

   bool isInIgnoreCallerPaths( const std::vector< walberla::regex > & regexes,
                               const std::string & callerPath, const int line ) const;

   inline void logError   ( const std::string & message, const std::string & callerPath, const int line );
   inline void logDevel   ( const std::string & message, const std::string & callerPath, const int line );
   inline void logResult  ( const std::string & message, const std::string & callerPath, const int line );
   inline void logWarning ( const std::string & message, const std::string & callerPath, const int line );
   inline void logInfo    ( const std::string & message, const std::string & callerPath, const int line );
   inline void logProgress( const std::string & message, const std::string & callerPath, const int line );
   inline void logDetail  ( const std::string & message, const std::string & callerPath, const int line );
   inline void logTracing ( const std::string & message, const std::string & callerPath, const int line );

   inline void logMessage( LogLevel logLevel, std::ostream & stream, const std::string & message );

   inline std::string getRankStamp() const;
   inline std::string resetLinebreaks( const std::string & message ) const;

   inline std::string getTimeStamp() const;
   static std::string getHeaderFooter( bool header );

   std::string createLog( const std::string & type, const std::string & message,
                          const std::string & callerPath, const int line ) const;



   LogLevel streamLogLevel_;
   LogLevel   fileLogLevel_;

   std::ofstream file_;

   uint_t processId_;
   uint_t numberOfProcesses_;

   double startTime_;
   bool showTimeStamp_;
   shared_ptr< CustomStamp > additionalStamp_;

   bool logCallerPath_;

   std::vector< walberla::regex > ignoreRegexes_;
   std::vector< walberla::regex > ignoreWarningRegexes_;
};



inline Logging::Logging() : singleton::Singleton<Logging>(),
   streamLogLevel_( WALBERLA_LOGLEVEL ), fileLogLevel_( WALBERLA_LOGLEVEL ),
   processId_( uint_c( mpi::MPIManager::instance()->worldRank() ) ),
   numberOfProcesses_( uint_c( mpi::MPIManager::instance()->numProcesses() ) ),
   startTime_( timing::WcPolicy::getTimestamp() ), showTimeStamp_( true ), logCallerPath_( false )
{
}



inline Logging::~Logging()
{
   stopLoggingToFile();
}



inline void Logging::stopLoggingToFile()
{
   if( file_.is_open() )
   {
      file_ << getHeaderFooter( false ) << std::endl;
      file_.close();
   }
}



inline void Logging::printHeaderOnStream()
{
   WALBERLA_ROOT_SECTION()
   {
      static bool headerWritten = false; // can be written only once!
      if( !headerWritten )
         std::cout << getHeaderFooter( true ) << std::endl;
      headerWritten = true;
   }
   WALBERLA_MPI_WORLD_BARRIER();
}



inline void Logging::printFooterOnStream()
{
   WALBERLA_MPI_WORLD_BARRIER();
   WALBERLA_ROOT_SECTION()
   {
      static bool footerWritten = false; // can be written only once!
      if( !footerWritten )
         std::cout << getHeaderFooter( false ) << std::endl;
      footerWritten = true;
   }
}



inline bool Logging::logDevel( const std::string & callerPath, const int line )
{
   return !isInIgnoreCallerPaths( ignoreRegexes_, callerPath, line );
}

inline bool Logging::logResult( const std::string & callerPath, const int line )
{
   return !isInIgnoreCallerPaths( ignoreRegexes_, callerPath, line );
}

inline bool Logging::logWarning( const std::string & callerPath, const int line )
{
   return !isInIgnoreCallerPaths( ignoreWarningRegexes_, callerPath, line );
}

#ifdef WALBERLA_LOGLEVEL_INFO
inline bool Logging::logInfo( const std::string & callerPath, const int line )
{
   return ( streamLogLevel_ >= INFO || ( fileLogLevel_ >= INFO && file_.is_open() ) ) &&
          !isInIgnoreCallerPaths( ignoreRegexes_, callerPath, line );
#else
inline bool Logging::logInfo( const std::string & /*callerPath*/, const int /*line*/ )
{
   return false;
#endif
}

#ifdef WALBERLA_LOGLEVEL_PROGRESS
inline bool Logging::logProgress( const std::string & callerPath, const int line )
{
   return ( streamLogLevel_ >= PROGRESS || ( fileLogLevel_ >= PROGRESS && file_.is_open() ) ) &&
          !isInIgnoreCallerPaths( ignoreRegexes_, callerPath, line );
#else
inline bool Logging::logProgress( const std::string & /*callerPath*/, const int /*line*/ )
{
   return false;
#endif
}

#ifdef WALBERLA_LOGLEVEL_DETAIL
inline bool Logging::logDetail( const std::string & callerPath, const int line )
{
   return ( streamLogLevel_ >= DETAIL || ( fileLogLevel_ >= DETAIL && file_.is_open() ) ) &&
          !isInIgnoreCallerPaths( ignoreRegexes_, callerPath, line );
#else
inline bool Logging::logDetail( const std::string & /*callerPath*/, const int /*line*/ )
{
   return false;
#endif
}

#ifdef WALBERLA_LOGLEVEL_TRACING
inline bool Logging::logTracing( const std::string & callerPath, const int line )
{
   return ( streamLogLevel_ >= TRACING || ( fileLogLevel_ >= TRACING && file_.is_open() ) ) &&
          !isInIgnoreCallerPaths( ignoreRegexes_, callerPath, line );
#else
inline bool Logging::logTracing( const std::string & /*callerPath*/, const int /*line*/ )
{
   return false;
#endif
}



inline void Logging::logError( const std::string & message, const std::string & callerPath, const int line )
{
#ifdef WALBERLA_THREAD_SAFE_LOGGING
#ifdef _OPENMP
   #pragma omp critical (logging)
   {
#endif
#endif
   logMessage( WARNING, std::cerr, createLog( ERROR_TAG, message, callerPath, line ) ); // WARNING has the highest priority and will always be logged!
#ifdef WALBERLA_THREAD_SAFE_LOGGING
#ifdef _OPENMP
   }
#endif
#endif
}

inline void Logging::logDevel( const std::string & message, const std::string & callerPath, const int line )
{
#ifdef WALBERLA_THREAD_SAFE_LOGGING
#ifdef _OPENMP
   #pragma omp critical (logging)
   {
#endif
#endif
   if( logDevel( callerPath, line ) )
      logMessage( WARNING, std::cout, createLog( DEVEL_TAG, message, callerPath, line ) ); // WARNING has the highest priority and will always be logged!
#ifdef WALBERLA_THREAD_SAFE_LOGGING
#ifdef _OPENMP
   }
#endif
#endif
}

inline void Logging::logResult( const std::string & message, const std::string & callerPath, const int line )
{
#ifdef WALBERLA_THREAD_SAFE_LOGGING
#ifdef _OPENMP
   #pragma omp critical (logging)
   {
#endif
#endif
   if( logResult( callerPath, line ) )
      logMessage( WARNING, std::cout, createLog( RESULT_TAG, message, callerPath, line ) ); // WARNING has the highest priority and will always be logged!
#ifdef WALBERLA_THREAD_SAFE_LOGGING
#ifdef _OPENMP
   }
#endif
#endif
}

inline void Logging::logWarning( const std::string & message, const std::string & callerPath, const int line )
{
#ifdef WALBERLA_THREAD_SAFE_LOGGING
#ifdef _OPENMP
   #pragma omp critical (logging)
   {
#endif
#endif
   if( logWarning( callerPath, line ) )
      logMessage( WARNING, std::cout, createLog( WARNING_TAG, message, callerPath, line ) );
#ifdef WALBERLA_THREAD_SAFE_LOGGING
#ifdef _OPENMP
   }
#endif
#endif
}

inline void Logging::logInfo( const std::string & message, const std::string & callerPath, const int line )
{
#ifdef WALBERLA_THREAD_SAFE_LOGGING
#ifdef _OPENMP
   #pragma omp critical (logging)
   {
#endif
#endif
   if( logInfo( callerPath, line ) )
      logMessage( INFO, std::cout, createLog( INFO_TAG, message, callerPath, line ) );
#ifdef WALBERLA_THREAD_SAFE_LOGGING
#ifdef _OPENMP
   }
#endif
#endif
}

inline void Logging::logProgress( const std::string & message, const std::string & callerPath, const int line )
{
#ifdef WALBERLA_THREAD_SAFE_LOGGING
#ifdef _OPENMP
   #pragma omp critical (logging)
   {
#endif
#endif
   if( logProgress( callerPath, line ) )
      logMessage( PROGRESS, std::cout, createLog( PROGRESS_TAG, message, callerPath, line ) );
#ifdef WALBERLA_THREAD_SAFE_LOGGING
#ifdef _OPENMP
   }
#endif
#endif
}

inline void Logging::logDetail( const std::string & message, const std::string & callerPath, const int line )
{
#ifdef WALBERLA_THREAD_SAFE_LOGGING
#ifdef _OPENMP
   #pragma omp critical (logging)
   {
#endif
#endif
   if( logDetail( callerPath, line ) )
      logMessage( DETAIL, std::cout, createLog( DETAIL_TAG, message, callerPath, line ) );
#ifdef WALBERLA_THREAD_SAFE_LOGGING
#ifdef _OPENMP
   }
#endif
#endif
}

inline void Logging::logTracing( const std::string & message, const std::string & callerPath, const int line )
{
#ifdef WALBERLA_THREAD_SAFE_LOGGING
#ifdef _OPENMP
   #pragma omp critical (logging)
   {
#endif
#endif
   if( logTracing( callerPath, line ) )
      logMessage( TRACING, std::cout, createLog( TRACING_TAG, message, callerPath, line ) );
#ifdef WALBERLA_THREAD_SAFE_LOGGING
#ifdef _OPENMP
   }
#endif
#endif
}



inline void Logging::logMessage( LogLevel logLevel, std::ostream & stream, const std::string & message )
{
   if( streamLogLevel_ >= logLevel )
      stream << message << std::endl;
   if( file_.is_open() && fileLogLevel_ >= logLevel )
      file_ <<  message << std::endl;
}



inline std::string Logging::getRankStamp() const
{
   std::stringstream rank;
   rank << std::setfill(' ') << std::setw( int_c( std::ceil( std::log10( real_c( numberOfProcesses_ ) ) ) ) ) << processId_;
   return std::string("[") + rank.str() + std::string("]");
}



inline std::string Logging::resetLinebreaks( const std::string & message ) const
{
   std::string newline = std::string("\n") + getRankStamp();
   newline.append( TAG_WIDTH + ( showTimeStamp_ ? TIMESTAMP_WIDTH : uint_c(0) )
                             + ( additionalStamp_ ? additionalStamp_->maxStampWidth() : uint_c(0) ) + uint_c(1), ' ' );
   return string_replace_all_copy( message, "\n", newline );
}



inline std::string Logging::getTimeStamp() const
{
   std::ostringstream oss;
   oss << "(" << std::fixed << std::setprecision( 3 ) << ( timing::WcPolicy::getTimestamp() - startTime_ ) << " sec)";
   return oss.str();
}



//======================================================================================================================
//
// UTILITY LOGGING MACROS
//
//======================================================================================================================

//////////////////////
// Logging Sections //
//////////////////////

#define WALBERLA_LOG_DEVEL_SECTION()            if( walberla::logging::Logging::instance()->logDevel   ( __FILE__, __LINE__ ) )
#define WALBERLA_LOG_DEVEL_ON_ROOT_SECTION()    if( walberla::logging::Logging::instance()->onRoot() && walberla::logging::Logging::instance()->logDevel   ( __FILE__, __LINE__ ) )

#define WALBERLA_LOG_RESULT_SECTION()           if( walberla::logging::Logging::instance()->logResult  ( __FILE__, __LINE__ ) )
#define WALBERLA_LOG_RESULT_ON_ROOT_SECTION()   if( walberla::logging::Logging::instance()->onRoot() && walberla::logging::Logging::instance()->logResult  ( __FILE__, __LINE__ ) )

#define WALBERLA_LOG_WARNING_SECTION()          if( walberla::logging::Logging::instance()->logWarning ( __FILE__, __LINE__ ) )
#define WALBERLA_LOG_WARNING_ON_ROOT_SECTION()  if( walberla::logging::Logging::instance()->onRoot() && walberla::logging::Logging::instance()->logWarning ( __FILE__, __LINE__ ) )

#ifdef WALBERLA_LOGLEVEL_INFO
#   define WALBERLA_LOG_INFO_SECTION()             if( walberla::logging::Logging::instance()->logInfo    ( __FILE__, __LINE__ ) )
#   define WALBERLA_LOG_INFO_ON_ROOT_SECTION()     if( walberla::logging::Logging::instance()->onRoot() && walberla::logging::Logging::instance()->logInfo    ( __FILE__, __LINE__ ) )
#else
#   define WALBERLA_LOG_INFO_SECTION()             if( false )
#   define WALBERLA_LOG_INFO_ON_ROOT_SECTION()     if( false )
#endif

#ifdef WALBERLA_LOGLEVEL_PROGRESS
#   define WALBERLA_LOG_PROGRESS_SECTION()         if( walberla::logging::Logging::instance()->logProgress( __FILE__, __LINE__ ) )
#   define WALBERLA_LOG_PROGRESS_ON_ROOT_SECTION() if( walberla::logging::Logging::instance()->onRoot() && walberla::logging::Logging::instance()->logProgress( __FILE__, __LINE__ ) )
#else
#   define WALBERLA_LOG_PROGRESS_SECTION()         if( false )
#   define WALBERLA_LOG_PROGRESS_ON_ROOT_SECTION() if( false )
#endif

#ifdef WALBERLA_LOGLEVEL_DETAIL
#   define WALBERLA_LOG_DETAIL_SECTION()           if( walberla::logging::Logging::instance()->logDetail  ( __FILE__, __LINE__ ) )
#   define WALBERLA_LOG_DETAIL_ON_ROOT_SECTION()   if( walberla::logging::Logging::instance()->onRoot() && walberla::logging::Logging::instance()->logDetail  ( __FILE__, __LINE__ ) )
#else
#   define WALBERLA_LOG_DETAIL_SECTION()           if( false )
#   define WALBERLA_LOG_DETAIL_ON_ROOT_SECTION()   if( false )
#endif

////////////////////
// Logging Macros //
////////////////////

namespace internal {
class LoggingMacroDispatcher
{
public:
   static inline void logDevel( const std::string & message, const std::string & callerPath, const int line )
   {
      walberla::logging::Logging::instance()->logDevel( message, callerPath, line );
   }
   static inline void logResult( const std::string & message, const std::string & callerPath, const int line )
   {
      walberla::logging::Logging::instance()->logResult( message, callerPath, line );
   }
   static inline void logWarning( const std::string & message, const std::string & callerPath, const int line )
   {
      walberla::logging::Logging::instance()->logWarning( message, callerPath, line );
   }
   static inline void logInfo( const std::string & message, const std::string & callerPath, const int line )
   {
      walberla::logging::Logging::instance()->logInfo( message, callerPath, line );
   }
   static inline void logProgress( const std::string & message, const std::string & callerPath, const int line )
   {
      walberla::logging::Logging::instance()->logProgress( message, callerPath, line );
   }
   static inline void logDetail( const std::string & message, const std::string & callerPath, const int line )
   {
      walberla::logging::Logging::instance()->logDetail( message, callerPath, line );
   }
};
} // namespace internal

///////////////////
// ERROR / ABORT //
///////////////////

// --> Abort.h (WALBERLA_ABORT)

///////////
// DEVEL //
///////////

#define WALBERLA_LOG_DEVEL(msg){\
   WALBERLA_LOG_DEVEL_SECTION(){\
      std::stringstream WALBERLA__LOGGING__OSS;\
      WALBERLA__LOGGING__OSS << msg;\
      walberla::logging::internal::LoggingMacroDispatcher::logDevel( WALBERLA__LOGGING__OSS.str(), __FILE__, __LINE__ );\
   }\
}
#define WALBERLA_LOG_DEVEL_ON_ROOT(msg){\
   if( walberla::logging::Logging::instance()->onRoot() ){\
   WALBERLA_LOG_DEVEL_SECTION(){\
      std::stringstream WALBERLA__LOGGING__OSS;\
      WALBERLA__LOGGING__OSS << msg;\
      walberla::logging::internal::LoggingMacroDispatcher::logDevel( WALBERLA__LOGGING__OSS.str(), __FILE__, __LINE__ );\
   }}\
}

#define WALBERLA_LOG_DEVEL_VAR(variable){\
   WALBERLA_LOG_DEVEL_SECTION(){\
   std::stringstream WALBERLA__LOGGING__OSS;\
   WALBERLA__LOGGING__OSS << #variable " = " << variable;\
   walberla::logging::internal::LoggingMacroDispatcher::logDevel( WALBERLA__LOGGING__OSS.str(), __FILE__, __LINE__ );\
}\
}
#define WALBERLA_LOG_DEVEL_VAR_ON_ROOT(variable){\
   if( walberla::logging::Logging::instance()->onRoot() ){\
   WALBERLA_LOG_DEVEL_SECTION(){\
   std::stringstream WALBERLA__LOGGING__OSS;\
   WALBERLA__LOGGING__OSS << #variable " = " << variable;\
   walberla::logging::internal::LoggingMacroDispatcher::logDevel( WALBERLA__LOGGING__OSS.str(), __FILE__, __LINE__ );\
}}\
}

////////////
// RESULT //
////////////

#define WALBERLA_LOG_RESULT(msg){\
   WALBERLA_LOG_RESULT_SECTION(){\
      std::stringstream WALBERLA__LOGGING__OSS;\
      WALBERLA__LOGGING__OSS << msg;\
      walberla::logging::internal::LoggingMacroDispatcher::logResult( WALBERLA__LOGGING__OSS.str(), __FILE__, __LINE__ );\
   }\
}
#define WALBERLA_LOG_RESULT_ON_ROOT(msg){\
   if( walberla::logging::Logging::instance()->onRoot() ){\
   WALBERLA_LOG_RESULT_SECTION(){\
      std::stringstream WALBERLA__LOGGING__OSS;\
      WALBERLA__LOGGING__OSS << msg;\
      walberla::logging::internal::LoggingMacroDispatcher::logResult( WALBERLA__LOGGING__OSS.str(), __FILE__, __LINE__ );\
   }}\
}

/////////////
// WARNING //
/////////////

#define WALBERLA_LOG_WARNING(msg){\
   WALBERLA_LOG_WARNING_SECTION(){\
      std::stringstream WALBERLA__LOGGING__OSS;\
      WALBERLA__LOGGING__OSS << msg;\
      walberla::logging::internal::LoggingMacroDispatcher::logWarning( WALBERLA__LOGGING__OSS.str(), __FILE__, __LINE__ );\
   }\
}
#define WALBERLA_LOG_WARNING_ON_ROOT(msg){\
   if( walberla::logging::Logging::instance()->onRoot() ){\
   WALBERLA_LOG_WARNING_SECTION(){\
      std::stringstream WALBERLA__LOGGING__OSS;\
      WALBERLA__LOGGING__OSS << msg;\
      walberla::logging::internal::LoggingMacroDispatcher::logWarning( WALBERLA__LOGGING__OSS.str(), __FILE__, __LINE__ );\
   }}\
}

//////////
// INFO //
//////////

#ifdef WALBERLA_LOGLEVEL_INFO
#define WALBERLA_LOG_INFO(msg){\
   WALBERLA_LOG_INFO_SECTION(){\
      std::stringstream WALBERLA__LOGGING__OSS;\
      WALBERLA__LOGGING__OSS << msg;\
      walberla::logging::internal::LoggingMacroDispatcher::logInfo( WALBERLA__LOGGING__OSS.str(), __FILE__, __LINE__ );\
   }\
}
#else
#define WALBERLA_LOG_INFO(msg) (void(0))
#endif
#ifdef WALBERLA_LOGLEVEL_INFO
#define WALBERLA_LOG_INFO_ON_ROOT(msg){\
   if( walberla::logging::Logging::instance()->onRoot() ){\
   WALBERLA_LOG_INFO_SECTION(){\
      std::stringstream WALBERLA__LOGGING__OSS;\
      WALBERLA__LOGGING__OSS << msg;\
      walberla::logging::internal::LoggingMacroDispatcher::logInfo( WALBERLA__LOGGING__OSS.str(), __FILE__, __LINE__ );\
   }}\
}
#else
#define WALBERLA_LOG_INFO_ON_ROOT(msg) (void(0))
#endif

//////////////
// PROGRESS //
//////////////

#ifdef WALBERLA_LOGLEVEL_PROGRESS
#define WALBERLA_LOG_PROGRESS(msg){\
   WALBERLA_LOG_PROGRESS_SECTION(){\
      std::stringstream WALBERLA__LOGGING__OSS;\
      WALBERLA__LOGGING__OSS << msg;\
      walberla::logging::internal::LoggingMacroDispatcher::logProgress( WALBERLA__LOGGING__OSS.str(), __FILE__, __LINE__ );\
   }\
}
#else
#define WALBERLA_LOG_PROGRESS(msg) (void(0))
#endif
#ifdef WALBERLA_LOGLEVEL_PROGRESS
#define WALBERLA_LOG_PROGRESS_ON_ROOT(msg){\
   if( walberla::logging::Logging::instance()->onRoot() ){\
   WALBERLA_LOG_PROGRESS_SECTION(){\
      std::stringstream WALBERLA__LOGGING__OSS;\
      WALBERLA__LOGGING__OSS << msg;\
      walberla::logging::internal::LoggingMacroDispatcher::logProgress( WALBERLA__LOGGING__OSS.str(), __FILE__, __LINE__ );\
   }}\
}
#else
#define WALBERLA_LOG_PROGRESS_ON_ROOT(msg) (void(0))
#endif

////////////
// DETAIL //
////////////

#ifdef WALBERLA_LOGLEVEL_DETAIL
#define WALBERLA_LOG_DETAIL(msg){\
   WALBERLA_LOG_DETAIL_SECTION(){\
      std::stringstream WALBERLA__LOGGING__OSS;\
      WALBERLA__LOGGING__OSS << msg;\
      walberla::logging::internal::LoggingMacroDispatcher::logDetail( WALBERLA__LOGGING__OSS.str(), __FILE__, __LINE__ );\
   }\
}
#else
#define WALBERLA_LOG_DETAIL(msg) (void(0))
#endif
#ifdef WALBERLA_LOGLEVEL_DETAIL
#define WALBERLA_LOG_DETAIL_ON_ROOT(msg){\
   if( walberla::logging::Logging::instance()->onRoot() ){\
   WALBERLA_LOG_DETAIL_SECTION(){\
      std::stringstream WALBERLA__LOGGING__OSS;\
      WALBERLA__LOGGING__OSS << msg;\
      walberla::logging::internal::LoggingMacroDispatcher::logDetail( WALBERLA__LOGGING__OSS.str(), __FILE__, __LINE__ );\
   }}\
}
#else
#define WALBERLA_LOG_DETAIL_ON_ROOT(msg) (void(0))
#endif

/////////////
// TRACING //
/////////////

// --> Tracing.h (WALBERLA_TRACE_IN)



} // namespace logging
} // namespace walberla
