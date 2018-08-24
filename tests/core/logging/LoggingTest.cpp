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
//! \file LoggingTest.cpp
//! \ingroup core
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#include "core/Abort.h"
#include "core/Environment.h"
#include "core/debug/TestSubsystem.h"
#include "core/logging/Logging.h"
#include "core/logging/Tracing.h"


using namespace walberla;


class MyStamp : public logging::Logging::CustomStamp
{
public:
   MyStamp() : step_(uint_c(0)) {}
   ~MyStamp() override = default;
   void step( uint_t s ) { step_ = s; }
   std::string stamp() override
   {
      std::ostringstream oss;
      oss << "[" << step_ << "]";
      return oss.str();
   }
   uint_t maxStampWidth() override { return 5; }
private:
   uint_t step_;
};


int fak( int n )
{
   WALBERLA_TRACE_IN;
   WALBERLA_LOG_DEVEL( "Processing: " << n );
   if( n == 0 )
      return 1;
   return n * fak( n - 1 );
}


int main( int argc, char** argv )
{
   debug::enterTestMode();

   Environment env( argc, argv );

   logging::Logging::instance()->includeLoggingToFile( "myLog.txt" );
   //logging::Logging::instance()->showTimeStamp( false );

   auto myStamp = make_shared< MyStamp >();
   logging::Logging::instance()->addCustomStamp( myStamp );
   myStamp->step(0);

   WALBERLA_LOG_DEVEL( "This is still in development!" );
   WALBERLA_LOG_RESULT( "This is a result!" );
   WALBERLA_LOG_WARNING( "This is a warning!" );
   WALBERLA_LOG_INFO( "This is just some information." );
   WALBERLA_LOG_PROGRESS( "I'm still here!" );
   WALBERLA_LOG_DETAIL( "These are some interesting details!" );

   logging::Logging::instance()->setLogLevel( logging::Logging::PROGRESS );
   myStamp->step(1);

   WALBERLA_LOG_DEVEL( "This is still in development!" );
   WALBERLA_LOG_RESULT( "This is a result!" );
   WALBERLA_LOG_WARNING( "This is a warning!" );
   WALBERLA_LOG_INFO( "This is just some information." );
   WALBERLA_LOG_PROGRESS( "I'm still here!" );
   WALBERLA_LOG_DETAIL( "These are some interesting details!" );

   logging::Logging::instance()->setLogLevel( logging::Logging::DETAIL );
   myStamp->step(2);

   WALBERLA_LOG_DEVEL( "This is still in development!" );
   WALBERLA_LOG_RESULT( "This is a result!" );
   WALBERLA_LOG_WARNING( "This is a warning!" );
   WALBERLA_LOG_INFO( "This is just some information." );
   WALBERLA_LOG_PROGRESS( "I'm still here!" );
   WALBERLA_LOG_DETAIL( "These are some interesting details!" );

   fak(4);

   //std::cin.get();

   logging::Logging::instance()->setFileLogLevel( logging::Logging::TRACING );
   myStamp->step(30);

   WALBERLA_LOG_DEVEL( "This is still in development!" );
   WALBERLA_LOG_RESULT( "This is a result!" );
   WALBERLA_LOG_WARNING( "This is a warning!" );
   WALBERLA_LOG_INFO( "This is just some information." );
   WALBERLA_LOG_PROGRESS( "I'm still here!" );
   WALBERLA_LOG_DETAIL( "These are some interesting details!" );

   fak(4);

   logging::Logging::instance()->setLogLevel( logging::Logging::WARNING );
   myStamp->step(40);

   WALBERLA_LOG_DEVEL( "This is still in development!" );
   WALBERLA_LOG_RESULT( "This is a result!" );
   WALBERLA_LOG_WARNING( "This is a warning!" );
   WALBERLA_LOG_INFO( "This is just some information." );
   WALBERLA_LOG_PROGRESS( "I'm still here!" );
   WALBERLA_LOG_DETAIL( "These are some interesting details!" );

   fak(4);

   myStamp->step(100);

   logging::Logging::instance()->setLogLevel( logging::Logging::INFO );
   bool runtimeErrorThrown = false;
   try
   {
      WALBERLA_LOG_INFO("The following Error is expected!");
      Abort::instance()->resetAbortFunction( &Abort::exceptionAbort );
      WALBERLA_CHECK(false); // This check is expected to fail
      // Should never be reached
      Abort::instance()->resetAbortFunction();
      WALBERLA_CHECK(false);
   }
   catch( const std::runtime_error & /*e*/ )
   {
      Abort::instance()->resetAbortFunction();
      runtimeErrorThrown = true;
   }
   WALBERLA_CHECK( runtimeErrorThrown );

   //runtimeErrorThrown = false;
   //try { WALBERLA_CHECK ( false ); }
   //catch( const std::runtime_error & /*e*/ ) { runtimeErrorThrown = true; }
   //WALBERLA_CHECK( runtimeErrorThrown );

   //WALBERLA_CHECK_EQUAL( 5, 23 );

   //WALBERLA_ASSERT_EQUAL( 5, 23, "I think 5 should be equal to 23!" );

   //WALBERLA_ASSERT_SECTION( 5 == 23 )
   //{
   //   WALBERLA_LOG_DEVEL( "Here we are ..." );
   //}

   //WALBERLA_ABORT( "Fatal error!" );

   return EXIT_SUCCESS;
}
