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
//! \file Shell.impl.h
//! \ingroup python_coupling
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

// Do not reorder includes - the include order is important

#include "PythonWrapper.h"

#ifdef WALBERLA_BUILD_WITH_PYTHON

#include "Shell.h"
#include "DictWrapper.h"
#include "Manager.h"

#include "core/logging/Logging.h"
#include "core/StringUtility.h"


namespace walberla {
namespace python_coupling {


   Shell::Shell( const std::string & prompt )
      : exposedVars_( new DictWrapper() )
   {
      Manager::instance()->triggerInitialization();

      using namespace boost::python;

      try {
         object main_module  = import("__main__");
         dict globals = extract<dict>( main_module.attr( "__dict__" ) );

         exec("import os \n"
              "import code \n"
              "try:\n"
              "   import readline \n"
              "   import rlcompleter \n"
              "   readline.parse_and_bind(\"tab: complete\") \n"
              "   histfile = os.path.join( os.path.expanduser('~'), '.waLBerla_history') \n"
              "   try:\n"
              "      readline.read_history_file(histfile)\n"
              "   except IOError:\n"
              "      pass\n"
              "except ImportError:\n"
              "   pass\n"
              "\n", globals );


         prompt1_ = prompt;
         for ( uint_t i=0; i < prompt.size(); ++i )
            prompt2_ += " ";

         prompt1_ += "> ";
         prompt2_ += "> ";

         exec("import sys", globals );
         exec( std::string( "sys.ps1 = '" + prompt1_ + "'" ).c_str(), globals );
         exec( std::string( "sys.ps2 = '" + prompt2_ + "'" ).c_str(), globals );
      }
      catch ( error_already_set & ) {
         PyErr_Print();
         WALBERLA_ABORT( "Error initializing Shell" );
      }
   }

   Shell::~Shell()
   {
      using namespace boost::python;

      object main_module  = import("__main__");
      dict globals = extract<dict>( main_module.attr( "__dict__" ) );
      exec("readline.write_history_file(histfile)", globals );
   }


   bool Shell::isCompleteCommand ( const std::string & code )
   {
      using namespace boost::python;
      object main_module  = import("__main__");
      dict globals = extract<dict>( main_module.attr( "__dict__" ) );

      boost::python::dict locals;
      locals["codeToTest"] = code;
      object compiledCode = eval( str( "code.compile_command(codeToTest)" ), globals, locals );
      return ( compiledCode != object() );
   }

   bool Shell::getCompleteCommand( std::string & result )
   {
      uint_t lineCounter = 0;

      while ( true )
      {
         char* line;

         if ( lineCounter == 0 )
            line = PyOS_Readline( stdin, stdout, (char*)prompt1_.c_str() );
         else
            line = PyOS_Readline( stdin, stdout, (char*)prompt2_.c_str() );

         if ( line == nullptr || *line == '\0' ) {  // interrupt or EOF
            result.clear();
            return false;
         }

         std::string strLine ( line );
         std::string strTrimmedLine( line );
         string_trim( strTrimmedLine );

         PyMem_Free( line );

         lineCounter++;
         result += strLine;

         bool commandComplete = isCompleteCommand( result );

         if ( lineCounter == 1 && commandComplete )
            return true;
         if ( strTrimmedLine.empty() && isCompleteCommand(result) ) // multiline commands have to end with empty line
            return true;
      }

      return false;
   }


   void Shell::operator() ()
   {
      using namespace boost::python;

      object main_module  = import("__main__");
      dict globals = extract<dict>( main_module.attr( "__dict__" ) );

      globals.update( exposedVars_->dict() );
      exec( "from waLBerla import *", globals );


      const int  MAX_LINE_LENGTH   = 1024;
      const char continueMarker    = 'c';
      const char stopMarker        = 's';

      WALBERLA_ROOT_SECTION()
      {
         std::string code;

         while ( true )
         {
            code.clear();

            try
            {
               if ( ! getCompleteCommand(code ) ) {
                  std::cout << "\n";
                  code.resize( MAX_LINE_LENGTH );
                  code[0] = stopMarker;
                  MPI_Bcast( (void*) code.c_str(), MAX_LINE_LENGTH, // send stop command
                             MPI_CHAR, 0, MPI_COMM_WORLD );
                  break;
               }
               else
               {
                  std::string codeToSend = continueMarker + code;

                  if ( codeToSend.size() >= uint_c( MAX_LINE_LENGTH ) )
                  {
                     WALBERLA_LOG_WARNING("Line length too big, only allowed " << MAX_LINE_LENGTH-1 << " characters" );
                     continue;
                  }
                  codeToSend.resize( MAX_LINE_LENGTH );
                  MPI_Bcast( (void*) codeToSend.c_str(), MAX_LINE_LENGTH, // send code snippet to other processes
                             MPI_CHAR, 0, MPI_COMM_WORLD );

                  PyRun_SimpleString( code.c_str() );
                  fflush( stderr );
                  WALBERLA_MPI_BARRIER();
               }
            }
            catch( boost::python::error_already_set & ) {
               PyErr_Print();
            }
         }
      }
      else
      {
         char * buffer = new char[MAX_LINE_LENGTH];

         while ( true)
         {

            MPI_Bcast( (void*) buffer, MAX_LINE_LENGTH, // send code snippet to other processes
                       MPI_CHAR, 0, MPI_COMM_WORLD );

            if ( *buffer == stopMarker )
               break;

            std::string code ( buffer+1 );

            PyRun_SimpleString( code.c_str() );
            fflush( stderr );
            WALBERLA_MPI_BARRIER();
         }


         delete [] buffer;
      }
   }


} // namespace python_coupling
} // namespace walberla


#else


#include "Shell.h"

namespace walberla {
namespace python_coupling {

   Shell::Shell( const std::string & )  {}
   Shell::~Shell() = default;
   void Shell::operator()() {}

} // namespace python_coupling
} // namespace walberla


#endif
