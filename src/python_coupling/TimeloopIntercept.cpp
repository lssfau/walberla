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
//! \file TimeloopIntercept.cpp
//! \ingroup python_coupling
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#include "TimeloopIntercept.h"
#include "DictWrapper.h"
#include "PythonCallback.h"
#include "Shell.h"

#include "core/logging/Logging.h"

#ifdef WALBERLA_BUILD_WITH_PYTHON



#if defined(_MSC_VER)

#else
#  include <csignal>
#endif


static bool signalTriggered = false;


namespace walberla {
namespace python_coupling {


#if defined(_MSC_VER)

void enableSignalHandler() {

}

#else

void customSignalHandler( int  ) {
   signalTriggered = true;
}

void enableSignalHandler() {
   std::signal (SIGUSR1, customSignalHandler );
}

#endif



TimeloopIntercept::TimeloopIntercept( const std::string & functionName, uint_t interval, bool enableSignalInterrupt )
         : callback_ ( walberla::make_shared<PythonCallback> ( functionName ) ),
           timestep_ ( 0 ),
           interval_ ( interval ),
           enableSignalInterrupt_ ( enableSignalInterrupt )
{
   if ( enableSignalInterrupt )
      enableSignalHandler();
}

TimeloopIntercept::TimeloopIntercept( const std::string & pythonFile, const std::string & functionName,
                                      uint_t interval, bool enableSignalInterrupt )
         : callback_ ( walberla::make_shared<PythonCallback> ( pythonFile, functionName ) ),
           timestep_ ( 0 ),
           interval_ ( interval ),
           enableSignalInterrupt_ ( enableSignalInterrupt )
{
   if ( enableSignalInterrupt )
      enableSignalHandler();
}


void TimeloopIntercept::operator() ()
{
   ++timestep_;
   if( interval_ > 0 &&  ( timestep_ % interval_ == 0 ) )
      (*callback_)();

   if ( enableSignalInterrupt_ && signalTriggered )
   {
      WALBERLA_LOG_INFO( "Interrupting Simulation at timestep " << timestep_ );
      signalTriggered = false;
      Shell shell;
      shell.data() = this->callback()->data();
      shell.run();
      WALBERLA_LOG_INFO( "Continue Simulation at timestep " << timestep_+1  << " ..." );
   }
}



} // namespace python_coupling
} // namespace walberla

#else

namespace walberla {
namespace python_coupling {

TimeloopIntercept::TimeloopIntercept( const std::string & functionName, uint_t interval, bool enableSignalInterrupt )
         : callback_ ( walberla::make_shared<PythonCallback> ( functionName ) ),
           timestep_ ( 0 ),
           interval_ ( interval ),
           enableSignalInterrupt_ ( enableSignalInterrupt )
{}

TimeloopIntercept::TimeloopIntercept( const std::string & pythonFile, const std::string & functionName,
                                      uint_t interval, bool enableSignalInterrupt )
         : callback_ ( walberla::make_shared<PythonCallback> ( pythonFile, functionName ) ),
           timestep_ ( 0 ),
           interval_ ( interval ),
           enableSignalInterrupt_ ( enableSignalInterrupt )
{}

void TimeloopIntercept::operator() (){}

} // namespace python_coupling
} // namespace walberla

#endif
