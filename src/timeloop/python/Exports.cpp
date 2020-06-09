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
//! \file Exports.cpp
//! \ingroup timeloop
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

// Do not reorder includes - the include order is important
#include "python_coupling/PythonWrapper.h"

#ifdef WALBERLA_BUILD_WITH_PYTHON

#include "python_coupling/Manager.h"
#include "python_coupling/helper/ModuleScope.h"

#include "timeloop/Timeloop.h"
#include "timeloop/SweepTimeloop.h"

using namespace boost::python;


namespace walberla {
namespace timeloop {

#ifdef WALBERLA_CXX_COMPILER_IS_GNU
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"
#endif
struct ITimeloopWrap : public ITimeloop, public wrapper<ITimeloop>
{
   void run() override                          {        this->get_override( "run" )();                     }
   void singleStep() override                   {        this->get_override( "singleStep" )();              }
   void stop() override                         {        this->get_override( "stop" )();                    }
   void synchronizedStop( bool s ) override     {        this->get_override( "synchronizedStop" )(s);       }
   void setCurrentTimeStep( uint_t ts) override {        this->get_override( "setCurrentTimeStep" )(ts);    }
   uint_t getCurrentTimeStep() const override   { return this->get_override( "getCurrentTimeStep" )();      }
   uint_t getNrOfTimeSteps()   const override   { return this->get_override( "getNrOfTimeSteps" )();        }
};

#ifdef WALBERLA_CXX_COMPILER_IS_GNU
#pragma GCC diagnostic pop
#endif

void exportModuleToPython()
{
   python_coupling::ModuleScope timeloopModule( "timeloop" );

   void ( Timeloop::*p_run1) ( const bool  ) = &Timeloop::run;
   void ( Timeloop::*p_run2) ( WcTimingPool & , const bool ) = &Timeloop::run;

   class_<ITimeloopWrap, boost::noncopyable > ("ITimeloop"  )
      .def( "getCurrentTimeStep", pure_virtual( &ITimeloop::getCurrentTimeStep  ) )
      .def( "getNrOfTimeSteps",   pure_virtual( &ITimeloop::getNrOfTimeSteps    ) )
      .def( "stop",               pure_virtual( &ITimeloop::stop                ) )
      .def( "synchronizedStop",   pure_virtual( &ITimeloop::synchronizedStop    ) )
      .def( "run",                pure_virtual( &ITimeloop::run                 ) )
      ;

   class_<Timeloop, bases<ITimeloop>, boost::noncopyable>( "CppTimeloop", no_init )
      .def( "run",                p_run1, args("logTimeStep") = true )
      .def( "run",                p_run2, ( args("timingPool"), args("logTimeStep") = true ) )
      ;
   class_< SweepTimeloop , bases<Timeloop>,  boost::noncopyable > ( "CppSweepTimeloop", no_init )
      .def( init< const shared_ptr<StructuredBlockStorage> & , uint_t  > () )
      ;

}


} // namespace timeloop
} // namespace walberla


#endif //WALBERLA_BUILD_WITH_PYTHON

