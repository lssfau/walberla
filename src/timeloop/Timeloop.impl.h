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
//! \file Timeloop.cpp
//! \ingroup timeloop
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#include "Timeloop.h"
#include "core/Abort.h"
#include "core/debug/Debug.h"
#include "core/uid/GlobalState.h"
#include "core/mpi/Reduce.h"
#include "core/perf_analysis/extern/likwid.h"

namespace walberla {
namespace timeloop {

template < typename TP >
void Timeloop<TP>::run( const bool logTimeStep )
{
   WALBERLA_LOG_PROGRESS( "Running timeloop for " << nrOfTimeSteps_ << " time steps" )
   while(curTimeStep_ < nrOfTimeSteps_) {
      singleStep( logTimeStep );
      if ( stop_ ) {
         stop_ = false;
         break;
      }
   }
   WALBERLA_LOG_PROGRESS( "Timeloop finished" )
}

template < typename TP >
void Timeloop<TP>::run(timing::TimingPool<TP> & tp, const bool logTimeStep )
{
   WALBERLA_LOG_PROGRESS( "Running timeloop for " << nrOfTimeSteps_ << " time steps" )

   while(curTimeStep_ < nrOfTimeSteps_) {
      singleStep( tp, logTimeStep );
      if ( stop_ ) {
         stop_ = false;
         break;
      }
   }

   WALBERLA_LOG_PROGRESS( "Timeloop finished" )
}

//*******************************************************************************************************************
/*! Stops the timeloop, has to be called on every process
*
*  While the timeloop is running and calling function or sweeps, a
*  called function can call this stop() function then the run() loop is stopped
*  before reaching nrOfTimeSteps
*/
//*******************************************************************************************************************

template < typename TP >
void Timeloop<TP>::stop()
{
   stop_ = true;
}


//*******************************************************************************************************************
/*! Similar to stop - useful if the stop signal is only known on a single process
*
* Typical scenario: One process gathers information about the simulation state to determine a
* stopping criterion ( convergence test ). If this process decides to stop this information has
* to be communicated, such that all processes stop - this is done using synchronizedStop()
*     -> If at least on process calls synchronizedStop(true) the timeloop is stopped
*/
//*******************************************************************************************************************

template < typename TP >
void Timeloop<TP>::synchronizedStop( bool stopVal )
{
   stop_ = stopVal;
   mpi::allReduceInplace( stop_, mpi::LOGICAL_OR );
}

template < typename TP >
void Timeloop<TP>::singleStep( const bool logTimeStep )
{
   LoggingStampManager const raii( make_shared<LoggingStamp>( *this ), logTimeStep );

   WALBERLA_LOG_PROGRESS( "Running time step " << curTimeStep_ )

   for(size_t i=0; i<beforeFunctions_.size(); ++i )
      executeSelectable( beforeFunctions_[i], uid::globalState(), "Pre-Timestep Function" );

   doTimeStep(uid::globalState());

   for(size_t i=0; i<afterFunctions_.size(); ++i )
      executeSelectable( afterFunctions_[i], uid::globalState(),"Post-Timestep Function" );

   ++curTimeStep_;
}

template < typename TP >
void Timeloop<TP>::singleStep( timing::TimingPool<TP> & tp, const bool logTimeStep )
{
   LoggingStampManager const raii( make_shared<LoggingStamp>( *this ), logTimeStep );

   WALBERLA_LOG_PROGRESS( "Running time step " << curTimeStep_ )

   for(size_t i=0; i<beforeFunctions_.size(); ++i )
      executeSelectable( beforeFunctions_[i], uid::globalState(), "Pre-Timestep Function", tp );

   doTimeStep(uid::globalState(), tp);

   for(size_t i=0; i<afterFunctions_.size(); ++i )
      executeSelectable( afterFunctions_[i], uid::globalState(),"Post-Timestep Function", tp );

   ++curTimeStep_;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////   Registering Functions   ///////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template < typename TP >
typename Timeloop<TP>::FctHandle
Timeloop<TP>::addFuncBeforeTimeStep(const VoidFctNoArguments& f, const std::string & id,
                                const Set<SUID> & r, const Set<SUID> & e )
{
    beforeFunctions_.emplace_back(f,r,e,id );
    return beforeFunctions_.size() - 1;
}

template < typename TP >
void Timeloop<TP>::addFuncBeforeTimeStep(const Timeloop::FctHandle & h,
                                     const VoidFctNoArguments& f, const std::string & id,
                                     const Set<SUID>&r, const Set<SUID> & e )
{
   WALBERLA_ASSERT_LESS( h, beforeFunctions_.size() ) //invalid FctHandle
   beforeFunctions_[h].add(f,r,e,id);
}


template < typename TP >
typename Timeloop<TP>::FctHandle
Timeloop<TP>::addFuncAfterTimeStep(const VoidFctNoArguments& f, const std::string & id,
                                      const Set<SUID> & r, const Set<SUID> & e )
{
    afterFunctions_.emplace_back(f,r,e,id );
    return afterFunctions_.size() - 1;
}

template < typename TP >
void Timeloop<TP>::addFuncAfterTimeStep(const Timeloop::FctHandle & h,
                                           const VoidFctNoArguments& f, const std::string & id,
                                           const Set<SUID>&r, const Set<SUID> & e )
{
   WALBERLA_ASSERT_LESS( h, afterFunctions_.size() ) //invalid FctHandle
   afterFunctions_[h].add(f,r,e,id);
}



template < typename TP >
void Timeloop<TP>::executeSelectable( const selectable::SetSelectableObject<VoidFctNoArguments,SUID> & selectable,
                                  const Set<SUID> & selector,
                                  const std::string & what )
{
   std::string objectName;
   const VoidFctNoArguments * exe = selectable.getUnique( selector,objectName );
   if( exe == nullptr )
      WALBERLA_ABORT( "Trying to selecting " << what << ": "
                      << "Multiple Matches found! Check your selector " << selector << std::endl
                      << "All registered objects: " << std::endl << selectable << std::endl )


   WALBERLA_LOG_PROGRESS("Running " << what << " \"" << objectName << "\"" )

   LIKWID_MARKER_START( objectName.c_str() );
   (*exe)();
   LIKWID_MARKER_STOP( objectName.c_str() );
}

template < typename TP >
void Timeloop<TP>::executeSelectable( const selectable::SetSelectableObject<VoidFctNoArguments,SUID> & selectable,
                                  const Set<SUID> & selector,
                                  const std::string & what,
                                  timing::TimingPool<TP> & timing )
{
   std::string objectName;
   const VoidFctNoArguments * exe = selectable.getUnique( selector, objectName );

   if( !exe)
      WALBERLA_ABORT( "Trying to select " << what << ": "
                      << "Multiple or no matches found! Check your selector " << selector << std::endl
                      << "All registered objects: " << std::endl << selectable << std::endl )

   WALBERLA_LOG_PROGRESS("Running " << what << " \"" << objectName << "\"" )

   timing[objectName].start();
   LIKWID_MARKER_START( objectName.c_str() );
   (*exe)();
   LIKWID_MARKER_STOP( objectName.c_str() );
   timing[objectName].end();
}



} // namespace timeloop
} // namespace walberla
