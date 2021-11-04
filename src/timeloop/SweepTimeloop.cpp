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
//! \file SweepTimeloop.cpp
//! \ingroup timeloop
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#include "SweepTimeloop.h"
#include "core/Abort.h"


namespace walberla {
namespace timeloop {


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////   Execution of Timeloop  ////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void SweepTimeloop::doTimeStep(const Set<SUID> &selectors)
{
   removeForDeletionMarkedSweeps();
   //iterate over all registered sweeps
   for( auto sweepIt = sweeps_.begin(); sweepIt != sweeps_.end(); ++sweepIt )
   {
      SweepAdder & s = * ( sweepIt->second );

      //select and execute before functions
      for( size_t j=0; j < s.beforeFuncs.size(); ++j )
         executeSelectable(s.beforeFuncs[j].selectableFunc_,selectors,"Pre-Sweep Function");

      // Loop over all blocks
      for( BlockStorage::iterator bi = blockStorage_.begin(); bi != blockStorage_.end(); ++bi )
      {
         if( s.sweep.empty() )
         {
            WALBERLA_ABORT("Selecting Sweep " << sweepIt->first << ": " <<
                           "No sweep has been registered! Did you only register a BeforeFunction or AfterFunction?" );
         }

         Sweep selectedSweep;
         size_t numSweeps = s.sweep.get(selectedSweep, selectors + bi->getState());

         // ensure that no more than one sweep has been added to a single SweepAdder object
         if (numSweeps == size_t(0)) {
            continue;
         } else {
            if (numSweeps > size_t(1)) {
               WALBERLA_ABORT("Only one sweep must be added to a single SweepAdder object. This error might be caused "
                              "by e.g. \"timeloop.add() << Sweep(A) << Sweep(B);\".")
            }
         }

         WALBERLA_LOG_PROGRESS_SECTION()
         {
            std::string sweepName;
            s.sweep.getUnique( selectors + bi->getState(), sweepName );
            WALBERLA_LOG_PROGRESS("Running sweep \"" << sweepName << "\" on block " << bi->getId() );
         }

         (selectedSweep.function_)( bi.get() );
      }

      // select and execute after functions
      for( size_t j=0; j < s.afterFuncs.size(); ++j )
         executeSelectable(s.afterFuncs[j].selectableFunc_,selectors,"Post-Sweep Function");
   }
}

void SweepTimeloop::doTimeStep(const Set<SUID> &selectors, WcTimingPool &timing)
{
   removeForDeletionMarkedSweeps();
   // On first run we extract all possible names of sweeps, independent of selectors
   // to register timer for them. Since all timingPools have to have the same entries on
   // all processes otherwise a timingPool.reduce() does not work
   // see test SweepTimeloopTimerReduction.cpp
   if ( firstRun_ || timing.empty() )
   {

      for( auto sweepIt = sweeps_.begin(); sweepIt != sweeps_.end(); ++sweepIt )
      {
         SweepAdder & s = * ( sweepIt->second );
         // loop over all possibilities in selectable object
         for( auto it = s.sweep.begin(); it != s.sweep.end(); ++it )
            timing.registerTimer( it.identifier() );
      }
      firstRun_ = false;
   }


   //iterate over all registered sweeps
   for( auto sweepIt = sweeps_.begin(); sweepIt != sweeps_.end(); ++sweepIt )
   {
      SweepAdder & s = * ( sweepIt->second );

      //select and execute before functions
      for( size_t j=0; j < s.beforeFuncs.size(); ++j )
         executeSelectable( s.beforeFuncs[j].selectableFunc_, selectors, "Pre-Sweep Function", timing );

      for( BlockStorage::iterator bi = blockStorage_.begin(); bi != blockStorage_.end(); ++bi )
      {
         Sweep selectedSweep;
         std::string sweepName;
         size_t numSweeps = s.sweep.get(selectedSweep, sweepName, selectors + bi->getState());

         // ensure that no more than one sweep has been added to a single SweepAdder object
         if (numSweeps == size_t(0)) {
            continue;
         } else {
            if (numSweeps > size_t(1)) {
               WALBERLA_ABORT("Only one sweep must be added to a single SweepAdder object. This error might be caused "
                              "by e.g. \"timeloop.add() << Sweep(A) << Sweep(B);\".")
            }
         }

         WALBERLA_LOG_PROGRESS("Running sweep \"" << sweepName << "\" on block " << bi->getId() );

         // loop over all blocks
         timing[sweepName].start();
         (selectedSweep.function_)( bi.get() );
         timing[sweepName].end();
      }

      // select and execute after functions
      for( size_t j=0; j < s.afterFuncs.size(); ++j )
         executeSelectable(s.afterFuncs[j].selectableFunc_,selectors,"Post-Sweep Function", timing );
   }
}



} // namespace timeloop
} // namespace walberla

