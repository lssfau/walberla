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
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//
//======================================================================================================================

#include "SweepTimeloop.h"
#include "core/Abort.h"


namespace walberla {
namespace timeloop {


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////   Execution of Timeloop  ////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template < typename TP >
void SweepTimeloop<TP>::doTimeStep(const Set<SUID> &selectors)
{
   removeForDeletionMarkedSweeps();
   //iterate over all registered sweeps
   for( auto sweepIt = sweeps_.begin(); sweepIt != sweeps_.end(); ++sweepIt )
   {
      SweepAdder & s = * ( sweepIt->second );

      //select and execute before functions
      for( size_t j=0; j < s.beforeFuncs.size(); ++j )
         this->executeSelectable(s.beforeFuncs[j].selectableFunc_, selectors, "Pre-Sweep Function");

      // Loop over all blocks
      for( BlockStorage::iterator bi = blockStorage_.begin(); bi != blockStorage_.end(); ++bi )
      {
         // ensure that at least one sweep has been registered (regardless of its selector)
         if( s.sweep.empty() )
         {
            WALBERLA_ABORT("Selecting Sweep " << sweepIt->first << ": " <<
                           "No sweep has been registered! Did you only register a BeforeFunction or AfterFunction?" )
         }

         // ensure that exactly one sweep has been registered that matches the specified selectors
         size_t const numSweeps = s.sweep.getNumberOfMatching(selectors + bi->getState());

         if (numSweeps == size_t(0)) {
            continue;
         } else {
            if (numSweeps > size_t(1)) {
               WALBERLA_ABORT("Only one sweep must be added to a single SweepAdder object. This error might be caused "
                              "by e.g. \"timeloop.add() << Sweep(A) << Sweep(B);\".")
            }
         }

         Sweep * selectedSweep = s.sweep.getUnique( selectors + bi->getState() );

         WALBERLA_LOG_PROGRESS_SECTION()
         {
            std::string sweepName;
            s.sweep.getUnique( selectors + bi->getState(), sweepName );
            WALBERLA_LOG_PROGRESS("Running sweep \"" << sweepName << "\" on block " << bi->getId() )
         }

         (selectedSweep->function_)( bi.get() );
      }

      // select and execute after functions
      for( size_t j=0; j < s.afterFuncs.size(); ++j )
         this->executeSelectable(s.afterFuncs[j].selectableFunc_, selectors, "Post-Sweep Function");
   }
}

template < typename TP >
void SweepTimeloop<TP>::doTimeStep(const Set<SUID> &selectors, timing::TimingPool<TP> &timing)
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
         this->executeSelectable( s.beforeFuncs[j].selectableFunc_, selectors, "Pre-Sweep Function", timing );

      for( BlockStorage::iterator bi = blockStorage_.begin(); bi != blockStorage_.end(); ++bi )
      {
         // ensure that at least one sweep has been registered (regardless of its selector)
         if( s.sweep.empty() )
         {
            WALBERLA_ABORT("Selecting Sweep " << sweepIt->first << ": " <<
                           "No sweep has been registered! Did you only register a BeforeFunction or AfterFunction?" )
         }

         // ensure that exactly one sweep has been registered that matches the specified selectors
         size_t const numSweeps = s.sweep.getNumberOfMatching(selectors + bi->getState());

         if (numSweeps == size_t(0)) {
            continue;
         } else {
            if (numSweeps > size_t(1)) {
               WALBERLA_ABORT("Only one sweep must be added to a single SweepAdder object. This error might be caused "
                              "by e.g. \"timeloop.add() << Sweep(A) << Sweep(B);\".")
            }
         }

         std::string sweepName;
         Sweep * selectedSweep = s.sweep.getUnique( selectors + bi->getState(), sweepName );

         WALBERLA_LOG_PROGRESS("Running sweep \"" << sweepName << "\" on block " << bi->getId() )

         // loop over all blocks
         timing[sweepName].start();
         (selectedSweep->function_)( bi.get() );
         timing[sweepName].end();
      }

      // select and execute after functions
      for( size_t j=0; j < s.afterFuncs.size(); ++j )
         this->executeSelectable(s.afterFuncs[j].selectableFunc_,selectors,"Post-Sweep Function", timing );
   }
}



} // namespace timeloop
} // namespace walberla

