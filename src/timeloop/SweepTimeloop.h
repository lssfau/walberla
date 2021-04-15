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
//! \file SweepTimeloop.h
//! \ingroup timeloop
//! \author Martin Bauer <martin.bauer@fau.de>
//! \brief Timeloop that runs a collection of sweeps.
//
//======================================================================================================================

#pragma once

#include "SelectableFunctionCreators.h"
#include "Timeloop.h"
#include "domain_decomposition/StructuredBlockStorage.h"

#include <functional>

#include <map>

namespace walberla {
namespace timeloop {


   //*******************************************************************************************************************
   /*!\brief Timeloop that runs a collection of sweeps.
    *
    * A Sweep is defined as an operation that processes a single block, i.e. a function "void sweepFunc ( IBlock * )"
    * A sweep is usually some numeric kernel, for example a LBM Stream-Collide step.
    * For all sweeps  Before- and AfterFunctions can be registered. This are void(void) functions that are run before
    * or after the sweep. Here one can register functions that have a close relation to the sweep, for example
    * communication steps that are necessary before the sweep can run.
    *
    * The SweepTimeloop supports the Selectable concept, so that the execution of sweeps can depend on the current
    * uid::globalState(), together with the current block state. This mechanism makes the registration of Sweeps somewhat
    * complex.
    *
    * \section sweepTimeloop_simple Simple Registration
    *
    * When registering sweeps (even without Selectable concept) the syntax may seem a little strange.
    * Lets do an example where we register simple C functions as sweeps:
      \code
      // The callback functions
      void mySweep          ( IBlock * ) { ... }
      void myBeforeFunction ()           { ... }
      void myAfterFunction  ()           { ... }

      SweepTimeloop timeloop;
      timeloop.add() << BeforeFunction( & myBeforeFunction, "MyBeforeFunction")
                     << Sweep         ( & mySweep,          "MySweep" )
                     << AfterFunction ( & myAfterFunction,  "MyAfterFunction );
     \endcode
    *
    * So we do not pass the parameters directly to add(), but add() returns an object,
    * where the different functions can be piped in.
    * The BeforeFunction, Sweep and AfterFunction take additional selector arguments, so
    * multiple of them can be piped into one sweep, as long a unique sweep can be determined using the selectors.
    *
    * Instead of registering plain C functions, it is usually a good idea to use functors as sweeps,
    * since they can have additional data. Be aware that the state of such functor sweeps is the same for
    * all blocks. To store data per block register it as BlockData or use SweepOnBlock:
    *
    * \section sweepTimeloop_onBlock Sweeps on Block
    *
    * Sometimes a sweep functor class needs to store data per block, not per process ( which is the default ).
    * The solution in this case would be, to put this data into a separate data structure and register that
    * as BlockData. The alternative is to register the sweep class itself as BlockData, whereas the actual
    * sweep function (which is the same for all blocks ) then only consists of fetching the sweep class and executing it.
    * This second solution is used by SweepOnBlock:
    *
    * The starting point is a class, that has members which should be placed block local:
      \code
      class MyLBMClass {
         public:
            MyLBMClass( real_t omega, PdfField * pdfField );
         private:
            Type someBlockLocalData;
      };
      \endcode
    * The second step is to write a so called creator class, which acts as the interface between the waLBerla
    * block data mechanism and your data. It basically is the initialization function for the blockdata:
      \code
      struct MyLBMClassCreator {
         MyLBMClassCreator( BlockDataID pdfFieldID ): pdfFieldID_(pdfFieldID) {}

         MyLBMClass * operator() ( IBlock * const block ) {
            return new MyLBMClass( block->getData<PdfField>( pdfFieldID_ ) );
         }

         private:
           BlockDataID pdfFieldID_
      }
      \endcode
    *
    * This is then registered like this:
    * \code
      timeloop.add() << SweepOnBlock<MyLBMClass>( MyLBMClassCreator(pdfFieldID) );
      \endcode
    *
    * \ingroup timeloop
    */
   //*******************************************************************************************************************
   class SweepTimeloop : public Timeloop
   {
   public:

      //****************************************************************************************************************
      /*!\name Constructor & Destructor */
      //@{

      SweepTimeloop( BlockStorage & blockStorage, uint_t nrOfTimeSteps )
         : Timeloop(nrOfTimeSteps), blockStorage_(blockStorage), nextId_(0),firstRun_(true)
      {}

      SweepTimeloop( const shared_ptr<StructuredBlockStorage> & structuredBlockStorage, uint_t nrOfTimeSteps )
         : Timeloop(nrOfTimeSteps), blockStorage_( structuredBlockStorage->getBlockStorage() ),
           nextId_(0), firstRun_(true)
      {}

      ~SweepTimeloop() override
      {
         for ( auto i = sweeps_.begin(); i != sweeps_.end(); ++i )
            delete i->second;
      }

      //@}
      //****************************************************************************************************************



      //****************************************************************************************************************
      /*!\name Registering of Sweeps and Before/After Functions */
      //@{

      SweepAdder & add() {
         ++nextId_;
         sweeps_[nextId_] = new SweepAdder ( blockStorage_,nextId_ );
         return *(  sweeps_[nextId_] );
      }

      void remove( SweepAdder & sweep ) {
         sweepsToDelete_.push_back( sweep.id_ );
      }

      //@}
      //****************************************************************************************************************

   protected:
      BlockStorage & blockStorage_;

      void removeForDeletionMarkedSweeps()
      {
         for( auto it = sweepsToDelete_.begin(); it != sweepsToDelete_.end(); ++it )
            sweeps_.erase( *it );
      }

      void doTimeStep(const Set<SUID> &selectors) override;
      void doTimeStep(const Set<SUID> &selectors, WcTimingPool &tp) override;

      uint_t nextId_;
      std::vector<uint_t> sweepsToDelete_;
      std::map< uint_t,SweepAdder* > sweeps_;

      bool firstRun_; ///< required to register timer in doTimeStep( selectors, timingPool)
   };


} // namespace timeloop
} // namespace walberla



//======================================================================================================================
//
//  EXPORT
//
//======================================================================================================================

namespace walberla {
   using timeloop::SweepTimeloop;

   using timeloop::Sweep;
   using timeloop::SweepOnBlock;
   using timeloop::BeforeFunction;
   using timeloop::AfterFunction;
}

