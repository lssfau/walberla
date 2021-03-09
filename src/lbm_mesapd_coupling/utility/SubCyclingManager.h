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
//! \file SubCyclingManager.h
//! \author Lukas Werner <lks.werner@fau.de>
//
//======================================================================================================================

#include "core/DataTypes.h"
#include "core/logging/Logging.h"
#include "core/timing/TimingPool.h"

namespace walberla {
namespace lbm_mesapd_coupling {

/*! \brief Handles the execution of subcycles in a timeloop to allow for finer grained time steps than the LBM ones.
*
* Supports registration of functions that are run before, during or after a subcycle.
* The SubCyclingManager itself has to be added to a parent (LBM) time loop, and will execute the functions registered
* for execution during the subcycling procedure numberOfSubCycles times.
*/

class SubCyclingManager {
public:
   using VoidVoidFunc = std::function<void ()>;

   /*! \name Construction & Destruction */
   //@{
   explicit SubCyclingManager(size_t numberOfSubCycles, shared_ptr<WcTimingPool> timingPool = nullptr);

   virtual ~SubCyclingManager() = default;
   //@}

   /*! \name Registration Functions */
   //@{
   using FuncHandle = size_t;

   FuncHandle addFuncBeforeSubCycles(const VoidVoidFunc &f, const std::string &identifier = "Other");
   FuncHandle addFuncDuringSubCycles(const VoidVoidFunc &f, const std::string &identifier = "Other");
   FuncHandle addFuncAfterSubCycles(const VoidVoidFunc &f, const std::string &identifier = "Other");
   //@}

   /*! \name Execution Functions */
   //@{
   void operator()();
   //@}

   /*! \name Bookkeeping Functions */
   //@{
   uint_t getCurrentTimeStep() const        { return currentTimeStep_;   }
   void setCurrentTimeStep(uint_t timestep) { currentTimeStep_ = timestep; }
   //@}

private:
   using IdentifiedFunc = std::pair<std::string, VoidVoidFunc>;

   inline void executeBeforeFunctions();
   inline void executeDuringFunctions();
   inline void executeAfterFunctions();

   inline void executeFunctions(std::vector<IdentifiedFunc>& functions);

   inline void startTiming(const std::string & name);
   inline void stopTiming(const std::string & name);

   size_t numberOfSubCycles_;
   shared_ptr<WcTimingPool> timingPool_;
   uint_t currentTimeStep_;

   std::vector<IdentifiedFunc> beforeFunctions_;
   std::vector<IdentifiedFunc> duringFunctions_;
   std::vector<IdentifiedFunc> afterFunctions_;
};

}
}