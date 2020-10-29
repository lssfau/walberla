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
//! \file SubCyclingManager.cpp
//! \author Lukas Werner <lks.werner@fau.de>
//
//======================================================================================================================


#include "SubCyclingManager.h"

#include <utility>

namespace walberla {
namespace lbm_mesapd_coupling {

SubCyclingManager::SubCyclingManager(size_t numberOfSubCycles, shared_ptr<WcTimingPool> timingPool)
   : numberOfSubCycles_(numberOfSubCycles), timingPool_(std::move(timingPool)), currentTimeStep_(0) {}

void SubCyclingManager::operator()() {
   executeBeforeFunctions();
   for (size_t s = 0; s < numberOfSubCycles_; ++s) {
      executeDuringFunctions();
   }
   executeAfterFunctions();
   ++currentTimeStep_;
}

SubCyclingManager::FuncHandle
SubCyclingManager::addFuncBeforeSubCycles(const VoidVoidFunc &f, const std::string &identifier) {
   IdentifiedFunc idFunc(identifier, f);
   beforeFunctions_.emplace_back(idFunc);
   return beforeFunctions_.size() - 1;
}

SubCyclingManager::FuncHandle
SubCyclingManager::addFuncDuringSubCycles(const VoidVoidFunc &f, const std::string &identifier) {
   IdentifiedFunc idFunc(identifier, f);
   duringFunctions_.emplace_back(idFunc);
   return duringFunctions_.size() - 1;
}

SubCyclingManager::FuncHandle
SubCyclingManager::addFuncAfterSubCycles(const VoidVoidFunc &f, const std::string &identifier) {
   IdentifiedFunc idFunc(identifier, f);
   afterFunctions_.emplace_back(idFunc);
   return afterFunctions_.size() - 1;
}

inline void SubCyclingManager::executeBeforeFunctions() {
   executeFunctions(beforeFunctions_);
}

inline void SubCyclingManager::executeDuringFunctions() {
   executeFunctions(duringFunctions_);
}

inline void SubCyclingManager::executeAfterFunctions() {
   executeFunctions(afterFunctions_);
}

inline void SubCyclingManager::executeFunctions(std::vector<IdentifiedFunc>& functions) {
   for (const auto &f: functions) {
      startTiming(f.first);
      f.second();
      stopTiming(f.first);
   }
}

inline void SubCyclingManager::startTiming(const std::string & name){
   if (timingPool_ != nullptr) {
      (*timingPool_)[name].start();
   }
}

inline void SubCyclingManager::stopTiming(const std::string & name){
   if (timingPool_ != nullptr) {
      (*timingPool_)[name].end();
   }
}

}
}