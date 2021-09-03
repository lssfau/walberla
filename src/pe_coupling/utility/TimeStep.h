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
//! \file TimeStep.h
//! \ingroup pe_coupling
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/timing/TimingTree.h"
#include "domain_decomposition/StructuredBlockStorage.h"
#include "pe/cr/ICR.h"

#include <functional>

namespace walberla {
namespace pe_coupling {

/*!\brief Carries out the the PE time steps, including sub iteration functionality.
 *
 * It executes \a numberOfSubIterations PE steps within one timestep of size \a timeStepSize.
 *
 * These PE sub iterations require, that the current external (e.g. hydrodynamic, gravitational, ...) forces and torques
 * acting on each particle remains unchanged. Thus, a map is set up internally that stores and re-sets these forces
 * and torques in each PE sub iteration.
 *
 * Additionally, a function \a forceEvaluationFunc can be given that allows to evaluate different forces before the PE
 * step is carried out. An example are particle-particle lubrication forces that have to be updated in each sub iteration.
 *
 */
class TimeStep
{
public:

   explicit TimeStep( const shared_ptr<StructuredBlockStorage> & blockStorage,
                      const BlockDataID & bodyStorageID,
                      pe::cr::ICR & collisionResponse,
                      const std::function<void (void)> & synchronizeFunc,
                      const real_t timeStepSize = real_t(1),
                      const uint_t numberOfSubIterations = uint_t(1),
                      const std::function<void (void)> & forceEvaluationFunc = [](){})
         : timeStepSize_( timeStepSize )
         , numberOfSubIterations_( ( numberOfSubIterations == 0 ) ? uint_t(1) : numberOfSubIterations )
         , blockStorage_( blockStorage )
         , bodyStorageID_( bodyStorageID )
         , collisionResponse_( collisionResponse )
         , synchronizeFunc_( synchronizeFunc )
         , forceEvaluationFunc_( forceEvaluationFunc )
   {}

   void operator()();

private:

   const real_t timeStepSize_;
   const uint_t numberOfSubIterations_;

   shared_ptr<StructuredBlockStorage> blockStorage_;
   const BlockDataID &  bodyStorageID_;

   pe::cr::ICR & collisionResponse_;
   std::function<void (void)> synchronizeFunc_;
   std::function<void (void)> forceEvaluationFunc_;

}; // class TimeStep



} // namespace pe_coupling
} // namespace walberla
