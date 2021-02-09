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
//! \file PlainIntegrator.h
//! \author Klaus Iglberger
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//! \brief Source file for the Plain Integrator solver
//
//======================================================================================================================


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include "pe/cr/PlainIntegrator.h"
#include "pe/ccd/ICCD.h"
#include "pe/fcd/IFCD.h"
#include "pe/bg/IBG.h"
#include "pe/rigidbody/BodyStorage.h"
#include "pe/rigidbody/RigidBody.h"

#include "core/logging/all.h"
#include "core/timing/TimingTree.h"

namespace walberla {
namespace pe {
namespace cr {

template< typename Integrator >
PlainIntegratorSolver<Integrator>::PlainIntegratorSolver( const Integrator & integrate,
                                                         const shared_ptr<BodyStorage>&      globalBodyStorage,
                                                         const shared_ptr<BlockStorage>&     blockStorage,
                                                         domain_decomposition::BlockDataID   storageID,
                                                         WcTimingTree*                       tt)
   : integrate_(integrate)
   , globalBodyStorage_(globalBodyStorage)
   , blockStorage_(blockStorage)
   , storageID_(storageID)
   , tt_(tt)
{

}

template< typename Integrator >
void PlainIntegratorSolver<Integrator>::timestep( const real_t dt )
{
   if (tt_!=nullptr) tt_->start("PlainIntegrator");
   if (tt_!=nullptr) tt_->start("Integrate Bodies");
   for (auto it = blockStorage_->begin(); it != blockStorage_->end(); ++it){
      IBlock & currentBlock = *it;
      Storage * storage = currentBlock.getData< Storage >( storageID_ );
      BodyStorage& localStorage = (*storage)[0];
      for (auto bd = localStorage.begin(); bd != localStorage.end(); ++bd){
         // Checking the state of the body
         WALBERLA_ASSERT( bd->checkInvariants(), "Invalid body state detected" );
         WALBERLA_ASSERT( !bd->hasSuperBody(), "Invalid superordinate body detected" );

         // Moving the body according to the acting forces (don't move a sleeping body)
         if( bd->isAwake() && !bd->hasInfiniteMass() ) {
            integrate_( bd.getBodyID(), dt, *this );
         }

         // Resetting the acting forces
         bd->resetForceAndTorque();

         // Checking the state of the body
         WALBERLA_ASSERT( bd->checkInvariants(), "Invalid body state detected" );
      }
   }
   if (tt_!=nullptr) tt_->stop("Integrate Bodies");
   if (tt_!=nullptr) tt_->stop("PlainIntegrator");
}

} // namespace cr
} // namespace pe
} // namespace walberla
