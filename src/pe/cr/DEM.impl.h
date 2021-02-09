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
//! \file DEM.impl.h
//! \author Klaus Iglberger
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//! \brief Source file for the DEM solver
//
//======================================================================================================================


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include "pe/cr/DEM.h"
#include "pe/ccd/ICCD.h"
#include "pe/fcd/IFCD.h"
#include "pe/bg/IBG.h"
#include "pe/rigidbody/BodyStorage.h"
#include "pe/rigidbody/RigidBody.h"
#include "pe/contact/Contact.h"
#include "pe/contact/ContactFunctions.h"
#include "pe/synchronization/SyncForces.h"

#include "core/logging/all.h"

namespace walberla {
namespace pe {
namespace cr {

template< typename Integrator, typename ContactResolver >
DEMSolver<Integrator,ContactResolver>::DEMSolver(
             const Integrator & integrate, const ContactResolver & resolveContact
           , const shared_ptr<BodyStorage>&    globalBodyStorage
           , const shared_ptr<BlockStorage>&   blockStorage
           , domain_decomposition::BlockDataID storageID
           , domain_decomposition::BlockDataID ccdID
           , domain_decomposition::BlockDataID fcdID
           , WcTimingTree*                     tt)
   : integrate_(integrate), resolveContact_(resolveContact)
   ,globalBodyStorage_(globalBodyStorage)
   , blockStorage_(blockStorage)
   , storageID_(storageID)
   , ccdID_(ccdID)
   , fcdID_(fcdID)
   , tt_(tt)
   , maxPenetration_(0)
   , numberOfContacts_(0)
   , numberOfContactsTreated_(0)
{

}

template< typename Integrator, typename ContactResolver >
void DEMSolver<Integrator,ContactResolver>::timestep( real_t dt )
{
   maxPenetration_          = real_c(0.0);
   numberOfContacts_        = 0;
   numberOfContactsTreated_ = 0;

   for (auto it = blockStorage_->begin(); it != blockStorage_->end(); ++it){
      IBlock & currentBlock = *it;

      ccd::ICCD* ccd = currentBlock.getData< ccd::ICCD >( ccdID_ );
      fcd::IFCD* fcd = currentBlock.getData< fcd::IFCD >( fcdID_ );
      ccd->generatePossibleContacts( tt_ );
      if (tt_ != nullptr) tt_->start("FCD");
      Contacts& cont = fcd->generateContacts( ccd->getPossibleContacts() );
      if (tt_ != nullptr) tt_->stop("FCD");

      for (auto cIt = cont.begin(); cIt != cont.end(); ++cIt){
         const real_t overlap( -cIt->getDistance() );
         if( overlap > maxPenetration_ )
            maxPenetration_ = overlap;
         if (shouldContactBeTreated( &(*cIt), currentBlock.getAABB() ))
         {
            ++numberOfContactsTreated_;
            resolveContact_( &(*cIt), dt);
         }
      }

      numberOfContacts_ += cont.size();

      cont.clear();
   }

//   if (numContacts > 0)
//      WALBERLA_LOG_DEVEL_ON_ROOT("#Contacts: " << numContacts << "." );

   if (tt_ != nullptr) tt_->start("ForceSync");
   reduceForces( *blockStorage_, storageID_, *globalBodyStorage_);
   if (tt_ != nullptr) tt_->stop("ForceSync");

   for (auto it = blockStorage_->begin(); it != blockStorage_->end(); ++it)
   {
      IBlock & currentBlock = *it;
      Storage * storage = currentBlock.getData< Storage >( storageID_ );
      BodyStorage& localStorage = (*storage)[0];
      BodyStorage& shadowStorage = (*storage)[1];

      // Updating the positions and velocities of all locally owned and global rigid bodies (but not shadow copies).
      WALBERLA_LOG_DETAIL( "Time integration starts...\n" );

      if (tt_ != nullptr) tt_->start("Integration");

      for( auto bodyIt = localStorage.begin(); bodyIt != localStorage.end(); ++bodyIt )
      {
         WALBERLA_LOG_DETAIL( "Time integration of body with system id " << bodyIt->getSystemID());// << "\n" << *bodyIt );

         // Checking the state of the body
         WALBERLA_ASSERT( bodyIt->checkInvariants(), "Invalid body state detected" );
         WALBERLA_ASSERT( !bodyIt->hasSuperBody(), "Invalid superordinate body detected" );
         
         // Moving the body according to the acting forces (don't move a sleeping body)
         if( bodyIt->isAwake() && !bodyIt->hasInfiniteMass() )
         {
            integrate_( bodyIt.getBodyID(), dt, *this );
         }
         
         // Resetting the acting forces
         bodyIt->resetForceAndTorque();
         
         // Checking the state of the rigid body
         WALBERLA_ASSERT( bodyIt->checkInvariants(), "Invalid body state detected" );

         // Resetting the acting forces
         bodyIt->resetForceAndTorque();
      }

      if (tt_ != nullptr) tt_->stop("Integration");

      // Reset forces of shadow copies
      for( auto& body : shadowStorage ) {
         body.resetForceAndTorque();
      }

   }
}

}  // namespace cr
} // namespace pe
}  // namespace walberla
