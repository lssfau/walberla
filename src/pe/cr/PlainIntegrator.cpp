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

PlainIntegrator::PlainIntegrator(const shared_ptr<BodyStorage>&      globalBodyStorage,
                                 const shared_ptr<BlockStorage>&     blockStorage,
                                 domain_decomposition::BlockDataID   storageID,
                                 WcTimingTree*                       tt)
   : globalBodyStorage_(globalBodyStorage)
   , blockStorage_(blockStorage)
   , storageID_(storageID)
   , tt_(tt)
{

}

void PlainIntegrator::timestep( const real_t dt )
{
   if (tt_!=NULL) tt_->start("PlainIntegrator");
   if (tt_!=NULL) tt_->start("Integrate Bodies");
   for (auto it = blockStorage_->begin(); it != blockStorage_->end(); ++it){
      IBlock & currentBlock = *it;
      Storage * storage = currentBlock.getData< Storage >( storageID_ );
      BodyStorage& localStorage = (*storage)[0];
      for (auto bd = localStorage.begin(); bd != localStorage.end(); ++bd){
         move(*bd, dt);
      }
   }
   if (tt_!=NULL) tt_->stop("Integrate Bodies");
   if (tt_!=NULL) tt_->stop("PlainIntegrator");
}

//=================================================================================================
//
//  SIMULATION FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Calculation of a time step of \a dt.
 *
 * \param dt Time step size.
 * \return void
 *
 * Calculating one single time step of size \a dt for the capsule. The global position, the
 * linear and angular velocity and the orientation of the capsule are changed depending on
 * the acting forces and the current velocities.
 */
void PlainIntegrator::move( BodyID id, real_t dt )
{
   // Checking the state of the body
   WALBERLA_ASSERT( id->checkInvariants(), "Invalid capsule state detected" );
   WALBERLA_ASSERT( !id->hasSuperBody(), "Invalid superordinate body detected" );

   // Resetting the contact node and removing all attached contacts
//      id->resetNode();
   id->clearContacts();

   // Moving the capsule according to the acting forces (don't move a sleeping capsule)
   if( id->isAwake() && !id->hasInfiniteMass() ) {
      // Calculating the linear acceleration by the equation
      //   force * m^(-1) + gravity
      const Vec3 vdot( id->getForce() * id->getInvMass() + getGlobalLinearAcceleration() );

      // Calculating the angular acceleration by the equation
      //   R * Iinv * R^T * torque
      const Vec3 wdot( id->getInvInertia() * id->getTorque() );

      // Updating the linear velocity
      id->setLinearVel(id->getLinearVel() + vdot * dt);

      // Updating the angular velocity
      id->setAngularVel( id->getAngularVel() + wdot * dt );

      // Calculating the translational displacement
      id->setPosition( id->getPosition() + id->getLinearVel() * dt );

      // Calculating the rotation angle
      const Vec3 phi( id->getAngularVel() * dt );

      // Calculating the new orientation
      id->rotate( Quat( phi, phi.length() ) );
      WALBERLA_ASSERT( realIsEqual( id->getRotation().getDeterminant(), real_c(1) ), "Corrupted rotation matrix determinant" );

      // Setting the axis-aligned bounding box
      id->calcBoundingBox();

      // Calculating the current motion of the capsule
      id->calcMotion();
   }

   // Resetting the acting forces
   id->resetForceAndTorque();

   // Checking the state of the capsule
   WALBERLA_ASSERT( id->checkInvariants(), "Invalid capsule state detected" );
}
//*************************************************************************************************

} // namespace cr
} // namespace pe
} // namespace walberla
