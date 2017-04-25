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
//! \file DEM.cpp
//! \ingroup pe
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

DEM::DEM(    const shared_ptr<BodyStorage>&    globalBodyStorage
           , const shared_ptr<BlockStorage>&   blockStorage
           , domain_decomposition::BlockDataID storageID
           , domain_decomposition::BlockDataID ccdID
           , domain_decomposition::BlockDataID fcdID
           , WcTimingTree*                     tt)
   : globalBodyStorage_(globalBodyStorage)
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

void DEM::timestep( real_t dt )
{
   maxPenetration_          = real_c(0.0);
   numberOfContacts_        = 0;
   numberOfContactsTreated_ = 0;

   for (auto it = blockStorage_->begin(); it != blockStorage_->end(); ++it){
      IBlock & currentBlock = *it;

      ccd::ICCD* ccd = currentBlock.getData< ccd::ICCD >( ccdID_ );
      fcd::IFCD* fcd = currentBlock.getData< fcd::IFCD >( fcdID_ );
      ccd->generatePossibleContacts( tt_ );
      if (tt_ != NULL) tt_->start("FCD");
      Contacts& cont = fcd->generateContacts( ccd->getPossibleContacts() );
      if (tt_ != NULL) tt_->stop("FCD");

      for (auto cIt = cont.begin(); cIt != cont.end(); ++cIt){
         const real_t overlap( -cIt->getDistance() );
         if( overlap > maxPenetration_ )
            maxPenetration_ = overlap;
         if (shouldContactBeTreated( &(*cIt), currentBlock.getAABB() ))
         {
            ++numberOfContactsTreated_;
            resolveContact( &(*cIt) );
         }
      }

      numberOfContacts_ += cont.size();

      cont.clear();
   }

//   if (numContacts > 0)
//      WALBERLA_LOG_DEVEL_ON_ROOT("#Contacts: " << numContacts << "." );

   if (tt_ != NULL) tt_->start("ForceSync");
   reduceForces( *blockStorage_, storageID_, *globalBodyStorage_);
   if (tt_ != NULL) tt_->stop("ForceSync");

   for (auto it = blockStorage_->begin(); it != blockStorage_->end(); ++it)
   {
      IBlock & currentBlock = *it;
      Storage * storage = currentBlock.getData< Storage >( storageID_ );
      BodyStorage& localStorage = (*storage)[0];
      BodyStorage& shadowStorage = (*storage)[1];

      // Updating the positions and velocities of all locally owned and global rigid bodies (but not shadow copies).
      WALBERLA_LOG_DETAIL( "Time integration starts...\n" );

      if (tt_ != NULL) tt_->start("Integration");

      for( auto bodyIt = localStorage.begin(); bodyIt != localStorage.end(); ++bodyIt )
      {
         WALBERLA_LOG_DETAIL( "Time integration of body with system id " << bodyIt->getSystemID());// << "\n" << *bodyIt );

         // Resetting the contact node and removing all attached contacts
      //      bodyIt->resetNode();
         bodyIt->clearContacts();

         move( *bodyIt, dt );

         // Resetting the acting forces
         bodyIt->resetForceAndTorque();
      }

      if (tt_ != NULL) tt_->stop("Integration");

      // Reset forces of shadow copies
      for( auto bodyIt = shadowStorage.begin(); bodyIt != shadowStorage.end(); ++bodyIt ) {
         bodyIt->clearContacts();
         bodyIt->resetForceAndTorque();
      }

   }
}

//=================================================================================================
//
//  SOLVER FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Resolves the given colliding contact.
 *
 * \param c The colliding contact.
 * \return The overlap of the contact.
 *
 * TODO
 */
void DEM::resolveContact( ContactID c ) const
{
   WALBERLA_LOG_DETAIL( "resolving contact: " << c->getID() );
   BodyID b1( c->getBody1()->getTopSuperBody() );
   BodyID b2( c->getBody2()->getTopSuperBody() );

   // Global position of contact
   const Vec3 gpos( c->getPosition() );

   // The absolute value of the penetration length
   real_t delta( -c->getDistance() );

   // Calculating the relative velocity in normal and tangential direction
   // The negative signs result from the different definition of the relative
   // normal velocity of the pe (see Contact::getType)
   const real_t relVelN( -c->getNormalRelVel() );
   const Vec3   relVel ( -c->getRelVel() );
   const Vec3   relVelT( relVel - ( relVelN * c->getNormal() ) );

   real_t fNabs( 0 );
   Vec3   fN;

   // Calculating the normal force based on the non-linear extended Hertz model
   // This force model is only applied in case of a sphere/sphere collision, since only
   // then an effective radius can be computed.
//   if( dem::forceModel == dem::hertz && c->hasEffectiveRadius() )
//   {
//      const real_t alpha   ( 1.5 );
//      const real_t beta    ( 0.5 );
//      const real_t Reff    ( c->getEffectiveRadius() );
//      const real_t Eeff    ( c->getEffectiveYoungModulus() );
//      const real_t k       ( ( real_c(4)/real_c(3) ) * Eeff * sqrt( Reff ) );

//      fNabs = k*std::pow( delta, alpha ) + c->getDampingN()*relVelN*std::pow( delta, beta );
//      if( fNabs < real_c(0) ) fNabs = real_c(0);
//      fN = fNabs * c->getNormal();
//   }

//   // Calculating the normal force based on a linear spring-dashpot force model
//   else
   {
      fNabs = getStiffness(c) * delta + getDampingN(c) * relVelN;
      if( fNabs < real_c(0) ) fNabs = real_c(0);
      fN = fNabs * c->getNormal();
   }

   // Calculating the tangential force based on the model by Haff and Werner
   const real_t fTabs( std::min( getDampingT(c) * relVelT.length(), getFriction(c) * fNabs ) );
   const Vec3   fT   ( fTabs * relVelT.getNormalizedOrZero() );

   // Add normal force at contact point
   b1->addForceAtPos(  fN, gpos );
   b2->addForceAtPos( -fN, gpos );

   // Add tangential force at contact point
   b1->addForceAtPos(  fT, gpos );
   b2->addForceAtPos( -fT, gpos );

   WALBERLA_LOG_DETAIL("nForce: " << fN << "\ntForce: " << fT);
}
//*************************************************************************************************

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
void DEM::move( BodyID id, real_t dt )
{
   // Checking the state of the body
   WALBERLA_ASSERT( id->checkInvariants(), "Invalid body state detected" );
   WALBERLA_ASSERT( !id->hasSuperBody(), "Invalid superordinate body detected" );

   // Moving the capsule according to the acting forces (don't move a sleeping capsule)
   if( id->isAwake() && !id->hasInfiniteMass() )
   {
      // Calculating the linear acceleration by the equation
      //   force * m^(-1) + gravity
      const Vec3 vdot( id->getForce() * id->getInvMass() + getGlobalLinearAcceleration() );

      // Calculating the angular acceleration by the equation
      //   R * Iinv * R^T * torque
      const Vec3 wdot( id->getInvInertia() * id->getTorque() );

      // Calculating the translational displacement
      id->setPosition( id->getPosition() + id->getLinearVel() * dt + 0.5 * vdot * dt * dt );

      // Calculating the rotation angle
      const Vec3 phi( id->getAngularVel() * dt + 0.5 * wdot * dt * dt);

      // Calculating the new orientation
      if (!floatIsEqual(phi.length(), 0))
         id->rotate( Quat( phi, phi.length() ) );
      WALBERLA_ASSERT_FLOAT_EQUAL( id->getRotation().getDeterminant(), real_c(1 ), "Corrupted rotation matrix determinant" );

      // Updating the linear velocity
      id->setLinearVel(id->getLinearVel() + vdot * dt);

      // Updating the angular velocity
      id->setAngularVel( id->getAngularVel() + wdot * dt );

      // Setting the axis-aligned bounding box
      id->calcBoundingBox();

      // Calculating the current motion of the capsule
      id->calcMotion();
   }

   // Resetting the acting forces
   id->resetForceAndTorque();

   // Checking the state of the rigid body
   WALBERLA_ASSERT( id->checkInvariants(), "Invalid capsule state detected" );
}
//*************************************************************************************************

}  // namespace cr
} // namespace pe
}  // namespace walberla
