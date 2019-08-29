#//======================================================================================================================
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
//! \file HCSITS.impl.h
//! \author Tobias Preclik
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//! \brief Source file for the hard contact solver
//
//======================================================================================================================

//*************************************************************************************************
// Includes
//*************************************************************************************************

#include "HCSITS.h"

#include "pe/bg/IBG.h"
#include "pe/ccd/ICCD.h"
#include "pe/fcd/IFCD.h"
#include "pe/contact/ContactFunctions.h"

#include "core/math/Constants.h"
#include "core/math/Limits.h"
#include "core/math/Shims.h"
#include "core/math/Utility.h"

#include "core/ConcatIterator.h"


namespace walberla {
namespace pe {
namespace cr {

inline void HardContactSemiImplicitTimesteppingSolvers::BodyCache::resize(const size_t n)
{
   v_.resize(n);
   w_.resize(n);
   dv_.resize(n);
   dw_.resize(n);
}

inline void HardContactSemiImplicitTimesteppingSolvers::ContactCache::resize(const size_t n)
{
   r1_.resize(n);
   r2_.resize(n);
   body1_.resize(n);
   body2_.resize(n);
   n_.resize(n);
   t_.resize(n);
   o_.resize(n);
   dist_.resize(n);
   mu_.resize(n);
   diag_nto_.resize(n);
   diag_nto_inv_.resize(n);
   diag_to_inv_.resize(n);
   diag_n_inv_.resize(n);
   p_.resize(n);
}

//=================================================================================================
//
//  CONSTRUCTOR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Constructor of the CollisionSystem class.
 */
inline HardContactSemiImplicitTimesteppingSolvers::HardContactSemiImplicitTimesteppingSolvers(
      const shared_ptr<BodyStorage>&    globalBodyStorage,
      const shared_ptr<BlockStorage>&   blockStorage,
      domain_decomposition::BlockDataID storageID,
      domain_decomposition::BlockDataID ccdID,
      domain_decomposition::BlockDataID fcdID,
      WcTimingTree*                     tt)
   : globalBodyStorage_(globalBodyStorage)
   , blockStorage_(blockStorage)
   , storageID_(storageID)
   , ccdID_(ccdID)
   , fcdID_(fcdID)
   , tt_(tt)
   , erp_              ( real_c(0.8) )
   , maxIterations_    ( 10 )
   , iteration_        ( 0 )
   , maxSubIterations_ ( 20 )
   , abortThreshold_   ( real_c(1e-7) )
   , relaxationModel_  ( InelasticFrictionlessContact )
   , overRelaxationParam_( real_c(1.0) )
   , relaxationParam_  ( real_c(0.75) )
   , maximumPenetration_ ( real_c(0.0) )
   , numContacts_      ( 0 )
   , numContactsTreated_( 0)
   , speedLimiterActive_( false )
   , speedLimitFactor_ ( real_c(1.0) )
   , requireSync_      ( false )
{
   // Logging the successful setup of the collision system
   WALBERLA_LOG_DETAIL( "Successfully initialized the collision system instance");
}
//*************************************************************************************************




//=================================================================================================
//
//  DESTRUCTOR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Destructor of the CollisionSystem class.
 */
inline HardContactSemiImplicitTimesteppingSolvers::~HardContactSemiImplicitTimesteppingSolvers()
{
   // Logging the successful destruction of the collision system
   WALBERLA_LOG_DETAIL( "Successfully destroyed the collision system instance");
}
//*************************************************************************************************




//=================================================================================================
//
//  SIMULATION FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Runs the simulation for one timestep.
 *
 * \param dt Size of the time step.
 * \return void
 *
 */
inline void HardContactSemiImplicitTimesteppingSolvers::timestep( const real_t dt )
{
   //WALBERLA_LOG_DETAIL( "New timestep!" );
   WALBERLA_ASSERT( !requireSync_, "Simulation requires synchronization before continuing." );

   const real_t dtinv( real_c(1) / dt );

   numContacts_        = 0;
   numContactsTreated_ = 0;
   maximumPenetration_ = 0;

   if (tt_ != NULL) tt_->start("Simulation Step");

   for (auto it = blockStorage_->begin(); it != blockStorage_->end(); ++it)
   {
      IBlock& currentBlock = *it;

      BodyCache&    bodyCache    = blockToBodyCache_[currentBlock.getId().getID()];
      ContactCache& contactCache = blockToContactCache_[currentBlock.getId().getID()];

      ccd::ICCD* ccd = currentBlock.getData< ccd::ICCD >( ccdID_ );
      fcd::IFCD* fcd = currentBlock.getData< fcd::IFCD >( fcdID_ );

      Storage * storage = currentBlock.getData< Storage >( storageID_ );
      BodyStorage& localStorage = (*storage)[0];
      BodyStorage& shadowStorage = (*storage)[1];

      // Detect all collisions
      WALBERLA_LOG_DETAIL( "Detecting contacts...");

      if (tt_ != NULL) tt_->start("Collision Detection");
      ccd->generatePossibleContacts( tt_ );
      if (tt_ != NULL) tt_->start("Fine");
      Contacts& contacts = fcd->generateContacts( ccd->getPossibleContacts() );
      if (tt_ != NULL) tt_->stop("Fine");

      if (tt_ != NULL) tt_->start("Filtering");
      // Filter out contacts
      size_t numContacts( contacts.size() );
      size_t numContactsMasked( 0 );

      std::vector<bool>   contactsMask_( numContacts );
      for( size_t i = 0; i < numContacts; ++i )
      {
         const ContactID c( &contacts[i] );
         if (shouldContactBeTreated(c, currentBlock.getAABB()))
         {
            contactsMask_[i] = true;
            ++numContactsMasked;
         } else
         {
            contactsMask_[i] = false;
         }
      }
      numContacts_        += numContacts;
      numContactsTreated_ += numContactsMasked;

      //      WALBERLA_LOG_DEVEL("contact filtering: " << numContactsMasked << "/" << contacts.size());
      if (tt_ != NULL) tt_->stop("Filtering");
      if (tt_ != NULL) tt_->stop("Collision Detection");

      if (tt_ != NULL) tt_->start("Collision Response Contact Caching");
      // Cache contact properties
      contactCache.resize( numContactsMasked );

      {

         size_t j = 0;
         for( size_t i = 0; i < numContacts; ++i )
         {
            if( !contactsMask_[i] )
               continue;

            const ContactID c( &contacts[i] );
            BodyID b1( c->getBody1()->getTopSuperBody() );
            BodyID b2( c->getBody2()->getTopSuperBody() );

            // the contact position is "in the middle"; r1 and r2 point to the same point but relative to body 1 and 2
            contactCache.body1_[j]    = b1;
            contactCache.body2_[j]    = b2;
            contactCache.r1_[j]       = c->getPosition() - b1->getPosition();
            contactCache.r2_[j]       = c->getPosition() - b2->getPosition();
            contactCache.n_[j]        = c->getNormal();   // points from body 2 towards body 1 and is normalized
            // construct vector perpendicular to the normal (cross product with cardinal basis vector where the 1 component is where the other vector has its component of smallest magnitude)
            if( std::fabs( contactCache.n_[j][0] ) < std::fabs( contactCache.n_[j][1] ) ) {
               if( std::fabs( contactCache.n_[j][0] ) < std::fabs( contactCache.n_[j][2] ) )
                  contactCache.t_[j] = Vec3( 0, contactCache.n_[j][2], -contactCache.n_[j][1] );   // = n x (1, 0, 0)^T
               else
                  contactCache.t_[j] = Vec3( contactCache.n_[j][1], -contactCache.n_[j][0], 0 );   // = n x (0, 0, 1)^T
            }
            else {
               if( std::fabs( contactCache.n_[j][1] ) < std::fabs( contactCache.n_[j][2] ) )
                  contactCache.t_[j] = Vec3( -contactCache.n_[j][2], 0, contactCache.n_[j][0] );   // = n x (0, 1, 0)^T
               else
                  contactCache.t_[j] = Vec3( contactCache.n_[j][1], -contactCache.n_[j][0], 0 );   // = n x (0, 0, 1)^T
            }
            normalize( contactCache.t_[j] );
            contactCache.o_[j]        = contactCache.n_[j] % contactCache.t_[j];
            Mat3 contactframe( contactCache.n_[j], contactCache.t_[j], contactCache.o_[j] );

            contactCache.dist_[j]  = c->getDistance();
            // If the distance is negative then penetration is present. This is an error and should be corrected.
            // Correcting the whole error is not recommended since due to the linearization the errors cannot
            // completely fixed anyway and the error reduction will introduce artificial restitution.
            // However, if the distance is positive then it is not about error correction but the distance that
            // can still be overcome without penetration and thus no error correction parameter should be applied.
            if( contactCache.dist_[j] < 0 ) {
               maximumPenetration_ = std::max( maximumPenetration_, -contactCache.dist_[j] );
               contactCache.dist_[j] *= erp_;
            }

            contactCache.mu_[j]       = getFriction(c);

            Mat3 diag    = -( math::skewSymCrossProduct(contactCache.r1_[j], math::skewSymCrossProduct(b1->getInvInertia(), contactCache.r1_[j]))
                              + math::skewSymCrossProduct(contactCache.r2_[j], math::skewSymCrossProduct(b2->getInvInertia(), contactCache.r2_[j])));
            diag[0]     += b1->getInvMass() + b2->getInvMass();
            diag[4]     += b1->getInvMass() + b2->getInvMass();
            diag[8]     += b1->getInvMass() + b2->getInvMass();
            diag         = contactframe.getTranspose() * diag * contactframe;

            // Diagonal block is know to be positive-definite and thus inverse always exists.
            contactCache.diag_nto_[j]      = diag;
            contactCache.diag_nto_inv_[j]  = diag.getInverse();
            contactCache.diag_n_inv_[j]    = math::inv(diag[0]);
            contactCache.diag_to_inv_[j]   = Mat2( diag[4], diag[5], diag[7], diag[8] ).getInverse();
            contactCache.p_[j] = Vec3();

            ++j;
         }
      }

      if (tt_ != NULL) tt_->stop("Collision Response Contact Caching");
      if (tt_ != NULL) tt_->start("Collision Response Body Caching");

      // Cache body properties (and time integrate v and w to the end of the time step by applying external forces, torques and accelerations)
      size_t numBodies( globalBodyStorage_->size() + localStorage.size() + shadowStorage.size() );
      bodyCache.resize( numBodies );

      size_t j = 0;
      // GLOBAL BODIES HAVE TO HAVE THE SAME INDEX WITHIN ALL BLOCKS!!!!
      for( auto body = globalBodyStorage_->begin(); body != globalBodyStorage_->end(); ++body, ++j ) {
         body->wake(); // BUGFIX: Force awaking of all bodies!
         body->index_ = j;
         WALBERLA_CHECK( body->hasInfiniteMass(), "Global bodies need to have infinite mass as they are not communicated!" );

         initializeVelocityCorrections( body.getBodyID(), bodyCache.dv_[j], bodyCache.dw_[j], dt ); // use applied external forces to calculate starting velocity

         bodyCache.v_[j] = body->getLinearVel();
         bodyCache.w_[j] = body->getAngularVel();
      }

      for( auto body = localStorage.begin(); body != localStorage.end(); ++body, ++j ) {
         body->wake(); // BUGFIX: Force awaking of all bodies!
         body->index_ = j;

         initializeVelocityCorrections( body.getBodyID(), bodyCache.dv_[j], bodyCache.dw_[j], dt ); // use applied external forces to calculate starting velocity

         if( body->isAwake() && !body->hasInfiniteMass() ) {
            bodyCache.v_[j] = body->getLinearVel() + getGlobalLinearAcceleration() * dt;
            bodyCache.w_[j] = body->getAngularVel() + dt * ( body->getInvInertia() * ( ( body->getInertia() * body->getAngularVel() ) % body->getAngularVel() ) );
         }
         else {
            bodyCache.v_[j] = body->getLinearVel();
            bodyCache.w_[j] = body->getAngularVel();
         }
      }

      for( auto body = shadowStorage.begin(); body != shadowStorage.end(); ++body, ++j )
      {
         body->wake(); // BUGFIX: Force awaking of all bodies!
         body->index_ = j;

         initializeVelocityCorrections( body.getBodyID(), bodyCache.dv_[j], bodyCache.dw_[j], dt );

         // Velocities of shadow copies will be initialized by velocity synchronization.
#ifndef NDEBUG
         bodyCache.v_[j] = Vec3( std::numeric_limits<real_t>::quiet_NaN(),
                                 std::numeric_limits<real_t>::quiet_NaN(),
                                 std::numeric_limits<real_t>::quiet_NaN() );
         bodyCache.w_[j] = Vec3( std::numeric_limits<real_t>::quiet_NaN(),
                                 std::numeric_limits<real_t>::quiet_NaN(),
                                 std::numeric_limits<real_t>::quiet_NaN() );
#endif
      }

      if (tt_ != NULL) tt_->stop("Collision Response Body Caching");
   }

   if (blockStorage_->size() == 0)
   {
      // create artificial block to handle global bodies even on processes where there are no blocks
      BodyCache&    bodyCache    = blockToBodyCache_[0];
      ContactCache& contactCache = blockToContactCache_[0];

      contactCache.resize( 0 );
      bodyCache.resize( globalBodyStorage_->size() );

      size_t j = 0;
      // GLOBAL BODIES HAVE TO HAVE THE SAME INDEX WITHIN ALL BLOCKS!!!!
      for( auto body = globalBodyStorage_->begin(); body != globalBodyStorage_->end(); ++body, ++j ) {
         body->wake(); // BUGFIX: Force awaking of all bodies!
         body->index_ = j;
         WALBERLA_CHECK( body->hasInfiniteMass(), "Global bodies need to have infinite mass as they are not communicated!" );

         initializeVelocityCorrections( body.getBodyID(), bodyCache.dv_[j], bodyCache.dw_[j], dt ); // use applied external forces to calculate starting velocity

         bodyCache.v_[j] = body->getLinearVel();
         bodyCache.w_[j] = body->getAngularVel();
      }
   }

   if (tt_ != NULL) tt_->start("Collision Response Resolution");
   const real_t rp = relaxationParam_;
   relaxationParam_ = real_c(1); // must be set to 1.0 such that dv and dw caused by external forces and torques are not falsely altered
   synchronizeVelocities( );
   relaxationParam_ = rp;

   // Iterate relaxation a constant number of times
   for( size_t it = 0; it < maxIterations_; ++it )
   {
      WALBERLA_LOG_DETAIL( "Iteration #" << it << "");

      real_t delta_max( 0 );

      iteration_ = it;

      if (tt_ != NULL) tt_->start("Collision Response Solving");
      for (auto blkIt = blockStorage_->begin(); blkIt != blockStorage_->end(); ++blkIt)
      {
         IBlock& currentBlock = *blkIt;

         BodyCache&    bodyCache    = blockToBodyCache_[currentBlock.getId().getID()];
         ContactCache& contactCache = blockToContactCache_[currentBlock.getId().getID()];

         switch( relaxationModel_ ) {
         case InelasticFrictionlessContact:
            delta_max = std::max( delta_max, relaxInelasticFrictionlessContacts( dtinv, contactCache, bodyCache ));
            break;

         case ApproximateInelasticCoulombContactByDecoupling:
            delta_max = relaxApproximateInelasticCoulombContactsByDecoupling( dtinv, contactCache, bodyCache );
            break;

         case ApproximateInelasticCoulombContactByOrthogonalProjections:
            delta_max = relaxInelasticCoulombContactsByOrthogonalProjections( dtinv, true, contactCache, bodyCache );
            break;

         case InelasticCoulombContactByDecoupling:
            delta_max = relaxInelasticCoulombContactsByDecoupling( dtinv, contactCache, bodyCache );
            break;

         case InelasticCoulombContactByOrthogonalProjections:
            delta_max = relaxInelasticCoulombContactsByOrthogonalProjections( dtinv, false, contactCache, bodyCache );
            break;

         case InelasticGeneralizedMaximumDissipationContact:
            delta_max = relaxInelasticGeneralizedMaximumDissipationContacts( dtinv, contactCache, bodyCache );
            break;

         default:
            throw std::runtime_error( "Unsupported relaxation model." );
         }
      }

      if (tt_ != NULL) tt_->stop("Collision Response Solving");

      synchronizeVelocities( );

      // Compute maximum impulse variation.
      // TODO:
      // - velocity variation would be better.
      // - do not reduce in every iteration



#if 0
      WALBERLA_MPI_SECTION()
      {
         if( tt_ != NULL ) tt_->start( "delta max reduction" );
         mpi::allReduceInplace(delta_max, mpi::MAX);
         if( tt_ != NULL ) tt_->stop( "delta max reduction" );
      }

      if( delta_max < abortThreshold_ ) {
         ++it;
         break;
      }
#else
      WALBERLA_UNUSED( delta_max );
#endif
   }

   if (tt_ != NULL) tt_->stop("Collision Response Resolution");
   if (tt_ != NULL) tt_->start("Collision Response Integration");

   // WARNING: Even though bodyCache.dv_[j] and bodyCache.dw_[j] _should_ be exactly 0 at all times for
   // bodies with infinite mass/inertia, this is not the case if the simulation breaks.
   // Infinite contact impulses can cause velocity corrections, which are not a number. If
   // added to the velocities of these bodies their positions on different
   // processes can get out of sync (some of them can get NaNs). To prevent these NaNs and
   // thus out of sync bodies, velocity corrections for bodies of infinite
   // mass/inertia are silently ignored.;
   {
      // choose arbitrary bodyChache. SHOULD BE ALL IN SYNC!
      BodyCache&    bodyCache    = blockToBodyCache_.begin()->second;
      size_t j = 0;
      for( auto body = globalBodyStorage_->begin(); body != globalBodyStorage_->end(); ++body, ++j )
      {
         WALBERLA_CHECK_EQUAL(j, body->index_);
         WALBERLA_LOG_DETAIL( "Integrating position of global body " << *body << " with velocity " << body->getLinearVel() << "" );
         if (body->hasInfiniteMass())
         {
            integratePositions( body.getBodyID(), bodyCache.v_[j], bodyCache.w_[j], dt );
         } else
         {
            integratePositions( body.getBodyID(), bodyCache.v_[j] + bodyCache.dv_[j], bodyCache.w_[j] + bodyCache.dw_[j], dt );
         }
         WALBERLA_LOG_DETAIL( "Result:\n" << *body << "");
      }
   }
   // Apply cached body properties (velocities) and time integrate positions
   for (auto it = blockStorage_->begin(); it != blockStorage_->end(); ++it)
   {
      IBlock& currentBlock = *it;

      BodyCache&    bodyCache    = blockToBodyCache_[currentBlock.getId().getID()];

      Storage * storage = currentBlock.getData< Storage >( storageID_ );
      BodyStorage& localStorage = (*storage)[0];
      BodyStorage& shadowStorage = (*storage)[1];

      for( auto body = ConcatIterator<BodyStorage::iterator>
           (localStorage.begin(),
            localStorage.end(),
            shadowStorage.begin(),
            shadowStorage.end()); body != ConcatIterator<BodyStorage::iterator>(); ++body )
      {
         if (!body->isCommunicating())
         {
            WALBERLA_CHECK(body->hasInfiniteMass(), "It is assumed that local non communicating bodies have infinite mass!");
            continue;
         }

         size_t j = body->index_;
         WALBERLA_LOG_DETAIL( "Integrating position of body with infinite mass " << *body << " with velocity " << bodyCache.v_[j] << "" );
         if( body->hasInfiniteMass() )
         {
            integratePositions( &(*body), bodyCache.v_[j], bodyCache.w_[j], dt );
         } else
         {
            integratePositions( &(*body), bodyCache.v_[j] + bodyCache.dv_[j], bodyCache.w_[j] + bodyCache.dw_[j], dt );
         }
         WALBERLA_LOG_DETAIL( "Result:\n" << *body << "" );
      }

      // NOTE: We might still need shadow copy updates if the user sets velocities or positions. Thus we might have to split up synchronize() again. It doesn't break anything if we still communicate the shadow copy updates though.
   }

   if (tt_ != NULL) tt_->stop("Collision Response Integration");

   blockToBodyCache_.clear();
   blockToContactCache_.clear();

   if (tt_ != NULL) tt_->stop("Simulation Step");
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Relaxes all contacts once. The contact model is for inelastic unilateral contacts without friction.
 *
 * \return The largest variation of contact impulses in the L-infinity norm.
 *
 * This function is to be called from resolveContacts(). Separating contacts are preferred over
 * persisting solutions if valid.
 */
inline real_t HardContactSemiImplicitTimesteppingSolvers::relaxInelasticFrictionlessContacts( real_t dtinv,
                                                                                              HardContactSemiImplicitTimesteppingSolvers::ContactCache& contactCache,
                                                                                              HardContactSemiImplicitTimesteppingSolvers::BodyCache& bodyCache )
{
   real_t delta_max( 0 );
   size_t numContactsMasked( contactCache.p_.size() );

   // Relax contacts
   for( size_t i = 0; i < numContactsMasked; ++i )
   {
      // Remove velocity corrections of this contact's reaction.
      bodyCache.dv_[contactCache.body1_[i]->index_] -= contactCache.body1_[i]->getInvMass() * contactCache.p_[i];
      bodyCache.dw_[contactCache.body1_[i]->index_] -= contactCache.body1_[i]->getInvInertia() * ( contactCache.r1_[i] % contactCache.p_[i] );
      bodyCache.dv_[contactCache.body2_[i]->index_] += contactCache.body2_[i]->getInvMass() * contactCache.p_[i];
      bodyCache.dw_[contactCache.body2_[i]->index_] += contactCache.body2_[i]->getInvInertia() * ( contactCache.r2_[i] % contactCache.p_[i] );

      // Calculate the relative contact VELOCITY in the global world frame (if no contact reaction is present at contact i)
      Vec3 gdot    ( ( bodyCache.v_[contactCache.body1_[i]->index_] + bodyCache.dv_[contactCache.body1_[i]->index_] ) -
            ( bodyCache.v_[contactCache.body2_[i]->index_] + bodyCache.dv_[contactCache.body2_[i]->index_] ) +
            ( bodyCache.w_[contactCache.body1_[i]->index_] + bodyCache.dw_[contactCache.body1_[i]->index_] ) % contactCache.r1_[i] -
            ( bodyCache.w_[contactCache.body2_[i]->index_] + bodyCache.dw_[contactCache.body2_[i]->index_] ) % contactCache.r2_[i] /* + diag_[i] * p */ );

      // Change from the global world frame to the contact frame
      Mat3 contactframe( contactCache.n_[i], contactCache.t_[i], contactCache.o_[i] );
      Vec3 gdot_nto( contactframe.getTranspose() * gdot );

      // The constraint in normal direction is actually a positional constraint but instead of g_n we use g_n/dt equivalently and call it gdot_n
      gdot_nto[0] += ( /* + trans( contactCache.n_[i] ) * ( contactCache.body1_[i]->getPosition() + contactCache.r1_[i] ) - ( contactCache.body2_[i]->getPosition() + contactCache.r2_[i] ) */ + contactCache.dist_[i] ) * dtinv;

      if( gdot_nto[0] >= 0 ) {
         // Contact is separating if no contact reaction is present at contact i.

         delta_max = std::max( delta_max, std::max( std::abs( contactCache.p_[i][0] ), std::max( std::abs( contactCache.p_[i][1] ), std::abs( contactCache.p_[i][2] ) ) ) );
         contactCache.p_[i] = Vec3();

         // No need to apply zero impulse.
      }
      else {
         // Contact is persisting.

         // Calculate the impulse necessary for a static contact expressed as components in the contact frame.
         Vec3 p_wf( contactCache.n_[i] * ( -contactCache.diag_n_inv_[i] * gdot_nto[0] ) );
         Vec3 dp( contactCache.p_[i] - p_wf );
         delta_max = std::max( delta_max, std::max( std::abs( dp[0] ), std::max( std::abs( dp[1] ), std::abs( dp[2] ) ) ) );

         contactCache.p_[i] = p_wf;

         // Apply impulse right away.
         bodyCache.dv_[contactCache.body1_[i]->index_] += contactCache.body1_[i]->getInvMass() * contactCache.p_[i];
         bodyCache.dw_[contactCache.body1_[i]->index_] += contactCache.body1_[i]->getInvInertia() * ( contactCache.r1_[i] % contactCache.p_[i] );
         bodyCache.dv_[contactCache.body2_[i]->index_] -= contactCache.body2_[i]->getInvMass() * contactCache.p_[i];
         bodyCache.dw_[contactCache.body2_[i]->index_] -= contactCache.body2_[i]->getInvInertia() * ( contactCache.r2_[i] % contactCache.p_[i] );
      }
   }

   return delta_max;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Relaxes all contacts once. The contact model is for inelastic unilateral contacts with approximate Coulomb friction.
 *
 * \return The largest variation of contact impulses in the L-infinity norm.
 *
 * This function is to be called from resolveContacts(). Separating contacts are preferred over
 * other solutions if valid. Static solutions are preferred over dynamic solutions. Dynamic
 * solutions are computed by decoupling the normal from the frictional components. That is
 * for a dynamic contact the normal component is relaxed first followed by the frictional
 * components. The determination of the frictional components does not perform any subiterations
 * and guarantees that the friction partially opposes slip.
 */
inline real_t HardContactSemiImplicitTimesteppingSolvers::relaxApproximateInelasticCoulombContactsByDecoupling( real_t dtinv,
                                                                                                                HardContactSemiImplicitTimesteppingSolvers::ContactCache& contactCache,
                                                                                                                HardContactSemiImplicitTimesteppingSolvers::BodyCache& bodyCache )
{
   real_t delta_max( 0 );
   size_t numContactsMasked( contactCache.p_.size() );

   // Relax contacts
   for( size_t i = 0; i < numContactsMasked; ++i ) {
      // Remove velocity corrections of this contact's reaction.
      bodyCache.dv_[contactCache.body1_[i]->index_] -= contactCache.body1_[i]->getInvMass() * contactCache.p_[i];
      bodyCache.dw_[contactCache.body1_[i]->index_] -= contactCache.body1_[i]->getInvInertia() * ( contactCache.r1_[i] % contactCache.p_[i] );
      bodyCache.dv_[contactCache.body2_[i]->index_] += contactCache.body2_[i]->getInvMass() * contactCache.p_[i];
      bodyCache.dw_[contactCache.body2_[i]->index_] += contactCache.body2_[i]->getInvInertia() * ( contactCache.r2_[i] % contactCache.p_[i] );

      // Calculate the relative contact velocity in the global world frame (if no contact reaction is present at contact i)
      Vec3 gdot    ( ( bodyCache.v_[contactCache.body1_[i]->index_] + bodyCache.dv_[contactCache.body1_[i]->index_] ) - ( bodyCache.v_[contactCache.body2_[i]->index_] + bodyCache.dv_[contactCache.body2_[i]->index_] ) + ( bodyCache.w_[contactCache.body1_[i]->index_] + bodyCache.dw_[contactCache.body1_[i]->index_] ) % contactCache.r1_[i] - ( bodyCache.w_[contactCache.body2_[i]->index_] + bodyCache.dw_[contactCache.body2_[i]->index_] ) % contactCache.r2_[i] /* + diag_[i] * p */ );

      // Change from the global world frame to the contact frame
      Mat3 contactframe( contactCache.n_[i], contactCache.t_[i], contactCache.o_[i] );
      Vec3 gdot_nto( contactframe.getTranspose() * gdot );

      //real_t gdot_n  ( trans( contactCache.n_[i] ) * gdot );  // The component of gdot along the contact normal n
      //Vec3 gdot_t  ( gdot - gdot_n * contactCache.n_[i] );  // The components of gdot tangential to the contact normal n
      //real_t g_n     ( gdot_n * dt /* + trans( contactCache.n_[i] ) * ( contactCache.body1_[i]->getPosition() + contactCache.r1_[i] ) - ( contactCache.body2_[i]->getPosition() + contactCache.r2_[i] ) */ + contactCache.dist_[i] );  // The gap in normal direction

      // The constraint in normal direction is actually a positional constraint but instead of g_n we use g_n/dt equivalently and call it gdot_n
      gdot_nto[0] += ( /* + trans( contactCache.n_[i] ) * ( contactCache.body1_[i]->getPosition() + contactCache.r1_[i] ) - ( contactCache.body2_[i]->getPosition() + contactCache.r2_[i] ) */ + contactCache.dist_[i] ) * dtinv;

      if( gdot_nto[0] >= 0 ) {
         // Contact is separating if no contact reaction is present at contact i.

         delta_max = std::max( delta_max, std::max( std::abs( contactCache.p_[i][0] ), std::max( std::abs( contactCache.p_[i][1] ), std::abs( contactCache.p_[i][2] ) ) ) );
         contactCache.p_[i] = Vec3();

         // No need to apply zero impulse.
      }
      else {
         // Contact is persisting (either static or dynamic).

         // Calculate the impulse necessary for a static contact expressed as components in the contact frame.
         Vec3 p_cf( -( contactCache.diag_nto_inv_[i] * gdot_nto ) );

         // Can p_cf[0] be negative even though -gdot_nto[0] > 0? Yes! Try:
         // A = [0.5 -0.1 +0.1; -0.1 0.5 -0.1; +0.1 -0.1 1];
         // b = [0.01 -1 -1]';
         // A\b    \approx [-0.19 -2.28 -1.21]'
         // eig(A) \approx [ 0.40  0.56  1.04]'

         real_t flimit( contactCache.mu_[i] * p_cf[0] );
         real_t fsq( p_cf[1] * p_cf[1] + p_cf[2] * p_cf[2] );
         if( fsq > flimit * flimit || p_cf[0] < 0 ) {
            // Contact cannot be static so it must be dynamic.
            // => Complementarity condition on normal reaction now turns into an equation since we know that the normal reaction is definitely not zero.

            // For simplicity we change to a simpler relaxation scheme here:
            // 1. Relax normal reaction with the tangential components equal to the previous values
            // 2. Relax tangential components with the newly relaxed normal reaction
            // Note: The better approach would be to solve the true 3x3 block problem!
            // Warning: Simply projecting the frictional components is wrong since then the normal action is no longer 0 and simulations break.

            // Add the action of the frictional reactions from the last iteration to the relative contact velocity in normal direction so we can relax it separately.
            // TODO This can be simplified:
            //p_cf = trans( contactframe ) * contactCache.p_[i];
            //p_cf[0] = 0;
            //p_[i] = contactframe * p_cf;
            Vec3 p_tmp = ( contactCache.t_[i] * contactCache.p_[i] ) * contactCache.t_[i] + ( contactCache.o_[i] * contactCache.p_[i] ) * contactCache.o_[i];

            //      |<-- This should vanish below since p_cf[0] = 0          -->|
            //gdot += ( contactCache.body1_[i]->getInvMass() + contactCache.body2_[i]->getInvMass() ) * p_tmp + ( contactCache.body1_[i]->getInvInertia() * ( contactCache.r1_[i] % p_tmp] ) ) % contactCache.r1_[i] + ( contactCache.body2_[i]->getInvInertia() * ( contactCache.r2_[i] % p_tmp ) ) % contactCache.r2_[i] /* + diag_[i] * p */;
            //real_t gdot_n = trans( contactCache.n_[i] ) * gdot;
            //gdot_n += ( /* + trans( contactCache.n_[i] ) * ( contactCache.body1_[i]->getPosition() + contactCache.r1_[i] ) - ( contactCache.body2_[i]->getPosition() + contactCache.r2_[i] ) */ + contactCache.dist_[i] ) * dtinv;

            real_t gdot_n = gdot_nto[0] + contactCache.n_[i] * ( ( contactCache.body1_[i]->getInvInertia() * ( contactCache.r1_[i] % p_tmp ) ) % contactCache.r1_[i] + ( contactCache.body2_[i]->getInvInertia() * ( contactCache.r2_[i] % p_tmp ) ) % contactCache.r2_[i] /* + diag_[i] * p */ );
            p_cf[0] = -( contactCache.diag_n_inv_[i] * gdot_n );

            // We cannot be sure that gdot_n <= 0 here and thus p_cf[0] >= 0 since we just modified it with the old values of the tangential reactions! => Project
            p_cf[0] = std::max( real_c( 0 ), p_cf[0] );

            // Now add the action of the normal reaction to the relative contact velocity in the tangential directions so we can relax the frictional components separately.
            p_tmp = contactCache.n_[i] * p_cf[0];
            Vec3 gdot2 = gdot + ( contactCache.body1_[i]->getInvInertia() * ( contactCache.r1_[i] % p_tmp ) ) % contactCache.r1_[i] + ( contactCache.body2_[i]->getInvInertia() * ( contactCache.r2_[i] % p_tmp ) ) % contactCache.r2_[i];
            Vec2 gdot_to;
            gdot_to[0] = contactCache.t_[i] * gdot2;
            gdot_to[1] = contactCache.o_[i] * gdot2;

            Vec2 ret = -( contactCache.diag_to_inv_[i] * gdot_to );
            p_cf[1] = ret[0];
            p_cf[2] = ret[1];

            flimit = contactCache.mu_[i] * p_cf[0];
            fsq = p_cf[1] * p_cf[1] + p_cf[2] * p_cf[2];
            if( fsq > flimit * flimit ) {
               const real_t f( flimit / std::sqrt( fsq ) );
               p_cf[1] *= f;
               p_cf[2] *= f;
            }
         }
         else {
            // Contact is static.
         }
         Vec3 p_wf( contactframe * p_cf );
         Vec3 dp( contactCache.p_[i] - p_wf );
         delta_max = std::max( delta_max, std::max( std::abs( dp[0] ), std::max( std::abs( dp[1] ), std::abs( dp[2] ) ) ) );

         contactCache.p_[i] = p_wf;

         // Apply impulse right away
         bodyCache.dv_[contactCache.body1_[i]->index_] += contactCache.body1_[i]->getInvMass() * contactCache.p_[i];
         bodyCache.dw_[contactCache.body1_[i]->index_] += contactCache.body1_[i]->getInvInertia() * ( contactCache.r1_[i] % contactCache.p_[i] );
         bodyCache.dv_[contactCache.body2_[i]->index_] -= contactCache.body2_[i]->getInvMass() * contactCache.p_[i];
         bodyCache.dw_[contactCache.body2_[i]->index_] -= contactCache.body2_[i]->getInvInertia() * ( contactCache.r2_[i] % contactCache.p_[i] );
      }

#if 0
      Vec3 gdot2   ( ( bodyCache.v_[contactCache.body1_[i]->index_] + bodyCache.dv_[contactCache.body1_[i]->index_] ) -
            ( bodyCache.v_[contactCache.body2_[i]->index_] + bodyCache.dv_[contactCache.body2_[i]->index_] ) +
            ( bodyCache.w_[contactCache.body1_[i]->index_] + bodyCache.dw_[contactCache.body1_[i]->index_] ) % contactCache.r1_[i] -
            ( bodyCache.w_[contactCache.body2_[i]->index_] + bodyCache.dw_contactCache.[contactCache.body2_[i]->index_] ) % contactCache.r2_[i] /* + diag_[i] * p */ );
      Vec3 gdot_nto2( contactframe.getTranspose() * gdot2 );
      WALBERLA_LOG_DETAIL( "gdot_n2 = " << gdot_nto2[0] );
      WALBERLA_LOG_DETAIL( "gdot_t2 = " << gdot_nto2[1] );
      WALBERLA_LOG_DETAIL( "gdot_o2 = " << gdot_nto2[2] );
      gdot_nto2[0] += ( /* + trans( contactCache.n_[i] ) * ( contactCache.body1_[i]->getPosition() + contactCache.r1_[i] ) - ( contactCache.body2_[i]->getPosition() + contactCache.r2_[i] ) */ + contactCache.dist_[i] ) * dtinv;
      WALBERLA_LOG_DETAIL( "gdot_n2' = " << gdot_nto2[0] );
#endif

      /*
       * compare DEM time-step with NSCD iteration:
       * - projections are the same
       * - velocities are the same if we use an explicit Euler discretization for the velocity time integration
       *
      f_cf[0] = -stiffness * contactCache.dist_ - damping_n * gdot_n = -[(stiffness * dt) * contactCache.dist_ * dtinv + damping_n * gdot_n] = -foo * (gdot_n + contactCache.dist_ * dtinv) where foo = stiffness * dt = damping_n;
      f_cf[1] = -damping_t * gdot_t                     = -damping_t * gdot_t;
      f_cf[2] = -damping_t * gdot_o                     = -damping_t * gdot_o;

      or: f_cf = -diag(foo, damping_t, damping_t) * gdot_nto   (since gdot_nto[0] is modified)
      vs. f_cf = -diaginv * gdot_nto in NSCD iteration

      => The NSCD iteration is more or less a DEM time step where we choose the stiffness and damping parameters such that penetration is non-existent after a time step and contacts are truly static (tangential rel. vel. is zero) unless the friction force hits its limit

      f_cf[0] = std::max( 0, f_cf[0] );

      flimit = contactCache.mu_ * f_cf[0];
      fsq = f_cf[1] * f_cf[1] + f_cf[2] * f_cf[2]
      if( fsq > flimit * flimit ) {
         f = flimit / sqrt( fsq );
         f_cf[1] *= f;
         f_cf[2] *= f;
      }

      f_wf = contactframe * f_cf;

      b1->addForceAtPos(  f_wf, gpos );
      b2->addForceAtPos( -f_wf, gpos );
      */

   }

   return delta_max;
}
//*************************************************************************************************


//*************************************************************************************************

/*!\brief Relaxes all contacts once. The contact model is for inelastic unilateral contacts with Coulomb friction.
 *
 * \return The largest variation of contact impulses in the L-infinity norm.
 *
 * This function is to be called from resolveContacts(). Separating contacts are preferred over
 * other solutions if valid. Static solutions are preferred over dynamic solutions. Dynamic
 * solutions are computed by decoupling the normal from the frictional components. That is
 * for a dynamic contact the normal component is relaxed first followed by the frictional
 * components. How much the frictional components directly oppose slip as required by the Coulomb
 * friction model depends on the number of subiterations performed. If no subiterations are
 * performed the friction is guaranteed to be at least partially dissipative.
 */
inline real_t HardContactSemiImplicitTimesteppingSolvers::relaxInelasticCoulombContactsByDecoupling( real_t dtinv,
                                                                                                     HardContactSemiImplicitTimesteppingSolvers::ContactCache& contactCache,
                                                                                                     HardContactSemiImplicitTimesteppingSolvers::BodyCache& bodyCache )
{
   real_t delta_max( 0 );
   size_t numContactsMasked( contactCache.p_.size() );

   // Relax contacts
   for( size_t i = 0; i < numContactsMasked; ++i ) {
      // Remove velocity corrections of this contact's reaction.
      bodyCache.dv_[contactCache.body1_[i]->index_] -= contactCache.body1_[i]->getInvMass() * contactCache.p_[i];
      bodyCache.dw_[contactCache.body1_[i]->index_] -= contactCache.body1_[i]->getInvInertia() * ( contactCache.r1_[i] % contactCache.p_[i] );
      bodyCache.dv_[contactCache.body2_[i]->index_] += contactCache.body2_[i]->getInvMass() * contactCache.p_[i];
      bodyCache.dw_[contactCache.body2_[i]->index_] += contactCache.body2_[i]->getInvInertia() * ( contactCache.r2_[i] % contactCache.p_[i] );

      // Calculate the relative contact velocity in the global world frame (if no contact reaction is present at contact i)
      Vec3 gdot    ( ( bodyCache.v_[contactCache.body1_[i]->index_] + bodyCache.dv_[contactCache.body1_[i]->index_] ) - ( bodyCache.v_[contactCache.body2_[i]->index_] + bodyCache.dv_[contactCache.body2_[i]->index_] ) + ( bodyCache.w_[contactCache.body1_[i]->index_] + bodyCache.dw_[contactCache.body1_[i]->index_] ) % contactCache.r1_[i] - ( bodyCache.w_[contactCache.body2_[i]->index_] + bodyCache.dw_[contactCache.body2_[i]->index_] ) % contactCache.r2_[i] /* + diag_[i] * p */ );

      // Change from the global world frame to the contact frame
      Mat3 contactframe( contactCache.n_[i], contactCache.t_[i], contactCache.o_[i] );
      Vec3 gdot_nto( contactframe.getTranspose() * gdot );

      //real_t gdot_n  ( trans( contactCache.n_[i] ) * gdot );  // The component of gdot along the contact normal n
      //Vec3 gdot_t  ( gdot - gdot_n * contactCache.n_[i] );  // The components of gdot tangential to the contact normal n
      //real_t g_n     ( gdot_n * dt /* + trans( contactCache.n_[i] ) * ( contactCache.body1_[i]->getPosition() + contactCache.r1_[i] ) - ( contactCache.body2_[i]->getPosition() + contactCache.r2_[i] ) */ + contactCache.dist_[i] );  // The gap in normal direction

      // The constraint in normal direction is actually a positional constraint but instead of g_n we use g_n/dt equivalently and call it gdot_n
      gdot_nto[0] += ( /* + trans( contactCache.n_[i] ) * ( contactCache.body1_[i]->getPosition() + contactCache.r1_[i] ) - ( contactCache.body2_[i]->getPosition() + contactCache.r2_[i] ) */ + contactCache.dist_[i] ) * dtinv;

      //WALBERLA_LOG_WARNING( "Contact #" << i << " is\nA = \n" << contactCache.diag_nto_[i] << "\nb = \n" << gdot_nto << "\nmu = " << contactCache.mu_[i] );

      if( gdot_nto[0] >= 0 ) {
         // Contact is separating if no contact reaction is present at contact i.

         delta_max = std::max( delta_max, std::max( std::abs( contactCache.p_[i][0] ), std::max( std::abs( contactCache.p_[i][1] ), std::abs( contactCache.p_[i][2] ) ) ) );
         contactCache.p_[i] = Vec3();
         //WALBERLA_LOG_WARNING( "Contact #" << i << " is separating." );

         // No need to apply zero impulse.
      }
      else {
         // Contact is persisting (either static or dynamic).

         // Calculate the impulse necessary for a static contact expressed as components in the contact frame.
         Vec3 p_cf( -( contactCache.diag_nto_inv_[i] * gdot_nto ) );

         // Can p_cf[0] be negative even though -gdot_nto[0] > 0? Yes! Try:
         // A = [0.5 -0.1 +0.1; -0.1 0.5 -0.1; +0.1 -0.1 1];
         // b = [0.01 -1 -1]';
         // A\b    \approx [-0.19 -2.28 -1.21]'
         // eig(A) \approx [ 0.40  0.56  1.04]'

         real_t flimit( contactCache.mu_[i] * p_cf[0] );
         real_t fsq( p_cf[1] * p_cf[1] + p_cf[2] * p_cf[2] );
         if( fsq > flimit * flimit || p_cf[0] < 0 ) {
            // Contact cannot be static so it must be dynamic.
            // => Complementarity condition on normal reaction now turns into an equation since we know that the normal reaction is definitely not zero.

            for (int j = 0; j < 20; ++j) {
               // For simplicity we change to a simpler relaxation scheme here:
               // 1. Relax normal reaction with the tangential components equal to the previous values
               // 2. Relax tangential components with the newly relaxed normal reaction
               // Note: The better approach would be to solve the true 3x3 block problem!
               // Warning: Simply projecting the frictional components is wrong since then the normal action is no longer 0 and simulations break.

               Vec3 gdotCorrected;
               real_t gdotCorrected_n;
               Vec2 gdotCorrected_to;

               // Calculate the relative contact velocity in the global world frame (if no normal contact reaction is present at contact i)
               p_cf[0] = 0;
               //                       |<- p_cf is orthogonal to the normal and drops out in next line ->|
               gdotCorrected   = /* ( contactCache.body1_[i]->getInvMass() + contactCache.body2_[i]->getInvMass() ) * p_cf  */ gdot + ( contactCache.body1_[i]->getInvInertia() * ( contactCache.r1_[i] % ( contactCache.t_[i] * p_cf[1] + contactCache.o_[i] * p_cf[2] ) ) ) % contactCache.r1_[i] + ( contactCache.body2_[i]->getInvInertia() * ( contactCache.r2_[i] % ( contactCache.t_[i] * p_cf[1] + contactCache.o_[i] * p_cf[2] ) ) ) % contactCache.r2_[i];
               gdotCorrected_n = contactCache.n_[i] * gdotCorrected + contactCache.dist_[i] * dtinv;

               // Relax normal component.
               p_cf[0] = std::max( real_c( 0 ), -( contactCache.diag_n_inv_[i] * gdotCorrected_n ) );

               // Calculate the relative contact velocity in the global world frame (if no frictional contact reaction is present at contact i)
               p_cf[1] = p_cf[2] = real_c( 0 );
               //                       |<- p_cf is orthogonal to the tangential plane and drops out   ->|
               gdotCorrected   = /* ( contactCache.body1_[i]->getInvMass() + contactCache.body2_[i]->getInvMass() ) * p_cf */ gdot + ( contactCache.body1_[i]->getInvInertia() * ( contactCache.r1_[i] % ( contactCache.n_[i] * p_cf[0] ) ) ) % contactCache.r1_[i] + ( contactCache.body2_[i]->getInvInertia() * ( contactCache.r2_[i] % ( contactCache.n_[i] * p_cf[0] ) ) ) % contactCache.r2_[i];
               gdotCorrected_to[0] = contactCache.t_[i] * gdotCorrected;
               gdotCorrected_to[1] = contactCache.o_[i] * gdotCorrected;

               // Relax frictional components.
               Vec2 ret = -( contactCache.diag_to_inv_[i] * gdotCorrected_to );
               p_cf[1] = ret[0];
               p_cf[2] = ret[1];

               flimit = contactCache.mu_[i] * p_cf[0];
               fsq = p_cf[1] * p_cf[1] + p_cf[2] * p_cf[2];
               if( fsq > flimit * flimit ) {
                  // 3.2.1 Decoupling
                  // \tilde{x}^0 = p_cf[1..2]

                  // Determine \tilde{A}
                  Mat2 diag_to( contactCache.diag_nto_[i](1, 1), contactCache.diag_nto_[i](1, 2), contactCache.diag_nto_[i](2, 1), contactCache.diag_nto_[i](2, 2) );

                        const real_t f( flimit / std::sqrt( fsq ) );
                  //p_cf[1] *= f;
                  //p_cf[2] *= f;

                  // Determine search interval for Golden Section Search
                  const real_t phi( real_c(0.5) * ( real_c(1) + std::sqrt( real_c( 5 ) ) ) );
                  real_t shift( std::atan2( -p_cf[2], p_cf[1] ) );
                  real_t acos_f( std::acos( f ) );

                  //WALBERLA_LOG_WARNING( acos_f << " " << shift );

                  real_t alpha_left( -acos_f - shift );
                  //Vec2 x_left( flimit * std::cos( alpha_left ), flimit * std::sin( alpha_left ) );
                  //real_t f_left( 0.5 * trans( x_left ) * ( diag_to * x_left ) - trans( x_left ) * ( -gdot_to ) );

                  real_t alpha_right( acos_f - shift );
                  //Vec2 x_right( flimit * std::cos( alpha_right ), flimit * std::sin( alpha_right ) );
                  //real_t f_right( 0.5 * trans( x_right ) * ( diag_to * x_right ) - trans( x_right ) * ( -gdot_to ) );

                  real_t alpha_mid( ( alpha_right + alpha_left * phi ) / ( 1 + phi ) );
                  Vec2 x_mid( flimit * std::cos( alpha_mid ), flimit * std::sin( alpha_mid ) );
                  real_t f_mid( real_c(0.5) * x_mid * ( diag_to * x_mid ) - x_mid * ( -gdotCorrected_to ) );

                  bool leftlarger = false;
                  for( size_t k = 0; k < maxSubIterations_; ++k ) {
                     real_t alpha_next( alpha_left + ( alpha_right - alpha_mid ) );
                     Vec2 x_next( flimit * std::cos( alpha_next ), flimit * std::sin( alpha_next ) );
                     real_t f_next( real_c(0.5) * x_next * ( diag_to * x_next ) - x_next * ( -gdotCorrected_to ) );
                     //WALBERLA_LOG_WARNING( "[(" << alpha_left << ", ?); (" << alpha_mid << ", " << f_mid << "); (" << alpha_right << ", ?)] <- (" << alpha_next << ", " << f_next << ")" );
                     //WALBERLA_LOG_WARNING( "left: " << alpha_mid - alpha_left << "  right: " << alpha_right - alpha_mid << "  ll: " << leftlarger );
                     //WALBERLA_ASSERT(leftlarger ? (alpha_mid - alpha_left > alpha_right - alpha_mid) : (alpha_mid - alpha_left < alpha_right - alpha_mid), "ll inconsistent!" );

                     if (leftlarger) {
                        // left interval larger
                        if( f_next < f_mid ) {
                           alpha_right = alpha_mid;
                           alpha_mid   = alpha_next;
                           x_mid       = x_next;
                           f_mid       = f_next;
                           leftlarger = true;
                        }
                        else {
                           alpha_left  = alpha_next;
                           leftlarger = false;
                        }
                     }
                     else {
                        // right interval larger
                        if( f_next < f_mid ) {
                           alpha_left = alpha_mid;
                           alpha_mid  = alpha_next;
                           x_mid      = x_next;
                           f_mid      = f_next;
                           leftlarger = false;
                        }
                        else {
                           alpha_right = alpha_next;
                           leftlarger = true;
                        }
                     }
                  }
                  //WALBERLA_LOG_WARNING( "dalpha = " << alpha_right - alpha_left );

                  p_cf[1] = x_mid[0];
                  p_cf[2] = x_mid[1];
               }
            }
            //WALBERLA_LOG_WARNING( "Contact #" << i << " is dynamic." );
         }
         else {
            // Contact is static.
            //WALBERLA_LOG_WARNING( "Contact #" << i << " is static." );
         }

         //WALBERLA_LOG_WARNING( "Contact reaction in contact frame: " << p_cf << "\n" << contactCache.diag_nto_[i]*p_cf + gdot_nto );
         Vec3 p_wf( contactframe * p_cf );
         Vec3 dp( contactCache.p_[i] - p_wf );
         delta_max = std::max( delta_max, std::max( std::abs( dp[0] ), std::max( std::abs( dp[1] ), std::abs( dp[2] ) ) ) );

         contactCache.p_[i] = p_wf;

         // Apply impulse right away
         bodyCache.dv_[contactCache.body1_[i]->index_] += contactCache.body1_[i]->getInvMass() * contactCache.p_[i];
         bodyCache.dw_[contactCache.body1_[i]->index_] += contactCache.body1_[i]->getInvInertia() * ( contactCache.r1_[i] % contactCache.p_[i] );
         bodyCache.dv_[contactCache.body2_[i]->index_] -= contactCache.body2_[i]->getInvMass() * contactCache.p_[i];
         bodyCache.dw_[contactCache.body2_[i]->index_] -= contactCache.body2_[i]->getInvInertia() * ( contactCache.r2_[i] % contactCache.p_[i] );
      }

#if 0
      Vec3 gdot2   ( ( bodyCache.v_[contactCache.body1_[i]->index_] + bodyCache.dv_[contactCache.body1_[i]->index_] ) -
            ( bodyCache.v_[contactCache.body2_[i]->index_] + bodyCache.dv_[contactCache.body2_[i]->index_] ) +
            ( bodyCache.w_[contactCache.body1_[i]->index_] + bodyCache.dw_[contactCache.body1_[i]->index_] ) % contactCache.r1_[i] -
            ( bodyCache.w_[contactCache.body2_[i]->index_] + bodyCache.dw_[contactCache.body2_[i]->index_] ) % contactCache.r2_[i] /* + diag_[i] * p */ );
      Vec3 gdot_nto2( contactframe.getTranspose() * gdot2 );
      WALBERLA_LOG_DETAIL( "gdot_n2 = " << gdot_nto2[0] );
      WALBERLA_LOG_DETAIL( "gdot_t2 = " << gdot_nto2[1] );
      WALBERLA_LOG_DETAIL( "gdot_o2 = " << gdot_nto2[2] );
   }
   gdot_nto2[0] += ( /* + trans( contactCache.n_[i] ) * ( contactCache.body1_[i]->getPosition() + contactCache.r1_[i] ) - ( contactCache.body2_[i]->getPosition() + contactCache.r2_[i] ) */ + contactCache.dist_[i] ) * dtinv;
   WALBERLA_LOG_DETAIL( "gdot_n2' = " << gdot_nto2[0] );
}
#endif

/*
       * compare DEM time-step with NSCD iteration:
       * - projections are the same
       * - velocities are the same if we use an explicit Euler discretization for the velocity time integration
       *
      f_cf[0] = -stiffness * contactCache.dist_ - damping_n * gdot_n = -[(stiffness * dt) * contactCache.dist_ * dtinv + damping_n * gdot_n] = -foo * (gdot_n + contactCache.dist_ * dtinv) where foo = stiffness * dt = damping_n;
      f_cf[1] = -damping_t * gdot_t                     = -damping_t * gdot_t;
      f_cf[2] = -damping_t * gdot_o                     = -damping_t * gdot_o;

      or: f_cf = -diag(foo, damping_t, damping_t) * gdot_nto   (since gdot_nto[0] is modified)
      vs. f_cf = -diaginv * gdot_nto in NSCD iteration

      => The NSCD iteration is more or less a DEM time step where we choose the stiffness and damping parameters such that penetration is non-existent after a time step and contacts are truly static (tangential rel. vel. is zero) unless the friction force hits its limit

      f_cf[0] = std::max( 0, f_cf[0] );

      flimit = contactCache.mu_ * f_cf[0];
      fsq = f_cf[1] * f_cf[1] + f_cf[2] * f_cf[2]
      if( fsq > flimit * flimit ) {
         f = flimit / sqrt( fsq );
         f_cf[1] *= f;
         f_cf[2] *= f;
      }

      f_wf = contactframe * f_cf;

      b1->addForceAtPos(  f_wf, gpos );
      b2->addForceAtPos( -f_wf, gpos );
      */
}

return delta_max;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Relaxes all contacts once. The contact model is for inelastic unilateral contacts with Coulomb friction.
 *
 * \param dtinv The inverse of the current time step.
 * \param approximate Use the approximate model showing bouncing.
 * \return The largest variation of contact impulses in the L-infinity norm.
 *
 * This function is to be called from resolveContacts(). The iterative method to solve the contact
 * problem is e.g. described in the article "A matrix-free cone complementarity approach for
 * solving large-scale, nonsmooth, rigid body dynamics" by A. Tasora and M. Anitescu in Computer
 * Methods in Applied Mechanics and Engineering (Volume 200, Issues 5–8, 15 January 2011,
 * Pages 439-453). The contact model is purely quadratic and convergence should be good but depends
 * on a parameter. The one-contact problem has a unique solution. The frictional reactions
 * for a dynamic contact converge to those that directly oppose slip. However, the contact is
 * not perfectly inelastic for dynamic contacts but bounces. These vertical motions tend to
 * go to zero for smaller time steps and can be interpreted as exaggerated vertical motions
 * coming from micro asperities (see "Optimization-based simulation of nonsmooth rigid multibody
 * dynamics" by M. Anitescu in Mathematical Programming (Volume 105, Issue 1, January 2006, Pages
 * 113-143). These motions can be prevented by a small change in the iteration proposed in "The
 * bipotential method: a constructive approach to design the complete contact law with friction and
 * improved numerical algorithms" by G. De Saxce and Z-Q. Feng in Mathematical and Computer
 * Modelling (Volume 28, Issue 4, 1998, Pages 225-245). Which iteration is used is controlled with
 * the approximate parameter.
 */
inline real_t HardContactSemiImplicitTimesteppingSolvers::relaxInelasticCoulombContactsByOrthogonalProjections( real_t dtinv,
                                                                                                                bool approximate,
                                                                                                                HardContactSemiImplicitTimesteppingSolvers::ContactCache& contactCache,
                                                                                                                HardContactSemiImplicitTimesteppingSolvers::BodyCache& bodyCache )
{
   real_t delta_max( 0 );
   size_t numContactsMasked( contactCache.p_.size() );

   // Relax contacts
   for( size_t i = 0; i < numContactsMasked; ++i ) {
      // Remove velocity corrections of this contact's reaction.
      bodyCache.dv_[contactCache.body1_[i]->index_] -= contactCache.body1_[i]->getInvMass() * contactCache.p_[i];
      bodyCache.dw_[contactCache.body1_[i]->index_] -= contactCache.body1_[i]->getInvInertia() * ( contactCache.r1_[i] % contactCache.p_[i] );
      bodyCache.dv_[contactCache.body2_[i]->index_] += contactCache.body2_[i]->getInvMass() * contactCache.p_[i];
      bodyCache.dw_[contactCache.body2_[i]->index_] += contactCache.body2_[i]->getInvInertia() * ( contactCache.r2_[i] % contactCache.p_[i] );

      // Calculate the relative contact velocity in the global world frame (as if no contact reaction was present at contact i)
      Vec3 gdot(   ( bodyCache.v_[contactCache.body1_[i]->index_] + bodyCache.dv_[contactCache.body1_[i]->index_] )
                 - ( bodyCache.v_[contactCache.body2_[i]->index_] + bodyCache.dv_[contactCache.body2_[i]->index_] )
                 + ( bodyCache.w_[contactCache.body1_[i]->index_] + bodyCache.dw_[contactCache.body1_[i]->index_] ) % contactCache.r1_[i]
                 - ( bodyCache.w_[contactCache.body2_[i]->index_] + bodyCache.dw_[contactCache.body2_[i]->index_] ) % contactCache.r2_[i] /* + diag_[i] * p */ );

      // Change from the global world frame to the contact frame
      Mat3 contactframe( contactCache.n_[i], contactCache.t_[i], contactCache.o_[i] );
      Vec3 gdot_nto( contactframe.getTranspose() * gdot );

      //real_t gdot_n( trans( contactCache.n_[i] ) * gdot );  // The component of gdot along the contact normal n
      //Vec3 gdot_t( gdot - gdot_n * contactCache.n_[i] );  // The components of gdot tangential to the contact normal n
      //real_t g_n( gdot_n * dt /* + trans( contactCache.n_[i] ) * ( contactCache.body1_[i]->getPosition() + contactCache.r1_[i] ) - ( contactCache.body2_[i]->getPosition() + contactCache.r2_[i] ) */ + contactCache.dist_[i] );  // The gap in normal direction

      // The constraint in normal direction is actually a positional constraint but instead of g_n we use g_n/dt equivalently and call it gdot_n
      gdot_nto[0] += ( /* + trans( contactCache.n_[i] ) * ( contactCache.body1_[i]->getPosition() + contactCache.r1_[i] )
                          - ( contactCache.body2_[i]->getPosition() + contactCache.r2_[i] ) */ + contactCache.dist_[i] ) * dtinv;

      // Tasora et al., 2011 starts here
      // -------------------------------------------------------------------------------------

      Vec3 p_cf( contactframe.getTranspose() * contactCache.p_[i] );
      if( approximate ) {
         // Calculate next iterate (Tasora et al., 2011).
         p_cf = p_cf - overRelaxationParam_ * ( contactCache.diag_nto_[i] * p_cf + gdot_nto );
      }
      else {
         // Calculate next iterate (De Saxce/Feng).
         Vec3 tmp( contactCache.diag_nto_[i] * p_cf + gdot_nto );
         tmp[0] += std::sqrt( math::sq( tmp[1] ) + math::sq( tmp[2] ) ) * contactCache.mu_[i];
         p_cf = p_cf - overRelaxationParam_ * tmp;
      }

      // Project on friction cone (Tasora et al., 2011; Fig. 2; Eq. (46)).
      real_t flimit( contactCache.mu_[i] * p_cf[0] );
      real_t fsq( p_cf[1] * p_cf[1] + p_cf[2] * p_cf[2] );

      if( p_cf[0] > 0 && fsq < flimit * flimit ) {
         // Unconstrained minimum is INSIDE THE UPPER CONE,
         // with slope n*mu
         // leading to a static contact and no projection is necessary.
      }
      else if( p_cf[0] < 0 && fsq < math::sq( p_cf[0] / contactCache.mu_[i] ) ) {
         // Unconstrained minimum is INSIDE THE LOWER CONE,
         // with slope n/mu
         // leading to a separating contact where no contact reaction is present,
         // i.e. the unconstrained minimum is projected to the tip of the lower cone.
         p_cf = Vec3();
      }
      else {
         // Unconstrained minimum is OUTSIDE THE CONE -> Project on cone surface, Eq. (45) in Tasora et al., 2011
         real_t f( std::sqrt( fsq ) ); // gamma_r
         p_cf[0] = ( f * contactCache.mu_[i] + p_cf[0] ) / ( math::sq( contactCache.mu_[i] ) + 1 ); // p_cf[0] = gamma_n
         
         real_t factor( contactCache.mu_[i] * p_cf[0] / f );
         p_cf[1] *= factor;
         p_cf[2] *= factor;
      }

      // Tasora et al., 2011 ends here
      // -------------------------------------------------------------------------------------

      Vec3 p_wf( contactframe * p_cf );
      Vec3 dp( contactCache.p_[i] - p_wf );
      delta_max = std::max( delta_max, std::max( std::abs( dp[0] ), std::max( std::abs( dp[1] ), std::abs( dp[2] ) ) ) );

      contactCache.p_[i] = p_wf;

      // Apply impulse right away
      bodyCache.dv_[contactCache.body1_[i]->index_] += contactCache.body1_[i]->getInvMass() * contactCache.p_[i];
      bodyCache.dw_[contactCache.body1_[i]->index_] += contactCache.body1_[i]->getInvInertia() * ( contactCache.r1_[i] % contactCache.p_[i] );
      bodyCache.dv_[contactCache.body2_[i]->index_] -= contactCache.body2_[i]->getInvMass() * contactCache.p_[i];
      bodyCache.dw_[contactCache.body2_[i]->index_] -= contactCache.body2_[i]->getInvInertia() * ( contactCache.r2_[i] % contactCache.p_[i] );
   }

   return delta_max;
}
//*************************************************************************************************


//*************************************************************************************************

/*!\brief Relaxes all contacts once. The contact model is for inelastic unilateral contacts with the generalized maximum dissipation principle for friction.
 *
 * \return The largest variation of contact impulses in the L-infinity norm.
 *
 * This function is to be called from resolveContacts(). Dynamic solutions are computed by
 * minimizing the kinetic energy along the intersection of the plane of maximum compression and
 * the friction cone.
 */
inline real_t HardContactSemiImplicitTimesteppingSolvers::relaxInelasticGeneralizedMaximumDissipationContacts( real_t dtinv,
                                                                                                               HardContactSemiImplicitTimesteppingSolvers::ContactCache& contactCache,
                                                                                                               HardContactSemiImplicitTimesteppingSolvers::BodyCache& bodyCache )
{
   real_t delta_max( 0 );
   size_t numContactsMasked( contactCache.p_.size() );

   // Relax contacts
   for( size_t i = 0; i < numContactsMasked; ++i ) {
      // Remove velocity corrections of this contact's reaction.
      bodyCache.dv_[contactCache.body1_[i]->index_] -= contactCache.body1_[i]->getInvMass() * contactCache.p_[i];
      bodyCache.dw_[contactCache.body1_[i]->index_] -= contactCache.body1_[i]->getInvInertia() * ( contactCache.r1_[i] % contactCache.p_[i] );
      bodyCache.dv_[contactCache.body2_[i]->index_] += contactCache.body2_[i]->getInvMass() * contactCache.p_[i];
      bodyCache.dw_[contactCache.body2_[i]->index_] += contactCache.body2_[i]->getInvInertia() * ( contactCache.r2_[i] % contactCache.p_[i] );

      // Calculate the relative contact velocity in the global world frame (if no contact reaction is present at contact i)
      Vec3 gdot    ( ( bodyCache.v_[contactCache.body1_[i]->index_] + bodyCache.dv_[contactCache.body1_[i]->index_] ) - ( bodyCache.v_[contactCache.body2_[i]->index_] + bodyCache.dv_[contactCache.body2_[i]->index_] ) + ( bodyCache.w_[contactCache.body1_[i]->index_] + bodyCache.dw_[contactCache.body1_[i]->index_] ) % contactCache.r1_[i] - ( bodyCache.w_[contactCache.body2_[i]->index_] + bodyCache.dw_[contactCache.body2_[i]->index_] ) % contactCache.r2_[i] /* + diag_[i] * p */ );

      // Change from the global world frame to the contact frame
      Mat3 contactframe( contactCache.n_[i], contactCache.t_[i], contactCache.o_[i] );
      Vec3 gdot_nto( contactframe.getTranspose() * gdot );

      // The constraint in normal direction is actually a positional constraint but instead of g_n we use g_n/dt equivalently and call it gdot_n
      gdot_nto[0] += ( /* + trans( contactCache.n_[i] ) * ( contactCache.body1_[i]->getPosition() + contactCache.r1_[i] ) - ( contactCache.body2_[i]->getPosition() + contactCache.r2_[i] ) */ + contactCache.dist_[i] ) * dtinv;

      //WALBERLA_LOG_WARNING( "Contact #" << i << " is\nA = \n" << contactCache.diag_nto_[i] << "\nb = \n" << gdot_nto << "\nmu = " << contactCache.mu_[i] );

      if( gdot_nto[0] >= 0 ) {
         // Contact is separating if no contact reaction is necessary without violating the penetration constraint.

         delta_max = std::max( delta_max, std::max( std::abs( contactCache.p_[i][0] ), std::max( std::abs( contactCache.p_[i][1] ), std::abs( contactCache.p_[i][2] ) ) ) );
         contactCache.p_[i] = Vec3();

         //WALBERLA_LOG_WARNING( "Contact #" << i << " is separating." );

         // No need to apply zero impulse.
      }
      else {
         // Contact is persisting (either static or dynamic).

         // Calculate the impulse necessary for a static contact expressed as components in the contact frame.
         Vec3 p_cf( -( contactCache.diag_nto_inv_[i] * gdot_nto ) );

         // Can p_cf[0] be negative even though -gdot_nto[0] > 0? Yes! Try:
         // A = [0.5 -0.1 +0.1; -0.1 0.5 -0.1; +0.1 -0.1 1];
         // b = [0.01 -1 -1]';
         // A\b    \approx [-0.19 -2.28 -1.21]'
         // eig(A) \approx [ 0.40  0.56  1.04]'

         real_t flimit( contactCache.mu_[i] * p_cf[0] );
         real_t fsq( p_cf[1] * p_cf[1] + p_cf[2] * p_cf[2] );
         if( fsq > flimit * flimit || p_cf[0] < 0 ) {
            // Contact cannot be static so it must be dynamic.
            // => Complementarity condition on normal reaction now turns into an equation since we know that the normal reaction is definitely not zero.

            // \breve{x}^0 = p_cf[1..2]

            // Eliminate normal component from 3x3 system: contactCache.diag_nto_[i]*p_cf + gdot_nto => \breve{A} \breve{x} - \breve{b}
            const real_t invA_nn( math::inv( contactCache.diag_nto_[i](0, 0) ) );
                                  const real_t offdiag( contactCache.diag_nto_[i](1, 2) - invA_nn * contactCache.diag_nto_[i](0, 1) * contactCache.diag_nto_[i](0, 2) );
                                  Mat2 A_breve( contactCache.diag_nto_[i](1, 1) - invA_nn *math::sq( contactCache.diag_nto_[i](0, 1) ), offdiag, offdiag, contactCache.diag_nto_[i](2, 2) - invA_nn *math::sq( contactCache.diag_nto_[i](0, 2) ) );
                                                                                                     Vec2 b_breve( -gdot_nto[1] + invA_nn * contactCache.diag_nto_[i](0, 1) * gdot_nto[0], -gdot_nto[2] + invA_nn * contactCache.diag_nto_[i](0, 2) * gdot_nto[0] );

                                                const real_t shiftI( std::atan2( -contactCache.diag_nto_[i](0, 2), contactCache.diag_nto_[i](0, 1) ) );
                                                                     const real_t shiftJ( std::atan2( -p_cf[2], p_cf[1] ) );
                                  const real_t a3( std::sqrt(math::sq( contactCache.diag_nto_[i](0, 1) ) +math::sq( contactCache.diag_nto_[i](0, 2) ) ) );
                                                                       const real_t fractionI( -contactCache.diag_nto_[i](0, 0) / ( contactCache.mu_[i] * a3 ) );
                                                                       const real_t fractionJ( std::min( invA_nn * contactCache.mu_[i] * ( ( -gdot_nto[0] ) / std::sqrt( fsq ) - a3 * std::cos( shiftI - shiftJ ) ), real_c( 1 ) ) );

                                                             // Search interval determination.
                                                             real_t alpha_left, alpha_right;
                                                   if( fractionJ < -1 ) {
                                                      // J is complete
                                                      const real_t angleI( std::acos( fractionI ) );
                                                      alpha_left = -angleI - shiftI;
                                                      alpha_right = +angleI - shiftI;
                                                      if( alpha_left < 0 ) {
                                                         alpha_left += 2 * math::pi;
                                                         alpha_right += 2 * math::pi;
                                                      }
                                                   }
                                                   else if( contactCache.diag_nto_[i](0, 0) > contactCache.mu_[i] * a3 ) {
               // I is complete
               const real_t angleJ( std::acos( fractionJ ) );
               alpha_left = -angleJ - shiftJ;
               alpha_right = +angleJ - shiftJ;
               if( alpha_left < 0 ) {
                  alpha_left += 2 * math::pi;
                  alpha_right += 2 * math::pi;
               }
            }
            else {
               // neither I nor J is complete
               const real_t angleJ( std::acos( fractionJ ) );
               real_t alpha1_left( -angleJ - shiftJ );
               real_t alpha1_right( +angleJ - shiftJ );
               if( alpha1_left < 0 ) {
                  alpha1_left += 2 * math::pi;
                  alpha1_right += 2 * math::pi;
               }
               const real_t angleI( std::acos( fractionI ) );
               real_t alpha2_left( -angleI - shiftI );
               real_t alpha2_right( +angleI - shiftI );
               if( alpha2_left < 0 ) {
                  alpha2_left += 2 * math::pi;
                  alpha2_right += 2 * math::pi;
               }

               // Swap intervals if second interval does not start right of the first interval.
               if( alpha1_left > alpha2_left ) {
                  std::swap( alpha1_left, alpha2_left );
                  std::swap( alpha1_right, alpha2_right );
               }

               if( alpha2_left > alpha1_right ) {
                  alpha2_right -= 2*math::pi;
                  if( alpha2_right > alpha1_right ) {
                     // [alpha1_left; alpha1_right] \subset [alpha2_left; alpha2_right]
                  }
                  else {
                     // [alpha2_left; alpha2_right] intersects the left end of [alpha1_left; alpha1_right]
                     alpha1_right = alpha2_right;
                  }
               }
               else {
                  alpha1_left = alpha2_left;
                  if( alpha2_right > alpha1_right ) {
                     // [alpha2_left; alpha2_right] intersects the right end of [alpha1_left; alpha1_right]
                  }
                  else {
                     // [alpha2_left; alpha2_right] \subset [alpha1_left; alpha1_right]
                     alpha1_right = alpha2_right;
                  }
               }

               alpha_left = alpha1_left;
               alpha_right = alpha1_right;
            }

            const real_t phi( real_c(0.5) * ( real_c(1) + std::sqrt( real_c( 5 ) ) ) );
                                                   real_t alpha_mid( ( alpha_right + alpha_left * phi ) / ( 1 + phi ) );
                                  Vec2 x_mid;
                  real_t f_mid;

            {
               real_t r_ub = contactCache.mu_[i] * ( -gdot_nto[0] ) / ( contactCache.diag_nto_[i](0, 0) + contactCache.mu_[i] * a3 * std::cos( alpha_mid + shiftI ) );
                     if( r_ub < 0 )
                     r_ub = math::Limits<real_t>::inf();
               x_mid = Vec2( r_ub * std::cos( alpha_mid ), r_ub * std::sin( alpha_mid ) );
               f_mid = real_c(0.5) * x_mid * ( A_breve * x_mid ) - x_mid * b_breve;
            }

            bool leftlarger = false;
            for( size_t k = 0; k < maxSubIterations_; ++k ) {
               real_t alpha_next( alpha_left + ( alpha_right - alpha_mid ) );
               real_t r_ub = contactCache.mu_[i] * ( -gdot_nto[0] ) / ( contactCache.diag_nto_[i](0, 0) + contactCache.mu_[i] * a3 * std::cos( alpha_next + shiftI ) );
                     if( r_ub < 0 )
                     r_ub = math::Limits<real_t>::inf();
               Vec2 x_next( r_ub * std::cos( alpha_next ), r_ub * std::sin( alpha_next ) );
               real_t f_next( real_c(0.5) * x_next * ( A_breve * x_next ) - x_next * b_breve );

               //WALBERLA_LOG_WARNING( "[(" << alpha_left << ", ?); (" << alpha_mid << ", " << f_mid << "); (" << alpha_right << ", ?)] <- (" << alpha_next << ", " << f_next << ")" );
               //WALBERLA_LOG_WARNING( "left: " << alpha_mid - alpha_left << "  right: " << alpha_right - alpha_mid << "  ll: " << leftlarger );
               //WALBERLA_ASSERT(leftlarger ? (alpha_mid - alpha_left > alpha_right - alpha_mid) : (alpha_mid - alpha_left < alpha_right - alpha_mid), "ll inconsistent!" );

               if (leftlarger) {
                  // left interval larger
                  if( f_next < f_mid ) {
                     alpha_right = alpha_mid;
                     alpha_mid   = alpha_next;
                     x_mid       = x_next;
                     f_mid       = f_next;
                     leftlarger = true;
                  }
                  else {
                     alpha_left  = alpha_next;
                     leftlarger = false;
                  }
               }
               else {
                  // right interval larger
                  if( f_next < f_mid ) {
                     alpha_left = alpha_mid;
                     alpha_mid  = alpha_next;
                     x_mid      = x_next;
                     f_mid      = f_next;
                     leftlarger = false;
                  }
                  else {
                     alpha_right = alpha_next;
                     leftlarger = true;
                  }
               }
            }
            //WALBERLA_LOG_DETAIL( "dalpha = " << alpha_right - alpha_left << "\n");
            {
               real_t alpha_init( std::atan2( p_cf[2], p_cf[1] ) );
               real_t r_ub = contactCache.mu_[i] * ( -gdot_nto[0] ) / ( contactCache.diag_nto_[i](0, 0) + contactCache.mu_[i] * a3 * std::cos( alpha_init + shiftI ) );
                     if( r_ub < 0 )
                     r_ub = math::Limits<real_t>::inf();
               Vec2 x_init( r_ub * std::cos( alpha_init ), r_ub * std::sin( alpha_init ) );
               real_t f_init( real_c(0.5) * x_init * ( A_breve * x_init ) - x_init * b_breve );

               if( f_init < f_mid )
               {
                  x_mid = x_init;
                  WALBERLA_LOG_DETAIL( "Replacing solution by primitive dissipative solution (" << f_init << " < " << f_mid << " at " << alpha_init << " vs. " << alpha_mid << ").\n");
               }
            }

            p_cf[0] = invA_nn * ( -gdot_nto[0] - contactCache.diag_nto_[i](0, 1) * x_mid[0] - contactCache.diag_nto_[i](0, 2) * x_mid[1] );
            p_cf[1] = x_mid[0];
            p_cf[2] = x_mid[1];
            //WALBERLA_LOG_DETAIL( "Contact #" << i << " is dynamic." );
         }
         else {
            // Contact is static.
            //WALBERLA_LOG_DETAIL( "Contact #" << i << " is static." );
         }
         Vec3 p_wf( contactframe * p_cf );
         Vec3 dp( contactCache.p_[i] - p_wf );
         delta_max = std::max( delta_max, std::max( std::abs( dp[0] ), std::max( std::abs( dp[1] ), std::abs( dp[2] ) ) ) );
         //WALBERLA_LOG_DETAIL( "Contact reaction in contact frame: " << p_cf << "\nContact action in contact frame: " << contactCache.diag_nto_[i]*p_cf + gdot_nto );

         contactCache.p_[i] = p_wf;

         // Apply impulse right away
         bodyCache.dv_[contactCache.body1_[i]->index_] += contactCache.body1_[i]->getInvMass() * contactCache.p_[i];
         bodyCache.dw_[contactCache.body1_[i]->index_] += contactCache.body1_[i]->getInvInertia() * ( contactCache.r1_[i] % contactCache.p_[i] );
         bodyCache.dv_[contactCache.body2_[i]->index_] -= contactCache.body2_[i]->getInvMass() * contactCache.p_[i];
         bodyCache.dw_[contactCache.body2_[i]->index_] -= contactCache.body2_[i]->getInvInertia() * ( contactCache.r2_[i] % contactCache.p_[i] );
      }

#if 0
      Vec3 gdot2   ( ( bodyCache.v_[contactCache.body1_[i]->index_] + bodyCache.dv_[contactCache.body1_[i]->index_] ) -
            ( bodyCache.v_[contactCache.body2_[i]->index_] + bodyCache.dv_[contactCache.body2_[i]->index_] ) +
            ( bodyCache.w_[contactCache.body1_[i]->index_] + bodyCache.dw_[contactCache.body1_[i]->index_] ) % contactCache.r1_[i] -
            ( bodyCache.w_[contactCache.body2_[i]->index_] + bodyCache.dw_[contactCache.body2_[i]->index_] ) % contactCache.r2_[i] /* + diag_[i] * p */ );
      Vec3 gdot_nto2( contactframe.getTranspose() * gdot2 );
      WALBERLA_LOG_DETAIL( "gdot_n2 = " << gdot_nto2[0] );
      WALBERLA_LOG_DETAIL( "gdot_t2 = " << gdot_nto2[1] );
      WALBERLA_LOG_DETAIL( "gdot_o2 = " << gdot_nto2[2] );

      gdot_nto2[0] += ( /* + trans( contactCache.n_[i] ) * ( contactCache.body1_[i]->getPosition() + contactCache.r1_[i] ) - ( contactCache.body2_[i]->getPosition() + contactCache.r2_[i] ) */ + contactCache.dist_[i] ) * dtinv;
      WALBERLA_LOG_DETAIL( "gdot_n2' = " << gdot_nto2[0] << "\n");
#endif

      /*
       * compare DEM time-step with NSCD iteration:
       * - projections are the same
       * - velocities are the same if we use an explicit Euler discretization for the velocity time integration
       *
      f_cf[0] = -stiffness * contactCache.dist_ - damping_n * gdot_n = -[(stiffness * dt) * contactCache.dist_ * dtinv + damping_n * gdot_n] = -foo * (gdot_n + contactCache.dist_ * dtinv) where foo = stiffness * dt = damping_n;
      f_cf[1] = -damping_t * gdot_t                     = -damping_t * gdot_t;
      f_cf[2] = -damping_t * gdot_o                     = -damping_t * gdot_o;

      or: f_cf = -diag(foo, damping_t, damping_t) * gdot_nto   (since gdot_nto[0] is modified)
      vs. f_cf = -diaginv * gdot_nto in NSCD iteration

      => The NSCD iteration is more or less a DEM time step where we choose the stiffness and damping parameters such that penetration is non-existent after a time step and contacts are truly static (tangential rel. vel. is zero) unless the friction force hits its limit

      f_cf[0] = std::max( 0, f_cf[0] );

      flimit = contactCache.mu_ * f_cf[0];
      fsq = f_cf[1] * f_cf[1] + f_cf[2] * f_cf[2]
      if( fsq > flimit * flimit ) {
         f = flimit / sqrt( fsq );
         f_cf[1] *= f;
         f_cf[2] *= f;
      }

      f_wf = contactframe * f_cf;

      b1->addForceAtPos(  f_wf, gpos );
      b2->addForceAtPos( -f_wf, gpos );
      */
   }

   return delta_max;
}
//*************************************************************************************************




//=================================================================================================
//
//  COMMUNICATION FUNCTIONS
//
//=================================================================================================

inline void HardContactSemiImplicitTimesteppingSolvers::parseVelocityCorrection(mpi::RecvBuffer& rb, BodyStorage& bodyStorage, BodyCache& bodyCache)
{
   //   if( tt_ != NULL ) tt_->start( "parseVelocityCorrection" );
   using namespace walberla::pe::communication;
   NotificationType notificationType;

   // receiving shadow copies [N], shadow copy updates [DN], (local) migrations [N], (remote) migrations [D], deletions [DN] and removal notifications [DN] from neighbors (N) and distant processes (D)
   unmarshal( rb, notificationType );

   switch( notificationType ) {
   case rigidBodyVelocityCorrectionNotification: {
      RigidBodyVelocityCorrectionNotification::Parameters objparam;
      unmarshal( rb, objparam );

      auto bodyIt = bodyStorage.find( objparam.sid_ );
      WALBERLA_ASSERT(bodyIt != bodyStorage.end(), "Body not found!");
      BodyID b( bodyIt.getBodyID() );
      bodyCache.v_[b->index_] += relaxationParam_*objparam.dv_;
      bodyCache.w_[b->index_] += relaxationParam_*objparam.dw_;

      WALBERLA_LOG_DETAIL( "Received rigid body velocity correction:\ndv = " << objparam.dv_ << "\ndw = " << objparam.dw_ << "\nv_total = " << b->getLinearVel() << "\nw_total = " << b->getAngularVel() << "\nbody = " << b << "\n");
      break;
   }
   default:
      throw std::runtime_error( "Received invalid notification type." );
   }
   //   if( tt_ != NULL ) tt_->stop( "parseVelocityCorrection" );
}

inline void HardContactSemiImplicitTimesteppingSolvers::parseVelocityCorrectionShadow(mpi::RecvBuffer& rb, BodyStorage& bodyStorage, BodyCache& bodyCache )
{
   //   if( tt_ != NULL ) tt_->start( "parseVelocityCorrectionShadow" );
   using namespace walberla::pe::communication;
   NotificationType notificationType;

   // receiving shadow copies [N], shadow copy updates [DN], (local) migrations [N], (remote) migrations [D], deletions [DN] and removal notifications [DN] from neighbors (N) and distant processes (D)
   unmarshal( rb, notificationType );

   switch( notificationType ) {
   case rigidBodyVelocityCorrectionNotification: {
      RigidBodyVelocityCorrectionNotification::Parameters objparam;
      unmarshal( rb, objparam );

      auto bodyIt = bodyStorage.find( objparam.sid_ );
      WALBERLA_ASSERT(bodyIt != bodyStorage.end(), "Body not found!");
      BodyID b( bodyIt.getBodyID() );
      WALBERLA_ASSERT( b->isRemote(), "Update notification must only concern shadow copies." );

      bodyCache.v_[b->index_] = objparam.dv_;
      bodyCache.w_[b->index_] = objparam.dw_;

      bodyCache.dv_[b->index_] = Vec3();
      bodyCache.dw_[b->index_] = Vec3();

      WALBERLA_LOG_DETAIL( "Received rigid body velocity correction ("<< b->index_ << "):\ndv = " << objparam.dv_ << "\ndw = " << objparam.dw_ << "\nv_total = " << b->getLinearVel() << "\nw_total = " << b->getAngularVel() << "\nbody = " << b << "\n");
      break;
   }
   default:
      throw std::runtime_error( "Received invalid notification type." );
   }

   //   if( tt_ != NULL ) tt_->stop( "parseVelocityCorrectionShadow" );
}


//*************************************************************************************************
/*!\brief TODO
 *
 * \return void
 */
inline void HardContactSemiImplicitTimesteppingSolvers::synchronizeVelocities( )
{
   if ((mpi::MPIManager::instance()->numProcesses() <= 1) && (blockStorage_->size() <= 1))
      return;

   if (tt_ != NULL) tt_->start("Velocity Sync");
   if (tt_ != NULL) tt_->start("Velocity Sync Correction Assembling");

   // Sending local force contributions of shadow copies to owner.
   WALBERLA_LOG_DETAIL( "Assembling of velocity correction message starts...\n" );

   using namespace walberla::pe::communication;
   //==========================================================
   // STEP1: Send velocities of shadow copies to owner
   //==========================================================

   mpi::BufferSystem syncVelBS( mpi::MPIManager::instance()->comm(),  256);
   std::set<mpi::MPIRank> recvRanks; // potential message senders
   for (auto blockIt = blockStorage_->begin(); blockIt != blockStorage_->end(); ++blockIt)
   {
      IBlock& currentBlock = *blockIt;
      Storage* storage  = currentBlock.uncheckedFastGetData< Storage >( storageID_ );
      BodyStorage& localStorage  = (*storage)[0];
      BodyStorage& shadowStorage = (*storage)[1];

      BodyCache&    bodyCache    = blockToBodyCache_[currentBlock.getId().getID()];

      const Owner me( int_c( mpi::MPIManager::instance()->rank() ), currentBlock.getId().getID() );

      size_t i = 0;
      for( auto body = shadowStorage.begin(); body != shadowStorage.end(); ++body )
      {
         i = body->index_;

         mpi::SendBuffer& sb = syncVelBS.sendBuffer( body->MPITrait.getOwner().rank_ );
         if (sb.isEmpty()) sb << walberla::uint8_c(0);

         if( bodyCache.dv_[i] == Vec3() && bodyCache.dw_[i] == Vec3() ) {
            // If we did not apply any corrections do not send anything.
            continue;
         }

         WALBERLA_LOG_DETAIL( "Sending velocity correction " << bodyCache.dv_[i] << ", " << bodyCache.dw_[i] << " of body " << body->getSystemID() << " to owner process " << body->MPITrait.getOwner() << ".");

         packNotificationWithoutSender(body->MPITrait.getOwner(), sb, RigidBodyVelocityCorrectionNotification( *body, bodyCache.dv_[i], bodyCache.dw_[i] ));
      }

      for( auto bodyIt = localStorage.begin(); bodyIt != localStorage.end(); ++bodyIt )
      {
         BodyID b(bodyIt.getBodyID());
         for (auto shadows = b->MPITrait.beginShadowOwners(); shadows != b->MPITrait.endShadowOwners(); ++shadows )
         {
            recvRanks.insert(shadows->rank_);
         }
      }
   }

   if (tt_ != NULL) tt_->stop("Velocity Sync Correction Assembling");

   //   for( ProcessIterator process = processstorage_.begin(); process != processstorage_.end(); ++process )
   //      sentVelocitiesSyncCorrections_.transfered( process->getSendBuffer().size() );

   if (tt_ != NULL) tt_->start("Velocity Sync Correction Communicate");

   WALBERLA_LOG_DETAIL( "Communication of velocity correction message starts..." );

   //   size_t sum = bs.size();
   //   mpi::reduceInplace(sum, mpi::SUM);
   //   WALBERLA_LOG_DEVEL_ON_ROOT("communication size: " << sum);
   syncVelBS.setReceiverInfo(recvRanks, true);
   syncVelBS.sendAll();

   if (tt_ != NULL) tt_->stop("Velocity Sync Correction Communicate");

   //   for( ProcessIterator process = processstorage_.begin(); process != processstorage_.end(); ++process )
   //      receivedVelocitiesSyncCorrections_.transfered( process->getRecvBuffer().size() );

   if (tt_ != NULL) tt_->start("Velocity Sync Correction Parsing");

   // Receiving force and torque contributions
   WALBERLA_LOG_DETAIL( "Parsing of velocity correction message starts...");

   for( auto it = syncVelBS.begin(); it != syncVelBS.end(); ++it )
   {
      //      if (tt_ != NULL) tt_->start("Inside Loop");
      walberla::uint8_t tmp;
      it.buffer() >> tmp;
      while( !it.buffer().isEmpty() )
      {
         //         IBlockID::IDType sender;
         IBlockID::IDType receiver;
         //         it.buffer() >> sender;
         it.buffer() >> receiver;
         auto blk = blockStorage_->getBlock(receiver);
         WALBERLA_CHECK(blk != NULL, receiver << " not on this process!");
         IBlock& block = *blk;
         Storage* storage  = block.uncheckedFastGetData< Storage >( storageID_ );
         BodyStorage& localStorage  = (*storage)[0];
         //         BodyStorage& shadowStorage = (*storage)[1];
         BodyCache&    bodyCache    = blockToBodyCache_[block.getId().getID()];
         parseVelocityCorrection(it.buffer(), localStorage, bodyCache);
      }
      //      if (tt_ != NULL) tt_->stop("Inside Loop");
   }

   if (tt_ != NULL) tt_->stop("Velocity Sync Correction Parsing");

   if (tt_ != NULL) tt_->start("Velocity Sync Update Assembling");
   WALBERLA_LOG_DETAIL( "Assembling of velocity update message starts...");

   //==========================================================
   // STEP2: Update velocities of remote shadow copies
   //==========================================================

   recvRanks.clear();

   for (auto blockIt = blockStorage_->begin(); blockIt != blockStorage_->end(); ++blockIt)
   {
      IBlock& block = *blockIt;
      Storage* storage  = block.uncheckedFastGetData< Storage >( storageID_ );
      BodyStorage& localStorage  = (*storage)[0];
      BodyStorage& shadowStorage = (*storage)[1];

      BodyCache&    bodyCache    = blockToBodyCache_[block.getId().getID()];

      const Owner me( int_c( mpi::MPIManager::instance()->rank() ), block.getId().getID() );

      size_t i = 0;
      for( auto body = localStorage.begin(); body != localStorage.end(); ++body ) {
         if( body->isGlobal() )
            continue;

         i = body->index_;

         bodyCache.v_[i] += relaxationParam_*bodyCache.dv_[i];
         bodyCache.w_[i] += relaxationParam_*bodyCache.dw_[i];
         bodyCache.dv_[i] = Vec3();
         bodyCache.dw_[i] = Vec3();

         for( auto shadow = body->MPITrait.beginShadowOwners(); shadow != body->MPITrait.endShadowOwners(); ++shadow ) {

            mpi::SendBuffer& sb = syncVelBS.sendBuffer( shadow->rank_ );
            if (sb.isEmpty()) sb << walberla::uint8_c(0);
            packNotificationWithoutSender(*shadow, sb, RigidBodyVelocityCorrectionNotification( *body, bodyCache.v_[i], bodyCache.w_[i] ));

            WALBERLA_LOG_DETAIL( "Sending velocity update " << bodyCache.v_[i] << ", " << bodyCache.w_[i] << " of body " << body->getSystemID() << " to process " << *shadow << " having a shadow copy.");
         }
      }

      for( auto bodyIt = shadowStorage.begin(); bodyIt != shadowStorage.end(); ++bodyIt )
      {
         BodyID b(bodyIt.getBodyID());
         recvRanks.insert(b->MPITrait.getOwner().rank_);
      }
   }

   if (tt_ != NULL) tt_->stop("Velocity Sync Update Assembling");

   //   for( ProcessIterator process = processstorage_.begin(); process != processstorage_.end(); ++process )
   //      sentVelocitiesSyncUpdates_.transfered( process->getSendBuffer().size() );

   if (tt_ != NULL) tt_->start("Velocity Sync Update Communincate");

   WALBERLA_LOG_DETAIL( "Communication of velocity update message starts...");

   syncVelBS.setReceiverInfo(recvRanks, true);
   syncVelBS.sendAll();

   if (tt_ != NULL) tt_->stop("Velocity Sync Update Communincate");

   //   for( ProcessIterator process = processstorage_.begin(); process != processstorage_.end(); ++process )
   //      receivedVelocitiesSyncUpdates_.transfered( process->getRecvBuffer().size() );

   if (tt_ != NULL) tt_->start("Velocity Sync Update Processing");

   // Receiving velocity updates
   WALBERLA_LOG_DETAIL( "Parsing of velocity update message starts...");

   for( auto it = syncVelBS.begin(); it != syncVelBS.end(); ++it )
   {
      //      if (tt_ != NULL) tt_->start("Inside Loop");
      walberla::uint8_t tmp;
      it.buffer() >> tmp;
      while( !it.buffer().isEmpty() )
      {
         //         IBlockID::IDType sender;
         IBlockID::IDType receiver;
         //         it.buffer() >> sender;
         it.buffer() >> receiver;
         auto blk = blockStorage_->getBlock(receiver);
         WALBERLA_CHECK(blk != NULL, receiver << " not on this process!");
         IBlock& block = *blk;
         Storage* storage  = block.uncheckedFastGetData< Storage >( storageID_ );
         //         BodyStorage& localStorage  = (*storage)[0];
         BodyStorage& shadowStorage = (*storage)[1];
         BodyCache&    bodyCache    = blockToBodyCache_[block.getId().getID()];
         parseVelocityCorrectionShadow(it.buffer(), shadowStorage, bodyCache);
      }
      //      if (tt_ != NULL) tt_->stop("Inside Loop");
   }

   if (tt_ != NULL) tt_->stop("Velocity Sync Update Processing");

   if (tt_ != NULL) tt_->start("Velocity Sync Globals");
   /*
   {
      size_t i;
      std::vector<real_t> reductionBuffer( globalBodyStorage_->size() * 6, real_c(0) );

      for (auto blockIt = blockStorage_->begin(); blockIt != blockStorage_->end(); ++blockIt)
      {
         IBlock& block = *blockIt;
         BodyCache&    bodyCache    = blockToBodyCache_[block.getId().getID()];
         i = 0;
         for( auto it = globalBodyStorage_->begin(); it != globalBodyStorage_->end(); ++it )
         {
            if (it->hasInfiniteMass()) continue;
            reductionBuffer[i++] += bodyCache.dv_[(*it)->index_][0];
            reductionBuffer[i++] += bodyCache.dv_[(*it)->index_][1];
            reductionBuffer[i++] += bodyCache.dv_[(*it)->index_][2];
            reductionBuffer[i++] += bodyCache.dw_[(*it)->index_][0];
            reductionBuffer[i++] += bodyCache.dw_[(*it)->index_][1];
            reductionBuffer[i++] += bodyCache.dw_[(*it)->index_][2];
         }
      }

      mpi::allReduceInplace(reductionBuffer, mpi::SUM);

      for (auto blockIt = blockStorage_->begin(); blockIt != blockStorage_->end(); ++blockIt)
      {
         IBlock& block = *blockIt;
         BodyCache&    bodyCache    = blockToBodyCache_[block.getId().getID()];
         i = 0;
         for( auto it = globalBodyStorage_->begin(); it != globalBodyStorage_->end(); ++it )
         {
            if (it->hasInfiniteMass()) continue;
            bodyCache.v_[(*it)->index_] += Vec3( reductionBuffer[i    ], reductionBuffer[i + 1], reductionBuffer[i + 2] );
            bodyCache.w_[(*it)->index_] += Vec3( reductionBuffer[i + 3], reductionBuffer[i + 4], reductionBuffer[i + 5] );
         }
      }
   }*/

   if (tt_ != NULL) tt_->stop("Velocity Sync Globals");
   if (tt_ != NULL) tt_->stop("Velocity Sync");
}
//*************************************************************************************************

//=================================================================================================
//
//  TIME-INTEGRATION FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Calculates the initial velocity corrections of a given body.
 *
 * \param body The body whose velocities to time integrate
 * \param dv On return the initial linear velocity correction.
 * \param w On return the initial angular velocity correction.
 * \param dt The time step size.
 * \return void
 *
 * Calculates the velocity corrections effected by external forces and torques in an explicit Euler
 * time integration of the velocities of the given body. For fixed objects the velocity corrections
 * are set to zero. External forces and torques are reset if indicated by the settings.
 */
inline void HardContactSemiImplicitTimesteppingSolvers::initializeVelocityCorrections( BodyID body, Vec3& dv, Vec3& dw, real_t dt ) const
{
   if( body->isAwake() ) {
      if( !body->hasInfiniteMass() )
      {
         dv = ( body->getInvMass() * dt ) * body->getForce();
         dw = dt * ( body->getInvInertia() * body->getTorque() );
      } else
      {
         dv = Vec3();
         dw = Vec3();
      }
   }
   else {
      WALBERLA_ASSERT( body->getLinearVel() == Vec3( 0, 0, 0 ) && body->getAngularVel() == Vec3( 0, 0, 0 ), "Sleeping body has non-zero velocities." );
      WALBERLA_ASSERT( body->getForce() == Vec3( 0, 0, 0 ) && body->getTorque() == Vec3( 0, 0, 0 ), "Sleeping body has non-zero forces or torques." );
   }

   body->resetForceAndTorque();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Time integration of the position and orientation of a given body.
 *
 * \param body The body whose position and orientation to time integrate
 * \param v The linear velocity to use for time integration of the position.
 * \param w The angular velocity to use for time integration of the orientation.
 * \param dt The time step size.
 * \return void
 *
 * Performs an Euler time integration of the positions of the given body. Velocities are damped if
 * indicated by the settings and stored back in the body properties. The bounding box is
 * recalculated and it is redetermined whether the body is awake or not. Also the data
 * structure tracking the contacts attached to the body are cleared and
 */
inline void HardContactSemiImplicitTimesteppingSolvers::integratePositions( BodyID body, Vec3 v, Vec3 w, real_t dt ) const
{
   if ( body->isFixed() )
      return;

   if( body->isAwake() )
   {
      if ( isSpeedLimiterActive() )
      {
         const auto speed = v.length();
         const auto aabb  = body->getAABB();
         const auto edge  = std::min(aabb.xSize(), std::min(aabb.ySize(), aabb.zSize()));
         if (speed * dt > edge * getSpeedLimitFactor() )
         {
            v = v * (edge * getSpeedLimitFactor() / dt / speed );
         }

         const real_t maxPhi = real_t(2) * math::pi * getSpeedLimitFactor();
         const real_t phi    = w.length() * dt;
         if (phi > maxPhi)
         {
            w = w / phi * maxPhi;
         }
      }

      // Calculating the translational displacement
      body->setPosition( body->getPosition() + v * dt );

      // Calculating the rotation angle
      const Vec3 phi( w * dt );

      // Calculating the new orientation
      const Quat dq( phi, phi.length() );
      WALBERLA_ASSERT(!math::isnan(dq), "dq: " << dq << " phi: " << phi << " w: " << w << " dt: " << dt << " angVel: " << body->getAngularVel());
      body->setOrientation( dq * body->getQuaternion() );

      // Storing the velocities back in the body properties
      body->setLinearVel( v );
      body->setAngularVel( w );

      // Setting the axis-aligned bounding box
      body->calcBoundingBox();

      // Calculating the current motion of the box
      body->calcMotion();
   }
}
//*************************************************************************************************

} // namespace cr
} // namespace pe

inline
void configure( const Config::BlockHandle& config, pe::cr::HCSITS& cr )
{
   using namespace pe;

   int HCSITSmaxIterations = config.getParameter<int>("HCSITSmaxIterations", 10);
   WALBERLA_LOG_INFO_ON_ROOT("HCSITSmaxIterations: " << HCSITSmaxIterations);

   real_t HCSITSRelaxationParameter = config.getParameter<real_t>("HCSITSRelaxationParameter", real_t(0.75) );
   WALBERLA_LOG_INFO_ON_ROOT("HCSITSRelaxationParameter: " << HCSITSRelaxationParameter);

   real_t HCSITSErrorReductionParameter = config.getParameter<real_t>("HCSITSErrorReductionParameter", real_t(0.8) );
   WALBERLA_LOG_INFO_ON_ROOT("HCSITSErrorReductionParameter: " << HCSITSErrorReductionParameter);

   std::string HCSITSRelaxationModelStr = config.getParameter<std::string>("HCSITSRelaxationModelStr", "ApproximateInelasticCoulombContactByDecoupling" );
   WALBERLA_LOG_INFO_ON_ROOT("HCSITSRelaxationModelStr: " << HCSITSRelaxationModelStr);

   cr::HCSITS::RelaxationModel HCSITSRelaxationModel;
   if (HCSITSRelaxationModelStr == "InelasticFrictionlessContact")
   {
      HCSITSRelaxationModel = cr::HCSITS::InelasticFrictionlessContact;
   } else if (HCSITSRelaxationModelStr == "ApproximateInelasticCoulombContactByDecoupling")
   {
      HCSITSRelaxationModel = cr::HCSITS::ApproximateInelasticCoulombContactByDecoupling;
   } else if (HCSITSRelaxationModelStr == "InelasticCoulombContactByDecoupling")
   {
      HCSITSRelaxationModel = cr::HCSITS::InelasticCoulombContactByDecoupling;
   } else if (HCSITSRelaxationModelStr == "InelasticGeneralizedMaximumDissipationContact")
   {
      HCSITSRelaxationModel = cr::HCSITS::InelasticGeneralizedMaximumDissipationContact;
   } else
   {
      WALBERLA_ABORT("Unknown HCSITSRelaxationModel: " << HCSITSRelaxationModelStr);
   }

   Vec3 globalLinearAcceleration = config.getParameter<Vec3>("globalLinearAcceleration", Vec3(0, 0, 0));
   WALBERLA_LOG_INFO_ON_ROOT("globalLinearAcceleration: " << globalLinearAcceleration);

   cr.setMaxIterations( uint_c(HCSITSmaxIterations) );
   cr.setRelaxationModel( HCSITSRelaxationModel );
   cr.setRelaxationParameter( HCSITSRelaxationParameter );
   cr.setErrorReductionParameter( HCSITSErrorReductionParameter );
   cr.setGlobalLinearAcceleration( globalLinearAcceleration );
}

} // namespace walberla

